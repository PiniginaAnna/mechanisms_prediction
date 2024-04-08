import re
from indigo import Indigo, IndigoException
from tqdm import tqdm
import gc


class Node:
    def __init__(self, molecules_smiles, rule_smarts=None, parent_id=None, level=0):
        self.molecules_smiles = molecules_smiles
        self.rules_smarts = [rule_smarts]
        self.parents = {parent_id}
        self.children = set()
        self.level = level
        self.is_product = False
        # TODO indigo change logic
        # self.indigo = None

    def create_specific_rules(self, indigo, molecules, rules_smarts, pseudoatoms_dict):
        def pseudoatom_replacement(rule_smarts, pseudoatom_name, pseudoatoms_list):
            specific_rules_smarts_part = []
            for pseudoatom in pseudoatoms_list:
                if indigo.substructureMatcher(molecules).match(indigo.loadSmarts(pseudoatom)):
                    all_atom = re.search(r'\[.+?\]', pseudoatom)[0][:-1]
                    atom_name = re.search(r'(\w|#)\w?', pseudoatom)[0]
                    other_group = re.search(r'\(.+\)', pseudoatom)
                    if other_group is None:
                        other_group = ''
                    else:
                        other_group = other_group[0]
                    specific_rule_smarts = re.sub(fr'\[{pseudoatom_name[:-2]}(.|..):(\d)\]',
                                                  fr'{all_atom}:\2]{other_group}', rule_smarts, count=1)
                    specific_rule_smarts = re.sub(fr'\[{pseudoatom_name[:-2]}(.|..):(\d)\]',
                                                  fr'[{atom_name}\1:\2]{other_group}', specific_rule_smarts, count=1)
                    specific_rules_smarts_part.append(specific_rule_smarts)
            return specific_rules_smarts_part

        def create_one_stage_rules(rules_smarts):
            nonlocal pseudoatoms_dict
            specific_rules_smarts = set()
            for rule_smarts in rules_smarts:
                for pseudoatom_name, pseudoatoms_list in pseudoatoms_dict.items():
                    if pseudoatom_name in rule_smarts.split('>>')[0]:
                        specific_rules_smarts = specific_rules_smarts.union(
                            pseudoatom_replacement(rule_smarts, pseudoatom_name, pseudoatoms_list))
                        break
                else:
                    specific_rules_smarts.add(rule_smarts)
            return list(specific_rules_smarts)

        # molecules = indigo.loadMolecule(self.molecules_smiles)
        specific_rules_smarts_first_stage = create_one_stage_rules(rules_smarts)
        specific_rules_smarts_second_stage = create_one_stage_rules(specific_rules_smarts_first_stage)

        return specific_rules_smarts_second_stage


class Graph:
    def __init__(self, root_molecules_smiles, products_smiles, limit_depth, limit_nodes_count, stop_when_product=False):
        self.nodes = {}
        self.current_level_nodes = {}
        self.new_level_nodes = {}
        self.root_molecules_smiles = root_molecules_smiles
        self.products_smiles = products_smiles
        self.limit_depth = limit_depth
        self.limit_nodes_count = limit_nodes_count
        self.stop_when_product = stop_when_product
        self.stop_flag = False

    def add_node(self, node_id, molecules_smiles, rule_smarts, parent_id, level):
        self.new_level_nodes[node_id] = Node(molecules_smiles, rule_smarts, parent_id, level)
        self.current_level_nodes[parent_id].children.add(node_id)

    def update_node(self, node_id, node, rule_smarts, parent_id):
        node.rules_smarts.append(rule_smarts)
        node.parents.add(parent_id)
        self.current_level_nodes[parent_id].children.add(node_id)

    @staticmethod
    def is_product(indigo, molecules_smiles, products_smiles):
        if '.' in products_smiles:
            all_products_smiles = products_smiles.split('.')
        else:
            all_products_smiles = [products_smiles]

        if '.' in molecules_smiles:
            all_molecules_smiles = molecules_smiles.split('.')
        else:
            all_molecules_smiles = [molecules_smiles]
        if all(any(indigo.exactMatch(indigo.loadMolecule(product_smiles), indigo.loadMolecule(molecule_smiles))
                   for molecule_smiles in all_molecules_smiles) for product_smiles in all_products_smiles):
            return True
        return False

    def grow_graph(self, rules_smarts, pseudoatoms_dict):
        self.current_level_nodes[0] = Node(self.root_molecules_smiles)
        reactants_charge_count = self.root_molecules_smiles.count('+') + self.root_molecules_smiles.count('-')
        n = 1

        for i in tqdm(range(self.limit_depth)):
            if self.stop_flag:
                self.nodes.update(self.current_level_nodes)
                return
            for node_id in list(self.current_level_nodes):
                node = self.current_level_nodes[node_id]
                if not node.is_product:
                    indigo = Indigo()
                    indigo.setOption("rpe-mode", "one-tube")
                    indigo.setOption("rpe-self-reaction", "1")

                    parent_id = node_id
                    parent_path_molecules = self.path_molecules(parent_id)
                    molecules = indigo.loadMolecule(node.molecules_smiles)
                    molecules.unfoldHydrogens()
                    molecules_atoms_count = molecules.countAtoms()
                    specific_rules_smarts = node.create_specific_rules(indigo, molecules, rules_smarts,
                                                                       pseudoatoms_dict)

                    monomers_table = indigo.createArray()
                    monomers_table.arrayAdd(indigo.createArray())
                    monomers_table.arrayAdd(indigo.createArray())
                    monomers_table.arrayAdd(indigo.createArray())
                    monomers_table.arrayAdd(indigo.createArray())
                    monomers_table.at(0).arrayAdd(molecules)

                    for specific_rule_smarts in specific_rules_smarts:
                        try:
                            rule = indigo.loadQueryReaction(specific_rule_smarts)
                            output_reactions = indigo.reactionProductEnumerate(rule, monomers_table)
                            for reaction in output_reactions.iterateArray():
                                if reaction.countReactants() == 1:  # miss 2 monomers
                                    reaction.unfoldHydrogens()
                                    for products in reaction.iterateProducts():
                                        if products.countAtoms() == molecules_atoms_count:  # miss added hydrogens
                                            products_smiles = products.canonicalSmiles()
                                            # miss loops
                                            if not any(products_smiles == path_molecules_smiles for
                                                       path_molecules_smiles in parent_path_molecules):
                                                # check penalty
                                                if products_smiles.count('+') + products_smiles.count('-') <= \
                                                        reactants_charge_count + 4:
                                                    # check charged carbons pairs
                                                    if not self.check_charged_atoms_pairs(products):
                                                        for j, new_node in self.new_level_nodes.items():
                                                            if products_smiles == new_node.molecules_smiles:
                                                                self.update_node(j, new_node, specific_rule_smarts,
                                                                                 parent_id)
                                                                break
                                                        else:
                                                            self.add_node(n, products_smiles, specific_rule_smarts,
                                                                          parent_id, level=i + 1)
                                                            is_product = self.is_product(indigo, products_smiles,
                                                                                         self.products_smiles)
                                                            self.new_level_nodes[n].is_product = is_product
                                                            if self.stop_when_product and is_product:
                                                                self.stop_flag = True
                                                            n += 1
                                                        # check limit nodes count
                                                        if n == self.limit_nodes_count:
                                                            del indigo
                                                            gc.collect()
                                                            self.nodes.update(self.current_level_nodes)
                                                            self.nodes.update(self.new_level_nodes)
                                                            return
                        except IndigoException as exc:
                            with open(f'error_rules_{self.__str__()}.txt', 'a') as file:
                                print(exc, file=file)
                                file.write('\n' + node.molecules_smiles + '\n')
                                file.write(specific_rule_smarts + '\n' + '\n')
                    del indigo
                    gc.collect()
            self.nodes.update(self.current_level_nodes)
            self.current_level_nodes = self.new_level_nodes
            self.new_level_nodes = {}
        self.nodes.update(self.current_level_nodes)

    @staticmethod
    def check_charged_atoms_pairs(molecule):
        for atom in molecule.iterateAtoms():
            if atom.charge() != 0:
                for neighbour_atom in atom.iterateNeighbors():
                    if neighbour_atom.charge() != 0:
                        if neighbour_atom.charge() == atom.charge():
                            if (atom.symbol() == 'N' or neighbour_atom.symbol() == 'N') and atom.charge() == 1:
                                continue
                            else:
                                return True
                        if atom.symbol() == 'C' and neighbour_atom.symbol() == 'C':
                            return True
        return False

    def paths_starts(self):
        paths_starts = []
        for i, node in self.nodes.items():
            if node.is_product:
                paths_starts.append(i)
        return paths_starts

    def level_nodes_count(self, level):
        n = 0
        for node in self.nodes.values():
            if node.level == level:
                n += 1
        return n

    def print_graph(self):
        for node_id, node in self.nodes.items():
            print(node.level, node.parents, node_id, node.children)

    def path_molecules(self, i):
        path_molecules = [self.current_level_nodes[i].molecules_smiles]
        new_parents = self.current_level_nodes[i].parents

        def add_parents(parents, path_molecules):
            new_parents = set()
            for parent in parents:
                if parent is None:
                    return path_molecules
                path_molecules.append(self.nodes[parent].molecules_smiles)
                new_parents.update(self.nodes[parent].parents)
            return add_parents(new_parents, path_molecules)

        return add_parents(new_parents, path_molecules)

import re
from indigo import Indigo, IndigoException
from tqdm import tqdm
import gc
# import ray
from memory_profiler import profile


class Node:
    def __init__(self, molecules_smiles, rule_smarts=None, parent_id=None, level=0, penalty=0):
        self.molecules_smiles = molecules_smiles
        self.rule_smarts = rule_smarts
        self.parent_id = parent_id
        self.level = level
        # TODO remove penalty
        self.penalty = penalty
        self.is_product = False
        self.children = []

    def create_specific_rules(self, rules_smarts, pseudoatoms_dict):
        def pseudoatom_replacement(rule_smarts, molecules_smiles, pseudoatom_name, pseudoatoms_list):
            specific_rules_smarts_part = []

            return specific_rules_smarts_part

        def create_one_stage_rules(rules_smarts):
            nonlocal pseudoatoms_dict
            specific_rules_smarts = set()
            for rule_smarts in rules_smarts:
                for pseudoatom_name, pseudoatoms_list in pseudoatoms_dict.items():
                    if pseudoatom_name in rule_smarts.split('>>')[0]:
                        specific_rules_smarts = specific_rules_smarts.union(
                            pseudoatom_replacement(rule_smarts, self.molecules_smiles, pseudoatom_name,
                                                   pseudoatoms_list))
                        break
                else:
                    specific_rules_smarts.add(rule_smarts)
            return list(specific_rules_smarts)

        specific_rules_smarts_first_stage = create_one_stage_rules(rules_smarts)
        specific_rules_smarts_second_stage = create_one_stage_rules(specific_rules_smarts_first_stage)
        
        return specific_rules_smarts_second_stage

# @profile
class Tree:
    def __init__(self, root_molecules_smiles, products_smiles, limit_depth, limit_nodes_count):
        self.nodes = {}
        self.root_molecules_smiles = root_molecules_smiles
        self.products_smiles = products_smiles
        self.limit_depth = limit_depth
        self.limit_nodes_count = limit_nodes_count

    def add_node(self, node_id, molecules_smiles, rule_smarts, parent_id, level, penalty):
        self.nodes[node_id] = Node(molecules_smiles, rule_smarts, parent_id, level, penalty)
        self.nodes[node_id].is_product = self.is_product(molecules_smiles, self.products_smiles)
        self.nodes[parent_id].children.append(node_id)

    @staticmethod
    def is_product(molecules_smiles, products_smiles):
        if '.' in products_smiles:
            all_products_smiles = products_smiles.split('.')
        else:
            all_products_smiles = [products_smiles]

        if '.' in molecules_smiles:
            all_molecules_smiles = molecules_smiles.split('.')
        else:
            all_molecules_smiles = [molecules_smiles]

        return False

    def grow_tree(self, rules_smarts, pseudoatoms_dict):
        self.nodes[0] = Node(self.root_molecules_smiles)
        reactants_charge_count = self.root_molecules_smiles.count('+') + self.root_molecules_smiles.count('-')

        for i in tqdm(range(self.limit_depth)):
            for node_id in list(self.nodes):
                node = self.nodes[node_id]
                if node.level == i and not node.is_product:
                    parent_id = node_id
                    import random
                    # Задаем размер массива в мегабайтах
                    array_size_mb = 50
                    # Вычисляем размер массива в байтах
                    array_size_bytes = array_size_mb * 1024 * 1024
                    # Создаем массив и заполняем его случайными значениями
                    array = bytearray(random.getrandbits(8) for _ in range(array_size_bytes))
                    self.add_node(len(self.nodes), 'cccccccccccc',
                                  'cccccccc>>ccccc', parent_id, level=i + 1,
                                  penalty=0)
                    # check limit nodes count
                    if len(self.nodes) == self.limit_nodes_count:
                        return
                    gc.collect()

    # def grow_tree_parallel(self, rules_smarts, pseudoatoms_dict, num_cpus):
    #     # TODO doesn't work
    #     @ray.remote(num_cpus=num_cpus)
    #     def process_node(node, node_id):
    #         indigo = Indigo()
    #         indigo.setOption("rpe-mode", "one-tube")
    #         indigo.setOption("rpe-self-reaction", "1")
    #
    #         molecules = indigo.loadMolecule(node.molecules_smiles)
    #         parent_id = node_id
    #         specific_rules_smarts = node.create_specific_rules(indigo, rules_smarts, pseudoatoms_dict)
    #         result_nodes = []
    #
    #         for specific_rule_smarts in specific_rules_smarts:
    #             try:
    #                 rule = indigo.loadQueryReaction(specific_rule_smarts)
    #                 monomers_table = indigo.createArray()
    #                 monomers_table.arrayAdd(indigo.createArray())
    #                 monomers_table.arrayAdd(indigo.createArray())
    #                 monomers_table.at(0).arrayAdd(molecules)
    #                 output_reactions = indigo.reactionProductEnumerate(rule, monomers_table)
    #                 for reaction in output_reactions.iterateArray():
    #                     if reaction.countReactants() == 1:  # miss 2 monomers
    #                         reaction.unfoldHydrogens()
    #                         for products in reaction.iterateProducts():
    #                             if products.countAtoms() == molecules.countAtoms():  # miss added hydrogens
    #                                 products_smiles = products.smiles()
    #                                 # miss loops
    #                                 if not any(indigo.exactMatch(products,
    #                                                              indigo.loadMolecule(path_molecules_smiles))
    #                                            for path_molecules_smiles in
    #                                            self.mechanism_path_molecules(parent_id)):
    #                                     if products_smiles.count('+') + products_smiles.count(
    #                                             '-') <= reactants_charge_count + 4:
    #                                         # check charged carbons pairs
    #                                         if not self.is_carbion_pair_charged(products):
    #                                             result_nodes.append((products_smiles, specific_rule_smarts, parent_id))
    #             except IndigoException as exc:
    #                 with open('error_rules.txt', 'a') as file:
    #                     print(exc, file=file)
    #                     file.write('\n' + node.molecules_smiles + '\n')
    #                     file.write(specific_rule_smarts + '\n' + '\n')
    #         gc.collect()
    #         return result_nodes
    #
    #     reactants_charge_count = self.nodes[0].molecules_smiles.count('+') + self.nodes[0].molecules_smiles.count('-')
    #
    #     for i in tqdm(range(self.limit_depth)):
    #         print(i)
    #         level_node_ids = [node_id for node_id in self.nodes.keys() if self.nodes[node_id].level == i]
    #         results = ray.get([process_node.remote(self.nodes[node_id], node_id) for node_id in level_node_ids])
    #         for result in results:
    #             for products_smiles, specific_rule_smarts, parent_id in result:
    #                 self.add_node(len(self.nodes), products_smiles,
    #                               specific_rule_smarts, parent_id, level=i + 1,
    #                               penalty=0)
    #                 # check limit nodes count
    #                 if len(self.nodes) == self.limit_nodes_count:
    #                     return

    # @staticmethod
    # def calculate_penalty(penalty, parent_molecule_smiles, child_molecule_smiles):
    #     parent_charge_count = parent_molecule_smiles.count('+') + parent_molecule_smiles.count('-')
    #     child_charge_count = child_molecule_smiles.count('+') + child_molecule_smiles.count('-')
    #     penalty += child_charge_count - parent_charge_count
    #     if penalty < 0:
    #         penalty = 0
    #     return penalty
    
    @staticmethod
    def is_carbion_pair_charged(molecule):
        for atom in molecule.iterateAtoms():
            if atom.charge() != 0 and atom.symbol() == 'C':
                for neighbour_atom in atom.iterateNeighbors():
                    if neighbour_atom.charge() != 0:
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

    def print_tree(self):
        for node_id, node in self.nodes.items():
            print(node.level, node.parent_id, node_id, node.children)

    def mechanism_path_molecules(self, i):
        mechanism_path_molecules = []
        while i is not None:
            mechanism_path_molecules.append(self.nodes[i].molecules_smiles)
            i = self.nodes[i].parent_id
        return mechanism_path_molecules


# if there is no charge in pseudoatom name it is necessary to add one symbol at the end
pseudoatoms_dict = {'Nu-': ['[O-;X1]', '[Cl-;X0]', '[N-;X2]', '[O-;X1]', '[C-;X1,X3]'],
                    'Nu:': ['[O+0,S+0;X1]', '[N+0;X3]', '[P+0;X3]', '[O+0;X2]'],
                    'E+': ['[C+;X3,X2]'],
                    'E:': ['C'],
                    'L:': ['[O+0;X2](P(=O)(Cl)Cl)', '[Cl+0;X1]', '[O+0;X2]([#1])', '[N+0;X3]([#1])([#1])'],
                    'L+': ['[O+;X3]([#1])([#1])'],
                    'L-': [],
                    'Bs-': ['[O-;X1]', '[C-;X3]'],
                    'Bs:': ['[O+0,X2]([#1])([#1])', '[N+0;X2,X3]', '[O+0;X2]'],
                    'LA:': [],
                    'LA-:': [],
                    'LB:': [],
                    'EWG:': [],
                    'EDG:': []}

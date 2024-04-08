# тестовые молекулы 10 первых + ограничения на дерево + правила 3 + 10 процесса

from mechanisms_tree import Tree
from pseudoatoms_dict import pseudoatoms_dict
from tqdm import tqdm
from pickle import dump
from multiprocessing import Pool
from CGRtools import RDFRead


def create_tree(reaction):
    with open('rules/rules3.smarts') as file:
        rules_smarts = file.read().split('\n')
    reaction_name = reaction.meta['r_c_1_freq']
    reaction.kekule()
    reaction.explicify_hydrogens()
    molecules = str(reaction).split('>')
    root_molecules_smiles = molecules[0]
    if molecules[1]:
        root_molecules_smiles = root_molecules_smiles + '.' + molecules[1]
    products_smiles = molecules[2]
    tree = Tree(root_molecules_smiles, products_smiles, limit_depth=15, limit_nodes_count=1000000,
                stop_when_product=True)
    tree.grow_tree(rules_smarts, pseudoatoms_dict)
    paths_starts = tree.paths_starts()
    return (reaction_name, tree, paths_starts)


with RDFRead('data/72_test_reactions_with_H_rc_0_with_reag.rdf', indexable=True) as file:
    reactions = file[:10]

with Pool(processes=10) as pool:
    for result in tqdm(pool.imap_unordered(create_tree, reactions, chunksize=1)):
        reaction_name, tree, paths_starts = result
        with open(f'results/trees/tree_freq_{reaction_name}.pkl', 'wb') as file:
            dump(tree, file)

        with open(f'results/trees/paths_starts_tree_freq_{reaction_name}.txt', 'w') as file:
            for j in paths_starts:
                file.write(str(j) + '\n')

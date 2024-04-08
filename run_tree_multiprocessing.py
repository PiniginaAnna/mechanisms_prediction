# тестовые молекулы 13-14 + старт 13 + ограничения на дерево + правила 3 + 2 процесса

import pandas as pd
from mechanisms_tree import Tree
from pseudoatoms_dict import pseudoatoms_dict
from tqdm import tqdm
from pickle import dump
from multiprocessing import Pool


def create_tree(reaction):
    with open('rules/rules3.smarts') as file:
        rules_smarts = file.read().split('\n')
    root_molecules_smiles = reaction[1]['SMILES_react']
    products_smiles = reaction[1]['SMILES_prod']
    tree = Tree(root_molecules_smiles, products_smiles, limit_depth=16, limit_nodes_count=1000000)
    tree.grow_tree(rules_smarts, pseudoatoms_dict)
    paths_starts = tree.paths_starts()
    return (tree, paths_starts)


test_reactions_df = pd.read_csv('test_mechanisms.csv', sep=';')

with Pool(processes=2) as pool:
    for n, result in tqdm(enumerate(pool.imap_unordered(create_tree, test_reactions_df.head(14).tail(2).iterrows(), chunksize=1), start=13)):
        tree, paths_starts = result
        with open(f'results/rv3_1000000/tree_{n}.pkl', 'wb') as file:
            dump(tree, file)

        with open(f'results/rv3_1000000/paths_starts_tree_{n}.txt', 'w') as file:
            for j in paths_starts:
                file.write(str(j) + '\n')

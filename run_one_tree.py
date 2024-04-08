# ограничения на дерево + правила 1 без С-Н

import pandas as pd
from mechanisms_tree import Tree
from pseudoatoms_dict import pseudoatoms_dict
from pickle import dump


# test_reactions_df = pd.read_csv('test_mechanisms.csv', sep=';')

root_molecules_smiles = '[H+].[H]ON=C(C([H])([H])C([H])([H])[H])C([H])([H])C([H])([H])[H]'
products_smiles = '[H]N(C(=O)C([H])([H])C([H])([H])[H])C([H])([H])C([H])([H])[H]'

with open('rules/rules3.smarts') as file:
    rules_smarts = file.read().split('\n')

tree = Tree(root_molecules_smiles, products_smiles, limit_depth=8, limit_nodes_count=1000)
tree.grow_tree(rules_smarts, pseudoatoms_dict)
paths_starts = tree.paths_starts()

with open(f'results/trees/tree.pkl', 'wb') as file:
    dump(tree, file)

with open(f'results/trees/paths_starts_tree.txt', 'w') as file:
    for j in paths_starts:
        file.write(str(j) + '\n')

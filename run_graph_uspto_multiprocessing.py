# реакции не все + ограничения на дерево + правила 3 + процессы

from mechanisms_graph import Graph
from pseudoatoms_dict import pseudoatoms_dict
from tqdm import tqdm
from pickle import dump
from multiprocessing import Pool
from CGRtools import RDFRead
from indigo import Indigo
ingido = Indigo()


def create_graph(reaction):
    with open('rules/rules3.smarts') as file:
        rules_smarts = file.read().split('\n')
    reaction_name = reaction.meta['r_c_1_freq']
    reaction.kekule()
    molecules = str(reaction).split('>')
    root_molecules_smiles = molecules[0]
    if molecules[1]:
        root_molecules_smiles = root_molecules_smiles + '.' + molecules[1]
    products_smiles = molecules[2]
    root_molecules_smiles = ingido.loadMolecule(root_molecules_smiles).canonicalSmiles()
    products_smiles = ingido.loadMolecule(products_smiles).canonicalSmiles()
    graph = Graph(root_molecules_smiles, products_smiles, limit_depth=11, limit_nodes_count=20000000,
                  stop_when_product=False)
    graph.grow_graph(rules_smarts, pseudoatoms_dict)
    paths_starts = graph.paths_starts()
    return (reaction_name, graph, paths_starts)


with RDFRead('data/72_test_reactions_with_H_rc_0_with_reag_(1).rdf', indexable=True) as file:
    reactions = file[:5]

with Pool(processes=5) as pool:
    for result in tqdm(pool.imap_unordered(create_graph, reactions, chunksize=1)):
        reaction_name, graph, paths_starts = result
        with open(f'results/graphs/graph_freq_{reaction_name}.pkl', 'wb') as file:
            dump(graph, file)

        with open(f'results/graphs/paths_starts_graph_freq_{reaction_name}.txt', 'w') as file:
            for j in paths_starts:
                file.write(str(j) + '\n')

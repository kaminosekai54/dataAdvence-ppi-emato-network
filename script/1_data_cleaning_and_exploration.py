################################################################################
# import of package
import pandas as pd
import scanpy as sc
from scipy import stats
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import requests
import json
from time import sleep
import os
from collections import Counter
import random


################################################################################

# data importation

#  we get the ppi dataset
ppiDf= pd.read_csv("/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/data/mippie_ppi_v1_0.tsv", sep = "\t")
# G = nx.from_numpy_array(ppiDf.values)
G=nx.Graph(name='full PPI Interaction Graph')
for (a,b) in zip(ppiDf.entrezA.tolist(), ppiDf.entrezB.tolist()):
    G.add_edge(a, b)

pos = nx.fruchterman_reingold_layout(G)
nx.draw(G, pos=pos, with_labels=False)
plt.savefig("/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/figures/full_ppi_network.png")

print(len(ppiDf))
print(len(ppiDf.entrezA.unique()))
print(len(ppiDf.entrezB.unique()))
print(ppiDf.head())
print(ppiDf.describe())


# now we import the scanpy dataset
# this dataset will be used  later
adata = sc.datasets.paul15()
# print(adata.obs_keys())
# print(adata.obs_names.tolist())
# print(adata.var_names.tolist())
#  if you want to have the csv that form the anotated data, uncomment the following line
# adata.write_csvs('/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/data/', )

#  after observation, we figure out that the gene ids used in the two different file are different, so we need to convert them

# first, we get the list of unique id to convert
nods = list(set(ppiDf.entrezA.unique().tolist()) | set(ppiDf.entrezB.unique().tolist()))


def chunkList(l, chunkSize):
    for i in range(len(l)):
        if i+ chunkSize < len(l): yield l[i:i+chunkSize]
        else: yield  l[i:]
    
        
    # ### not used anymore, but can steel be usefull
# convertENSIDToGeneID,
# this function will send a request to the biotools website, to convert the ENSID to a geneID
# this function is recursive to make it resillante,
# @param,
# @idToConvert, the ENSID that you want to convert,
# @nbTry, number of try the function called it self,
# @notFoundList, a list to add the potential not converted ensid, (can may be happend, so a hand made convertion will be requiered)
# @ return, a dictionary  on the shape {ensID: geneID}
def convertENSIDToGeneID(idToConvert, nbTry, notFoundList):
    try:
        url = "http://biotools.fr/mouse/ensembl_symbol_converter/"
        body = {'api': 1, 'id': idToConvert}
        r = requests.get(url, params=body)
        res= json.loads(r.text)
        return res
    except:
        # print(r.status_code)
        if nbTry < 5:
            sleep(nbTry+1)
            return scrap(idToConvert, nbTry+1, notFoundList)
        else:
            notFoundList.append(idToConvert)
            return {}



def convertENSIDList(listToConvert):
    geneSymbolDict = {}
    geneSymbolNotFound = []
    # pbar = tqdm(total=len(nods)) # We create a progress bar
    for idToConvert in listToConvert:
        geneSymbolDict   = {**geneSymbolDict , ** convertENSIDToGeneID(idToConvert, 0, geneSymbolNotFound )}
    # pbar.update(1)
    # print(len(geneSymbolDict))

    print(len(geneSymbolNotFound))
    print("sorry, but it seams you have to hand convert the following ensID : ", geneSymbolNotFound)
# pbar.close()
    print(geneSymbolDict)
    checkDf = pd.DataFrame.from_dict({"ens_id":geneSymbolDict.keys(), "gene_id":geneSymbolDict.values()})
    print(checkDf )
    checkDf.to_csv("/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/results/geneNameTranscription.csv", sep=",", index=False)
    return (checkDf, geneSymbolDict)

# geneTranscriptionDf = []
# geneSymbolDict ={}
# if not os.path.isfile('/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/results/geneNameTranscription.csv'): geneTranscriptionDf, geneSymbolDict   = convertENSIDList(nods)
# else: 
    # geneTranscriptionDf   = pd.read_csv("/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/results/geneNameTranscription.csv")
    # geneSymbolDict  = {k:v for (k,v) in zip(geneTranscriptionDf.ens_id.tolist(), geneTranscriptionDf.gene_id.tolist())}

newPPIDf = ppiDf.copy()
genIdToGeneSymbolDf = pd.read_csv("/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/data/mippie_proteins_v1_0.tsv", sep = "\t")
newPPIDf = newPPIDf.rename(columns={"entrezA":"A", "entrezB":"B"}).drop(["MIPPIE_score", "type","provider", "studies", "num_studies", "other_organisms", "num_organisms","exp_techniques","exp_scores"], axis=1)
newPPIDf.A = newPPIDf.A.apply(lambda x : str.lower(genIdToGeneSymbolDf[genIdToGeneSymbolDf.entrez == x].official_symbol.tolist()[0]))
newPPIDf.B = newPPIDf.B.apply(lambda x : str.lower(genIdToGeneSymbolDf[genIdToGeneSymbolDf.entrez == x].official_symbol.tolist()[0]))
newPPIDf = newPPIDf.drop_duplicates()
# newPPIDf = newPPIDf['A' != 'B']
newPPIDf= newPPIDf.query("A != B")


print(newPPIDf)
networkGene = [str.lower(id) for id in adata.var_names.unique().tolist()]
ppiGene = [str.lower(x) for x in set(newPPIDf.A.unique().tolist() + newPPIDf.B.unique().tolist())]
print(len(networkGene))
print(len(ppiGene))
commonGene = [id for id in networkGene if id in ppiGene]
notInPPI = [id for id in networkGene if id not in ppiGene]
# print(notInPPI[0:100])
print(len(notInPPI ))
print("nb gene common to the ppi ans single rna file : ", len(commonGene))

reducePPIIDDf = newPPIDf[newPPIDf.A.isin(commonGene) | newPPIDf.B.isin(commonGene)].drop_duplicates()
print(reducePPIIDDf )
reducePPIIDDf.to_csv("/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/results/reduce_emato_ppi_common.csv", index=False)
interractionCounter = Counter(reducePPIIDDf.A.tolist() + reducePPIIDDf.B.tolist())
# print(interractionCounter)
# print(interractionCounter.most_common())
countInterractionDf = pd.DataFrame.from_dict({"geneId":interractionCounter.keys(), "nbInterraction": interractionCounter.values()})
mostInterractiveDf = reducePPIIDDf[reducePPIIDDf.A.isin(countInterractionDf[countInterractionDf.nbInterraction > 20].geneId.tolist()) | reducePPIIDDf.B.isin(countInterractionDf[countInterractionDf.nbInterraction > 20].geneId.tolist())]
# print(mostInterractiveDf)
#  plot of the reduce network
if not os.path.isfile('/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/results/reduce_common_network.graphml'):
    G=nx.Graph(name='Protein Interaction Graph')
    for (a,b) in zip(reducePPIIDDf.A.tolist(), reducePPIIDDf.B.tolist()):
        G.add_edge(a, b)

    nx.draw(G)
    plt.savefig("/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/figures/reduce_network_only_common.png")
    nx.write_graphml(G, '/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/results/reduce_common_network.graphml')
    # nx.readwrite.write_gpickle
else : 
    G = nx.read_graphml('/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/results/reduce_common_network.graphml')
    pos = nx.fruchterman_reingold_layout(G)
    nx.draw(G, pos=pos, with_labels=False)
    plt.savefig("/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/figures/reduce_common_network.png")
    print("density of reduce network : ", nx.density(G))

# reducedNods = list(set(reducePPIIDDf.A.unique().tolist()) | set(reducePPIIDDf.B.unique().tolist()))
deg = nx.degree_centrality(G) # gives a dictionary!!!
# print(deg)
# print(G.degree())
degree_sequence = sorted([d for n, d in G.degree()], reverse=True) # gives array of degree values
# print(degree_sequence) 


# compute the density probability distribution

# Compute the degree distribution
degree_sequence = [d for n, d in G.degree()]

deggreDistribution = []
    
def bootStrapDensityDistribution(G, subGSize):
# Number of iterations
    n_iter = 10000
# Compute the bootstrap density
#list to store densities
    densities = []

    for i in range(n_iter):
        edges = random.sample(G.edges(), subGSize)
        subG = G.edge_subgraph(edges)
    # subG = G.subgraph(random.sample(G.nodes(), 25))

        densities.append(nx.density(subG))
    
# plot the density
    plt.hist(densities, bins=20)
    plt.xlabel('Density')
    plt.ylabel('Frequency')
    plt.title('Bootstrapped Density of the Graph')
    plt.show()
    plt.savefig("/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/figures/density_probability_distribution.png")


# modularity
# Compute the communities
communities = nx.community.greedy_modularity_communities(G)

# Compute the modularity
modularity = nx.community.modularity(G, communities)

print("modularity : ", modularity)

# 3. Draw network with node colors defined by degree
# plt.figure(figsize = (4, 4)) # set size of figure
# node_color = degree_sequence # assign node colors
# nx.draw(G, node_color = degree_sequence)
# plt.savefig("/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/figures/colored_nod_reduce_network_common.png")
################################################################################
# import of package
import pandas as pd
from scipy import stats
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import os
import community as community_louvain

################################################################################


path_df = pd.read_csv("/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/results/paga_path_erythrocytes.csv").drop('groups', axis=1)
# print(path_df)
def buildCorrelationMatrix(path_df, plotHitMap = False):
    corMatrixDf= pd.DataFrame(0, index=path_df.columns[1:], columns=path_df.columns[1:])
    num_cols = corMatrixDf.shape[1]
    corr_matrix = np.zeros((num_cols, num_cols))
    print(corr_matrix .shape)
    for i in range(num_cols):
        for j in range(i+1, num_cols):
            
    # print(corMatrixDf)
            corr, pval = stats.pearsonr(path_df.iloc[:,i], path_df.iloc[:,j])
            # Apply the Bonferroni correction to the p-value
            # pval_corrected = pval * num_cols
            # if pval_corrected < 0.05/ num_cols : corr_matrix[i,j] = pval_corrected
            corr_matrix[i,j] = abs(corr)
            
# remove the self loop
    np.fill_diagonal(corr_matrix, 0)
    corMatrixDf= pd.DataFrame(corr_matrix,  index=corMatrixDf.index, columns=corMatrixDf.columns)
    print("correlation matrix : ")
    print(corMatrixDf)

    corMatrixDf.to_csv("/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/results/corelationDf.csv")
    if plotHitMap:
#create a heatmap of the correlation matrix
        sns.heatmap(correlationMatrix, annot=True, cmap='coolwarm')
# show the plot
        plt.show()
        plt.savefig("/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/figures/correlation_matrix_hitMap.png")

def filterByDensity(corMatrixDf, desiredThreshold):
    tmpDf = corMatrixDf.copy()
    adj_matrix = tmpDf.values
    print(adj_matrix)
    tmpG = nx.from_numpy_array(adj_matrix)
    if nx.density(tmpG) <= desiredThreshold:
        return corMatrixDf, tmpG

    for threshold in np.arange(0, 1, 0.05):
        adj_matrix[adj_matrix < threshold] = 0
        tmpG = nx.from_numpy_array(adj_matrix)
        if nx.density(tmpG) <= desiredThreshold:
            print("correct density reached, with a correlation threshold = :", threshold)
            adj_matrix[adj_matrix >= threshold] = 1
            tmpG = nx.from_numpy_array(adj_matrix)
            tmpDf= pd.DataFrame(adj_matrix,  index=corMatrixDf.index, columns=corMatrixDf.columns)
            print(tmpDf)
            return tmpDf, tmpG

    # print("an error occured, impossible to correctly filter the matrix")
    # return tmpDf, tmpG



if not os.path.isfile("/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/results/corelationDf.csv") : buildCorrelationMatrix(path_df)
corMatrixDf = pd.read_csv("/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/results/corelationDf.csv")
corMatrixDf  = corMatrixDf .set_index(corMatrixDf.columns[0])
corMatrixDf= corMatrixDf.drop(corMatrixDf.columns[0], axis=0).drop(corMatrixDf.columns[0], axis=1)
# corMatrixDf = corMatrixDf.applymap(lambda x: 1 if x > 0 else x)
# adj_matrix = corMatrixDf.values
# G = nx.from_numpy_array(adj_matrix)

corMatrixDf, G = filterByDensity(corMatrixDf , 0.1)


print("density : ", nx.density(G))
# Add the column names as node labels
labels = {i:col for i, col in enumerate(corMatrixDf.columns)}
nx.relabel_nodes(G, labels, copy=False)
# Draw the initial graph
# pos = nx.fruchterman_reingold_layout(G, k=50)
nx.draw(G)
plt.show()
plt.savefig("/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/figures/trajectories_cluster_network.png")
#Determine the clusters using the Louvain method
partition = community_louvain.best_partition(G)
clusters = list(partition.values())

# print(partition)
community_sizes = {i: list(partition.values()).count(i) for i in set(partition.values())}
# print(community_sizes)
# Plot the network with node color representing the community
pos = nx.spring_layout(G)
nx.draw_networkx_nodes(G, pos, partition.keys(), node_size = [community_sizes[v]*100 for v in partition.values()])
nx.draw_networkx_edges(G, pos, alpha=0.5)
plt.show()


# Create a new graph with a color code by cluster
# color_map = [i for i in clusters]
# pos = nx.fruchterman_reingold_layout(G, k=50)

# nx.draw(G, pos=pos, node_color=color_map, with_labels=True)
# plt.show()
plt.savefig("/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/figures/colored_trajectories_cluster_network.png")

# superposition des graph
# G_PPI_common = nx.read_graphml('/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/results/reduce_common_network.graphml')
# pos2 = nx.fruchterman_reingold_layout(G_PPI_common)

# nx.draw_networkx_nodes(G_PPI_common, pos2, node_color='r')
# nx.draw_networkx_edges(G_PPI_common, pos2, edge_color='r')

# Draw the second graph
# pos2 = nx.spring_layout(G2)
# nx.draw_networkx_nodes(G, pos, node_color='b', alpha=0.5)
# nx.draw_networkx_edges(G, pos, edge_color='b', alpha=0.5)

# Add a legend
# plt.legend(['Graph 1', 'Graph 2'])
# plt.show()
# plt.savefig("/mnt/c/Users/culpi/Desktop/git-repo/dataAdvence-ppi-emato-network/figures/overlaped_network.png")

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
import os

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
# sc.logging.print_versions()
results_file = './write/paul15.h5ad'
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')  # low dpi (dots per inch) yields small inline figures
adata = sc.datasets.paul15()
# print(adata)
adata.X = adata.X.astype('float64')  # this is not required and results will be comparable without it

# Preprocessing and Visualization

sc.pp.recipe_zheng17(adata)

# Computes PCA coordinates, loadings and variance decomposition. Uses the implementation of scikit-learn [Pedregosa11].
sc.tl.pca(adata, svd_solver='arpack')



# Compute a neighborhood graph of observations [McInnes18].
sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
sc.tl.draw_graph(adata)
pl.show()
sc.pl.draw_graph(adata, color='paul15_clusters', legend_loc='on data')
pl.show()
# Denoising the graph
# To denoise the graph, we represent it in diffusion map space (and not in PCA space). Computing distances within a few diffusion components amounts to denoising the graph - we just take a few of the first spectral components. It’s very similar to denoising a data matrix using PCA. The approach has been used in a couple of papers, see e.g. Schiebinger et al. (2017) or Tabaka et al. (2018). It’s also related to the principles behind MAGIC Dijk et al. (2018).
sc.tl.diffmap(adata)
sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')
sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata, color='paul15_clusters', legend_loc='on data')


# Clustering and PAGA
# Note that today, we’d use sc.tl.leiden - here, we use sc.tl.louvain the sake of reproducing the paper results.
sc.tl.louvain(adata, resolution=1.0)

# For simple, coarse-grained visualization, compute the PAGA graph, a coarse-grained and simplified (abstracted) graph. Non-significant edges in the coarse- grained graph are thresholded away.
sc.tl.paga(adata, groups='louvain')
sc.pl.paga(adata, color=['louvain', 'Hba-a2', 'Elane', 'Irf8'])
sc.pl.paga(adata, color=['louvain', 'Itga2b', 'Prss34', 'Cma1'])
# Actually annotate the clusters — note that Cma1 is a Mast cell marker and only appears in a small fraction of the cells in the progenitor/stem cell cluster 8, see the single-cell resolved plot below.
print(adata.obs['louvain'].cat.categories)

adata.obs['louvain_anno'] = adata.obs['louvain']
adata.obs['louvain_anno'].cat.categories = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10/Ery', '11', '12',
       '13', '14', '15', '16/Stem', '17', '18', '19/Neu', '20/Mk', '21', '22/Baso', '23', '24/Mo']

    #    Let’s use the annotated clusters for PAGA.
sc.tl.paga(adata, groups='louvain_anno')
sc.pl.paga(adata, threshold=0.03, show=True)

# The following is just as well possible for a UMAP.
sc.tl.draw_graph(adata, init_pos='paga')
sc.pl.draw_graph(adata, color=['louvain_anno', 'Itga2b', 'Prss34', 'Cma1'], legend_loc='on data')

# Choose the colors of the clusters a bit more consistently.
pl.figure(figsize=(8, 2))
for i in range(28):
    pl.scatter(i, 1, c=sc.pl.palettes.zeileis_28[i], s=200)
pl.show()
zeileis_colors = np.array(sc.pl.palettes.zeileis_28)
new_colors = np.array(adata.uns['louvain_anno_colors'])
new_colors[[16]] = zeileis_colors[[12]]  # Stem colors / green
new_colors[[10, 17, 5, 3, 15, 6, 18, 13, 7, 12]] = zeileis_colors[[5, 5, 5, 5, 11, 11, 10, 9, 21, 21]]  # Ery colors / red
new_colors[[20, 8]] = zeileis_colors[[17, 16]]  # Mk early Ery colors / yellow
new_colors[[4, 0]] = zeileis_colors[[2, 8]]  # lymph progenitors / grey
new_colors[[22]] = zeileis_colors[[18]]  # Baso / turquoise
new_colors[[19, 14, 2]] = zeileis_colors[[6, 6, 6]]  # Neu / light blue
new_colors[[24, 9, 1, 11]] = zeileis_colors[[0, 0, 0, 0]]  # Mo / dark blue
new_colors[[21, 23]] = zeileis_colors[[25, 25]]  # outliers / grey
adata.uns['louvain_anno_colors'] = new_colors
# And add some white space to some cluster names. The layout shown here differs from the one in the paper, which can be found here. These differences, however, are only cosmetic. We had to change the layout as we moved from a randomized PCA and float32 to float64 precision.

sc.pl.paga_compare(
    adata, threshold=0.03, title='', right_margin=0.2, size=10, edge_width_scale=0.5,
    legend_fontsize=12, fontsize=12, frameon=False, edges=True, save=True)


    # Choose a root cell for diffusion pseudotime.
adata.uns['iroot'] = np.flatnonzero(adata.obs['louvain_anno']  == '16/Stem')[0]
sc.tl.dpt(adata)


# Select some of the marker gene names.
# reducePPIIDDf= pd.read_csv("/mnt/c/Users/alexis/Desktop/git-repo/dataAdvence-ppi-emato-network/results/reduce_emato_ppi_common.csv")
gene_names = sc.datasets.paul15().var_names.unique().tolist()
# list(set(reducePPIIDDf.A.unique().tolist() + reducePPIIDDf.B.unique().tolist()))
# ['Gata2', 'Gata1', 'Klf1', 'Epor', 'Hba-a2',  # erythroid
            #   'Elane', 'Cebpe', 'Gfi1',                    # neutrophil
            #   'Irf8', 'Csf1r', 'Ctsg']                     # monocyte

adata_raw = sc.datasets.paul15()
sc.pp.log1p(adata_raw)
sc.pp.scale(adata_raw)
adata.raw = adata_raw
sc.pl.draw_graph(adata, color=['louvain_anno', 'dpt_pseudotime'], legend_loc='on data')
paths = [('erythrocytes', [16, 12, 7, 13, 18, 6, 5, 10]),
         ('neutrophils', [16, 0, 4, 2, 14, 19]),
         ('monocytes', [16, 0, 4, 11, 1, 9, 24])]
adata.obs['distance'] = adata.obs['dpt_pseudotime']
adata.obs['clusters'] = adata.obs['louvain_anno']  # just a cosmetic change
adata.uns['clusters_colors'] = adata.uns['louvain_anno_colors']
# if not os.path.isdir("write") : os.mkdir("write")

_, axs = pl.subplots(ncols=3, figsize=(6, 2.5), gridspec_kw={'wspace': 0.05, 'left': 0.12})
pl.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2)
for ipath, (descr, path) in enumerate(paths):
    _, data = sc.pl.paga_path(
        adata, path, gene_names,
        show_node_names=False,
        ax=axs[ipath],
        ytick_fontsize=12,
        left_margin=0.15,
        n_avg=50,
        annotations=['distance'],
        show_yticks=True if ipath==0 else False,
        show_colorbar=False,
        color_map='Greys',
        groups_key='clusters',
        color_maps_annotations={'distance': 'viridis'},
        title='{} path'.format(descr),
        return_data=True,
        show=False)
    data.to_csv('/mnt/c/Users/alexis/Desktop/git-repo/dataAdvence-ppi-emato-network/results/paga_path_{}.csv'.format(descr))
pl.savefig('/mnt/c/Users/alexis/Desktop/git-repo/dataAdvence-ppi-emato-network/figures/paga_path_paul15.pdf')
pl.show()
print("all finish")

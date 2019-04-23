import palantir

# Plotting and miscellaneous imports
import os
import matplotlib
import matplotlib.pyplot as plt


palantir_dir = os.path.expanduser('./')
counts = palantir.io.from_csv(palantir_dir + 't_MGH54_mat.txt',delimiter=',')

##########
fig, ax = palantir.plot.plot_molecules_per_cell_and_gene(counts)
plt.show()
############

norm_df = palantir.preprocess.normalize_counts(counts)
norm_df = palantir.preprocess.log_transform(norm_df)

pca_projections, _ = palantir.utils.run_pca(norm_df)
dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=200)
ms_data = palantir.utils.determine_multiscale_space(dm_res)
tsne = palantir.utils.run_tsne(ms_data)

imp_df = palantir.utils.run_magic_imputation(norm_df, dm_res)



####################

min(tsne.ix[:,0])
max(tsne.ix[:,0])
min(tsne.ix[:,1])
max(tsne.ix[:,1])

fig, ax = palantir.plot.plot_tsne(tsne)
plt.show()

import numpy as np


tsne.index[np.where( (tsne.ix[:,0] > 5.5) * (tsne.ix[:,1] < -18) )]

tsne.index[np.where(tsne.index=='PDGFRA')]




fig, ax = palantir.plot.plot_tsne_by_cell_sizes(counts, tsne)
plt.show()

palantir.plot.plot_diffusion_components(tsne, dm_res)
plt.show()


palantir.plot.plot_gene_expression(norm_df, tsne, ['A1CF'])
plt.show()

palantir.plot.plot_gene_expression(imp_df, tsne, ['PDGFRA'])
plt.show()
palantir.plot.plot_gene_expression(norm_df, tsne, ['PDGFRA'])
plt.show()
##############################

cells=['MGH54_P7_G06']
palantir.plot.highlight_cells_on_tsne(tsne, cells)
plt.show()


start_cell='MGH54_P7_G06'
pr_res = palantir.core.run_palantir(ms_data, start_cell, num_waypoints=500)
palantir.plot.plot_palantir_results(pr_res, tsne)
plt.show()





genes = ['OLIG2','PDGFRA','ALDH1L1']
gene_trends = palantir.presults.compute_gene_trends( pr_res, imp_df.loc[:, genes])
palantir.plot.plot_gene_trends(gene_trends)
plt.show()

palantir.plot.plot_gene_trend_heatmaps(gene_trends)
plt.show()




#import rpy2
#import rpy2.robjects as RObjects
#from rpy2.robjects.packages import importr
#importr("name package", lib_loc = "/Library/Frameworks/R.framework/Versions/3.5/Resources/")

clusters = palantir.utils.determine_cell_clusters(pca_projections)
palantir.plot.plot_cell_clusters(tsne, clusters )
plt.show()


gene_trends = palantir.presults.compute_gene_trends(pr_res, 
                    imp_df.iloc[:, 0:1000], ['MGH54_P16_F01'])
trends = gene_trends['MGH54_P16_F01']['trends']
gene_clusters = palantir.presults.cluster_gene_trends(trends)
palantir.plot.plot_gene_trend_clusters(trends, gene_clusters)
plt.show()


###########
import pickle 
pickle.dump(xxx, open('XXX.obj','w'))

###########
xxx = pickle.load(open('XXX.obj'))



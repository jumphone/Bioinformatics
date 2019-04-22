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
dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=5)
ms_data = palantir.utils.determine_multiscale_space(dm_res)
tsne = palantir.utils.run_tsne(ms_data)
imp_df = palantir.utils.run_magic_imputation(norm_df, dm_res)



####################
fig, ax = palantir.plot.plot_tsne(tsne)
plt.show()

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




start_cell='MGH54_P2_H03'
pr_res = palantir.core.run_palantir(ms_data, start_cell, num_waypoints=200)
plt.show()


###########
import pickle 
pickle.dump(xxx, open('XXX.obj','w'))

###########
xxx = pickle.load(open('XXX.obj'))



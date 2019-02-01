import velocyto as vcy
import numpy as np

vlm = vcy.VelocytoLoom("picard_KA4FM.loom")
#vlm.normalize("S", size=True, log=True)
#vlm.dump_hdf5("my_velocyto_analysis")



vlm.filter_cells(bool_array=vlm.initial_Ucell_size > np.percentile(vlm.initial_Ucell_size, 0.5))
vlm.ca["ClusterName"]=['TMP'] * len(vlm.ca['CellID'])


vlm.set_clusters(vlm.ca["ClusterName"])
vlm.score_detection_levels(min_expr_counts=40, min_cells_express=30)
vlm.filter_genes(by_detection_levels=True)
vlm.score_cv_vs_mean(3000, plot=True, max_expr_avg=35)
vlm.filter_genes(by_cv_vs_mean=True)

vlm._normalize_S(relative_size=vlm.S.sum(0),
             target_size=vlm.S.sum(0).mean())
vlm._normalize_U(relative_size=vlm.U.sum(0),
             target_size=vlm.U.sum(0).mean())


vlm.perform_PCA()
vlm.knn_imputation(n_pca_dims=20, k=500, balanced=True, b_sight=2000, b_maxl=1500, n_jobs=16)


vlm.fit_gammas()
vlm.predict_U()
vlm.calculate_velocity()
vlm.calculate_shift(assumption="constant_velocity")
vlm.extrapolate_cell_at_t(delta_t=1.)

from sklearn.manifold import TSNE
bh_tsne = TSNE()
vlm.ts = bh_tsne.fit_transform(vlm.pcs[:, :25])

vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt", psc=1,
                             n_neighbors=2000, knn_random=True, sampled_fraction=0.5)
vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=True)

vlm.calculate_grid_arrows(smooth=0.8, steps=(40, 40), n_neighbors=300)

import matplotlib.pyplot as plt


plt.figure(None,(20,10))
#vlm.plot_grid_arrows(quiver_scale=0.6,
#                    scatter_kwargs_dict={"alpha":0.35, "lw":0.35, "edgecolor":"0.4", "s":38, "rasterized":True}, min_mass=24, angles='xy', scale_units='xy',
#                    headaxislength=2.75, headlength=5, headwidth=4.8, minlength=1.5,
#                    plot_random=True, scale_type="absolute")

vlm.plot_grid_arrows(quiver_scale=1.4,
                     scatter_kwargs_dict={"alpha":0.35, "lw":0.35, "edgecolor":"0.4", "s":38, "rasterized":True},
                     min_mass=5.5, angles='xy', scale_units='xy',
                     headaxislength=2.75, headlength=5, headwidth=4.8, minlength=1.5,
                     plot_random=False, scale_type="relative")


plt.savefig('foo.pdf')


import loompy

ds = loompy.connect("picard_KA4FM.loom")
S = np.array(ds["spliced"][:])
cellids = np.array(ds.ca.CellID)
genenames = np.array(ds.ra.Gene)
ds.close()


mask = np.in1d(cellids, vlm.ca["CellID"])
cellids = cellids[mask]
S = S[:, mask]

mask = (S.sum(1) > 20) & ((S>0).sum(1) > 15)
genenames = genenames[mask]
S = S[mask, :]

np.alltrue(cellids == vlm.ca["CellID"])

# Size normalize
S_sz = S / S.sum(0)
# Impute
Sx_sz = vcy.convolve_by_sparse_weights(S_sz, vlm.knn_smoothing_w)


plt.figure(None, (5,8))
gs = plt.GridSpec(3,2)
gene_list = ['Pdgfra', 'Cspg4', 
             'Olig1', 'Olig2',
              "Chd8", "Smarca4",]

for i, gene in enumerate(gene_list):
  
    plt.subplot(gs[i])
    
    colorandum = Sx_sz[np.where(genenames == gene)[0][0], :]   
    vcy.scatter_viz(vlm.ts[:,0], vlm.ts[:, 1], c=colorandum, cmap="magma_r", alpha=0.35, s=3, rasterized=True)
    plt.title(gene)
    plt.axis("off")
    
#plt.savefig("../figures/Haber_cellcycle_genes.pdf")

plt.savefig('gene.pdf')



plt.figure(None, (3,3))
gene='Pdgfra'
colorandum = Sx_sz[np.where(genenames == gene)[0][0], :]
vcy.scatter_viz(vlm.ts[:,0], vlm.ts[:, 1], c=colorandum, cmap="magma_r", alpha=0.35, s=3, rasterized=True)
plt.title(gene)
plt.axis("off")
plt.savefig('Pdgfra.pdf')


plt.figure(None, (3,3))
gene='Cspg4'
colorandum = Sx_sz[np.where(genenames == gene)[0][0], :]
vcy.scatter_viz(vlm.ts[:,0], vlm.ts[:, 1], c=colorandum, cmap="magma_r", alpha=0.35, s=3, rasterized=True)
plt.title(gene)
plt.axis("off")
plt.savefig('Cspg4.pdf')








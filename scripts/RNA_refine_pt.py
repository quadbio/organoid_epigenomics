import os
import h5py
import numpy as np
import scipy as sp
import anndata as ad
import pandas as pd
import scanpy as sc 
import scvelo as scv
import loompy as lp
import cellrank as cr

os.chdir('/home/jfleck/projects/cutntag/')

#### Read h5ad ####
rna = sc.read('data/RNA/RNA_all_v3.2cellrank.h5ad')
rna_meta = pd.read_csv('data/RNA/RNA_all_v2.3lines_meta.tsv', sep='\t')

cells_use = rna.obs.index.intersection(rna_meta.cell)
rna_use = rna[cells_use, :]
meta_use = rna_meta.set_index('cell').loc[cells_use, :]

#### Get ctx pseudotime ####
rna_ctx = rna_use[meta_use.ctx_traject, :]

sc.pp.highly_variable_genes(rna_ctx)
sc.pp.scale(rna_ctx)
sc.pp.pca(rna_ctx)
sc.pp.neighbors(rna_ctx, use_rep='X_pca', n_neighbors=100)
sc.tl.umap(rna_ctx)

p = sc.pl.umap(rna_ctx, color=['pseudotime_ranks', 'celltype_jf'], show=False)
p[0].figure.savefig('plots/RNA/ctx_rna_pt_umap.png')

scv.pp.moments(rna_ctx, n_neighbors=100, use_rep='X_css')
scv.tl.velocity(rna_ctx, mode='stochastic')
scv.tl.velocity_graph(rna_ctx, n_neighbors=100)

p = scv.pl.velocity_embedding_stream(rna_ctx, color='pseudotime_ranks', show=False)
p.figure.savefig('plots/RNA/ctx_rna_velo_pt_umap.png')

start_cell = rna_ctx.obs.index[rna_ctx.obs.clusters=='EB_2'][0]
end_clust = rna_ctx[rna_ctx.obs.clusters=='late_9', ]
end_cell = end_clust.obs.index[end_clust.obsm['X_umap'][:,1].argmax()]

start_idx = np.where(rna_ctx.obs.index==start_cell)[0][0]
end_idx = np.where(rna_ctx.obs.index==end_cell)[0][0]

# rna_ctx.obs['last_clust'] = (rna_ctx.obs.clusters=='late_9').astype(int)
# p = sc.pl.umap(rna_ctx, color='last_clust', show=False)
# p.figure.savefig('plots/RNA/ctx_rna_last_clust_umap.png')

scv.tl.velocity_pseudotime(rna_ctx, root_key=start_idx, end_key=end_idx)
scv.tl.velocity_pseudotime(rna_ctx)

rna_ctx.obs['pseudotime_ranks'] = rna_ctx.obs['velocity_pseudotime'].rank() / max(rna_ctx.obs['velocity_pseudotime'].rank())

p = scv.pl.velocity_embedding_stream(rna_ctx, color=['velocity_pseudotime', 'pseudotime_ranks'], show=False)
p[0].figure.savefig('plots/RNA/ctx_rna_velo_ctx_pt_umap.png')

ctx_pt_meta = rna_ctx.obs[['velocity_pseudotime', 'pseudotime_ranks']]
ctx_pt_meta.columns = ['ctx_vpt', 'ctx_vpt_ranks']
ctx_pt_meta.reset_index().to_csv('data/RNA/cellrank/RNA_ctx_vpt.tsv', sep='\t', index=False)

rna_ctx.write('data/RNA/RNA_ctx_v3.2vpt.h5ad')




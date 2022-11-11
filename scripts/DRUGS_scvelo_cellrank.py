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

#### Join all kallisto results ####
# all_samples = os.listdir('data/individual_samples/')
# RNA_samples = [s for s in all_samples
#     if ('scRNA' in s) and (('DMSO' in s) or ('A395' in s) or ('A485' in s) or ('drug' in s))
# ]
# 
# adata_list = []
# 
# for s in RNA_samples:
#     print(s)
# 
#     loom_path = f'data/individual_samples/{s}/{s}.loom'
#     ds = lp.connect(loom_path)
#     X = sp.sparse.csr_matrix(ds.layers[''].sparse().T)
#     X_mat = sp.sparse.csr_matrix(ds.layers[''].sparse().T)
#     X_spliced = sp.sparse.csr_matrix(ds.layers['spliced'].sparse().T)
#     X_unspliced = sp.sparse.csr_matrix(ds.layers['unspliced'].sparse().T)
# 
#     var = pd.DataFrame({
#         'Gene': ds.ra['Gene'],
#         'Accession': ds.ra['Accession']
#     }).set_index('Gene')
# 
#     obs = pd.DataFrame({
#         'CellID': ds.ca['CellID']
#     }).set_index('CellID')
# 
#     ds.close()
# 
#     adata = ad.AnnData(
#         X = X, layers={'spliced':X_spliced, 'unspliced':X_unspliced, 'matrix':X_mat}, obs=obs, var=var
#     )
#     adata.obs['sample'] = s
#     adata.var_names_make_unique()
#     adata_list.append(adata)
# 
# adata_combined = ad.concat(adata_list)
# 
# adata_combined.write('data/drugs/drugs_all_velo_v0raw.h5ad')


#### Read data ####

adata_combined = sc.read('data/drugs/drugs_all_velo_v0raw.h5ad')
drugs_d158_H3K27me3 = sc.read('data/drugs/drugs_d15_d18_A395_v1_v1.2lines.h5ad')

drugs_d158_A395 = drugs_d158_H3K27me3
drugs_d158_A395_cells = drugs_d158_A395.obs.index.intersection(adata_combined.obs.index)

adata_d158_A395 = adata_combined[drugs_d158_A395_cells, :].copy()
drugs_d158_A395_ = drugs_d158_A395[drugs_d158_A395_cells, :].copy()

adata_d158_A395.obs['RNA_snn_res.0.2'] = drugs_d158_A395_.obs['RNA_snn_res.0.2'].copy()
adata_d158_A395.obs['line'] = drugs_d158_A395_.obs['line'].copy()
adata_d158_A395.obs['inhib_annotation'] = drugs_d158_A395_.obs['inhib_annotation'].copy()
adata_d158_A395.obs['inhibitor_target'] = drugs_d158_A395_.obs['inhibitor_target'].copy()
adata_d158_A395.obs['orig.ident'] = drugs_d158_A395_.obs['orig.ident'].copy()

adata_d158_A395.obsm['X_harmony'] = drugs_d158_A395_.obsm['X_harmony']
adata_d158_A395.obsm['X_humap'] = drugs_d158_A395_.obsm['X_humap']


#### Preproc ####
sc.pp.log1p(adata_d158_A395)
adata_d158_A395.layers['lognorm'] = adata_d158_A395.X.copy()

sc.pp.highly_variable_genes(adata_d158_A395)
sc.pp.pca(adata_d158_A395)
sc.pp.neighbors(adata_d158_A395)
sc.tl.umap(adata_d158_A395)

p = sc.pl.scatter(adata_d158_A395, basis='humap', color='orig.ident', show=False)
p.figure.savefig('plots/drugs/velo/drugs_d158_A395_umap.png')

scv.pp.moments(adata_d158_A395, n_neighbors=20, use_rep='X_harmony')
scv.tl.velocity(adata_d158_A395, mode='stochastic')
scv.tl.velocity_graph(adata_d158_A395, n_neighbors=20, n_jobs=36)


# Plot velocities
p = scv.pl.velocity_embedding_stream(adata_d158_A395, basis='humap', color='orig.ident', show=False)
p.figure.savefig('plots/drugs/velo/drugs_d158_A395_velo_umap.png')

scv.tl.velocity_pseudotime(adata_d158_A395)

p = scv.pl.velocity_embedding_stream(adata_d158_A395, basis='humap', color='velocity_pseudotime', show=False)
p.figure.savefig('plots/drugs/velo/drugs_d158_A395_velo_pt_umap.png')

adata_d158_A395.write('data/drugs/drugs_d15_d18_A395_v2velo_v1.2lines.h5ad')


adata_d158_A395 = sc.read('data/drugs/drugs_d15_d18_A395_v2velo_v1.2lines.h5ad')

p = scv.pl.velocity_embedding_stream(adata_d158_A395, basis='humap', color='inhibitor_target', show=False)
p.figure.savefig('plots/drugs/velo/drugs_d158_A395_velo_treatment_umap.png')



















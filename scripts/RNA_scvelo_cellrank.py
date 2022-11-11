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


#### Join spliced and unspliced counts ####
all_samples = os.listdir('data/individual_samples/')
RNA_samples = [s for s in all_samples if ('scRNA' in s) and ('21d' not in s)]

adata_list = []

for s in RNA_samples:
    print(s)

    loom_path = f'data/individual_samples/{s}/{s}.loom'
    ds = lp.connect(loom_path)
    X = sp.sparse.csr_matrix(ds.layers[''].sparse().T)
    X_mat = sp.sparse.csr_matrix(ds.layers[''].sparse().T)
    X_spliced = sp.sparse.csr_matrix(ds.layers['spliced'].sparse().T)
    X_unspliced = sp.sparse.csr_matrix(ds.layers['unspliced'].sparse().T)

    var = pd.DataFrame({
        'Gene': ds.ra['Gene'],
        'Accession': ds.ra['Accession']
    }).set_index('Gene')

    obs = pd.DataFrame({
        'CellID': ds.ca['CellID']
    }).set_index('CellID')

    ds.close()

    adata = ad.AnnData(
        X = X, layers={'spliced':X_spliced, 'unspliced':X_unspliced, 'matrix':X_mat}, obs=obs, var=var
    )
    adata.obs['sample'] = s
    adata.var_names_make_unique()
    adata_list.append(adata)

adata_combined = ad.concat(adata_list)
adata_combined.obs.index = adata_combined.obs.index + '-1'

adata_combined.write('data/RNA/RNA_all_velo_v0raw.h5ad')


#### Subset and join with processed data ####
rna = sc.read('data/RNA/RNA_all_srt_v2.2matched.h5ad')
rna.layers['counts'] = rna.raw.X.copy()

common_cells = rna.obs.index.intersection(adata_combined.obs.index)
common_genes = rna.var.index.intersection(adata_combined.var.index)

rna_use = rna[common_cells, common_genes].copy()
spl_unspl_use = adata_combined[common_cells, common_genes].copy()

rna_use.layers['spliced'] = spl_unspl_use.layers['spliced'].copy()
rna_use.layers['unspliced'] = spl_unspl_use.layers['unspliced']
rna_use.layers['lognorm'] = rna_use.X.copy()

rna_use.write('data/RNA/RNA_all_v3velo.h5ad')


#### Preproc and run RNA velocity ####
rna_use = sc.read('data/RNA/RNA_all_v3velo.h5ad')

sc.pp.highly_variable_genes(rna_use)
sc.pp.neighbors(rna_use, use_rep='X_css')

# Plot to sanity check 
p = sc.pl.scatter(rna_use, basis='cssumap', show=False, color='celltype_jf')
p.figure.savefig('plots/RNA/rna_velo_umap.png')

scv.pp.moments(rna_use, n_neighbors=20, use_rep='X_css')
scv.tl.velocity(rna_use, mode='stochastic')
scv.tl.velocity_graph(rna_use, n_neighbors=20, n_jobs=40)

# Plot velocities
p = scv.pl.velocity_embedding_stream(rna_use, basis='cssumap', color='celltype_jf', show=False)
p.figure.savefig('plots/RNA/rna_velo_stream_umap.png')

rna_use.write('data/RNA/RNA_all_v3.1velo.h5ad')


#### Compute pseudotime ####
rna_use = sc.read('data/RNA/RNA_all_v3.1velo.h5ad')

start_cell = rna_use.obs.index[rna_use.obs.clusters=='EB_3'][0]
end_cell = rna_use.obs.index[rna_use.obs.clusters=='late_9'][0]

start_idx = np.where(rna_use.obs.index==start_cell)[0][0]
end_idx = np.where(rna_use.obs.index==end_cell)[0][0]

scv.tl.velocity_pseudotime(rna_use, root_key=start_idx, end_key=end_idx)

p = scv.pl.velocity_embedding_stream(rna_use, basis='cssumap', color='velocity_pseudotime', show=False)
p.figure.savefig('plots/RNA/rna_velo_pt_umap.png')

rna_use.obs['pseudotime_ranks'] = rna_use.obs['velocity_pseudotime'].rank() / max(rna_use.obs['velocity_pseudotime'].rank())

p = scv.pl.velocity_embedding_stream(rna_use, basis='cssumap', color='pseudotime_ranks', show=False)
p.figure.savefig('plots/RNA/rna_velo_pt_ranks_umap.png')


#### Annotate terminal states ####
rna_use.obs['terminal_states'] = np.nan
rna_use.obs['terminal_states'] = np.nan
rna_use.obs['terminal_states'][rna_use.obs.clusters=='late_9'] = 'ctx'
rna_use.obs['terminal_states'][rna_use.obs.clusters=='retina_21'] = 'dien'
rna_use.obs['terminal_states'][rna_use.obs.clusters=='late_5'] = 'rhom'
rna_use.obs['terminal_states'][rna_use.obs.clusters=='late_8'] = 'mesen'
rna_use.obs['terminal_states'][rna_use.obs.clusters=='retina_4'] = 'retina'
rna_use.obs['terminal_states'][rna_use.obs.clusters=='mo8_9'] = 'astro'
rna_use.obs['terminal_states'][rna_use.obs.clusters=='mo8_6'] = 'opc'
rna_use.obs['terminal_states'][rna_use.obs.clusters=='mid_11'] = 'nn'
rna_use.obs['terminal_states'][rna_use.obs.clusters=='mo8_0'] = 'chp'

rna_use.obs['terminal_states'] = np.nan
rna_use.obs['terminal_states'][rna_use.obs.celltype_jf=='ctx_ex'] = 'ctx'
rna_use.obs['terminal_states'][rna_use.obs.celltype_jf=='dien_ex'] = 'dien'
rna_use.obs['terminal_states'][rna_use.obs.celltype_jf=='rhom_ex'] = 'rhom'
rna_use.obs['terminal_states'][rna_use.obs.celltype_jf=='mesen_ex'] = 'mesen'
rna_use.obs['terminal_states'][rna_use.obs.celltype_jf=='RGC'] = 'retina'
rna_use.obs['terminal_states'][rna_use.obs.clusters=='mo8_9'] = 'astro'
rna_use.obs['terminal_states'][rna_use.obs.clusters=='mo8_6'] = 'opc'
rna_use.obs['terminal_states'][rna_use.obs.clusters=='mid_11'] = 'nn'
rna_use.obs['terminal_states'][rna_use.obs.clusters=='mo8_0'] = 'chp'

p = scv.pl.velocity_embedding_stream(rna_use, basis='cssumap', color='terminal_states', show=False)
p.figure.savefig('plots/RNA/rna_velo_terminal_umap.png')

rna_use.obs['terminal_states'] = rna_use.obs['terminal_states'].astype('category')


#### Run cellrank ####
terminal_states = rna_use.obs['terminal_states']

vk = cr.tl.kernels.VelocityKernel(rna_use)
vk.compute_transition_matrix()

pk = cr.tl.kernels.PalantirKernel(rna_use, time_key='velocity_pseudotime')
pk.compute_transition_matrix()

ck = cr.tl.kernels.ConnectivityKernel(rna_use)
ck.compute_transition_matrix()

combined_kernel = vk
combined_kernel = ck
combined_kernel = pk
combined_kernel = 0.3 * vk + 0.7 * pk 

g = cr.tl.estimators.GPCCA(combined_kernel)
g.set_terminal_states(terminal_states)
g.compute_absorption_probabilities(n_jobs=40, use_petsc=True)

rna_use.obs['to_ctx'] = rna_use.obsm['to_terminal_states']['ctx'].X
rna_use.obs['to_dien'] = rna_use.obsm['to_terminal_states']['dien'].X
rna_use.obs['to_rhom'] = rna_use.obsm['to_terminal_states']['rhom'].X
rna_use.obs['to_mesen'] = rna_use.obsm['to_terminal_states']['mesen'].X
rna_use.obs['to_retina'] = rna_use.obsm['to_terminal_states']['retina'].X
rna_use.obs['to_astro'] = rna_use.obsm['to_terminal_states']['astro'].X
rna_use.obs['to_opc'] = rna_use.obsm['to_terminal_states']['opc'].X
rna_use.obs['to_nn'] = rna_use.obsm['to_terminal_states']['nn'].X
rna_use.obs['to_chp'] = rna_use.obsm['to_terminal_states']['chp'].X

p = sc.pl.scatter(rna_use, basis='cssumap', color=['to_ctx','to_dien','to_rhom','to_mesen','to_retina','to_astro','to_opc','to_nn','to_chp'], show=False)
p[0].figure.savefig('plots/RNA/rna_velo_cellrank_probs_umap.png')


rna_use.obs['to_ctx_ranks'] = rna_use.obs['to_ctx'].rank() / rna_use.obs['to_ctx'].rank().max()
rna_use.obs['to_dien_ranks'] = rna_use.obs['to_dien'].rank() / rna_use.obs['to_dien'].rank().max()
rna_use.obs['to_rhom_ranks'] = rna_use.obs['to_rhom'].rank() / rna_use.obs['to_rhom'].rank().max()
rna_use.obs['to_mesen_ranks'] = rna_use.obs['to_mesen'].rank() / rna_use.obs['to_mesen'].rank().max()
rna_use.obs['to_retina_ranks'] = rna_use.obs['to_retina'].rank() / rna_use.obs['to_retina'].rank().max()
rna_use.obs['to_astro_ranks'] = rna_use.obs['to_astro'].rank() / rna_use.obs['to_astro'].rank().max()
rna_use.obs['to_opc_ranks'] = rna_use.obs['to_opc'].rank() / rna_use.obs['to_opc'].rank().max()
rna_use.obs['to_nn_ranks'] = rna_use.obs['to_nn'].rank() / rna_use.obs['to_nn'].rank().max()
rna_use.obs['to_chp_ranks'] = rna_use.obs['to_chp'].rank() / rna_use.obs['to_chp'].rank().max()

p = sc.pl.scatter(rna_use, basis='cssumap', color=['to_ctx_ranks','to_dien_ranks','to_rhom_ranks','to_mesen_ranks','to_retina_ranks','to_astro_ranks','to_opc_ranks','to_nn_ranks','to_chp_ranks'], show=False)
p[0].figure.savefig('plots/RNA/rna_velo_cellrank_probs_ranks_umap.png')

cellrank_meta = rna_use.obs[[
    'velocity_pseudotime', 'pseudotime_ranks',
    'to_ctx_ranks', 'to_dien_ranks', 'to_rhom_ranks', 'to_mesen_ranks', 'to_retina_ranks', 'to_astro_ranks', 
    'to_opc_ranks', 'to_nn_ranks', 'to_chp_ranks', 'to_ctx', 'to_dien', 'to_rhom',
    'to_mesen', 'to_retina', 'to_astro', 'to_opc', 'to_nn', 'to_chp']]
cellrank_meta.to_csv('data/RNA/cellrank/RNA_full_cellrank_probs.tsv', sep='\t')

trans_mat = combined_kernel._transition_matrix
sp.io.mmwrite('data/RNA/cellrank/velo_cr_transition.mtx', trans_mat)


#### PAGA with trans probs ####
rna_use.obs['clusters'] = rna_use.obs['clusters'].astype(str).astype('category')
scv.tl.paga(
    rna_use,
    groups='clusters',
    use_time_prior='velocity_pseudotime'
)

scv.pl.paga(rna_use, basis='cssumap')
con_mat = rna_use.uns['paga']['connectivities']
sp.io.mmwrite('data/RNA/cellrank/velo_paga_connectivities.mtx', con_mat)

rna_use.obs = rna_use.obs.drop('terminal_states_probs', axis=1)
rna_use.write('data/RNA/RNA_all_v3.2cellrank.h5ad')







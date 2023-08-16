import os
import numpy as np
import scanpy as sc
import scvelo as scv
import cellrank as cr
from scipy.io import mmwrite

os.chdir('/home/jfleck/projects/cutntag/')

adata_d158_A395 = sc.read('data/drugs/drugs_d15_d18_A395_v2velo_v1.2lines.h5ad')


sc.pp.neighbors(adata_d158_A395, use_rep='X_harmony', n_neighbors=20)
sc.tl.louvain(adata_d158_A395, neighbors_key='neighbors')


adata_d158_A395.obs['this'] = adata_d158_A395.obs['louvain'] == '14'
sc.pl.scatter(adata_d158_A395, basis='humap', color='this')


# 16 - psc
# 8 - nne
# 14 - neurons
# 1 - neural crest


adata_d158_A395.obs['terminal_states'] = np.nan
# adata_d158_A395.obs['terminal_states'][adata_d158_A395.obs.louvain=='16'] = 'psc'
adata_d158_A395.obs['terminal_states'][adata_d158_A395.obs.louvain=='8'] = 'nne'
adata_d158_A395.obs['terminal_states'][adata_d158_A395.obs.louvain=='14'] = 'neurons'
adata_d158_A395.obs['terminal_states'][adata_d158_A395.obs.louvain=='1'] = 'neural_crest'


scv.pl.velocity_embedding_stream(adata_d158_A395, basis='humap', color='terminal_states')
p = scv.pl.velocity_embedding_stream(adata_d158_A395, basis='humap', color='terminal_states', show=False)
p.figure.savefig('plots/drugs/drugs_d15_d18_terminal_velocity_umap.png')


#### Run cellrank ####
terminal_states = adata_d158_A395.obs['terminal_states']

vk = cr.tl.kernels.VelocityKernel(adata_d158_A395)
vk.compute_transition_matrix()

pk = cr.tl.kernels.PalantirKernel(adata_d158_A395, time_key='velocity_pseudotime')
pk.compute_transition_matrix()

ck = cr.tl.kernels.ConnectivityKernel(adata_d158_A395)
ck.compute_transition_matrix()

combined_kernel = vk
combined_kernel = ck
combined_kernel = pk
combined_kernel = 0.3 * vk + 0.7 * pk 

g = cr.tl.estimators.GPCCA(combined_kernel)
g.set_terminal_states(terminal_states)
g.compute_absorption_probabilities(n_jobs=40, use_petsc=True)


# adata_d158_A395.obs['to_psc'] = adata_d158_A395.obsm['to_terminal_states']['psc'].X
adata_d158_A395.obs['to_nne'] = adata_d158_A395.obsm['to_terminal_states']['nne'].X
adata_d158_A395.obs['to_neurons'] = adata_d158_A395.obsm['to_terminal_states']['neurons'].X
adata_d158_A395.obs['to_neural_crest'] = adata_d158_A395.obsm['to_terminal_states']['neural_crest'].X


# adata_d158_A395.obs['to_psc_ranks'] = adata_d158_A395.obs['to_psc'].rank() / adata_d158_A395.obs['to_psc'].rank().max()
adata_d158_A395.obs['to_nne_ranks'] = adata_d158_A395.obs['to_nne'].rank() / adata_d158_A395.obs['to_nne'].rank().max()
adata_d158_A395.obs['to_neurons_ranks'] = adata_d158_A395.obs['to_neurons'].rank() / adata_d158_A395.obs['to_neurons'].rank().max()
adata_d158_A395.obs['to_neural_crest_ranks'] = adata_d158_A395.obs['to_neural_crest'].rank() / adata_d158_A395.obs['to_neural_crest'].rank().max()


sc.pl.scatter(adata_d158_A395, basis='humap', color=['to_nne','to_neurons','to_neural_crest'])
# sc.pl.scatter(adata_d158_A395, basis='humap', color=['to_psc','to_nne','to_neurons','to_neural_crest'])
sc.pl.scatter(adata_d158_A395, basis='humap', color=['to_nne_ranks','to_neurons_ranks','to_neural_crest_ranks'])
sc.pl.scatter(adata_d158_A395, basis='humap', color=['to_psc_ranks','to_nne_ranks','to_neurons_ranks','to_neural_crest_ranks'])


p = sc.pl.scatter(adata_d158_A395, basis='humap', color=['to_psc_ranks','to_nne_ranks','to_neurons_ranks','to_neural_crest_ranks'], show=False)
p[0].figure.savefig('plots/drugs/drugs_d15_d18_terminal_probs_ranks_umap.png')


cellrank_meta = adata_d158_A395.obs[[
    'velocity_pseudotime',
    'to_psc_ranks','to_nne_ranks','to_neurons_ranks','to_neural_crest_ranks', 'to_psc','to_nne','to_neurons','to_neural_crest']]
cellrank_meta.to_csv('data/drugs/drugs_d15_d18_A395_cellrank_probs.tsv', sep='\t')

trans_mat = combined_kernel._transition_matrix
mmwrite('data/drugs/drugs_d15_d18_A395_velo_cr_transition.mtx', trans_mat)






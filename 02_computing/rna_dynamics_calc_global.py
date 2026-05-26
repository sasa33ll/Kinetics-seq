# %%
import os
import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
import warnings

# %%
adata= sc.read_h5ad('all_adata.h5ad')

# %%
t_label = 1.5
cluster_col = 'celltype1' 

# 1. adata
if type(adata.X) is np.ndarray:
    total_matrix = adata.X
else:
    total_matrix = adata.X.toarray()

if type(adata.layers['new']) is np.ndarray:
    new_matrix = adata.layers['new']
else:
    new_matrix = adata.layers['new'].toarray()
total_raw = pd.DataFrame(total_matrix, index=adata.obs_names, columns=adata.var_names).T
new_raw = pd.DataFrame(new_matrix, index=adata.obs_names, columns=adata.var_names).T

# 2. Library Size Normalization
lib_sizes = total_raw.sum(axis=0) 
lib_sizes[lib_sizes == 0] = 1 
scale_factors = 10000 / lib_sizes

total_norm = total_raw.multiply(scale_factors, axis=1)
new_norm = new_raw.multiply(scale_factors, axis=1)

# 3. k-NN
total_smooth_list = []
new_smooth_list = []
cell_order = []

for cluster in adata.obs[cluster_col].unique():
    cells_in_cluster = adata.obs[adata.obs[cluster_col] == cluster].index
    n_cells = len(cells_in_cluster)
    
    if n_cells <= 1:
        total_smooth_list.append(total_norm[cells_in_cluster])
        new_smooth_list.append(new_norm[cells_in_cluster])
        cell_order.extend(cells_in_cluster)
        continue
    ##
    sub_total = total_norm[cells_in_cluster]
    sub_new = new_norm[cells_in_cluster]
    ## PCA
    X_for_pca = np.log1p(sub_total.fillna(0)).T
    scaler = StandardScaler()
    X_for_pca_scaled = scaler.fit_transform(X_for_pca)
    curr_n_neighbors = min(10, n_cells)
    n_comps = min(30, n_cells - 1)
    pca = PCA(n_components=min(30, X_for_pca_scaled.shape[0]-1))
    X_pca = pca.fit_transform(X_for_pca_scaled)
    ##
    nbrs = NearestNeighbors(n_neighbors=curr_n_neighbors).fit(X_pca)
    indices = nbrs.kneighbors(X_pca, return_distance=False)
    ## 
    sub_total_smooth_np = np.array([sub_total.values[:, nn].mean(axis=1) for nn in indices]).T
    sub_new_smooth_np = np.array([sub_new.values[:, nn].mean(axis=1) for nn in indices]).T

    total_smooth_list.append(pd.DataFrame(sub_total_smooth_np, index=sub_total.index, columns=cells_in_cluster))
    new_smooth_list.append(pd.DataFrame(sub_new_smooth_np, index=sub_new.index, columns=cells_in_cluster))
    cell_order.extend(cells_in_cluster)

# 4. merge
total_smooth_unsorted = pd.concat(total_smooth_list, axis=1)
new_smooth_unsorted = pd.concat(new_smooth_list, axis=1)
total_smooth = total_smooth_unsorted.reindex(columns=adata.obs_names)
new_smooth = new_smooth_unsorted.reindex(columns=adata.obs_names)

# 5. Half-life

R = new_smooth / total_smooth
keep_mask = (total_smooth >= 3)
R_final = R.copy()
R_final[R_final <= 0] = np.nan
##
R_final[keep_mask & (R_final >= 1)] = 0.99
##
R_final[(~keep_mask) & (R_final >= 1)] = np.nan
R = R_final.copy()
##
half_life_raw = -t_label * np.log(2) / np.log(1 - R)

# 6. Gamma & Alpha
gamma = np.log(2) / half_life_raw
alpha = (new_smooth * gamma) / (1 - np.exp(-gamma * t_label))

true_zero_mask = (total_smooth <= 1e-6)
alpha[true_zero_mask] = 0.0

# 7. save
adata.layers['half_life'] = half_life_raw.T.values
adata.layers['alpha'] = alpha.T.values
adata.layers['gamma'] = gamma.T.values
# alpha_obs
alpha_matrix = adata.layers['alpha'].copy()
alpha_matrix[alpha_matrix < 0] = np.nan 
adata.obs['sum_alpha'] = np.nansum(alpha_matrix, axis=1) 
alpha_matrix_fill0 = np.nan_to_num(alpha_matrix, nan=0.0)
adata.obs['alpha'] = np.mean(alpha_matrix_fill0, axis=1)
# gamma_obs
gamma_matrix = adata.layers['gamma'].copy()
gamma_matrix[gamma_matrix <= 0] = np.nan 
gamma_matrix_fill0 = np.nan_to_num(gamma_matrix, nan=0.0)
adata.obs['gamma'] = np.mean(gamma_matrix_fill0, axis=1)


adata.write_h5ad('all_kinetics_adata.h5ad')



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
adata= sc.read_h5ad('cancercell_kinetics_adata.h5ad')

# %%
t_label = 1.5

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

n_cells = total_norm.shape[1]

X_for_pca = np.log1p(total_norm.fillna(0)).T
X_scaled = StandardScaler().fit_transform(X_for_pca)
X_pca = PCA(n_components=min(30, n_cells - 1)).fit_transform(X_scaled)

nbrs = NearestNeighbors(n_neighbors=min(10, n_cells)).fit(X_pca)
indices = nbrs.kneighbors(X_pca, return_distance=False)

total_smooth_np = np.array([total_norm.values[:, nn].mean(axis=1) for nn in indices]).T
new_smooth_np = np.array([new_norm.values[:, nn].mean(axis=1) for nn in indices]).T

total_smooth = pd.DataFrame(total_smooth_np, index=total_norm.index, columns=total_norm.columns)
new_smooth = pd.DataFrame(new_smooth_np, index=new_norm.index, columns=new_norm.columns)


# 4. Half-life

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

# 5. Gamma & Alpha
gamma = np.log(2) / half_life_raw
alpha = (new_smooth * gamma) / (1 - np.exp(-gamma * t_label))

true_zero_mask = (total_smooth <= 1e-6)
alpha[true_zero_mask] = 0.0

# 6. save
adata.layers['half_life'] = half_life_raw.T.values
adata.layers['total_smooth'] = total_smooth.T.values

# half_life
hl_matrix = adata.layers['half_life'].copy()
hl_matrix[hl_matrix <= 0] = np.nan

# var_aditional_strict
var_calc_threshold = 3
high_conf_mask = (total_smooth.T >= var_calc_threshold).values
hl_matrix_strict = np.where(high_conf_mask, hl_matrix, np.nan)
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)
    adata.var['half_life_gene_strict'] = np.nanmean(hl_matrix_strict, axis=0)

adata.write_h5ad('cancercell_kinetics_adata_add.h5ad')



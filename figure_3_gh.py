#%%
import sys
from scrna_utils import *
import os
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from sklearn.cluster import k_means
#%%
# Settings
os.chdir('/home2/s190548/work_mu/su_deng_rnaseq/new_dataset_for_rebuttal')
output_folder = 'results'
if not os.path.exists(output_folder):
    os.mkdir(output_folder)
sc.settings.figdir = output_folder
sns.set_theme(font='Arial', font_scale=1.5)
np.random.seed(1)
plt.rcParams["legend.markerscale"] = 2
sns.mpl.rcParams.update({'font.family':'Arial'})
#%%
# ==============================================================================
# ==============================================================================
# ==============================================================================

# Data preprocessing

genes = pd.read_csv('genes.csv').stack().dropna().unique()
pts_dfs = [x for x in os.listdir('cb_jiaotong_single') if 'exprs' in x]
cellmeta = pd.read_csv(
    os.path.join('cb_jiaotong_single', 'cell_meta.csv'), index_col=0
)

# Process p1
pt = pts_dfs[0]
df = pd.read_csv(
    os.path.join('cb_jiaotong_single', pt),
    sep = '\t')
df = df.set_index('Symbol')
df = df.iloc[:,1:].T
pt_id = pt.split('_')[1][1]
cell_id = cellmeta[
    (cellmeta['orig.ident'] == 'patient #{}'.format(pt_id)) & 
    (cellmeta.index.isin(df.index)) &
    (cellmeta.CellType == 'Epithelial cell')]
df = df.loc[cell_id.index]
adata_p1 = sc.AnnData(df)
adata_p1 = cal_qc(adata_p1)
adata_p1 = normalize_adata(adata_p1)
sc.pp.highly_variable_genes(adata_p1, min_mean=0.01, min_disp=0.5)
sc.pp.scale(adata_p1)
rps_genes = [x for x in adata_p1.var_names if x[:3] in ['RPS','RPL']]
mt_genes = [x for x in adata_p1.var_names if x[:3] == 'MT-']
adata_p1.var.loc[rps_genes+mt_genes,'highly_variable'] = False
sc.tl.pca(adata_p1,n_comps=100,random_state=0,use_highly_variable=True)

# Process p5
pt = pts_dfs[5]
df = pd.read_csv(
    os.path.join('cb_jiaotong_single', pt),
    sep = '\t')
df = df.set_index('Symbol')
df = df.iloc[:,1:].T
pt_id = pt.split('_')[1][1]
cell_id = cellmeta[
    (cellmeta['orig.ident'] == 'patient #{}'.format(pt_id)) & 
    (cellmeta.index.isin(df.index)) &
    (cellmeta.CellType == 'Epithelial cell')]
df = df.loc[cell_id.index]
adata_p5 = sc.AnnData(df)
adata_p5 = cal_qc(adata_p5)
adata_p5 = normalize_adata(adata_p5)
sc.pp.highly_variable_genes(adata_p5, min_mean=0.01, min_disp=0.5)
sc.pp.scale(adata_p5)
rps_genes = [x for x in adata_p5.var_names if x[:3] in ['RPS','RPL']]
mt_genes = [x for x in adata_p5.var_names if x[:3] == 'MT-']
adata_p5.var.loc[rps_genes+mt_genes,'highly_variable'] = False
sc.tl.pca(adata_p5,n_comps=20,random_state=0,use_highly_variable=True)
#%%
# ==============================================================================
# ==============================================================================
# ==============================================================================

# Figure 2g

# define gene lists
STEM = ['CD55','TACSTD2','KRT4','ATXN1']
AR_score = ['KLK3','PTGER4','ACSL3']
Jak_score = ['JAK1' ,'STAT1', 'IL6ST']
EMT_score = ['VIM','SNAI2', 'CDH11']
TP53 = ['TP53','RB1']
NEPC = ['ASCL1', 'SYP', 'CHGA','SPATS2L', 'ACTN4']
NEPC = [x for x in NEPC if x in set(adata_p1.var_names) & set(adata_p5.var_names)]

geneset = {
    'STEM' : STEM,
    'AR_score':AR_score,
    'Jak_score':Jak_score,
    'EMT_score':EMT_score,
    'TP53_score':TP53,
    'NEPC': NEPC,
    }

# PCA plots
for item in geneset.items():
    for adata, pt_name in zip([adata_p1, adata_p5],['P1','P5']):
        geneset_vector = adata.raw.to_adata().to_df()[item[1]].mean(axis=1)
        adata.obs[item[0]] = geneset_vector
        titles = [pt_name + '_' + x for x in ([item[0]] + item[1])]
        vmax = None
        if item[0] == 'AR_score':
            vmax= 1.5
        sc.pl.pca(
            adata, color = [item[0]] + item[1], s=100, ncols = 2,
            alpha = 0.75,
            show=False,
            title = titles, 
            save = 'fig_3g_' + pt_name + '_' + item[0],
            vmax = vmax
        )
#%%
# ==============================================================================
# ==============================================================================
# ==============================================================================

# Figure 2h

# Kmeans clustering to define TP53 high and low populations.
clusters = k_means(adata_p5.obsm['X_pca'],2)[1]
clusters = ['TP53/RB1-low' if x == 1 else 'TP53/RB1-high' for x in clusters]
adata_p5.obs['cluster'] = clusters
clusters = k_means(adata_p1.obsm['X_pca'],2)[1]
clusters = ['TP53/RB1-low' if x == 1 else 'TP53/RB1-high' for x in clusters]
adata_p1.obs['cluster'] = clusters

# Violin plots and T-tests

for geneset_key in geneset.keys():
    genes = geneset[geneset_key]
    p1_AR = adata_p1.raw.to_adata()[:,genes].to_df().mean(axis=1).to_frame(name=geneset_key)
    p1_AR['Patient'] = 'Patient#1 (CRPC-Adeno)'
    p1_AR['Cluster'] = adata_p1.obs.cluster
    p5_AR = adata_p5.raw.to_adata()[:,genes].to_df().mean(axis=1).to_frame(name=geneset_key)
    p5_AR['Patient'] = 'Patient#5 (CRPC-NE)'
    p5_AR['Cluster'] = adata_p5.obs.cluster
    plot_data = p1_AR.append(p5_AR)
    col = plot_data.columns[0]
    title = ''
    p1, p2, p3, p4 = plot_data.groupby(['Patient','Cluster'])
    pv13 = ttest_ind(p1[1][col].values, p3[1][col].values)[1]
    title += 'p1 TP53/RB1-high vs p5 TP53/RB1-high: {:.3e}\n'.format(pv13)
    pv12 = ttest_ind(p1[1][col].values, p2[1][col].values)[1]
    title += 'p1 TP53/RB1-high vs p1 TP53/RB1-low: {:.3e}\n'.format(pv12)
    pv24 = ttest_ind(p2[1][col].values, p4[1][col].values)[1]
    title += 'p1 TP53/RB1-low vs p5 TP53/RB1-low: {:.3e}\n'.format(pv24)
    pv34 = ttest_ind(p3[1][col].values, p4[1][col].values)[1]
    title += 'p5 TP53/RB1-high vs p5 TP53/RB1-low: {:.3e}\n'.format(pv34)

    _ = plt.figure(figsize=(8,6))
    sns.violinplot(
        data = plot_data, x = 'Patient', y = geneset_key, hue = 'Cluster', cut=0,
        inner='quartile', scale = 'width')
    plt.legend(bbox_to_anchor = (1, 0.5), loc = 'center left')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_folder + '/Figure_3h_{}_violin.pdf'.format(geneset_key))
#%%
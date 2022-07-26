#%%
import os
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from sklearn.preprocessing import minmax_scale
from collections import defaultdict
from cool_plots import plot_radar
from collections import defaultdict
from scutils import *

#%%
# Settings
os.chdir('/home2/s190548/work_mu/su_deng_rnaseq/singlecell')
output_path = './results'
if not os.path.exists(output_path):
    os.mkdir(output_path)
sc.settings.figdir = output_path
np.random.seed(1)
plt.rcParams["legend.markerscale"] = 2
#%%
# ==============================================================================
# ==============================================================================
# ==============================================================================

# preprocessing all data

adata = sc.read_10x_h5('./aggr/outs/count/filtered_feature_bc_matrix.h5')
cell_types = pd.read_csv('cell_types.csv', index_col=0)
for i, sample_name in enumerate([
    'sgNT-Veh','sgTP53/RB1-Veh','sgTP53/RB1/JAK1-Veh','sgNT-Enz',
    'sgTP53/RB1-Enz','sgTP53/RB1/JAK1-Enz']):
    batch_idx = [x for x in adata.obs_names if str(i+1) in x]
    adata.obs.loc[batch_idx,'Sample'] = sample_name

adata = cal_qc(adata)
adata.var_names_make_unique()
adata = mito_qc(
    adata, min_genes=100, max_genes=7000, percent_mito_cutoff=0.15)
adata = normalize_adata(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.01, min_disp=0.5)

sc.pp.scale(adata)
rps_genes = [x for x in adata.var_names if x[:3] in ['RPS','RPL']]
mt_genes = [x for x in adata.var_names if x[:3] == 'MT-']
adata.var.loc[rps_genes+mt_genes,'highly_variable'] = False
sc.tl.pca(adata,n_comps=50,random_state=0,use_highly_variable=True)
sc.pp.neighbors(adata, n_pcs=25, n_neighbors=15)
sc.tl.umap(adata, min_dist=0.5)

adata.obs['leiden'] = cell_types['leiden']
adata.obs['leiden_sub'] = cell_types['leiden_sub']

# Cell cycle scores
cell_cycle_genes = [x.strip() for x in open('regev_lab_cell_cycle_genes.txt')]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

# Save cell level annotations
adata.obs.to_csv('./results/Cell_annotations_full.csv')
#%%
# ==============================================================================
# ==============================================================================
# ==============================================================================

# Figure 6 A, B, C

conds = [
    ['sgNT-Veh','sgNT-Enz'],
    ['sgTP53/RB1-Veh','sgTP53/RB1-Enz'],
    ['sgTP53/RB1/JAK1-Veh','sgTP53/RB1/JAK1-Enz']
    ] # Evaluate each genotype with and without ENZ.
for c in ['leiden','Sample','phase']:
    for cond in conds:
        sc.pl.umap(
            adata[adata.obs['Sample'].isin(cond)], 
            color=c, 
            save='Figure_6_abc_{}_'.format(
                cond[0].split('-')[0].replace('/','_')
                ) + c, 
            show=False)

# Figure 6 E, G and H.
for c in ['leiden','Sample','phase']:
    sc.pl.umap(
        adata,
        color=c, save='Figure_6_egh' + c, show=False, alpha=1)

# Figure 6F
colors = pd.Series('',index = adata.obs.leiden_sub.astype(str).unique())
blues = sns.light_palette('Blue',4).as_hex()[1:3]
colors[['0','1']] = blues
greens = sns.light_palette('Green',5).as_hex()[1:4]
colors[['2-1','2-2','2-3']] = greens
reds = sns.light_palette('Red',5).as_hex()[1:4]
colors[['3-1','3-2','3-3']] = reds
purples = sns.light_palette('Purple',5).as_hex()[1:4]
colors[['4-1','4-2','4-3']] = purples
browns = sns.light_palette('Orange',4).as_hex()[1:3]
colors[['5-1','5-2']] = browns

_ = plt.figure(figsize=(12,12))
umap_x, umap_y = adata.obsm['X_umap'][:,0], adata.obsm['X_umap'][:,1]
sns.scatterplot(
    x = umap_x, y = umap_y, hue=adata.obs.leiden_sub.astype(str).values, 
    palette=colors.to_dict(), hue_order=colors.sort_index().index,
    linewidth=0, s=2)
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.legend(bbox_to_anchor = (1,0.5), loc = 'center left')
plt.savefig(output_path + '/Figure_6F.pdf', bbox_inches = 'tight')
plt.close()

# Figure S10 b to f
for gene in ['FKBP5', 'NKX3-1', 'PTGER4', 'HERC3', 'GNMT']:
    sc.pl.umap(
        adata,
        color=gene, save='Figure_S10_b_to_f{}.pdf'.format(gene), show=False, alpha=1)

# Figure S10 g to r
c4adata = adata[adata.obs.leiden == '4']
for gene in ['JAK1','JAK2','JAK3','TYK2', 'STAT1','STAT2', 'STAT3', 'IL6', 'IL6R', 'IL6ST']:
    sc.pl.umap(
        c4adata,
        color=gene, save='Figure_S10_g_to_r{}.pdf'.format(gene), show=False, alpha=1)
#%%
# ==============================================================================
# ==============================================================================
# ==============================================================================

# Figure 7

# Figure 7a
# Read genesets related to PCA.
pc_genes = pd.read_excel(
    'PCA scores alternative.xlsx')
ar_genes = sorted(pc_genes.iloc[:,3].dropna().unique())
basal_genes = sorted(pc_genes.iloc[:,0].dropna().unique())
luminal_genes = sorted(pc_genes.iloc[:,1].dropna().unique())
nepc_genes = sorted(pc_genes.iloc[:,2].dropna().unique())
emt_genes = sorted(pc_genes.iloc[:,5].dropna().unique())
stem_genes = sorted(pc_genes.iloc[:,4].dropna().unique())
geneset_dict = defaultdict()
for geneset, geneset_name in zip(
    [ar_genes, basal_genes ,luminal_genes ,nepc_genes ,emt_genes,stem_genes],
    ['AR score','BASAL','Luminal','NEPC','EMT', 'Stem']
    ):
    geneset_dict[geneset_name] = geneset

# Calculate gene scores.
for geneset_name, geneset in geneset_dict.items():
    sc.tl.score_genes(
        adata, geneset, score_name=geneset_name, ctrl_size=100)

# Plotting
gene_set_scores = adata.obs.loc[:,'AR score':'Stem']
gene_set_scores['Cluster'] = adata.obs.loc[gene_set_scores.index,'leiden'].values
gene_set_scores = gene_set_scores.groupby(['Cluster']).mean().fillna(0)
gene_set_scores.loc[:,:] = minmax_scale(gene_set_scores, axis=0)
g = sns.heatmap(
    gene_set_scores,cmap=sns.palettes.color_palette('YlGnBu_r', as_cmap=True))
plt.xlabel('')
plt.savefig(
    output_path + '/Figure 7a.pdf', bbox_inches='tight')

# Figure 6 B and C

cell_scores = adata.obs.loc[:,'BASAL':'Stem']
cell_scores = (cell_scores - cell_scores.min().min()) / np.ptp(cell_scores)
cluster_scores = cell_scores.groupby(adata.obs['leiden']).mean()
sample_scores = cell_scores.groupby(adata.obs['Sample']).mean()
cluster_scores.loc[:,:] = minmax_scale(cluster_scores, (0.05,1), axis=0)
sample_scores.loc[:,:] = minmax_scale(sample_scores, (0.05,1), axis=0)
plot_radar(cluster_scores, palette='tab10', figname=output_path + '/Figure_7b.pdf')
plot_radar(sample_scores, palette='tab20', figname=output_path + '/Figure_7c.pdf')

# Figure 6 d to i.
for c in ['AR score', 'BASAL', 'Luminal', 'NEPC', 'EMT', 'Stem']:
    plot_adata = adata.copy()
    # Coloring based on scaled score from 0.05 to 0.95 quantile.
    _scores = adata.obs[c].copy()
    _scores = np.clip(
        _scores, np.percentile(_scores, 5), np.percentile(_scores, 95))
    _scores = minmax_scale(_scores)
    plot_adata.obs[c] = _scores
    sc.pl.umap(
        plot_adata, color = c,save='Figure_7d_to_i' + c,
        show=False, legend_fontsize = 'large', color_map ='rainbow', s=25, alpha=0.75)
#%% 
# ==============================================================================
# ==============================================================================
# ==============================================================================

# Prep data for monocle

hvg = adata.var_names[adata.var.highly_variable]
adata.raw.to_adata().to_df()[hvg].round(2).to_csv(
    output_path + '/monocle_expr_raw_hvg.txt', sep='\t')
for col in adata.obs.columns[-8:]:
    norm_vals = np.clip(
        adata.obs[col], 
        a_min=np.percentile(adata.obs[col],0), 
        a_max=np.percentile(adata.obs[col],95))
    norm_vals = minmax_scale(norm_vals, (0,1))
    adata.obs[col] = norm_vals
adata.obs.to_csv(output_path + '/cell_meta.txt', sep = '\t')
# expr for just important geens
pc_genes = pd.read_excel('PCA scores alternative.xlsx')
pc_genes = pc_genes.stack().dropna().values
pc_genes = [x for x in set(pc_genes) if x in adata.var_names]
expr_pc_genes = adata.raw.to_adata().to_df()[pc_genes]
expr_pc_genes.round(2).to_csv(
    output_path + '/monocle_expr_pc_genes.txt', sep='\t')

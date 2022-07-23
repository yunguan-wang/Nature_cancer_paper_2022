#%%
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def cal_qc(adata):
    adata.var_names_make_unique()
    mito_genes = adata.var_names.str.startswith("mt-")
    mito_genes |= adata.var_names.str.startswith("MT-")
    all_genes = adata.var_names.values
    rps_genes = [
        x for x in all_genes if 
        ('rps' == x.lower()[:3]) |
        ('rpl' == x.lower()[:3]) |
        ('rps' == x.lower()[:3]) |
        ('rpl' == x.lower()[:3])]
    # for each cell compute fraction of counts in mito genes vs. all genes
    # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
    adata.obs["percent_mito"] = np.sum(
        adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1
    )
    adata.obs["percent_ribosomal"] = np.sum(
        adata[:, rps_genes].X, axis=1) / np.sum(adata.X, axis=1
    )
    # add the total counts per cell as observations-annotation to adata
    adata.obs["n_counts"] = adata.X.sum(axis=1)
    adata.obs["n_genes"] = np.sum(adata.X>0,axis=1)
    return adata

def plot_qc(adata, mito_high = None):
    fig = plt.figure(figsize=(18, 12))
    axes = fig.subplots(3, 2)
    axes = axes.ravel()
    sns.histplot(adata.obs["n_genes"], ax=axes[0], label="Number of genes")
    sns.histplot(adata.obs["n_counts"], ax=axes[1], label="Count depth")
    axes[0].set_title("Number of genes")
    axes[1].set_title("Count_depth")
    sns.histplot(adata.obs["percent_mito"], ax=axes[2])
    sns.histplot(adata.obs["percent_ribosomal"], ax=axes[3])
    if mito_high is None:
        norm = plt.Normalize(0, adata.obs["percent_mito"].max())
    else:
        norm = plt.Normalize(0, mito_high)
    sns.scatterplot(
        adata.obs["n_genes"],
        adata.obs["n_counts"],
        hue=adata.obs["percent_mito"],
        s=2,
        linewidth=0,
        alpha=0.6,
        ax=axes[4],
        legend=None,
        palette="gist_heat",
        hue_norm=norm,
    )
    sm = plt.cm.ScalarMappable(cmap="gist_heat", norm=norm)
    sm.set_array([])
    sns.scatterplot(
        adata.obs["percent_mito"],
        adata.obs["percent_ribosomal"],
        hue=adata.obs["percent_mito"],
        s=2,
        linewidth=0,
        alpha=0.6,
        ax=axes[5],
        legend=None,
        palette="gist_heat",
        hue_norm=norm,
    )
    axes[5].figure.colorbar(sm)
    plt.tight_layout()

def mito_qc(adata, min_genes=100, max_genes=4000, percent_mito_cutoff=0.5):
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_genes = max_genes)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs["percent_mito"] < percent_mito_cutoff, :]
    return adata
    # plt.savefig(pathout+'/QC.png')

def normalize_adata(adata):
    # sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    return adata
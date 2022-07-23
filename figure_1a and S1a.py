import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

os.chdir('/home2/s190548/work_mu/su_deng_rnaseq/results')
plt.rcParams["legend.markerscale"] = 2
plt.rcParams["font.size"] = 16
plt.rcParams["font.weight"] = 'bold'
for folder in ['run1','run2']:
    report = pd.DataFrame()
    for fn in [x for x in os.listdir(folder) if 'pathways.txt' in x]:
        _df = pd.read_csv(folder + '/' + fn, sep='\t', index_col=0)
        if 'shSox2' in fn:
            valid_pathways = _df[(_df.padj<.25)&(_df.NES<0)].pathway.values
            _df.NES = _df.NES
            label = 'Down'
        else:
            valid_pathways = _df[(_df.padj<.25)&(_df.NES>0)].pathway.values
            label = 'Up'
        _df = _df.set_index('pathway')[['NES','padj']]
        _df.columns = [fn.split(' ')[0] + '_' + label + '_' + x for x in _df.columns]
        report = pd.concat([report, _df.loc[valid_pathways]],axis=1)
    report = report[(report.iloc[:,[1,3,5]]<.25).sum(axis=1)>=2]
    report = report.iloc[:,[0,2,4]]
    report.columns = [x[:-4] for x in report.columns]
    report = report.fillna(0)
    report.index = [x[5:] for x in report.index]
    if folder == 'run1':
        g = sns.clustermap(
            report.iloc[:,[0,2,1]], figsize=(24,12), cmap='coolwarm', col_cluster=False)
    else:
        g = sns.clustermap(report, figsize=(24,12), cmap='coolwarm', center=0,col_cluster=False)
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=15)
    plt.savefig('Fig_1a_{}.pdf'.format(folder), bbox_inches='tight')

c_pos, c_neg = sns.diverging_palette(0,240,n=2).as_hex()
for folder in ['run1','run2']:
    report = pd.DataFrame()
    for fn in [x for x in os.listdir(folder) if 'pathways.txt' in x]:
        _df = pd.read_csv(folder + '/' + fn, sep='\t', index_col=0)
        valid_pathways = _df[_df.padj<.25]
        valid_pathways.pathway = [x[5:] for x in valid_pathways.pathway]
        valid_pathways.index = np.arange(len(valid_pathways))
        colors = np.array([c_pos if x >0 else c_neg for x in valid_pathways.NES])
        _ = plt.figure(figsize=(8,24))
        sns.barplot(
            data=valid_pathways, y='pathway', x = 'NES', palette=colors)
        plt.title(fn[:-4])
        jak_idx = valid_pathways[valid_pathways.pathway=='JAK_STAT_SIGNALING_PATHWAY'].index
        if jak_idx.shape[0] > 0:
            jak_idx = jak_idx[0]
            ax = plt.gca()
            ax.get_yticklabels()[jak_idx].set_weight("bold")
        plt.savefig(fn.replace('.txt','.pdf'), bbox_inches = 'tight')
        plt.close()
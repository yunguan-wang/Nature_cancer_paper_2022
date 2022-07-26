library('monocle3')
library('ggplot2')
library(devtools)
library('dplyr')
library('ggplot2')
library('circlize')
library('RColorBrewer')
library('ComplexHeatmap')
library(viridis)

clip = function(array_x, pmin, pmax) {
  pmin = quantile(array_x, pmin)
  pmax = quantile(array_x, pmax)
  array_x[array_x < pmin] = pmin
  array_x[array_x > pmax] = pmax
  return (array_x)
}

# Read processed data
setwd('A:/su_deng_rnaseq/singlecell/results/revision')
expr = read.table('monocle_expr_raw_hvg.txt',sep='\t', row.names = 1, header=T)
meta = read.table('cell_meta.txt', sep = '\t', row.names = 1, header=T)
meta$cell_cycle = meta$S_score + meta$G2M_score

########## Calculating pseudotime ##########
set.seed(1)
geneanno = data.frame(colnames(expr),row.names=colnames(expr))
colnames(geneanno)=c('gene_short_name')
cds = new_cell_data_set(
  as.matrix(t(expr)),
  cell_metadata = meta,
  gene_metadata = geneanno)
cds = preprocess_cds(
  cds, num_dim = 100, norm_method = 'none', verbose = T)
cds = reduce_dimension(
  cds, max_components = 50, reduction_method = 'UMAP', umap.min_dist = 0.5,
  umap.fast_sgd=TRUE, umap.n_neighbors=200, cores = 4)
cds = cluster_cells(cds, k=100)
cds = learn_graph(
  cds, 
  close_loop = T, 
  verbose = T, 
  use_partition = T,
  )
cds = order_cells(cds)

# Copy additional metadata to CDS object

meta$cell_cycle = meta$S_score + meta$G2M_score
meta$S_score = clip(meta$S_score, 0.05, 0.95)
meta$G2M_score = clip(meta$G2M_score, 0.05, 0.95)
cds = readRDS('monocle.Rdata')
for (col in colnames(meta)) {
  colData(cds)[,col] = meta[,col]
}
colData(cds)$leiden = paste("cluster", colData(cds)$leiden, sep='_')
c4_cells = row.names(meta)[meta$leiden=='4']
c4_cds = cds[,c4_cells]
########## custom plots ##########
for (i in 1:2) {
  if (i == 1) {
    plot_cds = cds
    prefix = ''
    n_clusters = 6
    cluster_names = meta$leiden
  } else {
    plot_cds = c4_cds
    prefix = 'C4_'
    n_clusters = 3
    cluster_names = colData(plot_cds)$leiden_sub
  }
  gl_table = read.table('../PCA scores alternative.csv',sep=',',header = T)
  gl_table = gl_table[,1:6]
  cluster_colors = setNames(
    brewer.pal(name="Set1", n=n_clusters), sort(unique(cluster_names))
  )
  pseudotime_colors = colorRamp2(
    seq(min(pseudotime(plot_cds)), max(pseudotime(plot_cds)), length = 9), rev(brewer.pal(name="YlGnBu", n=9)))
  topanno = HeatmapAnnotation(
    Pseudotime = pseudotime(plot_cds)[order(pseudotime(plot_cds))],
    Cluster = cluster_names[order(pseudotime(plot_cds))],
    col = list(Pseudotime=pseudotime_colors,Cluster=cluster_colors)
  )
  for (gl in colnames(gl_table)) {
    genes = unique(gl_table[gl][,1])
    genes = genes[genes %in% row.names(expr)]
    if (length(genes) == 0) {next}
    plot_data = expr[genes, order(pseudotime(plot_cds))]
    plot_data = plot_data[rowSums(plot_data)>0,]
    plot_data = t(apply(plot_data,1,function(x){smooth.spline(x,df=3)$y}))
    plot_data = t(apply(plot_data,1,function(x){(x-mean(x))/sd(x)}))
    psut = pseudotime(plot_cds)/max(pseudotime(plot_cds))
    pdf(
      paste(
        prefix, 'pseudotime_', gl, '_heatmap.pdf',sep=''), height=length(genes)/4)
    g = Heatmap(
      plot_data,
      col = colorRamp2(seq(from=-2,to=2,length=9),rev(brewer.pal(9, "RdBu"))),
      show_row_names = TRUE,
      show_column_names  = FALSE,
      row_title_rot = 0,
      top_annotation  = topanno,
      cluster_row_slices = FALSE,
      cluster_columns = FALSE,
      cluster_row = T)
    print(g)
    dev.off()
  }
}
# All scatter plots

cols = c('#5D58A6','#00B2E4','#8ACCA0','#F9A45F','#ED2024')

cat_vals = c("leiden_sub","leiden","Sample")
cont_vals = c(
  "pseudotime","S_score",'G2M_score','cell_cycle',"AR.score", "AR.Sub", 
  "BASAL", "Luminal", "NEPC", "EMT", "Stem", "JAK.activation")

for (val in cat_vals) {
  pdf(file=paste('Monocle_', val,'.pdf', sep=''))
  g = plot_cells(
    cds, color_cells_by = val, group_label_size = 5,
    label_cell_groups = F,
    label_groups_by_cluster = F,
    show_trajectory_graph = F,
    cell_size = 0.35)
  print(g)
  dev.off()
}

for (val in cont_vals) {
  pdf(file=paste('Monocle_', val,'.pdf', sep=''))
  if (val != 'S_score') {
    g = plot_cells(
      cds, color_cells_by = val, group_label_size = 5,
      label_cell_groups = F,
      label_groups_by_cluster = F,
      show_trajectory_graph = F,
      cell_size = 0.5) +
      scale_color_gradientn(colors = cols, oob = scales::squish, name=val)
  } else {
    g = plot_cells(
      cds, color_cells_by = val, group_label_size = 5,
      label_cell_groups = F,
      label_groups_by_cluster = F,
      show_trajectory_graph = F,
      cell_size = 0.5)
  }
  print(g)
  dev.off()
  if (val %in% c("leiden_sub","pseudotime","AR.score", "AR.Sub", "BASAL", 
                 "Luminal", "NEPC", "EMT", "Stem", "JAK.activation","S_score")) {
    pdf(file=paste('Monocle_C4_', val,'.pdf', sep=''))
    if (val == 'pseudotime'){
        v_min = quantile(pseudotime(c4_cds), 0.05)
        v_max = quantile(pseudotime(c4_cds), 0.95)
    } else {
        v_min = quantile(colData(c4_cds)[,val], 0.05)
        v_max = quantile(colData(c4_cds)[,val], 0.95)
    }
    if (val != 'S_score') {
    g = plot_cells(
      c4_cds, color_cells_by = val, group_label_size = 5,
      label_cell_groups = F,
      label_groups_by_cluster = F,
      show_trajectory_graph = F,
      cell_size = 1) + xlim(-3,1) + ylim(0,2) + 
      scale_color_gradientn(colors = cols,limits = c(v_min, v_max), oob = scales::squish, name=val)
    } else {
      g = plot_cells(
        c4_cds, color_cells_by = val, group_label_size = 5,
        label_cell_groups = F,
        label_groups_by_cluster = F,
        show_trajectory_graph = F,
        cell_size = 1) + xlim(-3,1) + ylim(0,2)}
    print(g)
    dev.off()
  }
}


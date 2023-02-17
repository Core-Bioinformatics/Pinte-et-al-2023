# Cellranger alignment for each sample using cellranger 6.1.1
# cellranger count --id=E3 --transcriptome=refdata-cellranger-GRCh38-3.1.0 --fastqs=raw/ --sample=SITTE3 --include-introns --localcores=10 --localmem=100

# Post-alignment analysis using Seurat

library(Seurat)

# define metadata information
metadata = data.frame('ID'=c('E3','F3','G3','H3'),
                      'sample'=c('NO KO','HY KO','NO WT','HY WT'),
                      'NO_HY'=c('NO','HY','NO','HY'),
                      'WT_KO'=c('KO','KO','WT','WT'))

# create Seurat object with all 4 samples
filelist = list.files()
for (file in filelist){
  sample.so = Read10X_h5(paste0(file,'/outs/filtered_feature_bc_matrix.h5'))
  sample.so = CreateSeuratObject(sample.so,project = "LV_Laetitia")
  sample.so$ID = file
  sample.so$sample = metadata[metadata$ID==file,]$sample
  sample.so$NO_HY = metadata[metadata$ID==file,]$NO_HY
  sample.so$WT_KO = metadata[metadata$ID==file,]$WT_KO
  if (file == filelist[[1]]){
    so = sample.so
  } else {
    so = merge(so, y = c(sample.so), project = "LV_Laetitia")
  }
}
sample.so=so

# get MT and RP counts and add percentages
mt.genes=grep("^MT-", rownames(so), value=FALSE, ignore.case=TRUE)
rp.genes=grep("^RP[SL]", rownames(so), value=FALSE, ignore.case=TRUE)
so[['percent.mt']] = PercentageFeatureSet(so, features=rownames(so)[mt.genes])
so[['percent.rp']] = PercentageFeatureSet(so, features=rownames(so)[rp.genes])

#-----------------------------------------------------------------------------------

# QC plots
VlnPlot(so,features=c('nCount_RNA','nFeature_RNA','percent.mt','percent.rp'),pt.size = 0,ncol=4)
plot(so@meta.data$nCount_RNA, 
     so@meta.data$nFeature_RNA, 
     xlab = 'n total counts', 
     ylab='n unique features', 
     col = rgb(red = 0, green = 0, blue = 0, alpha = 0.1), 
     main='Scatterplot of features-seq depth')
plot(so@meta.data$percent.mt, 
     so@meta.data$percent.rp, 
     xlab = 'MT %', 
     ylab='RP %', 
     col = rgb(red = 0, green = 0, blue = 0, alpha = 0.1), 
     pch = 16, 
     main='Scatterplot of MT-RP')
VlnPlot(so,features=c('nCount_RNA','nFeature_RNA','percent.mt','percent.rp'),
        pt.size = 0,ncol=4,group.by='sample')
VlnPlot(so,features=c('nCount_RNA','nFeature_RNA','percent.mt','percent.rp'),
        pt.size = 0,ncol=4,group.by='NO_HY')
VlnPlot(so,features=c('nCount_RNA','nFeature_RNA','percent.mt','percent.rp'),
        pt.size = 0,ncol=4,group.by='WT_KO')
barplot(table(so$ID))

#---------------------------------------------------------------------------------

# Filter out low quality cells
so=subset(so,nFeature_RNA>2500 & nCount_RNA < 75000 & percent.mt < 20 & percent.rp < 25 & nCount_RNA > 5000)

# Remove MT/RP genes and perform sctransform
so = so[-c(mt.genes, rp.genes),]
so = SCTransform(so,return.only.var.genes = F)

# Check cell cycle scores
so = CellCycleScoring(so, cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)

#---------------------------------------------------------------------------------

# ClustAssess feature stability analysis
library(ggplot2)
library(ClustAssess)
run.feature.stability <- function(seurat.object){
  n_repetitions = 30
  n_cores = 10
  # look at dim reduction and feature selection
  var_features = VariableFeatures(seurat.object)
  n_abundant = 3000
  most_abundant_genes = rownames(seurat.object)[order(Matrix::rowSums(GetAssayData(seurat.object, assay='SCT', slot='counts')), decreasing=TRUE)]
  # get feature sets for comparisons
  steps = seq(from = 500, to = 3000, by = 500)
  ma_hv_genes_intersection_sets = sapply(steps, function(x) intersect(most_abundant_genes[1:x], var_features[1:x]))
  ma_hv_genes_intersection = Reduce(union, ma_hv_genes_intersection_sets)
  ma_hv_steps = sapply(ma_hv_genes_intersection_sets, length)
  
  # perform evaluation
  pca_feature_stability_object = c(get_feature_stability(data_matrix = seurat.object@assays[["SCT"]]@scale.data,
                                                         feature_set = most_abundant_genes,
                                                         steps = steps,
                                                         n_repetitions = n_repetitions,
                                                         feature_type = "MA",
                                                         graph_reduction_type = "PCA",
                                                         npcs = 30,
                                                         min_dist = 0.3,
                                                         n_neighbors = 30,
                                                         metric = "cosine",
                                                         ncores = n_cores,
                                                         ecs_thresh = 1,
                                                         algorithm = 3),
                                   get_feature_stability(data_matrix = seurat.object@assays[["SCT"]]@scale.data,
                                                         feature_set = var_features,
                                                         steps = steps,
                                                         n_repetitions = n_repetitions,
                                                         feature_type = "HV",
                                                         graph_reduction_type = "PCA",
                                                         npcs = 30,
                                                         min_dist = 0.3,
                                                         n_neighbors = 30,
                                                         metric = "cosine",
                                                         ncores = n_cores,
                                                         ecs_thresh = 1,
                                                         algorithm = 3),
                                   get_feature_stability(data_matrix = seurat.object@assays[["SCT"]]@scale.data,
                                                         feature_set = ma_hv_genes_intersection,
                                                         steps = ma_hv_steps,
                                                         n_repetitions = n_repetitions,
                                                         feature_type = "MA_HV",
                                                         graph_reduction_type = "PCA",
                                                         npcs = 30,
                                                         min_dist = 0.3,
                                                         n_neighbors = 30,
                                                         metric = "cosine",
                                                         ncores = n_cores,
                                                         ecs_thresh = 1,
                                                         algorithm = 3))
  return(pca_feature_stability_object)
}
pca_feature_stability_object = run.feature.stability(so)
plot_feature_stability_boxplot(pca_feature_stability_object, text_size  = 2.5) +
                  theme(legend.position = c(1,0),
                        legend.justification = c(1,0))
plot_feature_stability_ecs_incremental(pca_feature_stability_object, dodge_width = 1, text_size = 2) +
                  theme(legend.position = c(1,0),
                        legend.justification = c(1,0))
plot_feature_stability_mb_facet(pca_feature_stability_object, text_size = 3)
plot_feature_stability_ecs_facet(pca_feature_stability_object)

# based on ClustAssess analysis, highly variable 2500 genes were selected both for all samples together and for normoxia and hypoxia samples separately
so = RunPCA(so,npcs = 30,approx = F,verbose = F,features = head(VariableFeatures(so,2500)))
so = RunUMAP(so,reduction = "pca",dims = 1:30,n.neighbors = 30,min.dist = 0.3,metric = "cosine",verbose = F)

# The following analysis was performed separately for only normoxia and only hypoxia samples
listofseurat = list('NO'=subset(so,NO_HY=='NO'),
                    'HY'=subset(so,NO_HY=='HY'),
                    'WT'=subset(so,WT_KO=='WT'))

# ClustAssess NN analysis
run.nn.stability <- function(so.obj){
  nn_conn_comps_object = c(get_nn_conn_comps(object = so.obj@reductions$pca@cell.embeddings,
                                             n_neigh_sequence = c(c(1,2,3,4), seq(from = 5, to = 30, by = 5)),
                                             n_repetitions = 30,
                                             graph_reduction_type = "UMAP",
                                             ncores = 10,
                                             min_dist = 0.3,
                                             n_neighbors = 30,
                                             metric = "cosine"),
                           get_nn_conn_comps(object = so.obj@assays[["SCT"]]@scale.data,
                                             n_neigh_sequence = c(c(1,2,3,4), seq(from = 5, to = 30, by = 5)),
                                             n_repetitions = 30,
                                             graph_reduction_type = "PCA",
                                             ncores = 10,
                                             nv = 30))
  nn_importance_object = mapply(c,
                                get_nn_importance(object = so.obj@assays[["SCT"]]@scale.data,
                                                  n_neigh_sequence = seq(from = 5, to = 30, by = 5),
                                                  n_repetitions = 30,
                                                  graph_reduction_type = "PCA",
                                                  ecs_thresh = 1,
                                                  ncores = 10,
                                                  algorithm = 1,
                                                  nv = 30),
                                get_nn_importance(object = so.obj@reductions$pca@cell.embeddings,
                                                  n_neigh_sequence = seq(from = 5, to = 30, by = 5),
                                                  n_repetitions = 30,
                                                  graph_reduction_type = "UMAP",
                                                  ecs_thresh = 1,
                                                  ncores = 10,
                                                  algorithm = 1,
                                                  min_dist = 0.3,
                                                  n_neighbors = 30,
                                                  metric = "cosine"),
                                SIMPLIFY = FALSE
  )
  print(plot_connected_comps_evolution(nn_conn_comps_object))
  print(plot_n_neigh_k_correspondence(nn_importance_object))
  print(plot_n_neigh_ecs(nn_importance_object))
}
run.nn.stability(listofseurat$NO)
run.nn.stability(listofseurat$HY)
# based on this analysis 30 was selected as appropriate for normoxia and 25 for hypoxia and SNN was more appropriate in both cases

# ClustAssess clustering stability analysis
library(patchwork)
run.clust.stability <- function(so.obj,k.param){
  adj_matrix = FindNeighbors(so.obj@reductions$umap@cell.embeddings, k.param = k.param, nn.method = "rann", verbose = F)$snn
  clustering_diff_obj = get_clustering_difference(graph_adjacency_matrix = adj_matrix,
                                                  resolution = c(seq(from = 0.01, to = 0.1, by = 0.01),seq(from=0.2,to=1,by=0.1)),
                                                  n_repetitions = 30,
                                                  ecs_thresh = 1,
                                                  ncores = 10,
                                                  algorithm = 1:3)
  plot_clustering_difference_boxplot(clustering_diff_obj)
  plot_clustering_difference_facet(clustering_diff_obj, so.obj@reductions$umap@cell.embeddings)
}
run.clust.stability(listofseurat$NO,30)
run.clust.stability(listofseurat$HY,25)
run.clust.stability(listofseurat$WT,30)
# based on this analysis SLM was selected as the most stable clustering approach

run.resolution.stability <- function(so.obj,k.param,clustering.method){
  resolution_gridsearch = get_resolution_importance(embedding = so.obj@reductions$umap@cell.embeddings,
                                                    resolution = c(seq(from = 0.01, to = 0.1, by = 0.01),seq(from=0.2, to=1, by = 0.1)),
                                                    n_neigh = k.param,
                                                    n_repetitions = 100,
                                                    clustering_method = clustering.method,
                                                    graph_type = 2,
                                                    ecs_thresh = 1,
                                                    ncores = 5)
  print(plot_k_resolution_corresp(resolution_gridsearch) +
          plot_annotation(title = "resolution - k correspondence with ecs threshold = 1"))
  print(plot_k_resolution_corresp(resolution_gridsearch, colour_information = "freq_k") +
          plot_annotation(title = "resolution - k correspondence with ecs threshold = 1"))
  print(plot_k_n_partitions(resolution_gridsearch) + 
          plot_annotation(title = "k - # partitions correspondence with ecs threshold = 1"))
  print(plot_k_n_partitions(resolution_gridsearch, colour_information = "frequency_partition") + 
          plot_annotation(title = "k - # partitions correspondence with ecs threshold = 1"))
}
run.resolution.stability(listofseurat$NO,30,3)
run.resolution.stability(listofseurat$HY,25,3)
run.resolution.stability(listofseurat$WT,30,3)
# the selected resolutions were 0.1 for NO, 0.05 for HY, 0.04 for WT
adj_matrix = FindNeighbors(listofseurat$NO@reductions$umap@cell.embeddings, k.param = 30, nn.method = "rann", verbose = F)$snn
clusters=FindClusters(adj_matrix,algorithm=3,resolution=0.1)
listofseurat[['NO']][['cluster_res_0.1']]=clusters
adj_matrix = FindNeighbors(listofseurat$HY@reductions$umap@cell.embeddings, k.param = 25, nn.method = "rann", verbose = F)$snn
clusters=FindClusters(adj_matrix,algorithm=3,resolution=0.05)
listofseurat[['HY']][['cluster_res_0.05']]=clusters
adj_matrix = FindNeighbors(listofseurat$WT@reductions$umap@cell.embeddings, k.param = 30, nn.method = "rann", verbose = F)$snn
clusters=FindClusters(adj_matrix,algorithm=3,resolution=0.05)
listofseurat[['WT']][['cluster_res_0.04']]=clusters

#----------------------------------------------------------------------------------

# cluster markers were identified using adjusted p-value <0.05 and |LFC|>0.25 and |LFC|>1 e.g.:
markers = FindMarkers(listofseurat$NO,group.by=cluster_res0.1,ident.1='0')
markers=markers[markers$p_val_adj<0.05,]
markers.lfc1 = markers[abs(markers$avg_log2FC)>1,]

# cell clusters were then labelled as follows:
listofseurat$NO$cell_type = ifelse(listofseurat$NO$cluster_res0.1=="0",'BC1',
                      ifelse(listofseurat$NO$cluster_res0.1=="1",'MCC2',
                             ifelse(listofseurat$NO$cluster_res0.1=="2",'MCC1',
                                    ifelse(listofseurat$NO$cluster_res0.1=="3",'MCC3',
                                           ifelse(listofseurat$NO$cluster_res0.1=="4",'BC2',
                                                  ifelse(listofseurat$NO$cluster_res0.1=="5",'Cycling BC',
                                                         ifelse(listofseurat$NO$cluster_res0.1=="6",'Deuterosomal',
                                                                'MCC4')))))))

listofseurat$HY$cell_type = ifelse(listofseurat$HY$cluster_res0.05=="0",'BC',
                      ifelse(listofseurat$HY$cluster_res0.05=="1",'MCC1',
                             ifelse(listofseurat$HY$cluster_res0.05=="2",'MCC2',
                                    ifelse(listofseurat$HY$cluster_res0.05=="4","Deuterosomal",  
                                           ifelse(listofseurat$HY$cluster_res0.05=="3",'Hypoxic BC','Cycling hypoxic BC')))))
# subclustering was performed through ClustAssess as above for only basal cells, i.e. BC1/BC2/Cycling BC for normoxia and BC/Hypoxic BC/Cycling hypoxic BC for hypoxia

#---------------------------------------------------------------------------------

# DE genes between WT and KO were identified in normoxia and hypoxia separately, overall and per cell cluster
markers = FindMarkers(listofseurat$NO,group.by = 'WT_KO',ident.1 = 'WT',ident.2 = 'KO')
markers = markers[markers$p_val_adj<0.05,]
markers.lfc1 = markers[abs(markers$avg_log2FC)>1,]
markers = FindMarkers(subset(listofseurat$NO,cell_type=='BC1'),group.by = 'WT_KO',ident.1 = 'WT',ident.2 = 'KO')
markers = markers[markers$p_val_adj<0.05,]
markers.lfc1 = markers[abs(markers$avg_log2FC)>1,]
# The same was performed between NO and HY in the wild type sample

#----------------------------------------------------------------------------------

# Enrichment analysis was performed on these DE genes using gprofiler2
library(magrittr)
library(gprofiler2)
features = read.csv('E3/outs/filtered_feature_bc_matrix/features.tsv',sep='\t',header = F)
selected.features = features[features$V2 %in% markers$gene_name,]
gprof = gprofiler2::gost(query = selected.features$V1,organism = 'hsapiens',
                         custom_bg = features[features$V2%in%rownames(so@assays$SCT),]$V1,
                         sources = c('GO','KEGG','REAC','TF','MIRNA'),
                         evcodes = T,
                         correction_method = 'fdr')
if (!is.null(gprof)){
  gprof$result <- gprof$result %>%
    dplyr::mutate(parents = sapply(.data$parents, toString),
                  intersection_names = sapply(.data$intersection, function(x){
                    ensids <- strsplit(x, split = ",")[[1]]
                    names <- selected.features$V2[match(ensids, selected.features$V1)]
                    paste(names, collapse = ",")
                  }))
}

#--------------------------------------------------------------------------------

# Pseudotime inference was performed using monocle3
library(SeuratWrappers)
library(monocle3)
run.pseudotime <- function(so.obj){
  # on all cells
  so.obj.cds <- as.cell_data_set(so.obj)
  so.obj.cds <- cluster_cells(cds = so.obj.cds, reduction_method = "UMAP")
  so.obj.cds <- learn_graph(so.obj.cds, use_partition = TRUE)
  so.obj.cds <- order_cells(so.obj.cds, 
                        reduction_method = "UMAP",
                        root_cells = names(subset(so.obj,TP63>0 & TOP2A>0 & MKI67>0 & KRT5>0 & FOXJ1==0 & SNTN==0)$orig.ident))
  print(plot_cells(
    cds = so.obj.cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE
  ))
  
  # on only wildtype
  so.obj.cds <- as.cell_data_set(subset(so.obj,WT_KO=='WT'))
  so.obj.cds <- cluster_cells(cds = so.obj.cds, reduction_method = "UMAP")
  so.obj.cds <- learn_graph(so.obj.cds, use_partition = TRUE)
  so.obj.WT = subset(so.obj,WT_KO=='WT')
  so.obj.cds <- order_cells(so.obj.cds, reduction_method = "UMAP",root_cells = names(subset(so.obj.WT,TP63>0 & TOP2A>0 & MKI67>0 & KRT5>0 & FOXJ1==0 & SNTN==0)$orig.ident))
  print(plot_cells(
    cds = so.obj.cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE
  ))

  # on only knockout
  so.obj.cds <- as.cell_data_set(subset(so.obj,WT_KO=='KO'))
  so.obj.cds <- cluster_cells(cds = so.obj.cds, reduction_method = "UMAP")
  so.obj.WT = subset(so.obj,WT_KO=='KO')
  so.obj.cds <- learn_graph(so.obj.cds, use_partition = TRUE)
  so.obj.cds <- order_cells(so.obj.cds, reduction_method = "UMAP",
                            root_cells = names(subset(so.obj.WT,TP63>0 & TOP2A>0 & MKI67>0 & KRT5>0 & FOXJ1==0 & SNTN==0)$orig.ident))
  print(plot_cells(
    cds = so.obj.cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE
  ))
}
run.pseudotime(listofseurat$NO)
run.pseudotime(listofseurat$HY)

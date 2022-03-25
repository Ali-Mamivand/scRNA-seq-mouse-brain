library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(AnnotationHub)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)

  ##################################################################
####-------------1.Load Data and Quality Control------------------####
  ##################################################################
#  Load E10 
E10 = Read10X(data.dir = 'C:/Users/ali/Desktop/git/scrna-seq/E10/')
E10=CreateSeuratObject(E10,project = 'E10')
sce_E10=as.SingleCellExperiment(E10)

#  Load E12 
E12 = Read10X(data.dir = 'C:/Users/ali/Desktop/git/scrna-seq/E12/') 
E12 = CreateSeuratObject(E12,project = 'E12')
sce_E12=as.SingleCellExperiment(E12)
#combine E10 and E12
sce_list=list(sce_E10,sce_E12)
# filter low quality cells
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
filter_cell= function(sce_object){
  chr.loc <- mapIds(ens.mm.v97, keys=rownames(sce_object),
                    keytype="SYMBOL", column="SEQNAME")
  is.mito.alt <- which(chr.loc=="MT")
  is.mito <- ((chr.loc)=="MT")
  qcstats <- perCellQCMetrics(sce_object,subsets= list(Mito=is.mito.alt))
  qc.lib <- isOutlier(qcstats$sum,log = TRUE, type = 'lower')
  qc.nexprs <- isOutlier(qcstats$detected,log = TRUE, type = 'lower')
  qc.mito <- qcstats$subsets_Mito_percent > 10
  discard <- qc.lib | qc.nexprs | qc.mito
  return(sce_object[,!discard])
}
sce_list=lapply(sce_list, filter_cell)
so_list=lapply(sce_list, as.Seurat)
  ##############################################################
####--------------2. data proccessing-------------------------####
  ##############################################################
so_list <- lapply(X = so_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = so_list)
# I prefer to use robust PCA(RPCA) as a more robust approach.
# you can conventional PCA as well by removing reduction argument in below
# code :
so_list <- lapply(X = so_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
so_anchors <- FindIntegrationAnchors(object.list = so_list, 
                                     anchor.features = features,reduction='rpca')
# this command creates an 'integrated' data assay
so.combined <- IntegrateData(anchorset = so_anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(so.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
so.combined <- ScaleData(so.combined, verbose = FALSE)
so.combined <- RunPCA(so.combined, npcs = 30, verbose = FALSE)
so.combined <- RunUMAP(so.combined, reduction = "pca", dims = 1:30)
so.combined <- FindNeighbors(so.combined, reduction = "pca", dims = 1:30)
so.combined <- FindClusters(so.combined, resolution = 0.5)
groups=so.combined$ident
so.combined <- AddMetaData(object = so.combined, metadata = groups,
                           col.name = "developmental_stage")
# Visualization
p1 <- DimPlot(so.combined, reduction = "umap", group.by = "developmental_stage")
p2 <- DimPlot(so.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
# saveRDS(so.combined,file = 'so_combine.RDS')
# so.combined = readRDS('so_combine.RDS')
DimPlot(so.combined, reduction = "umap", split.by = "developmental_stage")

# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(so.combined) <- "RNA"
nk.markers <- FindConservedMarkers(so.combined, ident.1 = 6, grouping.var = "developmental_stage", verbose = FALSE)
head(nk.markers)
goi=c('Sox2','Nes','Prrx2','Lhx2','Lhx9','Atoh1',
      'Gad1','Slc32a1','Ptf1a')
FeaturePlot(so.combined, features = goi,
            min.cutoff = "q9")


# find markers for every cluster compared to all remaining cells, report only the positive
# ones
so.combined.markers <- FindAllMarkers(so.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  ################################################################
####--------------3. Tarajectory Analysis(by monocle3)----------####
  ################################################################
counts=so.combined@assays$RNA@counts
cds <- new_cell_data_set(counts)
rowData(cds)$gene_short_name=rownames(cds)
cds <- preprocess_cds(cds, method = "PCA",num_dim = 30)
cds <- reduce_dimension(cds, preprocess_method = "PCA",
                        reduction_method = "UMAP")
reducedDims(cds)$UMAP[,1]= so.combined@reductions$umap[[,1]]
reducedDims(cds)$UMAP[,2]= so.combined@reductions$umap[[,2]]
# Running the clustering method. This is necessary to the construct the graph
cds <- cluster_cells(cds, reduction_method = "UMAP")
colData(cds)$clusters=so.combined@meta.data$seurat_clusters
# build trajectories in umap by monocle3:
cds = learn_graph(cds,use_partition = T)
# to find root cluster we plot stem cell makers:
(Stem_Cell_Markers=plot_cells(cds, genes= c('Sox2','Sox1','Prom1','Pax6'),
                              show_trajectory_graph=TRUE,
                              label_cell_groups=FALSE,
                              label_leaves=FALSE,
                              trajectory_graph_color = 'red',
                              trajectory_graph_segment_size = 0.5,
                              label_branch_points = FALSE))
cds = order_cells(cds,reduction_method = 'UMAP')
plot_cells(cds,reduction_method = 'UMAP',color_cells_by = 'pseudotime',
          label_cell_groups = F,label_branch_points = F,label_leaves = F)

# Draw UMAP by monocle3 clustering in monocle3:
plot_cells(cds, color_cells_by = "cluster", cell_size = 0.8,
           group_label_size = 5,label_branch_points = FALSE,
           label_roots = FALSE,label_leaves = FALSE)
# Draw UMAP by Seurat Clusters in monocle3:
(umap=plot_cells(cds, label_groups_by_cluster = FALSE, cell_size = 0.6,
                 color_cells_by = "clusters",label_roots = FALSE,
                 label_leaves = FALSE,show_trajectory_graph = FALSE,
                 group_label_size = 10,label_branch_points = FALSE,)+
    theme(text = element_text(size=20)))


# color umap by expression of specific genes:

(GABAergic_Markers=plot_cells(cds, genes= c('Gad1','Slc32a1'),
                              show_trajectory_graph=TRUE,
                              label_cell_groups=FALSE,
                              label_leaves=FALSE,
                              trajectory_graph_color = 'red',
                              trajectory_graph_segment_size = 0.5,
                              label_branch_points = FALSE))
pdf('GABAergic_Markers.pdf',width = 9,height = 5)
GABAergic_Markers
dev.off()
goi=c('Sox2','Nes','Prrx2','Lhx2','Lhx9','Atoh1',
      'Gad1','Slc32a1','Ptf1a')
(goi_plot=plot_cells(cds, genes= goi,
                     show_trajectory_graph=TRUE,
                     label_cell_groups=FALSE,
                     label_leaves=FALSE,
                     trajectory_graph_color = 'red',
                     trajectory_graph_segment_size = 0.3,
                     label_branch_points = FALSE,
                     label_roots = F))
pdf('Markers_umap.pdf',width = 13,height = 11)
goi_plot
dev.off()
# plot marker for each cluster:
plot_genes_by_group(cds,goi,group_cells_by = 'clusters',
                    ordering_type = 'maximal_on_diag')

# We find all the cells that are close to the starting point
cell_ids <- colnames(cds)[colData(cds)$clusters ==  "0"]
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
closest_vertex <- closest_vertex[cell_ids, ]
closest_vertex <- as.numeric(names(which.max(table(closest_vertex))))
mst <- principal_graph(cds)$UMAP
root_pr_nodes <- igraph::V(mst)$name[closest_vertex]

# We compute the trajectory
cds <- order_cells(cds, root_pr_nodes = root_pr_nodes)
plot_cells(cds, color_cells_by = "pseudotime")
goi= c('Sox2','Nes','Gad1','Lhx2')
cds_subset <- cds[rowData(cds)$gene_short_name %in% goi,]
pdf('genes_in_pseudotime.pdf',width = 9,height =9 )
plot_genes_in_pseudotime(cds_subset,  min_expr = NULL,
                         cell_size = 2,
                         nrow = NULL,
                         ncol = 1,
                         panel_order = NULL,
                         color_cells_by = "pseudotime",
                         trend_formula = "~ splines::ns(pseudotime, df=3)")+
  theme(text = element_text(size=20))
dev.off()
#------------------------------------------------
## also, you can perform pseudotime analysis in specific path:
cds_sub <- choose_graph_segments(cds)
cds_sub = preprocess_cds(cds_sub, method = 'PCA')
cds_sub = reduce_dimension(cds_sub,reduction_method = 'UMAP')
cds_sub = cluster_cells(cds_sub)
cds_sub <- learn_graph(cds_sub, use_partition = T)

#cds_sub <- order_cells(cds_sub,reduction_method = 'UMAP')
cell_ids <- colnames(cds_sub)[colData(cds_sub)$clusters ==  "0"]
closest_vertex <- cds_sub@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(cds_sub), ])
closest_vertex <- closest_vertex[cell_ids, ]
closest_vertex <- as.numeric(names(which.max(table(closest_vertex))))
mst <- principal_graph(cds_sub)$UMAP
root_pr_nodes <- igraph::V(mst)$name[closest_vertex]
cds_sub <- order_cells(cds_sub, root_pr_nodes = root_pr_nodes,
                       reduction_method = 'UMAP')
plot_cells(cds_sub, color_cells_by = "pseudotime",label_branch_points = F,
           label_leaves = F,label_roots = F,cell_size=0.8)
goi=c('Nes','Lhx9')
cds_subset <- cds_sub[rowData(cds_sub)$gene_short_name %in% goi,]

plot_genes_in_pseudotime(cds_subset,  min_expr = NULL,
                         cell_size = 2,
                         nrow = NULL,
                         ncol = 1,
                         panel_order = NULL,
                         color_cells_by = "pseudotime",
                         trend_formula = "~ splines::ns(pseudotime, df=3)")+
  theme(text = element_text(size=20))
#-----------------------------------------
# generate cluster markers using monocle3 
set.seed(219)
marker_test_res <- top_markers(cds, group_cells_by="cluster", 
                               reference_cells=1000, cores=6)
write.csv(marker_test_res,file='TopMarker_monocle_cluster.csv')
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=7)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

(markers=plot_genes_by_group(cds,
                             top_specific_marker_ids,
                             group_cells_by="cluster",
                             ordering_type="maximal_on_diag",
                             max.size=25)+
    theme(text = element_text(size=72,face = 'bold'),legend.key.height = unit(4,"line")))

pdf("top3Marker.pdf", width=52, height=45)
markers
dev.off()
##-----------------------------------------
#Finding modules of co-regulated genes
gene_module_df <- find_gene_modules(cds, resolution=1e-2,cores = 6)
gene_module_df$symbol = genes[(genes$Ens.ID %in% gene_module_df$id),2]
write.csv(gene_module_df,file = 'Madules.csv')
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=clusters(cds)[colnames(cds)])
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("clusters ", colnames(agg_mat))
pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=12,fontsize_row = 9)

# As you can see in below, we can color cells in UMAP by expressing genes in specific madules
# (madules 21 and 37 in this case):
(plot_cells_by_madules=plot_cells(cds, 
                                  genes=gene_module_df %>% filter(module %in% c(21,37)),
                                  show_trajectory_graph=FALSE))
pdf('plot_cells_by_madules.PDF',width =11,height = 6)
plot_cells_by_madules
dev.off()
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

sessionInfo()



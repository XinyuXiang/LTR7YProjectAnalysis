######################################################################
###         script for analysis of naive hESC scRNA-seq            ###
######################################################################

##### Setup #####
library(Seurat)
library(dplyr)
library(ggplot2)
library(sctransform)
rm(list = ls())
shhhh=gc() # perform garbage collection to free RAM

##### Input data #####
###gene###
scFINEmerge2_gene.data <- Read10X(data.dir = "./scFINEmerge2_gene/outs/filtered_feature_bc_matrix/")
# count distribution
hist(colSums(as.data.frame(scFINEmerge2_gene.data)), n=50, main = "scFINEmerge2_gene_count_distribution", xlab = "count")

###TE###
scFINEmerge2_TE.data <- Read10X(data.dir = "./1_cellranger/scFINEmerge2_TE/outs/filtered_feature_bc_matrix/")
# count distribution
hist(colSums(as.data.frame(scFINEmerge2_TE.data)), n=50, main = "scFINEmerge2_gene_count_distribution", xlab = "count")

##### Initialize Seurat object #####
scFINEmerge2_gene <- CreateSeuratObject(scFINEmerge2_gene.data, project = "scFINEmerge2_gene", min.cells = 3, min.features = 200)
scFINEmerge2_TE <- CreateSeuratObject(scFINEmerge2_TE.data, project = "scFINEmerge2_TE")


##### Standard pre-processing workflow #####
## QC and selecting cells for further analysis
# "percent.mt": The percentage of reads that map to the mitochondrial genome
# Low-quality / dying cells often exhibit extensive mitochondrial contamination
scFINEmerge2_gene[["percent.mt"]] <- PercentageFeatureSet(scFINEmerge2_gene, pattern = "^MT-")

### Visualize QC metrics as a violin plot
pdf("./VlnPlot_QC_scFINEmerge2_gene_QC.pdf")
VlnPlot(scFINEmerge2_gene, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
VlnPlot(scFINEmerge2_TE, features = c("nFeature_RNA", "nCount_RNA"), group.by = NULL,ncol = 2, pt.size = 0.01)
dev.off()

##### Remove unwanted cells from the dataset #####
scFINEmerge2_gene <- subset(scFINEmerge2_gene, subset = nFeature_RNA > 200 & percent.mt < 30)

##### Normalizing the data #####
# scale factor (10,000 by default), log-transforms the result, Normalized values are stored in [["RNA"]]@data
scFINEmerge2_gene <- NormalizeData(scFINEmerge2_gene, normalization.method = "LogNormalize", scale.factor = 10000)
scFINEmerge2_TE <- NormalizeData(object = scFINEmerge2_TE, normalization.method = "LogNormalize", scale.factor = 10000)

##### Identification of highly variable features #####
# calculate a subset of features that exhibit high cell-to-cell variation in the dataset
scFINEmerge2_gene <- FindVariableFeatures(scFINEmerge2_gene, selection.method = "vst", nfeatures = 2000)
# Identify the 50 most highly variable genes
top50 <- head(VariableFeatures(scFINEmerge2_gene), 50)
# Store the high variation gene in file
write.table(top50,"./HVGs_top50.txt")
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(scFINEmerge2_gene)
plot2 <- LabelPoints(plot = plot1, points = top10)
CombinePlots(plots = list(plot1, plot2))
dev.off()

scFINEmerge2_TE <- FindVariableFeatures(scFINEmerge2_TE, selection.method = "vst",nfeatures = 200)
# Identify the 50 most highly variable TEs
pdf("./dotplot_scFINEmerge2_TE_varTE.pdf", width = 15, height = 7)
top50 <- head(VariableFeatures(scFINEmerge2_TE), 50)
# Store the high variation TE in file
write.table(top50,"/HVTEs_top50.txt")
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(scFINEmerge2_TE)
plot2 <- LabelPoints(plot = plot1, points = top50, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()


##### Scale the data ######
# a standard pre-processing step prior to dimensional reduction techniques like PCA. 
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# The results of this are stored in [["RNA"]]@scale.data
all.genes <- rownames(scFINEmerge2_gene)
scFINEmerge2_gene <- ScaleData(scFINEmerge2_gene, features = all.genes)

all.TEs <- rownames(scFINEmerge2_TE)
scFINEmerge2_TE <- ScaleData(scFINEmerge2_TE, features = all.TEs)

##### Perform linear dimensional reduction #####
# perform PCA on the scaled data. 
# By default, only the previously determined variable features are used as input, but can be defined using features argument if you wish to choose a different subset.
scFINEmerge2_gene <- RunPCA(scFINEmerge2_gene, features = VariableFeatures(object = scFINEmerge2_gene))
VizDimLoadings(scFINEmerge2_gene, dims = 1:2, reduction = "pca")
DimPlot(scFINEmerge2_gene, reduction = "pca")

scFINEmerge2_TE <- RunPCA(scFINEmerge2_TE, features = VariableFeatures(object = scFINEmerge2_TE))
VizDimLoadings(scFINEmerge2_TE, dims = 1:2, reduction = "pca")
DimPlot(scFINEmerge2_TE, reduction = "pca")


##### Determine the ‘dimensionality’ of the dataset #####
ElbowPlot(scFINEmerge2_TE)


##### Cluster the cells #####
# Briefly, these methods embed cells in a graph structure - for example a K-nearest neighbor (KNN) graph, with edges drawn between cells with similar feature expression patterns, and then attempt to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’.
scFINEmerge2_gene <- FindNeighbors(scFINEmerge2_gene, dims = 1:20) # dims = previously defined dimensionality of the dataset (first 20PCs)
scFINEmerge2_gene <- FindClusters(scFINEmerge2_gene, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(scFINEmerge2_gene), 5)

# scFINEmerge2_TE <- FindNeighbors(scFINEmerge2_TE, dims = 1:20) 
# scFINEmerge2_TE <- FindClusters(scFINEmerge2_TE, resolution = 0.5)
# head(Idents(scFINEmerge2_TE), 5)


##### Save the object #####
saveRDS(scFINEmerge2_gene, file = "./scFINEmerge2_gene.rds")
saveRDS(scFINEmerge2_TE, file = "./scFINEmerge2_TE.rds")

##### Finding differentially expressed features #####
# find markers for every cluster compared to all remaining cells, report only the positive ones
scFINEmerge2_gene.markers <- FindAllMarkers(scFINEmerge2_gene, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scFINEmerge2_gene.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC) 
write.table(scFINEmerge2_gene.markers, "./list_cluster_markergenes_top30.txt")

scFINEmerge2_TE.markers <- FindAllMarkers(scFINEmerge2_TE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scFINEmerge2_TE.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.table(scFINEmerge2_TE.markers, "./list_cluster_markerTEs_top30.txt")

##### Set TE counts in gene #####
TE.count.data <- GetAssayData(object = scFINEmerge2_TE[["RNA"]], slot = "counts")
union_matrix2 = matrix(nrow=1180,ncol=ncol(scFINEmerge2_gene@assays$RNA))

for( i in 1:ncol(scFINEmerge2_gene@assays$RNA)){
  ifelse((colnames(scFINEmerge2_gene@assays$RNA)[i] %in% colnames(TE.count.data)),
         (TE.count.data.ingene = TE.count.data[,c(colnames(scFINEmerge2_gene@assays$RNA)[i])]),
         (TE.count.data.ingene = rep(0,1180)))
  union_matrix2[,i] = TE.count.data.ingene
}

colnames(union_matrix2) = colnames(scFINEmerge2_gene@assays$RNA)
row.names(union_matrix2) = row.names(scFINEmerge2_TE[["RNA"]])


scFINEmerge2_gene_TE = scFINEmerge2_gene
scFINEmerge2_gene_TE[["TEcounts"]] <- CreateAssayObject(counts = union_matrix2)
scFINEmerge2_gene_TE <- NormalizeData(scFINEmerge2_gene_TE, assay = "TEcounts", normalization.method = "LogNormalize", 
                                      scale.factor = 1000)
scFINEmerge2_gene_TE <- ScaleData(scFINEmerge2_gene_TE, assay = "TEcounts")
saveRDS(scFINEmerge2_gene_TE, file = "./scFINEmerge2_gene_TE.rds")


##### Feature plot #####
genes = read.table("genelist.txt",header = TRUE, stringsAsFactors = FALSE)
folder_featureplot = "./Featureplot/gene/"
folder_violinplot = "./Vlnplot/gene/"
folder_ridgeplot = "./Ridgeplot/gene/"

### featureplot of genes from genelist
pdf(paste0(folder_featureplot,"Featureplot.pdf"),width = 8,height = 8)
for (i in 1:nrow(genes)){
  if (genes[i,1] %in% (scFINEmerge2_gene@assays[["RNA"]]@counts@Dimnames[[1]])){
    p=FeaturePlot(scFINEmerge2_gene,features = genes[i,1])
    print(p)
  }
}
dev.off()




TEs = read.table("TElist.txt",header = TRUE, stringsAsFactors = FALSE)
folder_featureplot = "./Featureplot/TE/"
folder_violinplot = "./Vlnplot/TE/"
folder_ridgeplot = "./Ridgeplot/TE/"

### featureplot of TEs from TElist
pdf(paste0(folder_featureplot, "Featureplot_scFINEmerge2_TE_in_gene.pdf"),width = 8,height = 8)
for (i in 1:nrow(TEs)){
  if (TEs[i,1] %in% (scFINEmerge2_gene_TE@assays[["RNA"]]@counts@Dimnames[[1]])){
    p=FeaturePlot(scFINEmerge2_gene_TE, features = TEs[i,1])
    print(p)
  }
}
dev.off()



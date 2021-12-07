library(tidyverse)
library(readxl)
library(Seurat)

raw <- Read10X('data/', gene.column = 1)
UMI <- CreateSeuratObject(counts = raw, project = 'Data', min.cells = 1, min.features = 1)
UMI <- UMI[1:1000]

############################################################

#Create Ribo and Mito percent stats

UMI[['MitoPercent']] <- PercentageFeatureSet(UMI, pattern = "^MT-")
UMI[['RiboPercent']] <- PercentageFeatureSet(UMI, pattern = "^Rps-") + PercentageFeatureSet(UMI, pattern = "^Rpl")

UMI@meta.data <- UMI@meta.data %>% 
  rename(nCounts = nCount_RNA) %>% 
  rename(nGenes = nFeature_RNA)


#Graphs of counts content

jpeg("counts~genes.jpeg" , units="in", width=15, height=10, res=600)
UC_plot <- VlnPlot(UMI, features = c("nGenes", "nCounts"), ncol = 2)
UC_plot
dev.off()

jpeg("Ribo~Mito.jpeg" , units="in", width=15, height=10, res=600)
MR_plot <- VlnPlot(UMI, features = c("RiboPercent", "MitoPercent"), ncol = 2)
MR_plot
dev.off()


jpeg("counts~genes_QC.jpeg" , units="in", width=15, height=10, res=600)
CG_plot <- FeatureScatter(UMI, feature1 = "nCounts", feature2 = "nGenes")
CG_plot
dev.off()

####################################################################################

n_gen <- max(as.numeric(UMI@meta.data$nGenes))*0.75
cells_number <- length(Idents(UMI))


#####################################################################################

UMI <- NormalizeData(UMI, normalization.method = "LogNormalize", scale.factor = 1e6)


#######################################################################################

UMI <- FindVariableFeatures(UMI, selection.method = "vst", nfeatures = n_gen, binning.method = 'equal_frequency')

# Identify the 10 most highly variable genes

top20 <- head(VariableFeatures(UMI), 20)

plot1 <- VariableFeaturePlot(UMI)

jpeg("variable_genes.jpeg", units="in", width=10, height=7, res=600)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot2
dev.off()

#####################################################################################

all.genes <- rownames(UMI)
UMI <- ScaleData(UMI, features = all.genes)
UMI <- RunPCA(UMI, features = VariableFeatures(object = UMI))



################################

Elbow <- ElbowPlot(UMI, ndims = 50)

dims <- as.data.frame(Elbow$data$stdev)

#select the most variable reduction

{
  dim <- 1
  score <- c()
  element <- 0
  for (i in dims$`Elbow$data$stdev`) {
    element <- element + 1
    if (i-i*0.01 > dims$`Elbow$data$stdev`[element+1] & element < 50 | i-i*0.02 > dims$`Elbow$data$stdev`[element+2] & element < 49 | i-i*0.02 > dims$`Elbow$data$stdev`[element+3] & element < 48 | i-i*0.02 > dims$`Elbow$data$stdev`[element+4] & element < 47) {
      dim <- dim + 1
    } else 
      break
  }
  dim <- as.numeric(dim)
}


jpeg("Elbow.jpeg" , units="in", width=10, height=7, res=600)
Elbow + geom_vline(xintercept = dim, color = 'red')
dev.off()

#########################################################################################

UMI <- JackStraw(UMI, num.replicate = 10, dims = dim)
UMI <- ScoreJackStraw(UMI, dims = 1:dim)

#Select significient PCs
jc <- as.data.frame(UMI@reductions$pca@jackstraw@overall.p.values)
jc <- jc[jc$Score < 0.05,]
dim <- as.vector(jc$PC)

jpeg("JackStrawPlot.jpeg", units="in", width=10, height=7, res=600)
JackStrawPlot(UMI, dims = dim)
dev.off()

UMI <- FindNeighbors(UMI, dims = dim, reduction = 'pca')
UMI <- FindClusters(UMI, resolution = 0.5, n.start = 10, n.iter = 1000)


UMI <- RunUMAP(UMI, dims = dim, n.neighbors = 29, umap.method = "umap-learn")


width <- 15 + (length(unique(Idents(UMI))))/7

jpeg("UMAP.jpeg", units="in", width = width, height = 15, res=600)
DimPlot(UMI, reduction = "umap", raster = FALSE)
dev.off()

####################################################################################



#find markers for every cluster compared to all remaining cells, report only the positive ones

print('Searching for cluster marker genes')


UMI.markers <- FindAllMarkers(UMI, only.pos = TRUE, min.pct = 0.20, logfc.threshold = 0.25, test.use = 'MAST')

MAST_markers <- UMI.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.table(MAST_markers, file = "MAST_markers_clusters.csv", sep = ',')


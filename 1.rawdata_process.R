# 1.Loading packages

library(Seurat)
library(RCurl)
library(cowplot)
library(dplyr)
library(hdf5r) # install.packages('hdf5r')
options(future.globals.maxSize = 10000 * 1024^2)  # set allowed size to 2K MiB

# 2.Loadng ST-seq datasets

Control_1 <- Load10X_Spatial(
    data.dir = "1.raw_data/Control_1/outs/",
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "Control_1",
    filter.matrix = TRUE,
    to.upper = FALSE)

Control_2 <- Load10X_Spatial(
    data.dir = "1.raw_data/Control_2/outs/",
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "Control_2",
    filter.matrix = TRUE,
    to.upper = FALSE)

Control_3 <- Load10X_Spatial(
    data.dir = "1.raw_data/Control_3/outs/",
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "Control_3",
    filter.matrix = TRUE,
    to.upper = FALSE)

Control_4 <- Load10X_Spatial(
    data.dir = "1.raw_data/Control_4/outs/",
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "Control_4",
    filter.matrix = TRUE,
    to.upper = FALSE)

Control_5 <- Load10X_Spatial(
    data.dir = "1.raw_data/Control_5/outs/",
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "Control_5",
    filter.matrix = TRUE,
    to.upper = FALSE) # Mouse_Kidney_FFPE

Control_6 <- Load10X_Spatial(
    data.dir = "1.raw_data/Control_6/outs/",
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "Control_6",
    filter.matrix = TRUE,
    to.upper = FALSE)

# 3. Add metadata information

Control_1$sample <- "Control_1"
Control_2$sample <- "Control_2"
Control_3$sample <- "Control_3"
Control_4$sample <- "Control_4"
Control_5$sample <- "Control_5"
Control_6$sample <- "Control_6"

Control_1$type <- "Spatial"
Control_2$type <- "Spatial"
Control_3$type <- "Spatial"
Control_4$type <- "Spatial"
Control_5$type <- "FFPE"
Control_6$type <- "Spatial"

Control_1$source <- "This_dataset"
Control_2$source <- "This_dataset"
Control_3$source <- "This_dataset"
Control_4$source <- "This_dataset"
Control_5$source <- "10x Genomics"
Control_6$source <- "10x Genomics"

Control_1$slide <- "Sagittal"
Control_2$slide <- "Sagittal"
Control_3$slide <- "Sagittal"
Control_4$slide <- "Sagittal"
Control_5$slide <- "Sagittal"
Control_6$slide <- "Coronal"

# 4. Add percent_mito information

Control_1 <- PercentageFeatureSet(Control_1, "^mt-", col.name = "percent_mito")
Control_2 <- PercentageFeatureSet(Control_2, "^mt-", col.name = "percent_mito")
Control_3 <- PercentageFeatureSet(Control_3, "^mt-", col.name = "percent_mito")
Control_4 <- PercentageFeatureSet(Control_4, "^mt-", col.name = "percent_mito")
Control_5 <- PercentageFeatureSet(Control_5, "^mt-", col.name = "percent_mito")
Control_6 <- PercentageFeatureSet(Control_6, "^mt-", col.name = "percent_mito")

# 5. Merge datasets

merge_data <- merge(Control_1,c(Control_2,Control_3,Control_4,Control_5,Control_6))
merge_data
# An object of class Seurat
# 32285 features across 16993 samples within 1 assay 

table(merge_data$sample)
# Control_3     Control_4     Control_5     Control_6 Control_1 Control_2 
#   3213          2998          2991          3490       3259      2768

merge_data$sample = factor(merge_data$sample, levels = c("Control_1", "Control_2","Control_3", "Control_4","Control_5", "Control_6"))
merge_data$type = factor(merge_data$type, levels = c("Spatial", "FFPE"))
merge_data$source = factor(merge_data$source, levels = c("This_dataset", "10x Genomics"))
merge_data$slide = factor(merge_data$slide, levels = c("Sagittal","Coronal"))

merge_data2 <- merge(Control_1,c(Control_2,Control_3,Control_4))
merge_data2
saveRDS(merge_data2,"merge_data(4 samples).rds")

# 6. Datasets QC

pdf("merge_data_sample.pdf",height = 4,width = 10)
VlnPlot(merge_data, features = c("nFeature_Spatial","nCount_Spatial","percent_mito"), pt.size = 0,group.by = "sample",ncol = 3)
dev.off()

pdf("merge_data_sample2.pdf",height = 4,width = 10)
VlnPlot(merge_data2, features = c("nFeature_Spatial","nCount_Spatial","percent_mito"), pt.size = 0,group.by = "sample",ncol = 3)
dev.off()

p1 <- (SpatialFeaturePlot(merge_data,features = c("nFeature_Spatial"),alpha = 1,images = c("Control_1")) + color2) /
  (SpatialFeaturePlot(merge_data,features = c("nCount_Spatial"),alpha = 1,images = c("Control_1")) + color2) /
  (SpatialFeaturePlot(merge_data,features = c("percent_mito"),alpha = 1,images = c("Control_1")) + color2)

# p2 <- (SpatialFeaturePlot(ST,features = c("nFeature_Spatial"),alpha = 1,images = c("Control_5.4")) + color2) /
#   (SpatialFeaturePlot(merge_data,features = c("nCount_Spatial"),alpha = 1,images = c("Control_5.4")) + color2) /
#   (SpatialFeaturePlot(ST,features = c("percent_mito"),alpha = 1,images = c("Control_5.4")) + color2)
# 
# p3 <- (SpatialFeaturePlot(merge_data,features = c("nFeature_Spatial"),alpha = 1,images = c("Control_6.5")) + color2) /
#   (SpatialFeaturePlot(merge_data,features = c("nCount_Spatial"),alpha = 1,images = c("Control_6.5")) + color2) /
#   (SpatialFeaturePlot(merge_data,features = c("percent_mito"),alpha = 1,images = c("Control_6.5")) + color2)

p2 <- (SpatialFeaturePlot(ST2,features = c("nFeature_Spatial"),alpha = 1,images = c("Control_2.1")) + color2) /
  (SpatialFeaturePlot(ST2,features = c("nCount_Spatial"),alpha = 1,images = c("Control_2.1")) + color2) /
  (SpatialFeaturePlot(ST2,features = c("percent_mito"),alpha = 1,images = c("Control_2.1")) + color2)

p3 <- (SpatialFeaturePlot(ST2,features = c("nFeature_Spatial"),alpha = 1,images = c("Control_3.2")) + color2) /
  (SpatialFeaturePlot(ST2,features = c("nCount_Spatial"),alpha = 1,images = c("Control_3.2")) + color2) /
  (SpatialFeaturePlot(ST2,features = c("percent_mito"),alpha = 1,images = c("Control_3.2")) + color2)

p4 <- (SpatialFeaturePlot(ST2,features = c("nFeature_Spatial"),alpha = 1,images = c("Control_4.3")) + color2) /
  (SpatialFeaturePlot(ST2,features = c("nCount_Spatial"),alpha = 1,images = c("Control_4.3")) + color2) /
  (SpatialFeaturePlot(ST2,features = c("percent_mito"),alpha = 1,images = c("Control_4.3")) + color2)

pdf("QC-spatialFeature.pdf",height = 20,width = 15)
p1 | p2 | p3 | p4  
dev.off()

p1 | p2 | p3 # 10*10

# 7.Integrate Datasets 

split_data <- SplitObject(merge_data, split.by = "sample")
st.list = lapply(split_data, SCTransform, assay = "Spatial", method = "poisson")
st.features = SelectIntegrationFeatures(st.list, nfeatures = 3000, verbose = FALSE)
st.list <- PrepSCTIntegration(object.list = st.list, anchor.features = st.features, verbose = FALSE)
int.anchors <- FindIntegrationAnchors(object.list = st.list, normalization.method = "SCT", verbose = FALSE, anchor.features = st.features)
integrated_data <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)

# 8.Datasets reduction 
integrated_data <- RunPCA(object = integrated_data,verbose = FALSE)

# pdf("PCA_Elbowplot.pdf",height=10,width=10)
ElbowPlot(integrated_data,ndims = 50)
# dev.off()

integrated_data <- FindNeighbors(integrated_data,dim=1:40)
integrated_data <- FindClusters(integrated_data,resolution = 0.8)
integrated_data <- RunUMAP (integrated_data,reduction="pca", dims = 1:40)
integrated_data <- RunTSNE(integrated_data,dims = 1:40)

# 9.Dataset visulization

pdf("umap-cluster.pdf",height = 10,width = 10)
DimPlot(integrated_data,reduction = "umap" ,label = T,pt.size = 1) 
dev.off()

pdf("umap-sample.pdf",height = 10,width = 10)
DimPlot(integrated_data,reduction = "umap" ,group.by="sample",label = F,pt.size = 1)
dev.off()

pdf("umap-type.pdf",height = 10,width = 10)
DimPlot(integrated_data,reduction = "umap" ,group.by="type",label = T,pt.size = 1)
dev.off()

# pdf("tsne-cluster.pdf",height = 10,width = 10)
# DimPlot(integrated_data,reduction = "tsne" ,label = F,pt.size = 1)
# dev.off()
# 
# pdf("tsne-sample.pdf",height = 10,width = 10)
# DimPlot(integrated_data,reduction = "tsne" ,group.by="sample",label = F,pt.size = 1)
# dev.off()
# 
# pdf("tsne-type.pdf",height = 10,width = 10)
# DimPlot(integrated_data,reduction = "tsne" ,group.by="type",label = T,pt.size = 1)
# dev.off()

# 10.save dataset

saveRDS(integrated_data,"integrated_data.rds")

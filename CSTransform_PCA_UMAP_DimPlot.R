#CSTransform
#Violing plot stratification discussion https://github.com/satijalab/seurat/issues/1623

combined_merged <- readRDS("combined_merged_seurat_b4H.rds") #this is an object of merged combined controls and separately combined experimental variables 
combined_merged <- SCTransform(combined_merged, assay = "RNA", verbose = FALSE)#NB: time consuming 
#combined_merged <- FindVariableFeatures(combined_merged, selection.method = "vst", nfeatures = 2000) #SCT assay is comprised of multiple SCT models. To change the variable features, please set manually with VariableFeatures<-
#combined_merged <- ScaleData(combined_merged) for NormalizeData, not for SCTransform

#PCA (principal components analysis)
combined_merged <- RunPCA(combined_merged, features = VariableFeatures(combined_merged, assay = "SCT"))#somewhat time consuming

ElbowPlot(object = combined_merged, 
          ndims = 40)
#UMAP (Uniform Manifold Approximation and Projection for Dimension Reduction)
combined_merged <- RunUMAP(combined_merged, dims = 1:30) # 30 based on Elbow plot of SD vs PC, somewhat time consuming

#DimPlot(combined_merged, reduction = "umap", group.by = "Condition", label = FALSE)
DimPlot(combined_merged,
        reduction = "umap",
        group.by = "Condition",
        label = TRUE,
        cols = c("brown", "orange"),
        alpha = 0.7,
        pt.size = 0.6,
        label.size = 6.0,
        repel = TRUE)

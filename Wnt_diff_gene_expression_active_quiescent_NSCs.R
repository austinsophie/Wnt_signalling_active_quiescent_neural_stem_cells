#Required packages
library(Seurat)
library(tidyverse)

load("nsc.RData")
nsc.integrated$binary <- ifelse(nsc.integrated$quiescence == "active", yes = "active", no = "quiescent")
Idents(nsc.integrated) <- nsc.integrated$binary
DimPlot(nsc.integrated)
nsc.integrated <- SCTransform(nsc.integrated, return.only.var.genes = FALSE)

#this will test all genes as no thresholds
sig_genes_all <- FindMarkers(nsc.integrated, assay = "SCT", slot = "scale.data", ident.1 = "active",
                             ident.2 = "quiescent", min.pct = 0.1, min.diff.pct = 0.0, logfc.threshold = 0,
                             test.use = "t", min.cells.group = 0, min.cells.feature = 0) 

#change column 0 to column 1 with title genes                             
sig_genes_all<-tibble::rownames_to_column(sig_genes_all, "genes") 

#export FindMarkers results
write_csv(sig_genes_all, path = "/Users/austins/Desktop/sig_genes_all.csv", col_names = TRUE, append = FALSE)


#List of Wnt genes
wnt_genes2 <- c("Wnt7a", "Wnt7b", "Dkk3", "Sfrp1", "Fzd1", "Fzd2", "Fzd3", "Fzd8", "Fzd9", "Lrp5", "Lrp6","Tcf7l2","Dvl3", 
                "Apc", "Ctnnb1", "Gsk3b", "Axin2" )

#Heatmap of Wnt genes
DefaultAssay(nsc.integrated) <- "RNA"
nsc.integrated <- NormalizeData(nsc.integrated, assay = "RNA")
DoHeatmap(nsc.integrated, features = wnt_genes2, assay = "RNA", slot = "data") + scale_fill_gradientn(colors = c("blue", "white", "red"))




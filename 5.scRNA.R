rm(list=ls())
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(data.table)
  library(SingleR)
  library(celldex)
  library(harmony)
  library(clustree)
  library(ggplot2)
  library(openxlsx)
})
setwd("E:/sceQTL14_analysis/scRNA")
############################
## Read GSE181279 data
############################
sample_dirs <- list.dirs("GSE181279",recursive=FALSE)
sample_dirs
counts <- Read10X(data.dir=sample_dirs)
sc <- CreateSeuratObject(counts,project="GSE181279",min.cells=3,min.features=200)
############################
## Add metadata
############################
sample_names <- colnames(sc)
sc$Sample <- sub("^(.*?)_.*","\\1",sample_names)
sc$Group <- ifelse(grepl("GSM5494107|GSM5494110|GSM5494113",sc$Sample),"AD","Control")
############################
## QC
############################
sc[["percent.mt"]] <- PercentageFeatureSet(sc,pattern="^MT-")
sc[["percent.rb"]] <- PercentageFeatureSet(sc,pattern="^RP")
sc <- subset(sc, subset= nFeature_RNA>200 & nFeature_RNA<3000 & nCount_RNA>1000 & nCount_RNA<8000 & percent.mt<15)
saveRDS(sc,"GSE181279_QC.rds")
############################
## Normalization and HVG
############################
sc <- NormalizeData(sc,normalization.method="LogNormalize",scale.factor=10000)
sc <- FindVariableFeatures(sc,selection.method="vst",nfeatures=2000)
############################
## PCA
############################
sc <- ScaleData(sc)
sc <- RunPCA(sc,features=VariableFeatures(sc),npcs=40)
sc <- IntegrateLayers(
  object=sc,
  method=RPCAIntegration,
  orig.reduction="pca",
  new.reduction="integrated.rpca",
  verbose=FALSE)
############################
## Clustering
############################
pcSelect <- 15
sc <- FindNeighbors(sc,reduction="integrated.rpca",dims=1:pcSelect)
sc <- FindClusters(sc,resolution=1)
sc <- RunTSNE(sc,reduction="integrated.rpca",dims=1:pcSelect)
ggsave("TSNE_cluster.pdf",DimPlot(sc,label=TRUE),width=7,height=6)
############################
## Cell annotation
############################
ref <- HumanPrimaryCellAtlasData()
test <- GetAssayData(sc,layer="data")
pred <- SingleR(test=test,ref=ref,labels=ref$label.main)
sc$SingleR <- pred$labels
write.csv(data.frame(cell=colnames(sc),SingleR=sc$labels),"SingleR_annotation.csv",row.names=FALSE)
############################
## Marker based annotation
############################
markers <- FindAllMarkers(sc,only.pos=TRUE,min.pct=0.25,logfc.threshold=0.25)
write.csv(markers,"cluster_markers.csv",row.names=FALSE)

egene_file <- "coloc_PP_H4_gt_0.8.csv"
egene <- fread(egene_file)
genes <- unique(egene$Gene)
genes_exist <- intersect(genes,rownames(sc))
FeaturePlot(sc,features=genes_exist,reduction="umap")
ggsave("eGene_FeaturePlot.pdf",width=10,height=12)
p <- DotPlot(sc,
  features=genes_exist,
  group.by="SingleR")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1))
ggsave("eGene_DotPlot.pdf",p,width=12,height=6)
############################
## Differential expression
############################
Idents(sc) <- sc$Group
deg_results <- list()
celltypes <- unique(sc$SingleR)
for(ct in celltypes){cells <- WhichCells(sc,expression=SingleR==ct)
  sub <- subset(sc,cells=cells)
  if(length(unique(sub$Group))<2){next}
deg <- FindMarkers(
    sub,
    ident.1="AD",
    ident.2="Control",
    features=genes_exist,
    logfc.threshold=0,
    min.pct=0)
if(nrow(deg)>0){
    deg$Gene <- rownames(deg)
    deg$CellType <- ctdeg_results[[ct]] <- deg}}
deg_all <- bind_rows(deg_results)
fwrite(deg_all,"eGene_AD_vs_Control_scRNA.csv")
############################
##  Save object
############################
saveRDS(sc,"GSE181279_processed_Seurat.rds")
cat("Cells:",ncol(sc),"\n")
cat("Genes:",nrow(sc),"\n")


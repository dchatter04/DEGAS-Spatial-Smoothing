##############################################################
rm(list=ls())
library(Seurat)
library(hdf5r)
library(DEGAS)
library(dplyr)
library(ggplot2)
library(progressr)
library(progress)
library(ggrepel)
source("~/DEGAS_functions_for_smoothing.R")

################################# DEGAS model training##############

#Options
save_input = TRUE
save_output = TRUE
retrain = TRUE

# Loading bulk T2D data
Xpat = read.csv("~/bulkRNA_seq_Diabetes.csv",row.names=1)
Ypat = read.csv("~/patLab.csv",row.names=1)
Ypat = Ypat[colnames(Xpat),]


# Loading Xenium data
#######************ The Xenium T2D data is available in https://www.synapse.org/Synapse:syn68699752/files/ ***********
S0017155_ND = LoadXenium("~/output-XETG00126__0017155__ND16091_ND__20231010__204554") #this is ND_1 sample
S0017155_T2D = LoadXenium("~/output-XETG00126__0017155__ND16154_T2D__20231010__204554") $this is T2D_1 sample
S0017308_T2D = LoadXenium("~/output-XETG00126__0017308__ND16090_T2D__20231010__204554") #this is T2D_2 sample
S0017308_ND = LoadXenium("~/output-XETG00126__0017308__ND16091_ND__20231010__204554") # this is ND_2 sample
# Formatting Xenium data
nd1 = S0017155_ND@assays$Xenium$counts
rownames(nd1) = rownames(S0017155_ND)
colnames(nd1) = sub("-1","-ND1",colnames(S0017155_ND))
t2d1 = S0017155_T2D@assays$Xenium$counts
rownames(t2d1) = rownames(S0017155_T2D)
colnames(t2d1) = sub("-1","-T2D1",colnames(S0017155_T2D))
t2d2 = S0017308_T2D@assays$Xenium$counts
rownames(t2d2) = rownames(S0017308_T2D)
colnames(t2d2) = sub("-1","-T2D2",colnames(S0017308_T2D))
nd2 = S0017308_ND@assays$Xenium$counts
rownames(nd2) = rownames(S0017308_ND)
colnames(nd2) = sub("-1","-ND2",colnames(S0017308_ND))
Xst = cbind(nd1,t2d1,nd2,t2d2)

# Defining final matrices for training
feats1 = intersect(rownames(Xpat),rownames(Xst))
stDat = Xst[feats1,]
stDat = preprocessCounts(as.matrix(stDat))
stDat = as.matrix(stDat)
patDat = Xpat[feats1,]
patDat = preprocessCounts(as.matrix(patDat))
patLab = as.matrix(Ypat[rownames(patDat),])
rm(S0017155_ND,S0017155_T2D,S0017308_ND,S0017308_T2D,nd1,nd2,t2d1,t2d2,Xst)
if(save_input){
  saveRDS(stDat,file="~/Desktop/Xenium/stDat.rds")
}

## Training DEGAS model #########------------------------

model1 = runDEGASatlas(stDat,NULL,patDat,patLab,"~/Desktop/Xenium/tmp/","BlankClass","DenseNet",3,5,20000,123)
saveRDS(model1,file="~/Desktop/Xenium/model1.rds") #saving the DEGAS model
preds = predClassBag(model1,stDat,"pat")
saveRDS(preds,file="~/Desktop/Xenium/preds1.rds") #saving the predicted values


group = sub(".*[-]","",rownames(stDat))
mapper = c("ND","ND","T2D","T2D")
names(mapper) = c("ND1","ND2","T2D1","T2D2")
group2 = mapper[group]

# Loading UMAP coordinates
umap_nd1 = read.csv("~/output-XETG00126__0017155__ND16091_ND__20231010__204554/analysis/umap/gene_expression_2_components/projection.csv")
umap_nd1$Barcode = sub("-1","-ND1",umap_nd1$Barcode)
umap_t2d1 = read.csv("~/output-XETG00126__0017155__ND16154_T2D__20231010__204554/analysis/umap/gene_expression_2_components/projection.csv")
umap_t2d1$Barcode = sub("-1","-T2D1",umap_t2d1$Barcode)
umap_nd2 = read.csv("~/output-XETG00126__0017308__ND16090_T2D__20231010__204554/analysis/umap/gene_expression_2_components/projection.csv")
umap_nd2$Barcode = sub("-1","-T2D2",umap_nd2$Barcode)
umap_t2d2 = read.csv("~/output-XETG00126__0017308__ND16091_ND__20231010__204554/analysis/umap/gene_expression_2_components/projection.csv")
umap_t2d2$Barcode = sub("-1","-ND2",umap_t2d2$Barcode)
umap_coords = rbind(umap_nd1,umap_t2d1,umap_nd2,umap_t2d2)
rownames(umap_coords) = umap_coords$Barcode
umap_coords = umap_coords[rownames(stDat),]

# Loading cellular locations
locs_nd1 = read.csv("~/output-XETG00126__0017155__ND16091_ND__20231010__204554/cells.csv.gz")
locs_nd1$cell_id = sub("-1","-ND1",locs_nd1$cell_id)
locs_t2d1 = read.csv("~/output-XETG00126__0017155__ND16154_T2D__20231010__204554/cells.csv.gz")
locs_t2d1$cell_id = sub("-1","-T2D1",locs_t2d1$cell_id)
locs_nd2 = read.csv("~/output-XETG00126__0017308__ND16090_T2D__20231010__204554/cells.csv.gz")
locs_nd2$cell_id = sub("-1","-T2D2",locs_nd2$cell_id)
locs_t2d2 = read.csv("~/output-XETG00126__0017308__ND16091_ND__20231010__204554/cells.csv.gz")
locs_t2d2$cell_id = sub("-1","-ND2",locs_t2d2$cell_id)
locs_coords = rbind(locs_nd1,locs_t2d1,locs_nd2,locs_t2d2)
rownames(locs_coords) = locs_coords$cell_id
locs_coords = locs_coords[rownames(stDat),]

# Creating an overall dataframe for smoothing
df=data.frame(barcodes=umap_coords$Barcode,t2d_risk=preds[,1],umap1=umap_coords$UMAP.1,umap2=umap_coords$UMAP.2,locsx=locs_coords$x_centroid,locsy=locs_coords$y_centroid)
df$sample = sub(".*[-]","",df$barcodes)
rm(umap_nd1,umap_t2d1,umap_nd2,umap_t2d2,locs_nd1,locs_t2d1,locs_nd2,locs_t2d2)
gc()

df = na.omit(df)
saveRDS(df,file="~/dfv2.rds")


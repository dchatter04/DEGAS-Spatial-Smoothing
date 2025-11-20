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

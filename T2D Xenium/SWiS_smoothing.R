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


# Performing smoothing using sliding window
df = readRDS("~/Desktop/Xenium/dfv2.rds")
# Setting options
dat_ND1 = df[df$sample=="ND_1",]
sm_dat_ND1 = swSmoothAtlas(dat_ND1,1000,1000,500,100,50,5,"locsx","locsy","t2d_risk")
saveRDS(sm_dat_ND1,"~/swSmoothed_ND1.rds")

dat_ND2 = df[df$sample=="ND_2",]
sm_dat_ND2 = swSmoothAtlas(dat_ND2,1000,1000,500,100,50,5,"locsx","locsy","t2d_risk")
saveRDS(sm_dat_ND2,"~/swSmoothed_ND2.rds")

dat_T2D1 = df[df$sample=="T2D_1",]
sm_dat_T2D1 = swSmoothAtlas(dat_T2D1,1000,1000,500,100,50,5,"locsx","locsy","t2d_risk")
saveRDS(sm_dat_T2D1,"~/swSmoothed_T2D1.rds")

dat_T2D2 = df[df$sample=="T2D_2",]
sm_dat_T2D2 = swSmoothAtlas(dat_T2D2,1000,1000,500,100,50,5,"locsx","locsy","t2d_risk")
saveRDS(sm_dat_T2D2,"~/swSmoothed_T2D2.rds")


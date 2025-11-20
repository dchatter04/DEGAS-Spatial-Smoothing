rm(list = ls(all.names = TRUE))
library(devtools)
library(withr)
library(rprojroot)
library(tidyverse)
library(magrittr)
library(here)
library(pROC)
# If DEGAS hasn't been installed, do so using this command (after creating a r-packages file in your root directory)
#with_libpaths(new = '~/r-packages', install_github("tsteelejohnson91/DEGAS"))
library(DEGAS)
library(dplyr)
library(ggplot2)
library(ggpubr)


#required functions
################################################################################
# Differential expression

# This function negates in() and is a dependency of selectFeats()
`%notin%` <- Negate(`%in%`)
# This function uses a selection statistic to pick the most representative gene
# with duplicate gene ids.
# Note: This is most relevant for multiple probes to one gene or for converting
# between Ensembl gene IDs and HGNC gene symbols
selectFeats <- function(expression,features,selection_statistic)
{
  original_feats = row.names(expression)
  row.names(expression) = 1:dim(expression)[1]
  dup_feat_uniq = unique(features[duplicated(features)])
  message(length(dup_feat_uniq))
  if(length(dup_feat_uniq)>0){
    dup_feat = features[features %in% dup_feat_uniq]
    dup_feat_idx = which(features %in% dup_feat_uniq)
    rem_feat_idx = c()
    for(feat in dup_feat){
      feat_rowSums = apply(expression[dup_feat_idx[dup_feat==feat],],1,eval(parse(text=selection_statistic)))
      max_feat_idx = which(feat_rowSums==max(feat_rowSums))[1]
      rem_feat_idx = c(rem_feat_idx,as.numeric(names(feat_rowSums)[-max_feat_idx]))
    }
    expression = expression[-rem_feat_idx,]
    row.names(expression) = features[-rem_feat_idx]
    feat_df <- data.frame(new_feature=row.names(expression),original_feature=original_feats[-rem_feat_idx])
    row.names(feat_df) <- feat_df$new_feature
    return(list(expression,feat_df))
  }else{
    row.names(expression) = features
    feat_df <- data.frame(new_feature=row.names(expression),original_feature=original_feats[-rem_feat_idx])
    row.names(feat_df) <- feat_df$new_feature
    return(list(expression,feat_df))
  }
}

# Return non-parametric p-value
getPval = function(x,g){
  if(length(unique(g[!is.na(x)]))==2){
    tmp = wilcox.test(x~(g))
    return(tmp$p.value)
  }else{
    return(NA)
  }
}

# Return DGE for a single gene
getSparsityMeansLog2fcPval = function(x,g){
  sparsGpos = sum(x[g]>0,na.rm = TRUE)/sum(g[!is.na(x[g])])
  sparsGneg = sum(x[!g]>0,na.rm = TRUE)/sum(!g[!is.na(x[!g])])
  mnGpos = mean(x[g],na.rm=TRUE)
  mnGneg = mean(x[!g],na.rm=TRUE)
  log2FC = log2((mnGpos+1e-6)/(mnGneg+1e-6))
  out = data.frame(percG1 = sparsGpos, percG2 = sparsGneg,
                   meanG1 = mnGpos, meanG2 = mnGneg,
                   log2FC = log2FC, Pval = getPval(x,g))
  return(out)
}

# DGE one vs all for groups in g
FindClusterMarkers = function(X,g){
  MarkerTable = list()
  for(grp in unique(g)){
    #Pvals = apply(X,2,function(x) getPval(x,g==grp))
    #Log2FCs = apply(X,2,function(x) getLog2FC(x,g==grp))
    tmp = apply(X,2,function(x) getSparsityMeansLog2fcPval(x,g==grp))
    tmp = do.call(rbind,tmp)
    MarkerTable[[grp]] = data.frame(Gene = rownames(tmp),Cluster=rep(grp,dim(X)[2]),tmp)
  }
  MarkerTable = do.call(rbind,MarkerTable)
  return(MarkerTable)
}

JI <- function(list1,list2){
  return(length(intersect(list1,list2))/length(union(list1,list2)))
}

################################################################################
# Distance calculations and smoothing
# Return euclidean distance between two points
euclDist <- function(loc1,loc2){
  return(sqrt(sum((loc1 - loc2)^2)))
}

# Return a matrix of all pairwise distances
pairDist <- function(locs){
  N = dim(locs)[1]
  out = matrix(NA,N,N)
  for(i in 1:N){
    for(j in 1:i){
      out[i,j] = out[j,i] = euclDist(locs[i,],locs[j,])
    }
  }
  return(out)
}

knnSmoothAtlas <- function(sc_data,col_names,k,n_split){
  if(k < floor(0.5*n_split)){
    folds = splitKfoldCV(dim(sc_data)[1],floor(dim(sc_data)[1]/n_split))
    cor_list = list()
    i=0
    for (f in folds){
      i=i+1
      predsSmoothed = knnSmooth(as.numeric(sc_data[f,col_names[2]]),as.matrix(sc_data[f,c(col_names[3],col_names[4])]),k)
      #RespCor = toCorrCoeff(predsSmoothed)
      RespCor = predsSmoothed
      df = data.frame(cell_id = sc_data[f,col_names[1]], coord1 = sc_data[f,col_names[3]], coord2 = sc_data[f,col_names[4]], smoothed = RespCor)
      cor_list[[i]] = df
    }
    cor_list_df = do.call(rbind,cor_list)
    return(cor_list_df)
  }else{
    message(paste0("The number of neighbors (",k,") is greater than 50% of the sample (",n_split,")."))
  }
}


# Required Functions
knnSmoothAtlas2 <- function(sc_data,col_names,k,n_split){
  if(k < floor(0.5*n_split)){
    folds = splitKfoldCV(dim(sc_data)[1],floor(dim(sc_data)[1]/n_split))
    cor_list = list()
    i=0
    for (f in folds){
      i=i+1
      predsSmoothed = knnSmooth(as.numeric(sc_data[f,col_names[2]]),as.matrix(sc_data[f,c(col_names[3],col_names[4])]),k)
      #RespCor = toCorrCoeff(predsSmoothed)
      RespCor = predsSmoothed
      df = data.frame(cell_id = sc_data[f,col_names[1]], coord1 = sc_data[f,col_names[3]], coord2 = sc_data[f,col_names[4]], smoothed = RespCor)
      cor_list[[i]] = df
    }
    cor_list_df = do.call(rbind,cor_list)
    return(cor_list_df)
  }else{
    message(paste0("The number of neighbors (",k,") is greater than 50% of the sample (",n_split,")."))
  }
}


library(foreach)
library(doParallel)
library(progress)

#setup parallel backend to use many processors
cell_sets = cell_sets[sapply(cell_sets,function(x) dim(x)[1]>mincellsn)]
cores=detectCores()
cl <- makeCluster(cores[1]-4) #not to overload your computer
registerDoParallel(cl)

finalMatrix <- foreach(i=1:10, .combine=rbind) %dopar% {
  tempMatrix = knnSmooth(as.numeric(df[cell_sets[[i]]$barcode,"t2d_risk"]),as.matrix(df[cell_sets[[i]]$barcode,c("locsx","locsy")]),5) #calling a function
  
}
#stop cluster
stopCluster(cl)


tester <-function(a){
  pb <- progress_bar$new(total = 100)
  for (i in 1:100) {
    pb$tick()
    Sys.sleep(1 / 100)
  }
}
tester("a")


#make another fucntion swSmoothAtlasCOSMX, inctead of patch selection, use the FOV that are already provided
# replace the entire functions.

swSmoothAtlas <- function(dat,xsz,ysz,step,pad,mincellsn,kNN,colx,coly,colv){
  message("Start sliding window smoothing")
  xmin = min(dat[,colx])
  xmax = max(dat[,colx])
  ymin = min(dat[,coly])
  ymax = max(dat[,coly])
  xmarker = xmin
  ymarker = ymin
  
  # Retrieving barcodes for each patch
  message("Begin patch annotation")
  cell_sets = list()
  i=0
  finished = FALSE
  while(!finished){
    #get the barcodes for every FOV and put them here
    #each element in the list will be the barcodes for each FOV
    #set to keep=TRUE
    i = i+1
    if((xmarker+xsz)<xmax & (ymarker+ysz)<ymax){
      bc_tmp = dat$barcodes[dat[,colx]>xmarker&dat[,colx]<(xmarker+xsz)&dat[,coly]>ymarker&dat[,coly]<(ymarker+ysz)]
      keep_tmp = dat[dat[,colx]>xmarker&dat[,colx]<(xmarker+xsz)&dat[,coly]>ymarker&dat[,coly]<(ymarker+ysz),"locsx"]>(xmarker+pad) &
        dat[dat[,colx]>xmarker&dat[,colx]<(xmarker+xsz)&dat[,coly]>ymarker&dat[,coly]<(ymarker+ysz),"locsx"]<(xmarker+xsz-pad) &
        dat[dat[,colx]>xmarker&dat[,colx]<(xmarker+xsz)&dat[,coly]>ymarker&dat[,coly]<(ymarker+ysz),"locsy"]>(ymarker+pad) &
        dat[dat[,colx]>xmarker&dat[,colx]<(xmarker+xsz)&dat[,coly]>ymarker&dat[,coly]<(ymarker+ysz),"locsy"]<(ymarker+ysz+pad)
      cell_sets[[i]] = data.frame(barcode=bc_tmp,keep=keep_tmp) 
      xmarker = xmarker+step
    }else if((xmarker+xsz)>xmax & (ymarker+ysz)<ymax){
      bc_tmp = dat$barcodes[dat[,colx]>xmarker&dat[,colx]<(xmarker+xsz)&dat[,coly]>ymarker&dat[,coly]<(ymarker+ysz)]
      keep_tmp = dat[dat[,colx]>xmarker&dat[,colx]<(xmarker+xsz)&dat[,coly]>ymarker&dat[,coly]<(ymarker+ysz),"locsx"]>(xmarker+pad) &
        dat[dat[,colx]>xmarker&dat[,colx]<(xmarker+xsz)&dat[,coly]>ymarker&dat[,coly]<(ymarker+ysz),"locsx"]<(xmarker+xsz-pad) &
        dat[dat[,colx]>xmarker&dat[,colx]<(xmarker+xsz)&dat[,coly]>ymarker&dat[,coly]<(ymarker+ysz),"locsy"]>(ymarker+pad) &
        dat[dat[,colx]>xmarker&dat[,colx]<(xmarker+xsz)&dat[,coly]>ymarker&dat[,coly]<(ymarker+ysz),"locsy"]<(ymarker+ysz+pad)
      cell_sets[[i]] = data.frame(barcode=bc_tmp,keep=keep_tmp)
      xmarker = xmin
      ymarker = ymarker+step
    }else if((xmarker+xsz)<xmax & (ymarker+ysz)>ymax){
      bc_tmp = dat$barcodes[dat[,colx]>xmarker&dat[,colx]<(xmarker+xsz)&dat[,coly]>ymarker&dat[,coly]<(ymarker+ysz)]
      keep_tmp = dat[dat[,colx]>xmarker&dat[,colx]<(xmarker+xsz)&dat[,coly]>ymarker&dat[,coly]<(ymarker+ysz),"locsx"]>(xmarker+pad) &
        dat[dat[,colx]>xmarker&dat[,colx]<(xmarker+xsz)&dat[,coly]>ymarker&dat[,coly]<(ymarker+ysz),"locsx"]<(xmarker+xsz-pad) &
        dat[dat[,colx]>xmarker&dat[,colx]<(xmarker+xsz)&dat[,coly]>ymarker&dat[,coly]<(ymarker+ysz),"locsy"]>(ymarker+pad) &
        dat[dat[,colx]>xmarker&dat[,colx]<(xmarker+xsz)&dat[,coly]>ymarker&dat[,coly]<(ymarker+ysz),"locsy"]<(ymarker+ysz+pad)
      cell_sets[[i]] = data.frame(barcode=bc_tmp,keep=keep_tmp)
      xmarker = xmarker+step
    }else if((xmarker+xsz)>xmax & (ymarker+ysz)>ymax){
      bc_tmp = dat$barcodes[dat[,colx]>xmarker&dat[,colx]<(xmarker+xsz)&dat[,coly]>ymarker&dat[,coly]<(ymarker+ysz)]
      keep_tmp = dat[dat[,colx]>xmarker&dat[,colx]<(xmarker+xsz)&dat[,coly]>ymarker&dat[,coly]<(ymarker+ysz),"locsx"]>(xmarker+pad) &
        dat[dat[,colx]>xmarker&dat[,colx]<(xmarker+xsz)&dat[,coly]>ymarker&dat[,coly]<(ymarker+ysz),"locsx"]<(xmarker+xsz-pad) &
        dat[dat[,colx]>xmarker&dat[,colx]<(xmarker+xsz)&dat[,coly]>ymarker&dat[,coly]<(ymarker+ysz),"locsy"]>(ymarker+pad) &
        dat[dat[,colx]>xmarker&dat[,colx]<(xmarker+xsz)&dat[,coly]>ymarker&dat[,coly]<(ymarker+ysz),"locsy"]<(ymarker+ysz+pad)
      cell_sets[[i]] = data.frame(barcode=bc_tmp,keep=keep_tmp)
      finished = TRUE
    }else{
      message("ERROR: Encoding")
    }
  }
  len_ss = length(cell_sets)
  message(paste0("Patch Number: ",len_ss))
  message("Completed patch annotation")
  
  # Smoothing cells in each patch
  message("Begin patch smoothing")
  q21 = floor(quantile(1:len_ss,seq(from=0.0,to=1.0,by=0.05)))
  pb <- progress_bar$new(
    format = "[:bar] :percent",
    total = 21, clear = FALSE, width = 60)
  smoothed = list()
  j = 0
  for(cell_set in cell_sets){
    j=j+1
    if(j%in%q21){pb$tick()}
    Sys.sleep(2/100)
    if(dim(cell_set)[1]>mincellsn){
      smoothed_tmp = knnSmooth(as.numeric(dat[cell_set$barcode,colv]),as.matrix(dat[cell_set$barcode,c("locsx","locsy")]),kNN)
      cell_barcode_tmp = cell_set$barcode[cell_set$keep==TRUE]
      smoothed[[j]] = data.frame(barcode = cell_barcode_tmp, score = smoothed_tmp[cell_set$keep==TRUE])
    }
  }
  message("Completed patch smoothing")
  
  # Combining smoothed scores 
  message("Begin cell aggregation")
  smoothed = smoothed[sapply(smoothed,function(x) !is.null(dim(x)[1]))]
  smoothed = purrr::reduce(smoothed, full_join, by = "barcode")
  smoothed2 = smoothed
  rownames(smoothed2) = smoothed2$barcode
  smoothed2$barcode = NULL
  smoothed2$haz = rowMeans(smoothed2,na.rm=TRUE)
  out = data.frame(barcodes=rownames(smoothed2),smoothed=smoothed2$haz)
  message("Completed cell aggregation")
  return(out)
  message("Finished sliding window smoothing")
}
# set.seed(2)

meta_data_columns <- c("V1", "cell_id", "slide_ID_numeric", "Run_Tissue_name", "x_FOV_px", "y_FOV_px", "qcFlagsFOV", "cellType", "niche")

# What is the name of this experiment/sample?
sample_name <- "liver_coxmyx"

# Patient data
## RNA-seq
patDat <- data.table::fread('/N/project/degas_st/cosmyx/data/TCGA/patDat.csv', sep = ",") # Tumor and normal tissue
# patDat <- data.table::fread('/N/project/degas_st/cosmyx/data/TCGA/patDat_tumor.csv', sep = ",") # Only tumor tissue

## Clinical outcomes
patLab <- data.table::fread('/N/project/degas_st/cosmyx/data/Reference_liver/patLab.csv', sep = ",")

# Need non NA follow up time and status
patLab <- patLab %>% filter(!is.na(OS_time))

# only retain these same patients in patDat, and ordered the same way
patDat <- patDat %>% inner_join(patLab)

# ST data (cosmyx)
meta_df <- read_csv('/N/project/degas_st/CosMx_liver_data/cosmx_liver_metadata.csv') # To get the coordinates

st_data <- data.table::fread('/N/project/degas_st/cosmyx/data/processed/liver/cosmyx_liver_data.csv')


#plot(st_data$x_FOV_px,st_data$y_FOV_px)

colnames(st_data)

st_data <- st_data %>%
  mutate(FOV = str_replace(str_replace(cell_id, "c_[:digit:]_", ""), "_[:digit:].*$", ""), .before = cell_id)

common_cols <- intersect(meta_df %>% colnames, st_data %>% colnames)

st_data <- 
  left_join(
    st_data,
    meta_df %>% 
      select(
        all_of(common_cols), # the columns for joining data properly
    x_FOV_px, y_FOV_px)) # the actual real coordinates of the points to use for plotting

st_data <- st_data %>%
  relocate(contains('slide_mm'), .before = FOV)

# Find intersecting genes for feature selection
intersecting_genes <-
  intersect(colnames(st_data), colnames(patDat %>% select(-OS_status, -OS_time)))



# Feature selection (top 200 most variable genes)
genes <- patDat %>% select(all_of(intersecting_genes)) %>% apply(2, var) %>% sort(decreasing=TRUE) %>% names() %>% head(200)
#DEGAS.model = readRDS("/N/project/degas_st/cosmyx/degas_output/liver_coxmyx_200_2024-03-06.RDS")
DEGAS.model = readRDS("/N/project/degas_st/cosmyx/lihc_output_Nov14/lihc_allgenes_200_2025-11-15.RDS")

## Post processing smoothing

#must take both samples separately

#input sc data
#input the patient-level data
#run the prediction codes for each sample, call in the list of single cell data and then run the predictions for all
library(RColorBrewer)
unique(st_data$Run_Tissue_name)
st_data_normal = st_data[st_data$Run_Tissue_name=="NormalLiver",]
st_data_cancer=st_data[st_data$Run_Tissue_name=="CancerousLiver",]

st_data_list <- list(st_data_normal,st_data_cancer)
names(st_data_list) <- c("normal", "cancer")

# for( i in seq_along(st_data_list)){
#   st_data = st_data_list[[i]]

#   st_data2 = st_data%>%select(all_of(genes))
#   patDat2=patDat%>%select(all_of(genes))

#   scpatPreds = predClassBag(DEGAS.model,st_data2,"pat") #predClassBag(ccModel, Exp, scORpat)

#   model_output_scpat <-
#     tibble(
#       as_tibble(scpatPreds),
#       cell_id = st_data$cell_id,
#       fov = st_data$FOV,
#       coord1 = st_data$x_slide_mm,
#       coord2 = st_data$y_slide_mm,
#       sample = st_data$slide_ID_numeric,
#       cellType = st_data$cellType,
#       niche = st_data$niche)

#   colnames(model_output_scpat)[1] <- "Hazard"

#   model_output_scpat %<>% mutate(Haz_scaled = (Hazard - min(Hazard)) / (max(Hazard)- min(Hazard)), .before = 'Hazard')

#   write.csv(model_output_scpat,file=paste0("/N/project/degas_st/cosmyx/degas_Debolina/model_output_scpat_",names(st_data_list)[i],".csv"))
# }

model_output_scpat <- read.csv("/N/project/degas_st/cosmyx/lihc_output_Nov14/lihc_allgenes_200_2025-11-15.csv")

umap.df <- read.csv("/N/project/degas_st/cosmyx/degas_Debolina/umap_df.csv")
colnames(umap.df) <- c("cell_id","UMAP.1","UMAP.2")
# Create a new column based on row names
umap.df$sample <- ifelse(grepl("^c_1", umap.df$cell_id), "NormalLiver",
                         ifelse(grepl("^c_2", umap.df$cell_id), "CancerousLiver", NA))

#sum(is.na(umap.df$sample))
head(umap.df)

umap.df.normal <- umap.df[which(umap.df$sample=="NormalLiver"),] #umap embeddings for normal sample
umap.df.cancer <- umap.df[which(umap.df$sample=="CancerousLiver"),] #umap embeddings for cancer sample

meta_df <- read_csv('/N/project/degas_st/CosMx_liver_data/cosmx_liver_metadata.csv') # To get the spatial coordinates
#meta_df <- read_csv("/N/u/dchatter/Quartz/thindrives/OneDrive-IndianaUniversity/DEGAS/CosMx liver data/cosmx_liver_metadata.csv")
colnames(meta_df)
st_data <- data.table::fread('/N/project/degas_st/cosmyx/data/processed/liver/cosmyx_liver_data.csv')


#all.equal(meta_df$...1,meta_df$cell_id)
#separate this as normal and cancer
# Create a new column based on row names
meta_df$sample <- ifelse(grepl("^c_1", meta_df$cell_id), "NormalLiver",
                         ifelse(grepl("^c_2", meta_df$cell_id), "CancerousLiver", NA))
meta_df.normal <- meta_df[which(meta_df$sample=="NormalLiver"),] #umap embeddings for normal sample
meta_df.cancer <- meta_df[which(meta_df$sample=="CancerousLiver"),] #umap embeddings for cancer sample

# Sort both data frames by cell_id to ensure the rows align
meta_df.normal <- meta_df.normal[order(meta_df.normal$cell_id), ]
meta_df.cancer <- meta_df.cancer[order(meta_df.cancer$cell_id), ]

head(meta_df.normal)


st_data$sample <- ifelse(grepl("^c_1", st_data$cell_id), "NormalLiver",
                         ifelse(grepl("^c_2", st_data$cell_id), "CancerousLiver", NA))

st_data.normal <- st_data[which(st_data$sample=="NormalLiver"),] #umap embeddings for normal sample
st_data.cancer <- st_data[which(st_data$sample=="CancerousLiver"),] #umap embeddings for cancer sample

# Sort both data frames by cell_id to ensure the rows align
st_data.normal <- st_data.normal[order(st_data.normal$cell_id), ]
st_data.cancer <- st_data.cancer[order(st_data.cancer$cell_id), ]


#merge the UMAP coordiantes with the st_data
meta_df <- meta_df[order(meta_df$cell_id),]
meta_df <- merge(meta_df,umap.df,by="cell_id")
st_data.normal <- merge(st_data.normal,umap.df.normal, by="cell_id")
st_data.cancer <- merge(st_data.cancer,umap.df.cancer, by="cell_id")

head(st_data.normal)

vec <- meta_df$cell_id
# Split by underscore
split_vec <- strsplit(vec, "_")

# Convert to a data frame for easier viewing
df <- as.data.frame(do.call(rbind, split_vec))
colnames(df) <- c("prefix", "group", "fov", "id")
head(df)
meta_df_full <- data.frame(cbind(meta_df$cell_id,meta_df$Run_Tissue_name,meta_df$x_FOV_px,meta_df$y_FOV_px,meta_df$x_slide_mm,meta_df$y_slide_mm,meta_df$UMAP.1,meta_df$UMAP.2,df$group,df$fov,df$id))
colnames(meta_df_full) <- c("cell_id","Run_Tissue_name","x_FOV_px","y_FOV_px","x_slide_mm","y_slide_mm","UMAP.1","UMAP.2","group","fov","id")

str(meta_df_full)

#the normal sample-----

df_n <- meta_df_full[which(meta_df_full$group==1),]
# View result
head(df_n)
dim(df_n)
length(unique(df_n$fov))
length(unique(df_n$id))
# Find duplicated 'id' values across different 'subgroup' values
library(dplyr)

duplicates <- df_n %>%
  group_by(id) %>%
  filter(n_distinct(fov) >= 1) %>%
  arrange(id, fov)
dim(duplicates)
# View the results
#View(duplicates)


df_n_summary <- df_n %>%
  group_by(fov) %>%
  summarise(cell_count = n())
#View(df_n_summary)
summary(df_n_summary$cell_count)

# cancer sample----
df_c <- meta_df_full[which(meta_df_full$group==2),]
# View result
head(df_c)
dim(df_c)
length(unique(df_c$fov))
length(unique(df_c$id))
# Find duplicated 'id' values across different 'subgroup' values

duplicates <- df_c %>%
  group_by(id) %>%
  filter(n_distinct(fov) >= 1) %>%
  arrange(id, fov)
dim(duplicates)
# View the results
#View(duplicates)
df_c_summary <- df_c %>%
  group_by(fov) %>%
  summarise(cell_count = n())

df_c_summary
summary(df_c_summary$cell_count)


dim(meta_df.normal)
dim(normal_pred)
all.equal(meta_df.normal$cellType,normal_pred$cellType)

model_output_scpat[1:5,1:5]
###### Smoothing for the Normal sample
normal_pred <- model_output_scpat[model_output_scpat$sample==1,]
str(normal_pred)
normal_pred <- normal_pred[order(normal_pred$cell_id),]

#head(normal_pred)
##head(umap.df.normal)
sc_data_all <- merge(normal_pred,umap.df.normal,by="cell_id")
#head(sc_data_all)
#colnames(sc_data_all)

head(normal_pred)
# Sort both data frames by cell_id to ensure the rows align
normal_pred <- normal_pred[order(normal_pred$cell_id), ]
df_n <- df_n[order(df_n$cell_id), ]

normal.p <- merge(normal_pred,df_n, by="cell_id")
head(normal.p)

normal_pred$x_FOV_px <- df_n$x_FOV_px
normal_pred$y_FOV_px <- df_n$y_FOV_px
normal_pred$UMAP.1 <- df_n$UMAP.1
normal_pred$UMAP.2 <- df_n$UMAP.2
library(progress)


sc_data_all <- merge(normal_pred,umap.df.normal,by="cell_id")

head(sc_data_all)
dim(sc_data_all)
##### option 1: run the simple KNN k=20 and k =50 smoothing that was already done

#spatial smoothing only
sptial_sm <- knnSmoothAtlas2(sc_data_all,c("cell_id","Hazard","coord1","coord2"),50,10000)
rownames(sptial_sm) = sptial_sm$cell_id
sptial_sm$scaled_smoothed_haz= toCorrCoeff(sptial_sm$smoothed)
head(sptial_sm)

p1 <- ggplot(sptial_sm,aes(x=coord1,y=coord2,color=scaled_smoothed_haz)) + 
  geom_point(alpha=0.9, size=0.5) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(sptial_sm$scaled_smoothed_haz, na.rm = TRUE), 
                        name = "LIHC\nAssociation")+
  labs(
    x = "X coordinates",  # Replace with your desired label
    y = "Y coordinates"   # Replace with your desired label
  ) +
  theme_bw() +  # Sets white background
  ggtitle("Non-carcinoma")+  # Adds a title
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),         # Remove major gridlines
        panel.grid.minor = element_blank(),         # Remove minor gridlines 
        plot.title = element_text(size = 18, hjust=0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),  # Keep axis lines
        axis.ticks = element_line(color = "black"), # Keep axis ticks
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  )
#ggsave("/N/project/degas_st/cosmyx/degas_Debolina/Normal/spatial_smoothed_normal_Knn20_06262025.pdf",p1)
ggsave("/N/project/degas_st/cosmyx/degas_Debolina/Normal/spatial_smoothed_normal_Knn50_11132025.pdf",p1)


#### option 2: include he FOV based smoothign with it as well


smooth_by_fov <- function(df, kNN,locx,locy, mincells) {
  unique_fovs <- unique(df$fov)
  smoothed_list <- list()
  pb <- progress_bar$new(format = "[:bar] :percent", total = length(unique_fovs), clear = FALSE, width = 60)
  
  for (i in unique_fovs) {
    fov_id <- i
    pb$tick()
    
    fov_data <- df %>% filter(fov == fov_id)
    fov_data$x_FOV_px <- as.numeric(fov_data$x_FOV_px)
    fov_data$y_FOV_px <- as.numeric(fov_data$y_FOV_px)
    if (nrow(fov_data) >= mincells) {
      coords <- as.matrix(fov_data[, c("x_FOV_px", "y_FOV_px")])
      values <- fov_data$Hazard
      #values <- toCorrCoeff(fov_data$Hazard)
      smoothed_vals <- knnSmooth(values, coords, k = kNN)
      
      smoothed_list[[i]] <- fov_data %>%
        mutate(smoothed = smoothed_vals) %>%
        select(cell_id, fov, coord1, coord2, x_FOV_px, y_FOV_px, Hazard, smoothed)
      
      
    }else if (nrow(fov_data) > 0 & nrow(fov_data) < mincells) {
      # Not enough cells — keep raw predictions instead of smoothing
      message(paste("FOV", fov_id, "has fewer than", mincells, "cells. Using raw predictions."))
      
      smoothed_list[[i]] <- fov_data %>%
        mutate(smoothed = Hazard) %>%
        select(cell_id, fov, coord1, coord2, x_FOV_px, y_FOV_px, Hazard, smoothed)
      
      
    }
  }
  message("Completed FOV-based smoothing")
  
  # Combine results
  message("Begin aggregation")
  smoothed_combined <- do.call(rbind, smoothed_list)
  message("Finished smoothing and aggregation")
  
  return(smoothed_combined)
}

# run the smoothing
#if the number of cells < mincellsn, then do not perform smoothing. otherwise perform smoothing
result_normal <- smooth_by_fov(sc_data_all, kNN = 50, mincells = 51)
head(result_normal)
result_normal$scaled_FOV_smoothed <- toCorrCoeff(result_normal$smoothed)
saveRDS(result_normal,file="/N/project/degas_st/cosmyx/degas_Debolina/smoothed_normal_kNN50_11132025.RDS")
#plot the smoothed values

p1 <- ggplot(result_normal,aes(x=coord1,y=coord2,color=scaled_FOV_smoothed)) + 
  geom_point(alpha=0.9,size=0.5) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(result_normal$scaled_FOV_smoothed, na.rm = TRUE), 
                        name = "LIHC_risk")+
  ggtitle("Smoothed LIHC risk")
ggsave("/N/project/degas_st/cosmyx/degas_Debolina/Normal/scaled_FOV_smoothed_normal_Knn50_11132025.pdf",p1)

#outlier removed version also



# Remove outliers from result_normal$scaled_FOV_smoothed using the IQR method
# Step 1: Calculate bounds
x <- result_normal$scaled_FOV_smoothed
#x <- result_normal$smoothed

Q1 <- quantile(x, 0.25, na.rm = TRUE)
Q3 <- quantile(x, 0.75, na.rm = TRUE)
IQR_val <- Q3 - Q1

lower_bound <- Q1 - 1.5 * IQR_val
upper_bound <- Q3 + 1.5 * IQR_val

# Step 2: Identify non-outlier indices
non_outlier_idx <- which(x >= lower_bound & x <= upper_bound)

# Step 3: Filter the vector
result_normal$scaled_FOV_smoothed_no_outliers <- x
result_normal$scaled_FOV_smoothed_no_outliers[-non_outlier_idx] <- NA  # Set outliers to NA

# Optional: Create a new data frame with only non-outlier rows
result_normal_no_outliers <- result_normal[non_outlier_idx, ]

# View result
print(head(result_normal_no_outliers))
dim(result_normal_no_outliers)
head(result_normal_no_outliers)


p12 <- ggplot(result_normal_no_outliers,aes(x=coord1,y=coord2,color=scaled_FOV_smoothed_no_outliers)) + 
  geom_point(alpha=0.9, size=0.5) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(result_normal_no_outliers$scaled_FOV_smoothed_no_outliers, na.rm = TRUE), 
                        name = "LIHC\nAssociation")+
  labs(
    x = "X coordinates",  # Replace with your desired label
    y = "Y coordinates"   # Replace with your desired label
  ) +
  theme_bw() +  # Sets white background
  ggtitle("Non-carcinoma")+  # Adds a title
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),         # Remove major gridlines
        panel.grid.minor = element_blank(),         # Remove minor gridlines 
        plot.title = element_text(size = 18, hjust=0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),  # Keep axis lines
        axis.ticks = element_line(color = "black"), # Keep axis ticks
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  )

ggsave("/N/project/degas_st/cosmyx/degas_Debolina/Normal/FOV_smoothed_normal__no_outliers_Knn50_11132025.pdf",p12)




#################################################################################################

###### Smoothing for the Cancer sample----
cancer_pred <- read.csv("/N/project/degas_st/cosmyx/degas_Debolina/model_output_scpat_cancer.csv")
cancer_pred <- cancer_pred[order(cancer_pred$cell_id),]

sc_data_all <- merge(cancer_pred,umap.df.cancer,by="cell_id")
#head(sc_data_all)

head(cancer_pred)
cancer_pred$x_FOV_px <- df_c$x_FOV_px
cancer_pred$y_FOV_px <- df_c$y_FOV_px


sc_data_all_cancer <- merge(cancer_pred,umap.df.cancer,by="cell_id")
head(sc_data_all_cancer)
#spatial smoothing only

#spatial smoothing only
sptial_sm_cancer<- knnSmoothAtlas2(sc_data_all_cancer,c("cell_id","Hazard","coord1","coord2"),50,10000)
rownames(sptial_sm_cancer) = sptial_sm_cancer$cell_id
sptial_sm_cancer$scaled_smoothed_haz= toCorrCoeff(sptial_sm_cancer$smoothed)
head(sptial_sm_cancer)
saveRDS(sptial_sm_cancer,"/N/project/degas_st/cosmyx/degas_Debolina/Cancer/spatial_sm_cancer_knn50_11132025.rds" )
p1 <- ggplot(sptial_sm_cancer,aes(x=coord1,y=coord2,color=scaled_smoothed_haz)) + 
  geom_point(alpha=0.9, size=0.5) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(sptial_sm_cancer$scaled_smoothed_haz, na.rm = TRUE), 
                        name = "LIHC\nAssociation")+
  labs(
    x = "X coordinates",  # Replace with your desired label
    y = "Y coordinates"   # Replace with your desired label
  ) +
  theme_bw() +  # Sets white background
  ggtitle("Carcinoma")+  # Adds a title
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),         # Remove major gridlines
        panel.grid.minor = element_blank(),         # Remove minor gridlines 
        plot.title = element_text(size = 18, hjust=0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),  # Keep axis lines
        axis.ticks = element_line(color = "black"), # Keep axis ticks
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  )

ggsave("/N/project/degas_st/cosmyx/degas_Debolina/Cancer/spatial_smoothed_normal_Knn50_11132025.pdf",p1)




#FOV smoothing
#if the number of cells < mincellsn, then do not perform smoothing. otherwise perform smoothing
result_cancer <- smooth_by_fov(sc_data_all_cancer, kNN = 50, mincells = 51)
head(result_cancer)
result_cancer$scaled_FOV_smoothed <- toCorrCoeff(result_cancer$smoothed)
saveRDS(result_cancer,file="/N/project/degas_st/cosmyx/degas_Debolina/smoothed_cancer_kNN50_11132025.RDS")
#plot the smoothed values

p2 <- ggplot(result_cancer,aes(x=coord1,y=coord2,color=smoothed)) + 
  geom_point(alpha=0.9,size=0.5) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(result_cancer$smoothed, na.rm = TRUE), 
                        name = "LIHC_risk")+
  ggtitle("Smoothed LIHC risk")
ggsave("/N/project/degas_st/cosmyx/degas_Debolina/Cancer/FOV_smoothed_cancer_Knn50_11132025.pdf",p2)


p2 <- ggplot(result_cancer,aes(x=coord1,y=coord2,color=scaled_FOV_smoothed)) + 
  geom_point(alpha=0.9,size=0.5) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(result_cancer$scaled_FOV_smoothed, na.rm = TRUE), 
                        name = "LIHC_risk")+
  ggtitle("Smoothed LIHC risk")
ggsave("/N/project/degas_st/cosmyx/degas_Debolina/Cancer/scaled_FOV_smoothed_cancer_Knn50_11132025.pdf",p2)

x <- result_cancer$scaled_FOV_smoothed
#x <- result_cancer$smoothed

Q1 <- quantile(x, 0.25, na.rm = TRUE)
Q3 <- quantile(x, 0.75, na.rm = TRUE)
IQR_val <- Q3 - Q1

lower_bound <- Q1 - 1.5 * IQR_val
upper_bound <- Q3 + 1.5 * IQR_val

# Step 2: Identify non-outlier indices
non_outlier_idx <- which(x >= lower_bound & x <= upper_bound)

# Step 3: Filter the vector
result_cancer$scaled_FOV_smoothed_no_outliers <- x
result_cancer$scaled_FOV_smoothed_no_outliers[-non_outlier_idx] <- NA  # Set outliers to NA

# Optional: Create a new data frame with only non-outlier rows
result_cancer_no_outliers <- result_cancer[non_outlier_idx, ]

# View result
print(head(result_cancer_no_outliers))
dim(result_cancer_no_outliers)
head(result_cancer_no_outliers)


p12 <- ggplot(result_cancer_no_outliers,aes(x=coord1,y=coord2,color=scaled_FOV_smoothed_no_outliers)) + 
  geom_point(alpha=0.9, size=0.5) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(result_cancer_no_outliers$scaled_FOV_smoothed_no_outliers, na.rm = TRUE), 
                        name = "LIHC\nAssociation")+
  labs(
    x = "X coordinates",  # Replace with your desired label
    y = "Y coordinates"   # Replace with your desired label
  ) +
  theme_bw() +  # Sets white background
  ggtitle("Non-carcinoma")+  # Adds a title
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),         # Remove major gridlines
        panel.grid.minor = element_blank(),         # Remove minor gridlines 
        plot.title = element_text(size = 18, hjust=0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),  # Keep axis lines
        axis.ticks = element_line(color = "black"), # Keep axis ticks
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  )

ggsave("/N/project/degas_st/cosmyx/degas_Debolina/Cancer/FOV_smoothed_cancer__no_outliers_Knn50_11132025.pdf",p12)








##############################################################################################
#make the boxplots for Normal vs Cancer

result_normal=readRDS("/N/project/degas_st/cosmyx/degas_Debolina/smoothed_normal_kNN50_11132025.RDS")
result_cancer=readRDS("/N/project/degas_st/cosmyx/degas_Debolina/smoothed_cancer_kNN50_11132025.RDS")

head(result_normal)
head(result_cancer)

# Add group labels
result_normal$condition <- "Non-carcinoma"
result_cancer$condition <- "Carcinoma"

# Combine both data frames
combined <- rbind(result_normal, result_cancer)
combined$condition <- factor(combined$condition, levels = c("Non-carcinoma", "Carcinoma"))

head(combined)

library(ggpubr)

p3 <- ggboxplot(combined, x = "condition", y = "scaled_FOV_smoothed", 
          fill = "condition", palette = c("Non-carcinoma" = "blue", "Carcinoma" = "red")) +
  stat_compare_means(method = "t.test", label.y = max(combined$scaled_FOV_smoothe) * 1.05) +
  labs(y = "LIHC Association")

ggsave("/N/project/degas_st/cosmyx/degas_Debolina/Normal_vs_cancer_boxplot_raw.png",p3)
# Perform t-test
t_test_result <- t.test(Hazard ~ condition, data = combined)

p4 <- ggboxplot(combined, x = "condition", y = "smoothed", 
                fill = "condition", palette = c("normal" = "blue", "cancer" = "red")) +
  stat_compare_means(method = "t.test", label.y = max(combined$smoothed) * 1.05) +
  labs(y = "Smoothed Score", x = "Condition")

ggsave("/N/project/degas_st/cosmyx/degas_Debolina/Normal_vs_cancer_boxplot_smoothed_knn20.png",p4)
# Perform t-test
t_test_result <- t.test( scaled_FOV_smoothed_no_outliers ~ condition, data = combined)
if (t_test_result$p.value < 2.2e-16) {
  p_value_text <- "p < 2.2e-16"
} else {
  p_value_text <- paste0("p = ", formatC(t_test_result$p.value, format = "e", digits = 2))
}

combined$scaled_haz = toCorrCoeff(combined$Hazard)
combined$scaled_spatial_smoothed_haz = toCorrCoeff(combined$smoothed)
head(combined)
# Your plot with ggrepel
boxplot1 <- ggplot(combined, aes(x = condition, y = scaled_haz, fill = condition)) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),         # Remove major gridlines
        panel.grid.minor = element_blank(),         # Remove minor gridlines 
        plot.title = element_text(size = 18, hjust=0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),  # Keep axis lines
        axis.ticks = element_line(color = "black"), # Keep axis ticks
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  )+
  annotate("text", x = 1.5, y = max(combined$scaled_haz, na.rm = TRUE) * 1.05, 
           label = p_value_text, size = 5)+
  ylab("LIHC Association") +
  xlab("Sample") +
  scale_fill_manual(values = c("Non-carcinoma" = "blue", "Carcinoma" = "red")) 

ggsave("/N/project/degas_st/cosmyx/degas_Debolina/boxplot_scaled_raw_haz_06262025.pdf",boxplot1)

head(combined)

meta_df_sub <- meta_df[,c("cell_id","cellType")]
head(meta_df_sub)

combined_withcelltype <- merge(combined,meta_df_sub,by="cell_id")
head(combined_withcelltype)

library(ggplot2)
library(dplyr)

combined_withcelltype <- combined_withcelltype %>%
  mutate(scaled_haz = scaled_haz - median(scaled_haz, na.rm = TRUE))


# Compute mean scaled_haz per cellType and condition
mean_data <- combined_withcelltype %>%
  group_by(condition, cellType) %>%
  summarise(mean_scaled_haz = mean(scaled_haz, na.rm = TRUE), .groups = "drop")

# Optional: order cellType by mean value within each condition
mean_data <- mean_data %>%
  group_by(condition) %>%
  mutate(cellType = reorder(cellType, mean_scaled_haz)) %>%
  ungroup()

# Horizontal bar plot
barplot_means_horizontal <- ggplot(mean_data, aes(x = mean_scaled_haz, y = cellType)) +
  geom_col(fill = "gray40", width = 0.6) +
  facet_wrap(~ condition, scales = "free_y") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 16)
  ) +
  ylab(" ") +
  xlab("Average LIHC Association Scores")

# Save the plot
ggsave("/N/project/degas_st/cosmyx/degas_Debolina/barplot_scaledhaz_bycelltype_06262025.pdf",
       barplot_means_horizontal, width = 10, height = 6)




###############################

#spatial smoothing only
sptial_sm <- knnSmoothAtlas2(sc_data_all,c("cell_id","Hazard","coord1","coord2"),20,10000)
rownames(sptial_sm) = sptial_sm$cell_id

sptial_sm$scaled_smoothed_haz= toCorrCoeff(sptial_sm$smoothed)


p1 <- ggplot(sptial_sm,aes(x=coord1,y=coord2,color=scaled_smoothed_haz)) + 
  geom_point(alpha=0.9, size=0.5) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(sptial_sm$scaled_smoothed_haz, na.rm = TRUE), 
                        name = "LIHC\nAssociation")+
  labs(
    x = "X coordinates",  # Replace with your desired label
    y = "Y coordinates"   # Replace with your desired label
  ) +
  theme_bw() +  # Sets white background
  ggtitle("Non-carcinoma")+  # Adds a title
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),         # Remove major gridlines
        panel.grid.minor = element_blank(),         # Remove minor gridlines 
        plot.title = element_text(size = 18, hjust=0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),  # Keep axis lines
        axis.ticks = element_line(color = "black"), # Keep axis ticks
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  )
ggsave("/N/project/degas_st/cosmyx/degas_Debolina/Normal/spatial_smoothed_normal_Knn20_06262025.pdf",p1)

head(sptial_sm)

normal_UMAP_sm = readRDS("/N/project/degas_st/cosmyx/degas_Debolina/Normal/Normal_sm_20_10000.rds")
colnames(normal_UMAP_sm)
normal_UMAP_smoothed <- normal_UMAP_sm[,c("cell_id","shit")]
colnames(normal_UMAP_smoothed) <- c("cell_id","UMAP_sm_pred")

normal_pred_v2 <- merge(normal_pred,normal_UMAP_smoothed,by="cell_id")
colnames(normal_pred_v2)

smooth_by_fov <- function(df, kNN,locx,locy, mincells) {
  unique_fovs <- unique(df$fov)
  smoothed_list <- list()
  pb <- progress_bar$new(format = "[:bar] :percent", total = length(unique_fovs), clear = FALSE, width = 60)

  for (i in seq_along(unique_fovs)) {
    fov_id <- unique_fovs[i]
    pb$tick()

    fov_data <- df %>% filter(fov == fov_id)
    fov_data$x_FOV_px <- as.numeric(fov_data$x_FOV_px)
    fov_data$y_FOV_px <- as.numeric(fov_data$y_FOV_px)
    if (nrow(fov_data) >= mincells) {
      coords <- as.matrix(fov_data[, c("x_FOV_px", "y_FOV_px")])
      values <- fov_data$Hazard

      smoothed_vals <- knnSmooth(values, coords, k = kNN)

      smoothed_list[[i]] <- fov_data %>%
        mutate(smoothed = smoothed_vals) %>%
        select(cell_id, fov, coord1, coord2, x_FOV_px, y_FOV_px, Hazard, smoothed)


    }else if (nrow(fov_data) > 0 & nrow(fov_data) < mincells) {
    # Not enough cells — keep raw predictions instead of smoothing
    message(paste("FOV", fov_id, "has fewer than", mincells, "cells. Using raw predictions."))
    
    smoothed_list[[i]] <- fov_data %>%
        mutate(smoothed = Hazard) %>%
        select(cell_id, fov, coord1, coord2, x_FOV_px, y_FOV_px, Hazard, smoothed)
 

    }
  }
 message("Completed FOV-based smoothing")
  
  # Combine results
  message("Begin aggregation")
  smoothed_combined <- do.call(rbind, smoothed_list)
  message("Finished smoothing and aggregation")
  
  return(smoothed_combined)
}

# run the smoothing
#if the number of cells < mincellsn, then do not perform smoothing. otherwise perform smoothing

result_normal <- smooth_by_fov(sc_data_all, kNN = 20, mincells = 21)
head(result_normal)

saveRDS(result_normal,file="/N/project/degas_st/cosmyx/degas_Debolina/FOV_smoothed_normal_kNN20_06262025.RDS")
#plot the smoothed values
result_normal$scaled_FOV_smoothed = toCorrCoeff(result_normal$smoothed)
summary(result_normal$scaled_FOV_smoothed)
boxplot((result_normal$scaled_FOV_smoothed))
p11 <- ggplot(result_normal,aes(x=coord1,y=coord2,color=scaled_FOV_smoothed)) + 
  geom_point(alpha=0.9, size=0.5) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(result_normal$scaled_FOV_smoothed, na.rm = TRUE), 
                        name = "LIHC\nAssociation")+
  labs(
    x = "X coordinates",  # Replace with your desired label
    y = "Y coordinates"   # Replace with your desired label
  ) +
  theme_bw() +  # Sets white background
  ggtitle("Non-carcinoma")+  # Adds a title
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),         # Remove major gridlines
        panel.grid.minor = element_blank(),         # Remove minor gridlines 
        plot.title = element_text(size = 18, hjust=0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),  # Keep axis lines
        axis.ticks = element_line(color = "black"), # Keep axis ticks
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  )
ggsave("/N/project/degas_st/cosmyx/degas_Debolina/Normal/FOV_smoothed_normal_Knn20_06262025.png",p11)


# Remove outliers from result_normal$scaled_FOV_smoothed using the IQR method

# Step 1: Calculate bounds
x <- result_normal$scaled_FOV_smoothed
Q1 <- quantile(x, 0.25, na.rm = TRUE)
Q3 <- quantile(x, 0.75, na.rm = TRUE)
IQR_val <- Q3 - Q1

lower_bound <- Q1 - 1.5 * IQR_val
upper_bound <- Q3 + 1.5 * IQR_val

# Step 2: Identify non-outlier indices
non_outlier_idx <- which(x >= lower_bound & x <= upper_bound)

# Step 3: Filter the vector
result_normal$scaled_FOV_smoothed_no_outliers <- x
result_normal$scaled_FOV_smoothed_no_outliers[-non_outlier_idx] <- NA  # Set outliers to NA

# Optional: Create a new data frame with only non-outlier rows
result_normal_no_outliers <- result_normal[non_outlier_idx, ]

# View result
print(head(result_normal_no_outliers))
dim(result_normal_no_outliers)



p12 <- ggplot(result_normal_no_outliers,aes(x=coord1,y=coord2,color=scaled_FOV_smoothed)) + 
  geom_point(alpha=0.9, size=0.5) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(result_normal_no_outliers$scaled_FOV_smoothed, na.rm = TRUE), 
                        name = "LIHC\nAssociation")+
  labs(
    x = "X coordinates",  # Replace with your desired label
    y = "Y coordinates"   # Replace with your desired label
  ) +
  theme_bw() +  # Sets white background
  ggtitle("Non-carcinoma")+  # Adds a title
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),         # Remove major gridlines
        panel.grid.minor = element_blank(),         # Remove minor gridlines 
        plot.title = element_text(size = 18, hjust=0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),  # Keep axis lines
        axis.ticks = element_line(color = "black"), # Keep axis ticks
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  )
ggsave("/N/project/degas_st/cosmyx/degas_Debolina/Normal/FOV_smoothed_normal__no_outliers_Knn20_06262025.pdf",p12)








##cancer

result_cancer <- smooth_by_fov(sc_data_all_cancer, kNN = 20, mincells = 21)
saveRDS(result_cancer,file="/N/project/degas_st/cosmyx/degas_Debolina/FOV_smoothed_cancer_kNN20_06262025.RDS")
#plot the smoothed values
result_cancer$scaled_FOV_smoothed = toCorrCoeff(result_cancer$smoothed)
p22 <- ggplot(result_cancer,aes(x=coord1,y=coord2,color=scaled_FOV_smoothed)) + 
  geom_point(alpha=0.9, size=0.5) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(result_cancer$scaled_FOV_smoothed, na.rm = TRUE), 
                        name = "LIHC\nAssociation")+
  labs(
    x = "X coordinates",  # Replace with your desired label
    y = "Y coordinates"   # Replace with your desired label
  ) +
  theme_bw() +  # Sets white background
  ggtitle("Carcinoma")+  # Adds a title
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),         # Remove major gridlines
        panel.grid.minor = element_blank(),         # Remove minor gridlines 
        plot.title = element_text(size = 18, hjust=0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),  # Keep axis lines
        axis.ticks = element_line(color = "black"), # Keep axis ticks
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  )
ggsave("/N/project/degas_st/cosmyx/degas_Debolina/Cancer/FOV_smoothed_cancer_Knn20_06262025.png",p22)

# Step 1: Calculate bounds
x <- result_cancer$scaled_FOV_smoothed
Q1 <- quantile(x, 0.25, na.rm = TRUE)
Q3 <- quantile(x, 0.75, na.rm = TRUE)
IQR_val <- Q3 - Q1

lower_bound <- Q1 - 1.5 * IQR_val
upper_bound <- Q3 + 1.5 * IQR_val

# Step 2: Identify non-outlier indices
non_outlier_idx <- which(x >= lower_bound & x <= upper_bound)

# Step 3: Filter the vector
result_cancer$scaled_FOV_smoothed_no_outliers <- x
result_cancer$scaled_FOV_smoothed_no_outliers[-non_outlier_idx] <- NA  # Set outliers to NA

# Optional: Create a new data frame with only non-outlier rows
result_cancer_no_outliers <- result_cancer[non_outlier_idx, ]

# View result
print(head(result_cancer_no_outliers))
dim(result_cancer_no_outliers)


p22 <- ggplot(result_cancer_no_outliers,aes(x=coord1,y=coord2,color=scaled_FOV_smoothed)) + 
  geom_point(alpha=0.9, size=0.5) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(result_cancer_no_outliers$scaled_FOV_smoothed, na.rm = TRUE), 
                        name = "LIHC\nAssociation")+
  labs(
    x = "X coordinates",  # Replace with your desired label
    y = "Y coordinates"   # Replace with your desired label
  ) +
  theme_bw() +  # Sets white background
  ggtitle("Carcinoma")+  # Adds a title
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),         # Remove major gridlines
        panel.grid.minor = element_blank(),         # Remove minor gridlines 
        plot.title = element_text(size = 18, hjust=0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black"),  # Keep axis lines
        axis.ticks = element_line(color = "black"), # Keep axis ticks
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  )
ggsave("/N/project/degas_st/cosmyx/degas_Debolina/Cancer/FOV_smoothed_cancer__no_outliers_Knn20_06262025.png",p22)




library(pROC)
model_output_patpat <- read.csv("/N/project/degas_st/cosmyx/degas_Debolina/model_output_patpat.csv")
head(model_output_patpat)
  # AUC in patients
roc_result <- roc(model_output_patpat$OS_status,model_output_patpat$Hazard)

pdf(file=paste0("/N/project/degas_st/cosmyx/degas_Debolina/ROC_curve_06262025.pdf"))
plot(roc_result, col = "blue", lwd = 3, legacy.axes = TRUE,
     main = "", xlab = "", ylab = "", axes = FALSE, frame.plot = FALSE)
# Add custom axes
axis(1, cex.axis = 1.4)
axis(2, cex.axis = 1.4, las=1)
mtext("1 - Specificity", side = 1, line = 3, cex = 1.4, font = 1)
mtext("Sensitivity", side = 2, line = 3, cex = 1.4, font = 1)
title(main = "ROC Curve", cex.main = 1.4, font.main = 1)

# Add AUC text
text(0.5, 0.1, labels = paste("AUC =", round(as.numeric(roc_result$auc), 4)), col = "red", cex = 1.4, font = 2)

dev.off()
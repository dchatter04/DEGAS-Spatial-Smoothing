########################################################################
# MAGIC Lab library v1.0
# New functions (3/6/2025)

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


################################################################################
# Normalization, standardization, and transformation
#z-score normalization
normFunc <- function(x){return((x-mean(x, na.rm = T))/(sd(x, na.rm = T)+1e-3))}

# scaling from 0-1
scaleFunc <- function(x){return((x- min(x, na.rm = T)) /(max(x, na.rm = T)-min(x, na.rm = T)+1e-3))}

# Preprocess count data
normalizeScale <-function(X){
  return(t(apply(t(apply(as.matrix(t(X)),1,normFunc)),1,scaleFunc)))
}

# Sigmoid activation function
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

# Log sum exp transformation (for softmax)
logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

# Softmax activation function
softmax <- function (X) {
  return(t(apply(X,1,function(x) exp(x - logsumexp(x)))))
}

# Onehot labels to list of label names
fromOneHot <- function(labMat){
  return(apply(labMat,1,function(x) colnames(labMat)[which(x==1)]))
}

# List of label names to a onehot matrix with labels as column names
toOneHot <- function(labels){
  labs = unique(labels)
  out = matrix(0,length(labels),length(labs))
  colnames(out) = labs
  row.names(out) = row.names(labels)
  for(i in 1:length(labels)){
    out[i,labels[i]] = 1
  }
  return(out)
}

# Convert matrix of output weights to max value for each row (row max = 1 and not row max = 0)
probtoOneHot <- function(probMat){
  idx = apply(probMat,1,function(x) which(x==max(x)))
  probMat = probMat*0
  for (i in 1:length(idx)){
    probMat[i,idx[i]] = 1
  }
  return(probMat)
}

# Quantile normalization
quantNorm <- function(df,center='median',rescale=TRUE,rescale_mult=1e4){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  meds = eval(parse(text=paste0("apply(df_final, 2, ",center,")")))
  df_final = sweep(df_final, 2, meds, '-')
  if(rescale){
    df_final[df_final>0] = eval(parse(text=paste0("log2(",rescale_mult,"*df_final[df_final>0])")))
    df_final[df_final<0] = eval(parse(text=paste0("-log2(-",rescale_mult,"*df_final[df_final<0])")))
    meds = eval(parse(text=paste0("apply(df_final, 2, ",center,")")))
    df_final = sweep(df_final, 2, meds, '-')
  }
  return(df_final)
}

# This function simultaneously normalizes the rows and columns
# of a matrix.
rowcolNorm <- function(mat,f=1.0,e=1E-6){
  while(!(all(abs(rowSums(mat)-1)<e) & all(abs(colSums(mat)-1)<e))){
    mat = f*t(t(mat)/colSums(mat))
    mat = f*mat/rowSums(mat)
  }
  return(mat)
}

# Convert DEGAS output [0,1] to an association [-1,1]
toCorrCoeff <- function(probs){
  k=dim(probs)[2]
  if(k<2 || is.null(k)){k=2}
  l=2
  return(2*((probs-1/k)/(l-l/k) + 1/l)-1)
}

################################################################################
# Sampling functions
# Split up samples for k fold cross validation
splitKfoldCV <- function(N,k){
  if(k<3){
    stop("Please use 3 or more folds")
  }
  Idx = as.numeric(sample(1:N,N,replace=FALSE))
  sz = rep(floor(N/k),k)
  rem = N-sum(sz)
  if(rem>0){
    cntr=0
    for(i in 1:rem){
      if(cntr==k){
        cntr = 1
      }else{
        cntr = cntr+1
      }
      sz[cntr] = sz[cntr]+1
    }
  }
  cntr = 0
  grpIdx = list()
  for (i in 1:k){
    grpIdx[[i]] = Idx[(cntr+1):(cntr+sz[i])]
    cntr = cntr + sz[i]
  }
  return(grpIdx)
}

# return random sample (of s) which is evenly distributed across sample groups (g)
# where each group has n samples.
# Note: If a group has less than n samples, then all samples in that group are used.
evenSamp <- function(s,g,n){
  groups = unique(g)
  out = list()
  for (group in groups){
    if(sum(g==group)>=n){
      out[[group]] = sample(s[g==group],n,replace=FALSE)
    }else if(sum(g==group)>0){
      out[[group]] = sample(s[g==group],sum(g==group),replace=FALSE)
    }else{
      #Adding nothing
    }
  }
  out = unlist(out)
  return(out)
}

# kNN smoothing for smaller datasets. This is a dependency of knnSmoothAtlas()
knnSmooth <- function(probs,locs,k=5){
  out = probs
  dists = pairDist(locs)
  if(class(probs)[1]=="numeric"){
    N = length(probs)
    for(i in 1:N){
      idx = order(dists[i,])
      out[i] = mean(probs[idx[1:k]],na.rm=TRUE)
    }
  }else{
    N = dim(locs)[1]
    for(i in 1:N){
      idx = order(dists[i,])
      out[i,] = colMeans(probs[idx[1:k],],na.rm=TRUE)
    }
  }
  return(out)
}


# Required Functions
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

### the sliding window smoothing function
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
  
  
  while (!finished) {
    i = i + 1
    
    xstart = xmarker
    xend = xmarker + xsz
    ystart = ymarker
    yend = ymarker + ysz
    
    # Extract current patch
    sub_dat = dat[dat[, colx] > xstart & dat[, colx] < xend &
                    dat[, coly] > ystart & dat[, coly] < yend, ]
    
    # Skip empty patches
    if (nrow(sub_dat) > 0) {
      keep_tmp = sub_dat$locsx > (xstart + pad) & sub_dat$locsx < (xend - pad) &
        sub_dat$locsy > (ystart + pad) & sub_dat$locsy < (yend - pad)
      
      bc_tmp = sub_dat$barcodes
      cell_sets[[i]] = data.frame(barcode = bc_tmp, keep = keep_tmp)
    } else {
      message(sprintf("Patch %d is empty, skipping.", i))
    }
    
    # Advance the scanning window
    if ((xend + step) < xmax) {
      xmarker = xmarker + step
    } else if ((yend + step) < ymax) {
      xmarker = xmin
      ymarker = ymarker + step
    } else {
      finished = TRUE
    }
  }
  len_ss = length(cell_sets)
  message(paste0("Patch Number: ", len_ss))
  message("Completed patch annotation")
  
  # Smoothing cells in each patch
  message("Begin patch smoothing")
  q21 = floor(quantile(1:len_ss, seq(from = 0.0, to = 1.0, by = 0.05)))
  pb <- progress_bar$new(
    format = "[:bar] :percent",
    total = 21, clear = FALSE, width = 60)
  
  #remove NULL entries from cell_sets
  cell_sets <- cell_sets[!sapply(cell_sets, is.null)]
  
  #getting smoothed values
  smoothed = list()
  j = 0
  
  for (cell_set in cell_sets) {
    j = j + 1
    if (j %in% q21) pb$tick()
    Sys.sleep(2/100)
    
    # Skip if patch is empty or no TRUE in 'keep'
    if (nrow(cell_set) == 0 || sum(cell_set$keep, na.rm = TRUE) < mincellsn) {
      next
    }
    
    # Run smoothing only on 'kept' cells
    kept_indices = which(cell_set$keep)
    if (length(kept_indices) == 0) next
    
    barcodes_kept = cell_set$barcode[kept_indices]
    dat.kept = dat[which(dat$barcodes %in%barcodes_kept),]
    values = as.numeric(dat.kept[, colv])
    coords = as.matrix(dat.kept[, c("locsx", "locsy")])
    
    smoothed_tmp = knnSmooth(values, coords, kNN)
    
    smoothed[[length(smoothed) + 1]] = data.frame(
      barcode = barcodes_kept,
      score = smoothed_tmp
    )
  }
  message("Completed patch smoothing")
  
  # Combining smoothed scores
  message("Begin cell aggregation")
  # Filter out empty elements
  smoothed = smoothed[sapply(smoothed, function(x) nrow(x) > 0)]
  
  # Combine safely
  if (length(smoothed) > 0) {
    smoothed_combined = purrr::reduce(smoothed, full_join, by = "barcode")
    rownames(smoothed_combined) = smoothed_combined$barcode
    smoothed_combined$barcode = NULL
    
    smoothed_combined$haz = rowMeans(smoothed_combined, na.rm = TRUE)
    
    out = data.frame(barcodes = rownames(smoothed_combined),
                     smoothed = smoothed_combined$haz)
  } else {
    out = data.frame(barcodes = character(), smoothed = numeric())
  }
  
  message("Completed cell aggregation")
  return(out)
  message("Finished sliding window smoothing")
}




#####setup parallel backend to use many processors

library(foreach)
library(doParallel)
cell_sets = cell_sets[sapply(cell_sets,function(x) dim(x)[1]>mincellsn)
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

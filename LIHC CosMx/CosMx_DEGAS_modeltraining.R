
# r/4.4.1
# python/3.10.10

library(devtools)
library(withr)
library(rprojroot)
library(tidyverse)
library(magrittr)
library(here)
# If DEGAS hasn't been installed, do so using this command (after creating a r-packages file in your root directory)
#with_libpaths(new = '~/r-packages', install_github("tsteelejohnson91/DEGAS"))
#library(DEGAS)
library(DEGAS, lib.loc = "/N/project/degas_st/DEGAS_package")

set.seed(2)

meta_data_columns <- c("V1", "cell_id", "slide_ID_numeric", "Run_Tissue_name", "x_FOV_px", "y_FOV_px", "qcFlagsFOV", "cellType", "niche")

# What is the name of this experiment/sample?
sample_name <- "lihc_200genes"

# Patient data
## RNA-seq
patDat <- data.table::fread('/N/project/degas_st/cosmyx/data/TCGA/patDat.csv', sep = ",") # Tumor and normal tissue

## Clinical outcomes
patLab <- data.table::fread('/N/project/degas_st/cosmyx/data/Reference_liver/patLab.csv', sep = ",")

# Need non NA follow up time and status
patLab <- patLab %>% filter(!is.na(OS_time))

# only retain these same patients in patDat, and ordered the same way
patDat <- patDat %>% inner_join(patLab)

# ST data (cosmyx)
meta_df <- read_csv('/N/project/degas_st/CosMx_liver_data/cosmx_liver_metadata.csv') # To get the coordinates

st_data <- data.table::fread('/N/project/degas_st/cosmyx/data/processed/liver/cosmyx_liver_data.csv')

st_data <- st_data %>%
  mutate(FOV = str_replace(str_replace(cell_id, "c_[:digit:]_", ""), "_[:digit:].*$", ""), .before = cell_id)

common_cols <- intersect(meta_df %>% colnames, st_data %>% colnames)

st_data <- 
  left_join(
    st_data,
    meta_df %>% 
      select(
        all_of(common_cols), # the columns for joining data properly
        x_slide_mm, y_slide_mm)) # the actual real coordinates of the points to use for plotting

st_data <- st_data %>%
  relocate(contains('slide_mm'), .before = FOV)

# Find intersecting genes for feature selection
intersecting_genes <-
  intersect(colnames(st_data), colnames(patDat %>% select(-OS_status, -OS_time)))



# Feature selection (top 200 most variable genes)
genes <- patDat %>% select(all_of(intersecting_genes)) %>% apply(2, var) %>% sort(decreasing=TRUE) %>% names() %>% head(200)


# Run the DEGAS model training for given  scSRT gene-expression data, patient label, and patient gene-expresiion

runDEGASBlankCox <- function(stDat, patDat, patLab, genes, iter) {
  
  path.data = '/N/project/degas_st/cosmyx/data/DEGAS/' #give the path where the data are stored
  path.result = '/N/project/degas_st/cosmyx/lihc_output_Nov14/' #this is your output directory
  initDEGAS()
  DEGAS.toolsPath <<- '/N/project/degas_st/DEGAS_package/DEGAS/tools/' #this is where DEGAS tools and required functions can be called from.  
  DEGAS.pyloc <<- '/N/soft/rhel8/deeplearning/Python-3.10.10/bin/python3.10' # overwrite the DEGAS path with your python location
  tmpDir = paste0(path.result, 'tmp/') #this is where Activation functions, biaseses etc. will be stored after model training
  

  set_seed_term(2)
  ccoxModel_PRAD =
    runCCMTLBag(
      scExp = stDat %>% select(all_of(genes)),
      scLab = NULL, #stLab %>% pull(cluster) %>% toOneHot(),
      patExp = patDat %>% select(all_of(genes)),
      patLab = patDat %>% select(OS_time, OS_status) %>% apply(2, as.numeric),
      tmpDir = tmpDir,
      'BlankCox','DenseNet',3,5)

  saveRDS(ccoxModel_PRAD, file = paste0('/N/project/degas_st/cosmyx/lihc_output_Nov14/',iter, "_",  Sys.Date(), '.RDS'))

  model_output <-
    tibble(
      as_tibble(predClassBag(ccoxModel_PRAD,
                             stDat %>% select(all_of(genes)),
                             'pat')),
      cell_id = stDat$cell_id,
      fov = stDat$FOV,
      coord1 = stDat$x_slide_mm,
      coord2 = stDat$y_slide_mm,
      sample = stDat$slide_ID_numeric,
      cellType = stDat$cellType,
      niche = stDat$niche)

  colnames(model_output)[1] <- "Hazard"

  model_output %<>% mutate(Haz_scaled = (Hazard - min(Hazard)) / max(Hazard), .before = 'Hazard')

  return(model_output)
}


i=length(genes)

  output <-
    runDEGASBlankCox(stDat = st_data,
                     patDat = patDat,
                     patLab = patLab,
                     genes = genes[1:i],
                     iter = i)
  iter = i

  med <- output %>% pull(Hazard) %>% median

  write_csv(output, file = paste0('/N/project/degas_st/cosmyx/lihc_output_Nov14/',iter, "_",  Sys.Date(), '.csv'))



# ########## new DEGAS----
# if there are too many cells and system runs out of memory, use the subsampled version of the 
# model training function

# path.result = '/N/project/degas_st/cosmyx/lihc_output_Nov14/'
# tmpDir = paste0(path.result, 'tmp/')

# initDEGAS <- function(){
#   #DEGAS.pyloc <<- "python3"
#   #DEGAS.toolsPath <<- paste0(.libPaths()[1],"/DEGAS/tools/")
#   DEGAS.pyloc <<- '/N/soft/rhel8/deeplearning/Python-3.10.10/bin/python3.10'
#  DEGAS.toolsPath <<- '/N/project/degas_st/DEGAS_package/DEGAS/tools/'
#   # overwrite the DEGAS path with your python location (module load deeplearning/2.9.1)
#   DEGAS.train_steps <<- 2000
#   DEGAS.scbatch_sz <<- 200
#   DEGAS.patbatch_sz <<- 50
#   DEGAS.hidden_feats <<- 50
#   DEGAS.do_prc <<- 0.5
#   DEGAS.lambda1 <<- 3.0
#   DEGAS.lambda2 <<- 3.0
#   DEGAS.lambda3 <<- 3.0
#   DEGAS.seed <<- "NULL"
# }
# # Training DEGAS model

#   model1 = runDEGASatlas(stDat,NULL,patDat,patLab,tmpDir,"BlankCox","DenseNet",3,5,20000,123)
#   saveRDS(model1,file="/N/project/degas_st/cosmyx/lihc_output_Nov14/runDEGASatlas_model1_Nov14.rds")
#   preds = predClassBag(model1,stDat,"pat")
#   saveRDS(preds,file="/N/project/degas_st/cosmyx/lihc_output_Nov14/runDEGASatlas_preds1_Nov14.rds")


##### Do this after you run the bulkRNAseq data cleanup code and save the patDat.csv and patLab.csv files

#r/4.4.1
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
library(DEGAS, lib.loc = "~/DEGAS_package") #this is the folder that has DEGAS installed

set.seed(2)

meta_data_columns <- c("V1", "cell_id", "slide_ID_numeric", "Run_Tissue_name", "x_FOV_px", "y_FOV_px", "qcFlagsFOV", "cellType", "niche")

# What is the name of this experiment/sample?
sample_name <- "lihc_200genes"

# Patient data
## RNA-seq
patDat <- data.table::fread('~/patDat.csv', sep = ",") # Tumor and normal tissue

## Clinical outcomes
patLab <- data.table::fread('~/patLab.csv', sep = ",")

# Need non NA follow up time and status
patLab <- patLab %>% filter(!is.na(OS_time))

# only retain these same patients in patDat, and ordered the same way
patDat <- patDat %>% inner_join(patLab)

################## # ST data (CosMx)########
#this is the metadata associated with the CosMx data obtained from https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/human-liver-rna-ffpe-dataset/###
meta_df <- read_csv('~/cosmx_liver_metadata.csv') # To get the coordinates

st_data <- data.table::fread('~/cosmyx_liver_countdata.csv')

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


# Initialize function to run DEGAS iteratively across increasing number of feature sets.
runDEGASBlankCox <- function(stDat, patDat, patLab, genes, iter) {
  
  path.data = '/N/project/degas_st/cosmyx/data/DEGAS/'
  path.result = '/N/project/degas_st/cosmyx/lihc_output_Nov14/'
  initDEGAS()
  DEGAS.toolsPath <<- '/N/project/degas_st/DEGAS_package/DEGAS/tools/'
  # overwrite the DEGAS path with your python location
  DEGAS.pyloc <<- '/N/soft/rhel8/deeplearning/Python-3.10.10/bin/python3.10'
  tmpDir = paste0(path.result, 'tmp/')
 
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

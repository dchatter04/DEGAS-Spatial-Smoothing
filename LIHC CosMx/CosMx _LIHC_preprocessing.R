library(tidyverse)

st_data <- data.table::fread('/N/project/degas_st/cosmyx/cosmx_liver_metadata.csv')

meta_data_columns <- c("V1", "cell_id", "slide_ID_numeric", "Run_Tissue_name", "x_FOV_px", "y_FOV_px", "qcFlagsFOV", "cellType", "niche")

meta_data_df <- st_data %>% select(all_of(meta_data_columns))
rna_df <- st_data %>% select(-all_of(meta_data_columns))

st_data_out <- cbind(
        meta_data_df,
        preprocessCounts(t(rna_df) %>% apply(2, as.numeric))    
)
colnames(st_data_out)[(length(meta_data_columns)+1): length(st_data_out)] <- colnames(rna_df)

data.table::fwrite(st_data_out, '/N/project/degas_st/cosmyx/data/processed/liver/cosmyx_liver_countdata.csv')

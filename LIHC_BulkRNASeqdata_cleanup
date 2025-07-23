# rm(list = ls(all.names = TRUE))

library(tidyverse)
library(DEGAS)
####### ******* the LIHC reference data was downloaded from http://firebrowse.org/?cohort=LIHC&download_dialog=true *********************

patDat <- read.delim('~/gdac.broadinstitute.org_LIHC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/LIHC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt', header = FALSE)

patDat <- patDat %>% t() %>% as_tibble()

colnames(patDat) <- slice(patDat, 1)

patDat <- slice(patDat, -1)

patDat <- patDat %>% select(-c(`SLC35E2|728661`,`SLC35E2|9906`))

colnames(patDat) <- gsub(pattern = "\\|.*",replacement = "", colnames(patDat))

cols_to_delete <- stringr::str_detect(colnames(patDat), pattern = "\\?") %>% which() #%>% sapply(isTRUE)

patDat <- patDat %>% select(-cols_to_delete)

patDat <- patDat %>% mutate(across(.cols = 3:ncol(patDat),.fns = as.numeric))

patDat <- cbind(patDat[1:2], t(preprocessCounts(patDat[-c(1,2)]))) %>% as_tibble()

patDat_tumor <- patDat %>% filter(!str_detect(`Hybridization REF`, "11A")) # 11A means normal, 01A means tumor.

# Here is subsetting of the file name. I want to make one object with only cancer tissues
patDat <- patDat %>% mutate(`Hybridization REF` = str_replace_all(toupper(`Hybridization REF`),  pattern = "-[0-9][0-9]A.*", replacement = "")) %>% select(-gene_id)

patDat_tumor <- patDat_tumor %>% mutate(`Hybridization REF` = str_replace_all(toupper(`Hybridization REF`),  pattern = "-[0-9][0-9]A.*", replacement = "")) %>% select(-gene_id)

colnames(patDat)[1] <- "PATIENT_ID"

colnames(patDat_tumor)[1] <- "PATIENT_ID"


patDat <- patDat %>% mutate(PATIENT_ID = tolower(PATIENT_ID))


patDat_tumor <- patDat_tumor %>% mutate(PATIENT_ID = tolower(PATIENT_ID))

# Here's the RNA-seq gene expression file

data.table::fwrite(patDat, file = '~/patDat.csv', sep = ",")


# This is the clinical data with survival information
clin_data <- read.delim('~/gdac.broadinstitute.org_LIHC.Clinical_Pick_Tier1.Level_4.2016012800.0.0/LIHC.clin.merged.picked.txt', header = FALSE) %>% t() %>% as.data.frame()

colnames(clin_data) <- slice(clin_data, 1)
clin_data <- clin_data[-1,]
clin_data <- clin_data %>% rename(PATIENT_ID = `Hybridization REF`)

clin_data %>% select(PATIENT_ID, days_to_last_known_alive, days_to_last_followup, days_to_death, vital_status) %>% view()

clin_data %>%
  mutate(OS_status = vital_status,
         OS_time = if_else(is.na(days_to_last_followup),
                           days_to_death,
                           days_to_last_followup)) %>%
  select(PATIENT_ID, OS_time, OS_status) %>%
  data.table::fwrite('~/patLab.csv')

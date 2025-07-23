library(tidyverse)
library(DEGAS)

# Read in data
patDat <- data.table::fread('~/T2D_bulkRNAseqdat.csv')

# Fix column names
patnames <- colnames(patDat)[-1]

patLab <- tibble(patnames, status = ifelse(str_detect(patnames, 'ND'), 'no_disease', 't2d'),
    t2d = ifelse(str_detect(patnames, 'T2D'), 1L, 0L),
    no_disease = ifelse(str_detect(patnames, 'ND'), 1L, 0L))

# Output the cleaned dataset (Clinical data)
patLab %>%
    select(patnames, t2d, no_disease) %>%
    write_csv('~/T2D_patLab.csv')


# Output the cleaned dataset (Transcriptomic data)
gene_names <- patDat[,1] %>% pull()

patDat <- patDat %>% select(-genes) %>%
    apply(2, as.numeric) %>%
    t() %>% as.data.frame() %>% rownames_to_column('patient') %>% tibble()

colnames(patDat) <- append(c('patient'), gene_names)

patDat <- cbind(patDat %>% select(patient), t(preprocessCounts(patDat[-1])))

patDat %>% data.table::fwrite('~/T2D_patDat.csv')

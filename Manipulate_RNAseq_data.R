# script to manipulate RNAseq data August 1
# source: 
# setwd("/Users/mokira/TOOLS/Microarray_RMA")

# load the packages 
install.packages("dplyr")
library(dplyr)
library(tidyverse)
library(GEOquery)

# read the data
dat <- read.csv(file = "GSE183947_fpkm.csv")

# get metadata
gse <- getGEO(GEO = 'GSE183947', GSEMatrix = TRUE)
gse

# get the phenotype data
metadata <- pData(phenoData(gse[[1]]))

# subset the coloumns of metadata required
# first parameter is the data and second is the coloumn number required

metadata.modified <- metadata %>%
  select(1,10,11,17) %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue:", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis:", "", metastasis))

# reshaping the data
# in gather function key is the class of colnames and value are the values
# and -gene represent, not to touch that particular coloumn

dat.long <- dat %>%
  rename(gene = X) %>%
  gather(key = 'samples', value = 'FPKM', -gene)

# join the dataframes  = dat.long + metadata.modified

dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("samples" = "description"))

# explore the data
# compare the expression of BRCA1 and BRCA2 gene in tumour vs normal sample

dat.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  group_by(gene, tissue) %>%
  summarise(mean_FPKM = mean(FPKM),
            median_FPKM = median(FPKM)) %>%
  arrange(mean_FPKM)

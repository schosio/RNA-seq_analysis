# script to manipulate RNAseq data
# source of data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154169
# setwd("/Users/mokira/MPLAB/RNA-seq_analysis/GSE154169/")

# load the packages 
library(dplyr)
library(tidyverse)
library(GEOquery)

# read the data
dat <- read.csv(file = "GSE154169_normalised_CPM.txt", sep = "\t")

# get metadata
gse <- getGEO(GEO = 'GSE154169', GSEMatrix = TRUE)
gse

# get the phenotype data
metadata <- pData(phenoData(gse[[1]]))

# subset the coloumns of metadata required
metadata.modified <- metadata %>%
  rename(genotype = characteristics_ch1.1) %>%
  rename(treatment = characteristics_ch1.2) %>%
  select(1,11,12) %>%
  mutate(genotype = gsub("genotype/variation:", "", genotype)) %>%
  mutate(treatment = gsub("treatment:", "", treatment)) %>%
  mutate(title = gsub(".*rep1", "", title)) %>%
  mutate(title = gsub(".*rep2", "", title)) %>%
  mutate(title = gsub(".*rep3", "", title)) %>%
  mutate(title = gsub(" [", "",title, fixed = T)) %>%
  mutate(title = gsub("]", "", title, fixed = T))
  
# reshaping the data
# in gather function key is the class of colnames and value are the values
# and -gene represent, not to touch that particular coloumn

dat.long <- dat %>%
  rename(gene = Gene) %>%
  gather(key = 'samples', value = 'CPM', -gene)

# join the dataframes  = dat.long + metadata.modified

dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("samples" = "title"))


# explore the data

dat.long %>%
  filter( gene == 'aac') %>%
  group_by(gene, genotype) %>%
  summarise(mean_CPM = mean(CPM),
            median_CPM = median(CPM))


# Plotting and visualization

library(ggplot2)

# 1. barplot

dat.long %>%
  filter(gene == 'whiB3') %>%
  ggplot(., aes(x = samples, y = CPM, fill = genotype, shape = treatment)) +
  geom_col()

#2. Density
dat.long %>%
  filter(gene == 'whiB3') %>%
  ggplot(., aes(x = CPM, fill = genotype)) +
  geom_density(alpha = 0.3)

#3. heatmap
genes.of.interest <- c('whiB3', 'whiB1', 'whiB2', 'whiB4', 'whiB5', 'whiB6',
                       'whib7')

dat.long %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = samples, y = gene, fill = CPM)) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red')


# 4. Boxplot
dat.long %>%
  filter(gene == 'whiB3') %>%
  ggplot(., aes(x = genotype, y = CPM)) +
  geom_boxplot()
  #violin plot
  #geom_violin()

# 5. Line graph

dat.long %>%
  filter(gene == 'whiB3') %>%
  ggplot(., aes(x = samples, y = CPM, colour = genotype, shape = treatment,
                group = interaction(genotype, treatment))) +
  geom_point(size = 3) +
  # connecting the points; group part in aes is neccessary for this
  geom_line() +
  # Changing the labels
  labs(title = "whiB3 expression", colour = "Genotype", shape = "Treatment",
       x = "Samples") +
  # Removal of x lables and ticks
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
  

# 5. Line graph

dat.long %>%
  filter(gene == 'ansA') %>%
  ggplot(., aes(x = samples, y = CPM, colour = genotype, shape = treatment,
                group = interaction(genotype, treatment))) +
  geom_point(size = 3) +
  # connecting the points; group part in aes is neccessary for this
  geom_line() +
  # Changing the labels
  labs(title = "whiB6 expression", colour = "Genotype", shape = "Treatment",
       x = "Samples") +
  # Removal of x lables and ticks
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

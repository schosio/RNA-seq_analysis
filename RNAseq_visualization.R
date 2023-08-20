# script to visualize RNAseq data
# source: https://www.youtube.com/watch?v=RukuTtiY4Sg 
# Barplot, density plot, boxplot, scatterplot, heatmap
# setwd("/Users/mokira/PC/")

# load the libraries
library(tidyverse)
library(ggplot2)

#data

# 1. barplot

dat.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = samples, y = FPKM, fill = tissue)) +
  geom_col()

#2. Density
dat.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = FPKM, fill = tissue)) +
  geom_density(alpha = 0.3)

# 3. Boxplot
dat.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = metastasis, y = FPKM)) +
  #geom_boxplot()
  # violin plot
  geom_violin()

# 4. Scatterplot
dat.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  # to plot one gene on x and other on y, data needs to be transformed
  spread(key = gene, value = FPKM) %>%
  ggplot(., aes(x = BRCA1, y = BRCA2)) +
  geom_point()

# 5. Heatmap
genes.of.interest <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')

pdf("heatmap_save2.pdf", width = 10, height = 8)
dat.long %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = samples, y = gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red')

dev.off()
# save the plot

ggsave(p, filename = 'heatmap_save1.png', height = 8, width = 10 )
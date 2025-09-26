################################################################################
# Script name: figure_4.R
# Purpose: Generate Figure 4 – visualization of mesohabitat contributions 
#                             on the SOM grid
# Description: This script visualizes the influence of mesohabitats on 
#              self-organizing maps (SOMs), showing:
#              (a) functional feeding groups,
#              (b) trophic groups,
#              (c) habit groups.
# Input: data/som_cluster_data.csv
################################################################################

# The required packages to run the following script-----------------------------
library(kohonen) # for the kohonen SOM
library(aweSOM)
library(SOMbrero) # this script used the SOMbrero package
library(GGally)
library(ggpubr) # arranging the plots made by ggplot2
library(cowplot) # arranging the plots made under SOMbrero
library(indicspecies) # for indicator taxa
library(psych) # for the whole summary including the mean, median and Q1 and Q2
library(factoextra)
library(cluster)
library(tidyverse)
library(rstatix)


# Importing data for analysis ------------------------------------------------
som_clus_data <- read_csv("data/som_cluster_data.csv")


#########################################################################
####### creating stack bar plot for all the cluster with respect to RPBS
########################################################################


# Stacked + percent
pretty_palette2 <- c("#10E205","#E7B800","#FC4E07")

stack_ffg <- ggplot(data = som_clus_data,
                    aes(x = cluster_ffg,  fill = habitat)) +
  geom_bar(position = "fill", alpha = 0.6)+
  ggplot2::scale_fill_manual(values = pretty_palette2)+
  theme_classic()+
  theme(legend.position="top", axis.title.x=element_blank())

stack_tg <- ggplot(data = som_clus_data,
                   aes(x = cluster_tg,  fill = habitat)) +
  geom_bar(position = "fill", alpha = 0.6)+
  ggplot2::scale_fill_manual(values = pretty_palette2)+
  theme_classic()+
  theme(legend.position="top", axis.title.x=element_blank())

stack_hg <- ggplot(data = som_clus_data,
                   aes(x = cluster_hg,  fill = habitat)) +
  geom_bar(position = "fill", alpha = 0.6)+
  ggplot2::scale_fill_manual(values = pretty_palette2)+
  theme_classic()+
  theme(legend.position="top", axis.title.x=element_blank())

ggarrange(stack_ffg,stack_tg,stack_hg,
          ncol = 3, common.legend = TRUE,
          legend="bottom", labels="auto")
ggsave("Fig_add_var.pdf", width = 340, height = 180, units = "mm", dpi = 300)
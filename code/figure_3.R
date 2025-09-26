################################################################################
# Script name: figure_3.R
# Purpose: Generate Figure 3 – distribution patterns of functional feeding,
#          trophic, and habit groups using self-organizing maps (SOMs)
# Description: This script visualizes macroinvertebrate community patterns
#              with self-organizing maps (SOMs). It shows:
#              (a–e) functional feeding groups,
#              (f–j) trophic groups,
#              (k–o) habit groups,
#              based on estimated abundance values.
# Input: data/som.csv
#        data/cluster_env_data_som.csv
#
# This script uses all mesohabitat data (rapid, pool, bench).
# Before running, rows and columns containing only zeros or NA values 
# were removed from the dataset.
# All variables have been arranged in columns for consistency.
################################################################################

################################################################################
#### The main analysis starts here #############################################
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


# Importing data for SOM analysis (This is the master file)----Step-1-----------
som.data <- read.csv("data/som.csv", stringsAsFactors = TRUE)

# Creating summary of the data for the paper------------------------------------
# Importing the data after running the SOM cluster
# this step should run after the SOM
cluster_som <- read.csv("data/cluster_env_data_som.csv")
table_env <- describeBy(cluster_som, 
                        group=cluster_som$cluster, quant =c(.25,.75))
table_env
#-------------------------------------------------------------------------------

# Checking for the autocorrelation among env and bio variables####Step-2########
som_ffg_corr <- som.data[,c(116:123)]
som_tg_corr <- som.data[,c(124:135)]
som_hg_corr <- som.data[,c(136:140)]
som_div_corr <- som.data[,c(141:145)]
som_bio_cor <- som.data[,c(22:115)]

som.env_cor <- som.data[,c(5:21)]

# plotting the results of the correlations########################Step3#########

# plot of the species variables autocorrelation with spearman method
corr_1 <- ggcorr(som_ffg_corr, method = c("pairwise", "spearman"),
                 nbreaks = 6,
                 label = TRUE,
                 label_size = 4,
                 color = "grey50",
                 size = 4)
corr_2 <- ggcorr(som_tg_corr, method = c("pairwise", "spearman"),
                 nbreaks = 6,
                 label = TRUE,
                 label_size = 4,
                 color = "grey50",
                 size = 4)
corr_3 <- ggcorr(som_hg_corr, method = c("pairwise", "spearman"),
                 nbreaks = 6,
                 label = TRUE,
                 label_size = 4,
                 color = "grey50",
                 size = 4)
corr_4 <- ggcorr(som_div_corr, method = c("pairwise", "spearman"),
                 nbreaks = 6,
                 label = TRUE,
                 label_size = 4,
                 color = "grey50",
                 size = 4)
# plot of the soil environment variables autocorrelation########################
# with pearson method
corr_env <- ggcorr(som.env_cor, method = c("pairwise", "pearson"),
                   nbreaks = 6,
                   label = TRUE,
                   label_size = 2,
                   color = "grey50",
                   size = 4)



# Final filtering and making ready data for the SOM analysis##########Step-4####
som.unique <- som.data[,c(116:123)]
som.unique1 <- som.data[,c(124:135)]
som.unique2 <- som.data[,c(136:140)]
som.unique3 <- som.data[,c(141:145)]
som.unique4 <- som.data[,c(22:32,34:36,38:59,
                           61:69,71:90,92:108,
                           112:114,115)]

library(rcompanion)
plotNormalHistogram(som.unique, prob = FALSE,
                    main = "Normal Distribution overlay on Histogram",
                    length = 1000 )
plotNormalHistogram(som.unique1, prob = FALSE,
                    main = "Normal Distribution overlay on Histogram",
                    length = 1000 )
plotNormalHistogram(som.unique2, prob = FALSE,
                    main = "Normal Distribution overlay on Histogram",
                    length = 1000 )
plotNormalHistogram(som.unique3, prob = FALSE,
                    main = "Normal Distribution overlay on Histogram",
                    length = 1000 )



# Log (x + 1) Transormation on the data#########################Step-5##########
som.trans <- as.matrix(log(som.unique+1))
som.trans1 <- as.matrix(log(som.unique1+1))
som.trans2 <- as.matrix(log(som.unique2+1))
som.trans3 <- as.matrix(log(som.unique3+1))
som.trans4 <- as.matrix(log(som.unique4+1))
plotNormalHistogram(som.trans3, prob = FALSE,
                    main = "Normal Distribution overlay on Histogram",
                    length = 1000 )

# SOM training data - FFG's ###################################Step-6###########
set.seed(1105)
# Run the SOM algorithm with 10 intermediate backups and 500 iterations
my.som <- trainSOM(x.data=som.trans, type='numeric',
                   topo='hexagonal',
                   dimension=c(6,8),
                   affectation='standard', 
                   dist.type='euclidean', maxit=500,
                   scaling='unitvar',
                   init.proto='random', nb.save=10,
                   radius.type='gaussian', eps0=1)


# Quality of the SOM-----------------------------------Step-7-------------------
quality(my.som, quality.type = "all")

# The energy evolves seen byt the following plot################################
plot_e <- plot(my.som, what="energy")+
  ggplot2::ggtitle("")

# Clustering distribution
# Hitmap visualization
plot(my.som, what = "obs",
     type = "hitmap")
#+theme_transparent()+
ggplot2::scale_fill_grey(start = 0.9, end = 0.3)

# umatrix#######################################################################
plot_u <- plot(my.som, what = "prototypes",
               type = "umatrix")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  guides(fill = guide_legend(element_blank()))+
  theme_transparent()

# Distance between prototype
plot_pro <- plot(my.som, what = "prototypes",
                 type = "smooth.dist")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("")

# Building super-classes from the resulting SOM-------------Step-8--------------
plot(superClass(my.som))
my.sc <- superClass(my.som, k = 5)
plot(my.sc)
summary(my.sc) # significant with the cluster
som_clust <- plot(my.sc, plot.var = F, h = 0.4) 


plot(my.sc, show.names = F, type = "grid") +
  ggplot2::scale_fill_manual(values = pretty_palette)

plot_hit <-plot(my.sc, what = "obs",
                type = "hitmap",
                show.names = FALSE)+
  theme_bw() # Number of observation

#-Clustering interpretation-----------------------------Step-9------------------------ 
class(som.data$habitat)
levels(som.data$habitat)
plot(my.som, what = "add", type = "pie", variable = som.data$habitat) +
  scale_fill_brewer(type = "qual") + 
  guides(fill = guide_legend(title = "Mesohabitat"))

plot(my.som, what = "add", type = "names", variable = som.data$stn)

plot(my.som, what = "prototypes", type = "poly.dist")

plot(my.som, what = "add", type = "names", variable = som.data$habitat)
plot(my.som, what = "add", type = "words", variable = som.data[,5:6])
plot(my.som, what = "add", type = "words", variable = som.data[,14:21])
plot(my.som, what = "add", type = "color", variable = som.data$S_TP)
plot(my.som, what = "prototypes", type = "color", variable = 1)

###--Ploting------------------------
# SOM Grid ### This plot should be done after cluster finalization #############
#pretty_palette<- c('#7fc97f', '#f0027f', '#fdc086', "#ffff99", "#386cb0") # , '', '#9467bd', ''
#pretty_palette<- c('#66c2a5', '#fc8d62', '#a6d854', "#e78ac3", "#8da0cb")
pretty_palette<- c('#F8766d', '#7cae00', '#d39200', "#00bfc4", "#c77cff")
##################################################################
# Labeling all the cells in the SOM ##############################
##################################################################
label_1 <- "I"
label_2 <- "II"
label_3 <- "III"
label_4 <- "IV"
label_5 <- "V"

##########################################################################
som_grid <- plot(my.sc, type = "grid", show.names=F) + 
  ggplot2::ggtitle("SOM")+ 
  #ggplot2::scale_fill_grey(start = 0.9, end = 0.3) +
  ggplot2::scale_fill_manual(values = pretty_palette) +
  ggplot2::scale_alpha_manual(values=c(0.01,0.01,0.01,0.01,0.01,0.01))+
  theme_transparent()+
  annotate ("text",size= 16,fontface = "bold",
            x = 2,y = 1.8,label = label_1)+
  annotate ("text",size= 16, fontface = "bold",
            x = 1.5, y = 6.1,label = label_2)+
  annotate ("text",size= 16, fontface = "bold",
            x = 5.5, y = 6.1,label = label_5)+
  annotate ("text",size= 16,fontface = "bold",
            x = 4.0, y = 3.5,label = label_4)+
  annotate ("text",size= 16,fontface = "bold",
            x = 5.0, y = 1.7,label = label_3)+
  theme(legend.position="none")

### ploting ffg's on the prototype

FC <- plot(my.sc,show.names = FALSE, what = "prototypes",
           type = "color", variable = "FC")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey10")+
  ggplot2::ggtitle("Filter-collectors", subtitle="(22.56***)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        legend.title=element_blank())+
  theme(legend.key.height= unit(1.2, 'cm'),
        legend.key.width= unit(0.4, 'cm'))

GC <- plot(my.sc,show.names = FALSE, what = "prototypes",
           type = "color", variable = "GC")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey10")+
  ggplot2::ggtitle("Gatherer-collectors", subtitle="(19.17***)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        legend.title=element_blank())+
  theme(legend.key.height= unit(1.2, 'cm'),
        legend.key.width= unit(0.4, 'cm'))

P <- plot(my.sc,show.names = FALSE, what = "prototypes",
          type = "color", variable = "P")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey10")+
  ggplot2::ggtitle("Predators", subtitle="(14.69***)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        legend.title=element_blank())+
  theme(legend.key.height= unit(1.2, 'cm'),
        legend.key.width= unit(0.4, 'cm'))

SC <- plot(my.sc,show.names = FALSE, what = "prototypes",
           type = "color", variable = "SC")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey10")+
  ggplot2::ggtitle("Scrapers", subtitle="(33.89***)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        legend.title=element_blank())+
  theme(legend.key.height= unit(1.2, 'cm'),
        legend.key.width= unit(0.4, 'cm'))

MF <- plot(my.sc,show.names = FALSE, what = "prototypes",
           type = "color", variable = "MF")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey10")+
  ggplot2::ggtitle("Multi-FFG's", subtitle="(22.07***)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        legend.title=element_blank())+
  theme(legend.key.height= unit(1.2, 'cm'),
        legend.key.width= unit(0.4, 'cm'))

################################################################################
###### Habit groups ##########################################################
################################################################################
# SOM training data - HG's ###################################Step-6###########
set.seed(1225)
# Run the SOM algorithm with 10 intermediate backups and 500 iterations
my.som1 <- trainSOM(x.data=som.trans1, type='numeric',
                    topo='hexagonal',
                    dimension=c(6,8),
                    affectation='standard', 
                    dist.type='euclidean', maxit=500,
                    scaling='unitvar',
                    init.proto='random', nb.save=10,
                    radius.type='gaussian', eps0=1)


# Quality of the SOM-----------------------------------Step-7-------------------
quality(my.som1, quality.type = "all")

# The energy evolves seen byt the following plot################################
plot_e1 <- plot(my.som1, what="energy")+
  ggplot2::ggtitle("")

# Clustering distribution
# Hitmap visualization
plot(my.som1, what = "obs",
     type = "hitmap")
#+theme_transparent()+
ggplot2::scale_fill_grey(start = 0.9, end = 0.3)

# umatrix#######################################################################
plot_u1 <- plot(my.som1, what = "prototypes",
                type = "umatrix")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  guides(fill = guide_legend(element_blank()))+
  theme_transparent()

# Distance between prototype
plot_pro1 <- plot(my.som1, what = "prototypes",
                  type = "smooth.dist")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("")

# Building super-classes from the resulting SOM-------------Step-8--------------
plot(superClass(my.som1))
my.sc1 <- superClass(my.som1, k = 4)
plot(my.sc1)
summary(my.sc1) # significant with the cluster
som_clust1<- plot(my.sc1, plot.var = F, h = 0.4) # cluster
plot(my.sc1, show.names = F, type = "grid") 

plot_hit1 <- plot(my.sc1, what = "obs",
                  type = "hitmap",
                  show.names = FALSE)+
  theme_bw() # Number of observation

#-Clustering interpretation-----------------------------Step-9------------------------ 
class(som.data$habitat)
levels(som.data$habitat)
plot(my.som1, what = "add", type = "pie", variable = som.data$habitat) +
  scale_fill_brewer(type = "qual") + 
  guides(fill = guide_legend(title = "Mesohabitat"))

plot(my.som1, what = "add", type = "names", variable = som.data$stn)

plot(my.som1, what = "prototypes", type = "poly.dist")

plot(my.som1, what = "add", type = "names", variable = som.data$habitat)
plot(my.som1, what = "add", type = "words", variable = som.data[,5:6])
plot(my.som1, what = "add", type = "color", variable = som.data$S_TP)
plot(my.som1, what = "prototypes", type = "color", variable = 1)

###--Ploting------------------------
# SOM Grid ### This plot should be done after cluster finalization #############
pretty_palette<- c('#7fc97f', '#f0027f', '#fdc086', "#ffff99", "#386cb0") # , '', '#9467bd', ''
pretty_palette<- c('#66c2a5', '#fc8d62', '#a6d854', "#e78ac3", "#8da0cb")
##################################################################
# Labeling all the cells in the SOM ##############################
##################################################################
label_1 <- "I"
label_2 <- "II"
label_3 <- "III"
label_4 <- "IV"


##########################################################################
som_grid1<- plot(my.sc1, type = "grid", show.names=F) + 
  ggplot2::ggtitle("SOM")+ 
  #ggplot2::scale_fill_grey(start = 0.9, end = 0.3) +
  ggplot2::scale_fill_manual(values = pretty_palette) +
  ggplot2::scale_alpha_manual(values=c(0.01,0.01,0.01,0.01,0.01,0.01))+
  theme_transparent()+
  annotate ("text",size= 16,fontface = "bold",
            x = 2,y = 1.8,label = label_1)+
  annotate ("text",size= 16, fontface = "bold",
            x = 1.5, y = 6.1,label = label_2)+
  annotate ("text",size= 16,fontface = "bold",
            x = 5.5, y = 6.1,label = label_3)+
  annotate ("text",size= 16,fontface = "bold",
            x = 6.0, y = 3.5,label = label_4)+
  theme(legend.position="none")

### ploting ffg's on the prototype

CB <- plot(my.sc1,show.names = FALSE, what = "prototypes",
           type = "color", variable = "CB")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey10")+
  ggplot2::ggtitle("Climbers", subtitle="(22.56***)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        legend.title=element_blank())+
  theme(legend.key.height= unit(1.2, 'cm'),
        legend.key.width= unit(0.4, 'cm'))

CN <- plot(my.sc1,show.names = FALSE, what = "prototypes",
           type = "color", variable = "CN")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey10")+
  ggplot2::ggtitle("Clingers", subtitle="(36.64)***")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        legend.title=element_blank())+
  theme(legend.key.height= unit(1.2, 'cm'),
        legend.key.width= unit(0.4, 'cm'))

SP <- plot(my.sc1,show.names = FALSE, what = "prototypes",
           type = "color", variable = "SP")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey10")+
  ggplot2::ggtitle("Sprawlers", subtitle="(43.22)***")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        legend.title=element_blank())+
  theme(legend.key.height= unit(1.2, 'cm'),
        legend.key.width= unit(0.4, 'cm'))

CB.CN <- plot(my.sc1,show.names = FALSE, what = "prototypes",
              type = "color", variable = "CB.CN")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey10")+
  ggplot2::ggtitle("Climber-Clingers", subtitle="(3.51)*")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        legend.title=element_blank())+
  theme(legend.key.height= unit(1.2, 'cm'),
        legend.key.width= unit(0.4, 'cm'))

SW.CB <- plot(my.sc1,show.names = FALSE, what = "prototypes",
              type = "color", variable = "SW.CB")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey10")+
  ggplot2::ggtitle("Swimmers-Climbers", subtitle="(2.89)*")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        legend.title=element_blank())+
  theme(legend.key.height= unit(1.2, 'cm'),
        legend.key.width= unit(0.4, 'cm'))

SW.CN <- plot(my.sc1,show.names = FALSE, what = "prototypes",
              type = "color", variable = "SW.CN")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey10")+
  ggplot2::ggtitle("Swimmers-Clingers", subtitle="(21.42)***")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        legend.title=element_blank())+
  theme(legend.key.height= unit(1.2, 'cm'),
        legend.key.width= unit(0.4, 'cm'))

SW.CW <- plot(my.sc1,show.names = FALSE, what = "prototypes",
              type = "color", variable = "SW.CW")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey10")+
  ggplot2::ggtitle("Swimmers-Crawlers", subtitle="(6.14)***")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        legend.title=element_blank())+
  theme(legend.key.height= unit(1.2, 'cm'),
        legend.key.width= unit(0.4, 'cm'))

MH <- plot(my.sc1,show.names = FALSE, what = "prototypes",
           type = "color", variable = "MH")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey10")+
  ggplot2::ggtitle("Multi-habits", subtitle="(20.61)***")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        legend.title=element_blank())+
  theme(legend.key.height= unit(1.2, 'cm'),
        legend.key.width= unit(0.4, 'cm'))
################################################################################
###### Trophic groups ##########################################################
################################################################################
# SOM training data - TG's ###################################Step-6###########
set.seed(1335)
# Run the SOM algorithm with 10 intermediate backups and 500 iterations
my.som2 <- trainSOM(x.data=som.trans2, type='numeric',
                    topo='hexagonal',
                    dimension=c(6,8),
                    affectation='standard', 
                    dist.type='euclidean', maxit=520,
                    scaling='unitvar',
                    init.proto='random', nb.save=10,
                    radius.type='gaussian', eps0=1)


# Quality of the SOM-----------------------------------Step-7-------------------
quality(my.som2, quality.type = "all")

# The energy evolves seen byt the following plot################################
plot_e2 <- plot(my.som2, what="energy")+
  ggplot2::ggtitle("")

# Clustering distribution
# Hitmap visualization
plot(my.som2, what = "obs",
     type = "hitmap")
#+theme_transparent()+
ggplot2::scale_fill_grey(start = 0.9, end = 0.3)

# umatrix#######################################################################
plot_u2 <- plot(my.som2, what = "prototypes",
                type = "umatrix")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  guides(fill = guide_legend(element_blank()))+
  theme_transparent()

# Distance between prototype
plot_pro2 <- plot(my.som2, what = "prototypes",
                  type = "smooth.dist")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("")

# Building super-classes from the resulting SOM-------------Step-8--------------
plot(superClass(my.som2))
my.sc2 <- superClass(my.som2, k = 5)
plot(my.sc2)
summary(my.sc2) # significant with the cluster
som_clust2 <- plot(my.sc2, plot.var = F, h = 0.4) # cluster
plot(my.sc2, show.names = F, type = "grid") 

plot_hit2<- plot(my.sc2, what = "obs",
                 type = "hitmap",
                 show.names = FALSE)+
  theme_bw() # Number of observation



###--Ploting------------------------
# SOM Grid ### This plot should be done after cluster finalization #############
pretty_palette<- c('#7fc97f', '#f0027f', '#fdc086', "#ffff99", "#386cb0") # , '', '#9467bd', ''
pretty_palette<- c('#66c2a5', '#fc8d62', '#a6d854', "#e78ac3", "#8da0cb")
##################################################################
# Labeling all the cells in the SOM ##############################
##################################################################
label_1 <- "I"
label_2 <- "II"
label_3 <- "III"
label_4 <- "IV"
label_5 <- "V"

##########################################################################
som_grid2 <-plot(my.sc2, type = "grid", show.names=F) + 
  ggplot2::ggtitle("SOM")+ 
  #ggplot2::scale_fill_grey(start = 0.9, end = 0.3) +
  ggplot2::scale_fill_manual(values = pretty_palette) +
  ggplot2::scale_alpha_manual(values=c(0.01,0.01,0.01,0.01,0.01,0.01))+
  theme_transparent()+
  annotate ("text",size= 16,fontface = "bold",
            x = 2.5,y = 0.9,label = label_1)+
  annotate ("text",size= 16, fontface = "bold",
            x = 2.0, y = 3.5,label = label_2)+
  annotate ("text",size= 16,fontface = "bold",
            x = 5.5, y = 6.0,label = label_5)+
  annotate ("text",size= 16,fontface = "bold",
            x = 6.0, y = 1.7,label = label_4)+
  annotate ("text",size= 16,fontface = "bold",
            x = 2.5, y = 6.1,label = label_3)+
  theme(legend.position="none")

### ploting ffg's on the prototype

CA <- plot(my.sc2,show.names = FALSE, what = "prototypes",
           type = "color", variable = "CA")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey10")+
  ggplot2::ggtitle("Carnivores", subtitle="(11.34)***")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        legend.title=element_blank())+
  theme(legend.key.height= unit(1.2, 'cm'),
        legend.key.width= unit(0.4, 'cm'))

DE <- plot(my.sc2,show.names = FALSE, what = "prototypes",
           type = "color", variable = "DE")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey10")+
  ggplot2::ggtitle("Detritivores", subtitle="(16.44)***")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        legend.title=element_blank())+
  theme(legend.key.height= unit(1.2, 'cm'),
        legend.key.width= unit(0.4, 'cm'))

DH <- plot(my.sc2,show.names = FALSE, what = "prototypes",
           type = "color", variable = "DH")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey10")+
  ggplot2::ggtitle("Detri-Herbivores", subtitle="(6.03)***")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        legend.title=element_blank())+
  theme(legend.key.height= unit(1.2, 'cm'),
        legend.key.width= unit(0.4, 'cm'))

HE <- plot(my.sc2,show.names = FALSE, what = "prototypes",
           type = "color", variable = "HE")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey10")+
  ggplot2::ggtitle("Herbivores", subtitle="(19.17)***")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        legend.title=element_blank())+
  theme(legend.key.height= unit(1.2, 'cm'),
        legend.key.width= unit(0.4, 'cm'))

MT <- plot(my.sc2,show.names = FALSE, what = "prototypes",
           type = "color", variable = "MT")+
  ggplot2::scale_fill_continuous(low = "white", high = "grey10")+
  ggplot2::ggtitle("Multi-trophic", subtitle="(88.46)***")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        legend.title=element_blank())+
  theme(legend.key.height= unit(1.2, 'cm'),
        legend.key.width= unit(0.4, 'cm'))


################################################################################
#  Grouping and plotting #######################################################
################################################################################


plot_ffg <- ggarrange(FC, GC, P, SC, MF,
                      ncol =5, nrow = 1,
                      labels =c("a","b","c","d","e"))

comb_ffg <- annotate_figure(plot_ffg, 
                            left = text_grob("Functional feeding groups",
                                                       color = "black", 
                                             face = "bold", size = 12,
                                                       rot = 90))

plot_tg <- ggarrange(CA, DE, DH, HE, MT,
                     ncol =5, nrow = 1,
                     labels =c("f","g","h","i","j"))

comb_tg <- annotate_figure(plot_tg, 
                           left = text_grob("Trophic groups",
                                                     color = "black", 
                                            face = "bold", size = 12,
                                                     rot = 90))

plot_hg <- ggarrange(CB, CN, SP, CB.CN, MH,
                     ncol =5, nrow = 1,
                     labels =c("k","l","m","n","o"))

comb_hg <- annotate_figure(plot_hg, 
                           left = text_grob("Habit groups",
                                                     color = "black", 
                                            face = "bold", size = 12,
                                                     rot = 90))

ggarrange(comb_ffg,comb_tg,comb_hg, nrow = 3)


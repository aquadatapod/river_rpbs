################################################################################
# Script name: figure5_SOM_taxa_clusters.R
# Purpose: Generate Figure 5 – trained SOM showing individual MBI taxa 
#          represented by categorized sub-clusters
# Description: This script visualizes self-organizing maps (SOMs) of MBI taxa. 
#              Each sub-cluster is labeled with:
#              - F values derived from ANOVA
#              - Stars indicating significance levels
# Input: data/som.csv
#        data/cluster_env_data_som.csv
################################################################################

###############################################################################
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
library(rstatix)


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

som.env_cor <- som.data[,c(5:21)]

# plotting the results of the correlations########################Step3#########

# plot of the species variables autocorrelation with spearman method
ggcorr(som_ffg_corr, method = c("pairwise", "spearman"),
       nbreaks = 6,
       label = TRUE,
       label_size = 4,
       color = "grey50",
       size = 4)
ggcorr(som_tg_corr, method = c("pairwise", "spearman"),
       nbreaks = 6,
       label = TRUE,
       label_size = 4,
       color = "grey50",
       size = 4)
ggcorr(som_hg_corr, method = c("pairwise", "spearman"),
       nbreaks = 6,
       label = TRUE,
       label_size = 4,
       color = "grey50",
       size = 4)
ggcorr(som_div_corr, method = c("pairwise", "spearman"),
       nbreaks = 6,
       label = TRUE,
       label_size = 4,
       color = "grey50",
       size = 4)
# plot of the soil environment variables autocorrelation########################
# with pearson method
ggcorr(som.env_cor, method = c("pairwise", "pearson"),
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
plotNormalHistogram(som.trans3, prob = FALSE,
                    main = "Normal Distribution overlay on Histogram",
                    length = 1000 )

# SOM training data############################################Step-6###########
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
plot(my.som, what="energy")+
  ggplot2::ggtitle("")

# Clustering distribution
# Hitmap visualization
plot(my.som, what = "obs",
     type = "hitmap")
#+theme_transparent()+
ggplot2::scale_fill_grey(start = 0.9, end = 0.3)

# umatrix#######################################################################
plot(my.som, what = "prototypes",
     type = "umatrix")
#+
ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  guides(fill = guide_legend(element_blank()))+
  theme_transparent()

# Distance between prototype
plot(my.som, what = "prototypes",
     type = "smooth.dist")
#+
ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("")


# SOM Grid ### This plot should be done after cluster finalization #############
pretty_palette<- c('#bddd22', '#8BD3E6', '#FF6D6A', "#EFBE7D") # , '', '#9467bd', ''

##################################################################
# Labeling all the cells in the SOM ##############################
##################################################################
label_1 <- "I"
label_2 <- "II"
label_3 <- "III"
label_4 <- "IV"
label_5 <- "V"

C1 <- "JR1 \n JP2 \n CR4 \n JR3 " 
C1_2 <- "\n BR1 \n BR9 \n CR3 \n JR2 \n CR6 "
C1_3 <- "CR5 \n BR6 \n FR1 \n JP3 "
C2 <- "BR5 \n BR7"
C3 <- "CR2 \n BR4"
C4 <- "CR1 \n CB8"
C5 <- "WR2"
C6 <- "CP3 \n BP1 \n BP3 \n WP1"
C8 <- " CP5 \n CP7 \n CP8 \n WP2 "
C8_1 <- "BP2 \n BP9 \n CP2 \n FR5 \n FP9 "
C8_2 <- "CP4 \n BP6 \n CP6 \n JP9"
C9 <- "WR1 \n WR3"
C10 <- "FR3"
C13 <- "FR4"
C14 <- "BP7"
C15 <- "JP8"
C16 <- "WB1"
C17 <- "CR7 \n BR8"
C18 <- "FR2"
C19 <- "JR6"
C21 <- "WR4"
C22 <- "JR9 \n FP3"
C24 <- "JB9 \n JP7 \n CP1 \n CB1"
C24_1 <- "CB2 \n FP5 \n JB8"
C25 <- "CR8"
C28 <- "JR5"
C29 <- "WR6"
C32 <- "WP4 \n FP4 \n WP3 \n CB3 "
C32_1 <- "BP8 \n FB8 \n FP8 \n Wb5"
C33 <- "JP1 \n BR3 \n CR9"
C35 <- "FR8 \n BB9 \n BB8 \n FR7 \n JR8"
C36 <- "FR9 \n BB2 \n BB4 \n JR7"
C37 <- "WB4 \n FB2 \n JB2 \n BR2"
C37_1 <- "BB3 \n CB9 \n FB9 \n BB1"
C37_2 <- "FB3\n FP2 \n WB8 \n BP4"
C38 <- "BB6 \n JB1"
C39 <- "WB7 \n CB5 \n FB1 \n BB7 "
C39_1 <- "CB4 \n CB6 \n FB4"
C40 <- "WB3 \n FP7 \n FP6 \n CP9"
C40_1 <- "CB7 \n FB5 \n FB6 \n FR6"

##########################################################################
p5 <- plot(my.sc, type = "grid", show.names=F) + 
  ggplot2::ggtitle("SOM")+ 
  ggplot2::scale_fill_grey(start = 0.9, end = 0.3) +
  #ggplot2::scale_fill_manual(values = pretty_palette) +
  ggplot2::scale_alpha_manual(values=c(0.01,0.01,0.01,0.01,0.01,0.01))+
  theme_transparent()+
  annotate ("text",size= 15,fontface = "bold",
            x = 2,y = 1.8,label = label_1)+
  annotate ("text",size= 15, fontface = "bold",
            x = 2, y = 3.5,label = label_2)+
  annotate ("text",size= 15, fontface = "bold",
            x = 1.5, y = 6.1,label = label_3)+
  annotate ("text",size= 15,fontface = "bold",
            x = 4.0, color= "white",y = 5.2,label = label_4)+
  annotate ("text",size= 15,fontface = "bold",
            x = 4.5, color= "white",y = 2.6,label = label_5)+
  
  annotate ("text",size= 3, x = 1.2, y = 0.9, label = C1)+
  annotate ("text",size= 3, x = 1.5, y = 0.9, label = C1_2)+
  annotate ("text",size= 3, x = 1.8, y = 0.9, label = C1_3)+
  annotate ("text",size= 3, x = 1.0, y = 1.7, label = C2)+
  annotate ("text",size= 3, x = 1.5, y = 2.6, label = C3)+
  annotate ("text",size= 3, x = 1.0, y = 3.5, label = C4)+
  annotate ("text",size= 3, x = 1.5, y = 4.4, label = C5)+
  annotate ("text",size= 3, x = 1.0, y = 5.2, label = C6)+
  annotate ("text",size= 3, x = 0.7, y = 7.0, label = C8)+
  annotate ("text",size= 3, x = 1.0, y = 7.0, label = C8_1)+
  annotate ("text",size= 3, x = 1.3, y = 7.0, label = C8_2)+
  annotate ("text",size= 3, x = 2.5, y = 0.9, label = C9)+
  annotate ("text",size= 3, x = 2.0, y = 1.5, label = C10)+
  annotate ("text",size= 3, x = 2.5, y = 4.4, label = C13)+
  annotate ("text",size= 3, x = 2.0, y = 5.2, label = C14)+
  annotate ("text",size= 3, x = 2.5, y = 6.1, label = C15, color ="white")+
  annotate ("text",size= 3, x = 2.0, y = 6.9, label = C16)+
  annotate ("text",size= 3, x = 3.5, y = 0.9, label = C17)+
  annotate ("text",size= 3, x = 3.0, y = 1.7, label = C18)+
  annotate ("text",size= 3, x = 2.5, y = 2.6, label = C19)+
  annotate ("text",size= 3, x = 3.5, y = 4.4, label = C21, color ="white")+
  annotate ("text",size= 3, x = 3.0, y = 5.2, label = C22, color ="white")+
  annotate ("text",size= 3, x = 2.9, y = 6.9, label = C24, color ="white")+
  annotate ("text",size= 3, x = 3.2, y = 6.9, label = C24_1, color ="white")+
  annotate ("text",size= 3, x = 4.5, y = 0.9, label = C25)+
  annotate ("text",size= 3, x = 4.0, y = 3.5, label = C28, color ="white")+
  annotate ("text",size= 3, x = 4.5, y = 4.4, label = C29, color ="white")+
  annotate ("text",size= 3, x = 3.9, y = 6.9, label = C32, color ="white")+
  annotate ("text",size= 3, x = 4.2, y = 6.9, label = C32_1, color ="white")+
  annotate ("text",size= 3, x = 5.5, y = 0.9, label = C33)+
  annotate ("text",size= 3, x = 5.5, y = 2.6, label = C35, color ="white")+
  annotate ("text",size= 3, x = 5.0, y = 3.5, label = C36, color ="white")+
  annotate ("text",size= 3, x = 5.2, y = 4.3, label = C37, color ="white")+
  annotate ("text",size= 3, x = 5.5, y = 4.3, label = C37_1, color ="white")+
  annotate ("text",size= 3, x = 5.8, y = 4.3, label = C37_2, color ="white")+
  annotate ("text",size= 3, x = 5.0, y = 5.2, label = C38, color ="white")+
  annotate ("text",size= 3, x = 5.4, y = 6.1, label = C39, color ="white")+
  annotate ("text",size= 3, x = 5.7, y = 6.1, label = C39_1, color ="white")+
  annotate ("text",size= 3, x = 4.9, y = 6.9, label = C40, color ="white")+
  annotate ("text",size= 3, x = 5.2, y = 6.9, label = C40_1, color ="white")+
  theme(legend.position="none")

# Ploting together for the paper ###############################################
first_plot = plot_grid(p5,  labels = c('a'), label_size = 25)
second_plot = plot_grid(p1,p3,p4,p8,
                        labels = c('b','c','d', 'e'), label_size = 25)

plot_grid(first_plot,
          second_plot,
          labels=c('', '', '','',''), ncol=2)

ggsave("Fig2.pdf", width = 40, height = 25, units = "cm", dpi = 300)


# Quality of the SOM-----------------------------------Step-7-------------------
quality(my.som, quality.type = "all")

# Summary of the SOM
summary(my.som)

# Building super-classes from the resulting SOM-------------Step-8--------------
plot(superClass(my.som))
my.sc <- superClass(my.som, k = 5)
plot(my.sc)

summary(my.sc) # significant with the cluster
a <-  plot(my.sc, plot.var = F, main= "SOM cluster",
           sub = "", cex= 0.5,
           xlab = "", ylab = "Euclidean distance")

plot(my.sc, show.names = T, type = "grid") # SOM grid
plot(my.sc, what = "obs",
     type = "hitmap",
     show.names = FALSE)+
  theme_bw() # Number of observation

#Indicator taxa -------------------------------Step-last---------------------------
#spe.all <- read.csv("data/rivermacro_spe.csv", stringsAsFactors = TRUE)
groups = c(rep(1,40), rep(2,35), rep(3,35)) # we used all species record
groups

matrix = som.data[,22:115]
indval = multipatt(matrix, groups, duleg=F, max.order = 3,
                   func = "IndVal.g", control = how(nperm = 999))
summary(indval, indvalcomp=TRUE)
summary(indval)

#-------------------------------------------------Step-9------------------------ 
# Clustering interpretation
class(som.data$habitat)
levels(som.data$habitat)
plot(my.som, what = "add", type = "pie", variable = som.data$habitat) +
  scale_fill_brewer(type = "qual") + 
  guides(fill = guide_legend(title = "Mesohabitat"))

plot(my.som, what = "add", type = "names", variable = som.data$stn)

plot(my.som, what = "prototypes", type = "poly.dist")


plot(my.som, what = "add", type = "names", variable = som.data$habitat)
plot(my.som, what = "add", type = "words", variable = som.data[,5:6])
plot(my.som, what = "add", type = "color", variable = som.data$S_TP)
plot(my.som, what = "prototypes", type = "color", variable = 1)

plot(my.som, what = "prototypes",
     variable = "FC")

# Species ploting according to the cluster

# And some others identify the super clusters with titles:

# Cluster I-------------------------Step-10-------------------------------------
plot(my.sc,show.names = FALSE, what = "prototypes",
     type = "color", variable = "P")

#+
ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Alotanypus sp.", subtitle="(3.0*)")+
  theme(plot.title = element_text(size = 10))+
  theme(legend.position="none")+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5,face="italic"),
        plot.subtitle=element_text(hjust=0.5))

b<-  plot(my.sc,show.names = FALSE, what = "prototypes",
          type = "color", variable = "trifsciatus")+
  #guides(fill = guide_legend(element_blank()))+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Cricotopus trifasciatus", subtitle="(19.75***)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

c<-   plot(my.sc,show.names = FALSE, what = "prototypes",
           type = "color", variable = "euiefferiella")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Eufiefferiella sp.", subtitle="(4.73**)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

d<-   plot(my.sc,show.names = FALSE, what = "prototypes",
           type = "color", variable = "halocladius")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Halocladius sp.", subtitle="(10.42***)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

e<-  plot(my.sc,show.names = FALSE, what = "prototypes",
          type = "color", variable = "pentaneurini")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Pentaneurini spp.", subtitle="(6.39***)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

f<- plot(my.sc,show.names = FALSE, what = "prototypes",
         type = "color", variable = "trichocladius")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Trichocladius sp.", subtitle="(3.45*)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

g<-  plot(my.sc,show.names = FALSE, what = "prototypes",
          type = "color", variable = "baetis")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Baetis sp.", subtitle="(64.29***)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

h<-  plot(my.sc,show.names = FALSE, what = "prototypes",
          type = "color", variable = "electrogena")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Electrogena sp.", subtitle="(6.74***)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

i<- plot(my.sc,show.names = FALSE, what = "prototypes",
         type = "color", variable = "serratella")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Serratella sp.", subtitle="(3.81**)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

j<-  plot(my.sc,show.names = FALSE, what = "prototypes",
          type = "color", variable = "ephemerella")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Ephemerella spp.", subtitle="(12.16***)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

k<-  plot(my.sc,show.names = FALSE, what = "prototypes",
          type = "color", variable = "Hydropsychae")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Hydropsyche spp.", subtitle="(106.74***)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

l<-  plot(my.sc,show.names = FALSE, what = "prototypes",
          type = "color", variable = "polyplectropus")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Polyplectropus sp.", subtitle="(6.51***)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

m<-  plot(my.sc,show.names = FALSE, what = "prototypes",
          type = "color", variable = "simulidae")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Simulidae", subtitle="(6.64***)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))

n<-  plot(my.sc,show.names = FALSE, what = "prototypes",
          type = "color", variable = "limnichidae")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("limnichidae", subtitle="(2.47*)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))



clustre_I <- plot_grid(a,b,c,d,e,f,g,
                       labels =c('','','','',
                                 '','','','',
                                 '','','','',''),
                       ncol=7, nrow =1)

clustre_I2 <- plot_grid(h,i,j,k,l,m,n,
                        l,m, labels =c('','','','',
                                       '','','','',
                                       '','','','',''),
                        ncol=7, nrow =1)


# now add the title
title <- ggdraw() + 
  draw_label(
    "I",
    fontface = 'bold',
    x = 0,
    hjust = 0,
    size = 24,
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )


fig3a <- plot_grid(
  title,
  clustre_I,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

# Cluster II ------------------------------------------------------------------- 

o <- plot(my.sc,show.names = FALSE, what = "prototypes",
          type = "color", variable = "psectrocladius")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Psectrocladius spp.", subtitle="(5.86***)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

p<-  plot(my.sc,show.names = FALSE, what = "prototypes",
          type = "color", variable = "paf")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("P. flavus", subtitle="(3.47*)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

q <-  plot(my.sc,show.names = FALSE, what = "prototypes",
           type = "color", variable = "sordidellus")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("P. sordidellus", subtitle="(5.33***)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

r<- plot(my.sc,show.names = FALSE, what = "prototypes",
         type = "color", variable = "chematopsyche")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Cheumatopsyche sp.", subtitle="(3.85**)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

s<-  plot(my.sc,show.names = FALSE, what = "prototypes",
          type = "color", variable = "metrichia")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Metrichia sp.", subtitle="(4.96**)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))


t<-  plot(my.sc,show.names = FALSE, what = "prototypes",
          type = "color", variable = "ecdyonurus")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Ecdyonurus sp.", subtitle="(3.78**)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

u<- plot(my.sc,show.names = FALSE, what = "prototypes",
         type = "color", variable = "limnichidae")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Limnichidae", subtitle="(3.78**)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

clustre_II <- plot_grid(o,p,q,r,s,t,u,
                        labels =c('','','','','',''),
                        ncol=7, nrow =1)

# now add the title
title2 <- ggdraw() + 
  draw_label(
    "II",
    fontface = 'bold',
    x = 0,
    hjust = 0,
    size = 24,
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

fig3b <- plot_grid(
  title2,
  clustre_II,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

# Cluster III ------------------------------------------------------------------

v<- plot(my.sc,show.names = FALSE, what = "prototypes",
         type = "color", variable = "cryptochironomus")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Cryptochironomus sp.", subtitle="(5.43***)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

w<- plot(my.sc,show.names = FALSE, what = "prototypes",
         type = "color", variable = "microtendipes")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Microtendipes sp.", subtitle="(2.74*)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

x<- plot(my.sc,show.names = FALSE, what = "prototypes",
         type = "color", variable = "polypedilum")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Polypedilum sp.", subtitle="(10.82***)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

y<- plot(my.sc,show.names = FALSE, what = "prototypes",
         type = "color", variable = "ephemera")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Ephemera sp.", subtitle="(12.17***)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

z<- plot(my.sc,show.names = FALSE, what = "prototypes",
         type = "color", variable = "tubifex")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Tubifex sp.", subtitle="(2.71*)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

a1<- plot(my.sc,show.names = FALSE, what = "prototypes",
          type = "color", variable = "aphelocheiridae")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Aphelocheirus sp.", subtitle="(3.22*)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))


clustre_III <- plot_grid(v,w,x,y,z,a1, 
                         labels =c('','','',
                                   '','','',''),
                         ncol=7, nrow =1)


# now add the title
title3 <- ggdraw() + 
  draw_label(
    "III",
    fontface = 'bold',
    x = 0,
    hjust = 0,
    size = 24,
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

fig3c <- plot_grid(
  title3,
  clustre_III,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

# Cluster IV -------------------------------------------------------------------

b1<- plot(my.sc,show.names = FALSE, what = "prototypes",
          type = "color", variable = "chironomus")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Chironomus sp.", subtitle="(3.00*)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5, face="italic"),
        plot.subtitle=element_text(hjust=0.5))

c1<- plot(my.sc,show.names = FALSE, what = "prototypes",
          type = "color", variable = "ephydridae")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Ephydridae", subtitle="(2.60*)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))

d1<- plot(my.sc,show.names = FALSE, what = "prototypes",
          type = "color", variable = "dolichopodidae")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Dolichopodidae", subtitle="(2.60*)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))

clustre_IV <- plot_grid(b1,c1,d1,
                        labels =c('',''),
                        nrow =1, ncol=4)

# now add the title
title4 <- ggdraw() + 
  draw_label(
    "IV",
    fontface = 'bold',
    x = 0,
    hjust = 0,
    size = 24,
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

fig3d <- plot_grid(
  title4,
  clustre_IV,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

# Cluster V

e1<-plot(my.sc,show.names = FALSE, what = "prototypes",
         type = "color", variable = "cardiocladius")+
  guides(fill = FALSE)+
  ggplot2::scale_fill_continuous(low = "white", high = "grey20")+
  ggplot2::ggtitle("Cardiocladius sp.", subtitle="(4.57**)")+
  theme(plot.title = element_text(size = 10))+
  theme_void()+
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))

clustre_V <- plot_grid(e1, cluster_p,
                       labels =c(''), 
                       nrow =1, ncol=3)

# now add the title
title5 <- ggdraw() + 
  draw_label(
    "V",
    fontface = 'bold',
    x = 0,
    hjust = 0,
    size = 24,
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

fig3e <- plot_grid(
  title5,
  clustre_V,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)



#----------------------------------------
cluster_p <- plot(my.sc, type = "grid", show.names=F) + 
  ggplot2::ggtitle("SOM")+ 
  ggplot2::scale_fill_grey(start = 0.9, end = 0.3) +
  #ggplot2::scale_fill_manual(values = pretty_palette) +
  ggplot2::scale_alpha_manual(values=c(0.01,0.01,0.01,0.01,0.01,0.01))+
  theme_transparent()+
  annotate ("text",size= 8,fontface = "bold",
            x = 2,y = 1.8,label = label_1)+
  annotate ("text",size= 8, fontface = "bold",
            x = 2, y = 3.5,label = label_2)+
  annotate ("text",size= 8, fontface = "bold",
            x = 1.5, y = 6.1,label = label_3)+
  annotate ("text",size= 8,fontface = "bold",
            x = 4.0, color= "white",y = 5.2,label = label_4)+
  annotate ("text",size= 8,fontface = "bold",
            x = 4.5, color= "white",y = 2.6,label = label_5)+
  theme(legend.position="none")


fig3de <- plot_grid(fig3d, fig3e,
                    nrow = 1, 
                    # rel_heights values control vertical title margins
                    rel_heights = c(0.1, 1)
)

# final Figure 3 and exporting--------------------------------------------------

plot_grid(fig3a, clustre_I2, fig3b, fig3c, fig3de,
          ncol=1)
ggsave("Figure_5.pdf", width = 300, height = 340, units = "mm", dpi = 300)

# ------------------------------------------------------------------------------


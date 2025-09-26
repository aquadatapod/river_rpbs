
################################################################################
# Script name: figure_2.R 
# Purpose: Generate Figure 2 – community structure and diversity patterns
# 
# Description: This script analyzes macroinvertebrate data from warm temperate rivers.
#              It runs Tukey post-hoc tests, ANOSIM for community differences,
#              and species accumulation curves, and combines them into Figure 2.
################################################################################


# Load packages ----------------------------------------------------------------
library(vegan)
library(plyr)
library(stats)
library(ape)
library(tidyverse)
library(car)
library(indicspecies) # for indicator taxa
library(wPerm)
attach(sem)

# loading data -----------------------------------------------------------------
tuk_spe <- read.csv("data/sem_all_river.csv", stringsAsFactors = T)

#--post hoc test on alpha diversity---------------------------------------------
# comparision between habitat

summary(spe.aov<-aov(richness~habitat,data=tuk_spe))
tukey <- TukeyHSD(spe.aov)
tukey
summary(tukey)
plot(tukey, col = c('#46a3fa','#e4d735','#e28320'))




# loading data for multivariate stats (other than Tukey-test)-------------------

# Species community whole data frame (fish abundance): "DoubsSpe.csv"
spe <- read.csv("data/rivermacro_spe.csv", stringsAsFactors = TRUE)
spe <- spe[-c(43,44,48,101,61,63,108,112,113,72,
              78,118,40,119,82,88,84,122,42,
              71,120,121,81,49,73),] 


# Environmental data frame: "DoubsEnv.csv"
env <- read.csv("data/riversed_env.csv", stringsAsFactors = TRUE)

env2 <- env[,-c(3,5,6,10)]

# selected species data
spe.data <- read.csv("data/spe.indicator.csv", stringsAsFactors = T)
spe.data<- spe.data[,-c(28,29)]


# knowing structure of the data ------------------------------------------------

dim(spe.data)
dim(env)

### Species distibution, all species confounded
(ab <- table(unlist(spe.data)))
barplot(ab, las=1, xlab="Abundance class", ylab="Frequency", col=grey(7:0/7))

### Number of absences
sum(spe.data==0)

### Proportion of zeros in the community data set
sum(spe.data==0)/(nrow(spe.data)*ncol(spe.data))

pairs(env, main="Bivariate Plots of the Environmental Data" )

# We can visually look for correlations between variables:
heatmap(abs(cor(env[, -c(1,2)])), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))
# Transformation species and env data ------------------------------------------
# Hellinger transformation
spe.hel<-decostand(spe[, -c(1,2)], method="hellinger")
spe.hel2<-decostand(spe.data[, -c(1,2)], method="hellinger")

spe_x <- (spe.data[, -c(1,2)]+1)
spe.log<-decostand(spe_x, method="log")

# Standardization of the 16 environmental data
env.z <- decostand(env[, -c(1,2)], method="standardize")
round(apply(env.z, 2, mean), 1) # the data are now centered (means~0)
apply(env.z, 2, sd)   # the data are now scaled (standard deviations=1)


### Plot the results
#### quick plots scaling 1 and 2
# grouping the factors habitat
levels(env$habitat) <- c("Rapid","Pool","Bench")
hab <- env$habitat
bg <- c("#e28320","#e4d735","#46a3fa") # 3nice colors for our habitats


# species accumulation curve-------------------------------------
sac <- specaccum(spe[, -c(1,2)])
sac2 <- specaccum(spe[, -c(1,2)],"random")
plot(sac, ci.type="polygon", ci.col="yellow")

#subset each habitat into its own df
spe %>% filter(habitat == "Rapid") -> rap
spe%>% filter(habitat == "Pool") -> pl
spe %>% filter(habitat == "Bench") -> bl

#calc species accumulation curve for each habitat
curve_rapid = specaccum(rap[, 3:40], method = "random")
curve_pool = specaccum(pl[, 3:29], method = "random")
curve_bench = specaccum(bl[, 3:28], method = "random")

# Step 1: Call the pdf command to start the plot
pdf(file = "plots/a.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
# Step 2: Create the plot with R code
plot(sac, ci.type="poly", col="grey", lwd=2, #plot curve_all first
     ci.lty=0, ci.col="lightgrey")
plot(curve_rapid, add = T, #then plot the rest
     ci.type="poly", col="orange",
     lwd=2, ci.lty=0, ci.col="#e28320") #col is COLOUR setting, so change it to something else if you want
plot(curve_pool, add = TRUE,
     ci.type="poly", col="yellow",
     lwd=2, ci.lty=0, ci.col="#e4d735")
plot(curve_bench, add = TRUE,
     ci.type="poly", col="blue",
     lwd=2, ci.lty=0, ci.col="#46a3fa")
legend("bottomright", legend=levels(hab), bty="n", col="gray32", pch=22, cex=1.5, pt.bg=bg)



# PERMANOVA---------------------------------------------------------------------
perm.oneway.anova(richness, habitat, R = 999)
perm.oneway.anova(richness, river, R = 999)
perm.oneway.anova(richness, site, R = 999)

# ANOSIM------------------------------------------------------------------------
#spe.dist$habitat = spe$habitat
grouping <- spe %>% dplyr::pull(habitat)
col.group1 <- c(rep("#e28320", times=40), rep("#e4d735",  times=28), rep("#46a3fa", times=29))
#grouping <- spe[["habitat"]]
spe.dist <- vegdist(spe.hel,"bray")
#spe.ano <- with(env, anosim(spe.dist, grouping))
#plot(spe.ano)
spe.ano <- anosim(spe.dist, grouping)
summary(spe.ano)

pdf(file = "plots/anosim.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches

plot(spe.ano, col = c('#bbbbbb','#46a3fa','#e4d735','#e28320'))

# Step 3: Run dev.off() to create the file!
dev.off()
#spe.ano2 <- with(env, anosim(spe.dist, habitat))
#plot(spe.ano2, col= env$habitat)
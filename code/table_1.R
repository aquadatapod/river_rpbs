
# This analysis is for the indicator taxa analysis 
# (Main manuscript Table 1)

library(vegan)
library(tidyverse)
library(indicspecies) # for indicator taxa



# Importing data for SOM analysis (This is the master file)----Step-1-----------
som.data <- read.csv("som.csv", stringsAsFactors = TRUE)


#Indicator taxa -------------------------------Step-last------------------------

groups = c(rep(1,40), rep(2,35), rep(3,35)) # we used all species record
groups

matrix = som.data[,22:115]
indval = multipatt(matrix, groups, duleg=F, max.order = 3,
                   func = "IndVal.g", control = how(nperm = 999))
summary(indval, indvalcomp=TRUE)
summary(indval)


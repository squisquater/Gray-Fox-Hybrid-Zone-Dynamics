# Load packages -- the order here is important because some pkg functions overwrite others.
#install.packages('ENMeval')
#install.packages('rasterVis')
#install.packages('rasterExtras')

#remotes::install_github("rspatial/rspatial")
#devtools::install_git('https://github.com/rsh249/rasterExtras')
#devtools::install_github("jamiemkass/ENMeval")

library(ENMeval)
library(raster)
library(dplyr)
library(rasterVis)
library(dismo)
library(ggplot2)


setwd("~/Desktop/ENMeval")

# Set a random seed in order to be able to reproduce this analysis.
set.seed(48)

#########################################################################
### You can search online databases like GBIF using the spocc package ###
#########################################################################

#uc <- spocc::occ('Urocyon cinereoargenteus', 'gbif', limit=500, has_coords=TRUE)
#uc.occs <- as.data.frame(uc$gbif$data$Urocyon_cinereoargenteus[,2:3])

#write.csv(uc.occs, file = "uc.occs.500.csv")

#uc.occs <- read.csv('uc.occs.1200.csv')

# Removing occurrences that have the same coordinates is good practice to avoid pseudoreplication.
#uc.occs <- uc.occs[!duplicated(uc.occs),]

#occs <- uc.occs

##################################################
##### I am going to use inaturalist instead ######
##################################################

#install.packages("rinat")
library(rinat)

inat_GF <- get_inat_obs(query = "Urocyon cinereoargenteus", quality = "research", geo = TRUE, maxresults = 10000)

#Double check the DF and remove any lines that have non gray fox that may have gotten through the search. 
inat_GF <- inat_GF[-c(209,3388,6780),]

## adding a consistent species ID column
inat_GF$spec <- "Urocyon cinereoargenteus"

## spatially thin the data
library(spThin)

## 25km Thinning parameter ###
## you may want to test out different thinning parameters ##
## This also writes an output file with thinned results as a .csv ##
inat_GF_thin25km <-
  thin( loc.data = inat_GF, 
        lat.col = "latitude", long.col = "longitude", 
        spec.col = "spec", 
        thin.par = 25, reps = 10, 
        locs.thinned.list.return = TRUE, 
        write.files = TRUE, 
        out.dir = "~/Desktop/ENMeval",
        write.log.file = FALSE)

inat_GF_thin25km_optimal <- read.csv("thinned_data_thin1_new_new.csv")

#### Note the above dataset can also be found in this folder for reproducing the ENM models


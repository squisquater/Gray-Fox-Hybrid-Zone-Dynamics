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


###########################
## Download Climate Data ##
###########################

#remotes::install_github("rspatial/geodata")
library(geodata)
wc = worldclim_global(var='bio', res = 5, path=tempdir())

### I would recommend downlaoding the CHELSA version which is more compatible 
###  with other paleoclim data. 
### If you use the worldclim dataset (above) some of the variables are on 
### opposite scales and you will end up with nonsense models when you predict 
### paleoclim models

#remotes::install_github("joeroe/rpaleoclim")
library("rpaleoclim")
library("terra") # For plotting

##### STEP 1: Download Data #####
### See Paleoclim website for info on different datasets ###
### http://www.paleoclim.org ###

##Set Extent
GrayFoxRange <- c(-125, -60, 0, 50)

paleoclim("cur", "5m", region = GrayFoxRange, cache_path = "~/Desktop/ENMeval")
paleoclim("lgm", "5m", region = GrayFoxRange, cache_path = "~/Desktop/ENMeval")
paleoclim("lig", "5m", region = GrayFoxRange, cache_path = "~/Desktop/ENMeval")


##### STEP 2: Load Data #####
#Current (1979 – 2013): Anthropocene, v1.2b**
envs.current <- load_paleoclim("CHELSA_cur_V1_2B_r5m.zip", as = c("terra", "raster"))

#Pleistocene: Last Glacial Maximum (ca. 21 ka), v1.2b**, NCAR CCSM4
envs.LGM <- load_paleoclim("chelsa_LGM_v1_2B_r5m.zip", as = c("terra", "raster"))

#Last Interglacial (ca. 130 ka), v1.
envs.LIG <- load_paleoclim("LIG_v1_5m.zip", as = c("terra", "raster"))


### Crop to extent
envs.current = crop(envs.current, GrayFoxRange)
envs.LGM = crop(envs.LGM, GrayFoxRange)
envs.LIG = crop(envs.LIG, GrayFoxRange)


## Choose predictor variables ##
predvars = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
preds = envs.current[[predvars]]

## Run model using maxnet model##
eval.19bioclim = ENMevaluate(occ=inat_GF_thin25km_optimal[,c('longitude', 'latitude')], env = preds, algorithm = 'maxnet', partitions ='randomkfold', parallel=TRUE, numCores = 12, fc=c("L", "Q", "P", "T", "H", "LQ","LQP", "LQH", "LQHP", "LQHPT"), RMvalues=seq(0.5, 2, 0.5), rasterPreds=NULL)

eval.19bioclim.results <- eval.19bioclim@results
write.csv(eval.19bioclim.results, "eval.19bioclim.results.csv")
bestmod = which(eval.19bioclim@results$AICc==min(eval.19bioclim@results$AICc))
eval.19bioclim@results[bestmod,]

# make prediction
pr = predict(preds, eval.19bioclim@models[[bestmod]], type = 'cloglog', na.rm = TRUE)
pr_df = as.data.frame(pr, xy=T)

#heatmap of current ENM
plot.current <- ggplot() +
  geom_raster(data = pr_df, aes(x = x, y = y, fill = lyr1)) +
  geom_point(data=inat_GF_thin25km_optimal, aes(x=longitude, y=latitude), col='white', cex=0.03) +
  coord_quickmap() +
  theme_classic() + 
  scale_fill_gradientn(colours=viridis::viridis(99, option = "magma", begin = 0, end = 1),
                       na.value = "black") +
  labs(title = "Current Ecological Niche Model",
    subtitle = "Anthropocene: (1979 – 2013)",
    x = "Longitude",
    y = "Latitude",
    colour = "Gears"
  )

plot.current

### Use the generated model to predict occurence during different time periods ###

  ###Last Glacial Maximum###
envs.LGM.predvars = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
envs.LGM.preds = envs.LGM[[envs.LGM.predvars]]

envs.LGM.pr = predict(envs.LGM, eval.19bioclim@models[[bestmod]], type = 'cloglog', na.rm = TRUE)
envs.LGM.pr_df = as.data.frame(envs.LGM.pr, xy=T)

plot.LGM <- ggplot() +
  geom_raster(data = envs.LGM.pr_df, aes(x = x, y = y, fill = lyr1)) +
  coord_quickmap() +
  theme_classic() + 
  scale_fill_gradientn(colours=viridis::viridis(99, option = "magma", begin = 0, end = 1),
                       na.value = "black") +
  labs(title = "Last Glacial Maximum Ecological Niche Model",
       subtitle = "ca. 21kya",
       x = "Longitude",
       y = "Latitude",
       colour = "Gears"
  )

plot.LGM 


###Last Interglacial###
envs.LIG.predvars = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
envs.LIG.preds = envs.LIG[[envs.LIG.predvars]]

envs.LIG.pr = predict(envs.LIG, eval.19bioclim@models[[bestmod]], type = 'cloglog', na.rm = TRUE)
envs.LIG.pr_df = as.data.frame(envs.LIG.pr, xy=T)

plot.LIG <- ggplot() +
  geom_raster(data = envs.LIG.pr_df, aes(x = x, y = y, fill = lyr1)) +
  coord_quickmap() +
  theme_classic() + 
  scale_fill_gradientn(colours=viridis::viridis(99, option = "magma", begin = 0, end = 1),
                       na.value = "black") +
  labs(title = "Last Interglacial Ecological Niche Model",
       subtitle = "ca. 130kya",
       x = "Longitude",
       y = "Latitude",
       colour = "Gears"
  )

plot.LIG 


##### Below you can find other datasets that you can download for predicting occurence in the past

### See Paleoclim website for info on different datasets ###
### http://www.paleoclim.org ###
paleoclim("mis19", "5m", region = GrayFoxRange, cache_path = "~/Desktop/ENMeval") 
paleoclim("hs1", "5m", region = GrayFoxRange, cache_path = "~/Desktop/ENMeval")
paleoclim("lh", "5m", region = GrayFoxRange, cache_path = "~/Desktop/ENMeval")
paleoclim("mh", "5m", region = GrayFoxRange, cache_path = "~/Desktop/ENMeval")
paleoclim("lig", "5m", region = GrayFoxRange, cache_path = "~/Desktop/ENMeval")
paleoclim("mpwp", "5m", region = GrayFoxRange, cache_path = "~/Desktop/ENMeval")

##Load the data as a rasterstack

#Pliocene: mid-Pliocene warm period (3.205 Ma), v1.0*
envs.mPWP <- load_paleoclim("mPWP_v1_r5m.zip", as = c("terra", "raster"))

#Pleistocene: mid-Holocene, Northgrippian (8.326-4.2 ka)
envs.MH <- load_paleoclim("MH_v1_5m.zip", as = c("terra", "raster"))

#Pleistocene: MIS19 (ca. 787 ka), v1.0*
envs.MIS19 <- load_paleoclim("MIS19_v1_r5m.zip", as = c("terra", "raster"))

#Pleistocene: Heinrich Stadial 1 (17.0-14.7 ka), v1.0
envs.HS1 <- load_paleoclim("HS1_v1_5m.zip", as = c("terra", "raster"))

#Pleistocene: late-Holocene, Meghalayan (4.2-0.3 ka), v1.0
envs.LH <- load_paleoclim("LH_v1_5m.zip", as = c("terra", "raster"))

#Last Interglacial (ca. 130 ka), v1.
envs.LIG <- load_paleoclim("LIG_v1_5m.zip", as = c("terra", "raster"))


#You can also set a threshold and create a categorical layer for habitat suitability
pr_df$Suitability <- ifelse(pr_df$lyr1 > 0.5,1,0)


threshold.current <- ggplot(data = pr_df, aes(x, y)) +
  geom_tile(aes(fill = Suitability)) + 
  scale_fill_gradient(low="gray", high="red") +
  theme_classic()

threshold.current


### Plot Fossil Data ###

#Read in associated fossil data for this time period
setwd("~/Desktop/GrayFox-ME-ManuscriptRevision/FigureUpdates/ENM")
list.files()
fossil.data <- read.csv("Urocyon_Cineroargenteus_PostLGPFossils_PBDB.FAUNMAP.csv")

plot.fossil.PreLGP <- ggplot() +
  geom_raster(data = pr_df, aes(x = x, y = y, fill = lyr1)) +
  geom_point(data=fossil.data, aes(x=longitude, y=latitude), col='black', cex=2, alpha=0.6) +
  coord_quickmap() +
  theme_classic() +
  scale_fill_gradient(low="darkgray",high="darkgray")




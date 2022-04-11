### This script has the code for running the combine models (eastern and western gray foxes).
### For running East and West separately see HZAR_East.R and HZAR_West.R. 

### Note that if you want information on AICc values, cline widths, and cline locations you will need to
### record these as you run the scripst or re-write the script to output to a summary file. 
### Otherwise the second set will overwrite the first set. PDFs of the resulting top cline model will be 
### output separately if you run the entire script at once. 

### Need to run a few different sets of models here
##### 1) Using the 50% Ancestry spline
##### 2) Using the modified ancestry spline (based on AHMM timing)

### For each set of models I will test 5 different tail shapes.
### All pmin and pmax will be set to 0 and 1 respectively.

library('hzar')

setwd("C:/Users/Sophie/Desktop/GrayFox/GBS/")

#################################
#####EBK 0.5 Ancestry Spline#####
#################################

## Read in Data
Nuclear <- read.csv("GrayFox_HZAR_Combined.csv", header=T, stringsAsFactors=F) 

## First create vector with distance based intervals
intervals <- seq(from = -2120, to = 1400, by =20)

## Add a new column to the matrix that assigns a location bin to the data
Nuclear$D.int <- findInterval(Nuclear$km, intervals, rightmost.closed = TRUE) 

## Quick plot to check data
plot(Nuclear$D.int, Nuclear$E.anc)

#####################################################################
### Make a distance file and summary matrix unique to this dataset###
#####################################################################

## Calculate the mean and Standard Error of each bin
unique <- unique(Nuclear$D.int)

## Make a list of locations
LocIDs <- Nuclear$D.int

p <- Nuclear$E.anc

## Read in distances of "sampling locations"
d <- read.csv("GrayFox_distances_ForHZAR_CombinedData_0.5splinemodel.csv", header=T, stringsAsFactors=F)

## Output vector of distances
d$distance -> dists

## Get population means
sapply(unique(LocIDs),function(x){mean(p[which(LocIDs %in% x)])},USE.NAMES=F) %>% round(4) -> means
means

## Get population variances
sapply(unique(LocIDs),function(x){var(p[which(LocIDs %in% x)])},USE.NAMES=F) -> vars

## Get population counts
sapply(unique(LocIDs),function(x){length(p[which(LocIDs %in% x)])},USE.NAMES=F) -> counts

## Create data frame for input into Hzar
data.frame(LocIDs=unique(LocIDs),dists,means,vars,counts) -> m

## quick plot
plot(m$dists, m$means)

################
##Model Clines##
################

mydata <- hzar.doMolecularData1DPops(
  distance=m$dists,
  pObs=m$means,
  nEff=m$counts,
  siteID=m$LocIDs)

clineModel_none <- hzar.makeCline1DFreq(data=mydata,tails="none")
clineModel_none <- hzar.model.addCenterRange(clineModel_none, -500,500)
clineModel_none <- hzar.model.addMaxWidth(meta.model=clineModel_none,maxValue=1000)

clineModel_right <- hzar.makeCline1DFreq(data=mydata,tails="right")
clineModel_right<- hzar.model.addCenterRange(clineModel_right, -500,500)
clineModel_right <- hzar.model.addMaxWidth(meta.model=clineModel_right,maxValue=1000)

clineModel_left <- hzar.makeCline1DFreq(data=mydata,tails="left")
clineModel_left <- hzar.model.addCenterRange(clineModel_left, -500,500)
clineModel_left <- hzar.model.addMaxWidth(meta.model=clineModel_left,maxValue=1000)

clineModel_both <- hzar.makeCline1DFreq(data=mydata,tails="both")
clineModel_both <- hzar.model.addCenterRange(clineModel_both, -500,500)
clineModel_both <- hzar.model.addMaxWidth(meta.model=clineModel_both,maxValue=1000)

clineModel_mirror <- hzar.makeCline1DFreq(data=mydata,tails="mirror")
clineModel_mirror <- hzar.model.addCenterRange(clineModel_both, -500,500)
clineModel_mirror <- hzar.model.addMaxWidth(meta.model=clineModel_mirror,maxValue=1000)

## Set initial pMin and pMax

0 -> hzar.meta.init(clineModel_none)$pMin
1 -> hzar.meta.init(clineModel_none)$pMax

0 -> hzar.meta.init(clineModel_right)$pMin
1 -> hzar.meta.init(clineModel_right)$pMax

0 -> hzar.meta.init(clineModel_left)$pMin
1 -> hzar.meta.init(clineModel_left)$pMax

0 -> hzar.meta.init(clineModel_both)$pMin
1 -> hzar.meta.init(clineModel_both)$pMax

0 -> hzar.meta.init(clineModel_mirror)$pMin
1 -> hzar.meta.init(clineModel_mirror)$pMax

hzar.meta.fix(clineModel_none)$pMin <- TRUE
hzar.meta.fix(clineModel_none)$pMax <- TRUE

hzar.meta.fix(clineModel_right)$pMin <- TRUE
hzar.meta.fix(clineModel_right)$pMax <- TRUE

hzar.meta.fix(clineModel_left)$pMin <- TRUE
hzar.meta.fix(clineModel_left)$pMax <- TRUE

hzar.meta.fix(clineModel_both)$pMin <- TRUE
hzar.meta.fix(clineModel_both)$pMax <- TRUE

hzar.meta.fix(clineModel_mirror)$pMin <- TRUE
hzar.meta.fix(clineModel_mirror)$pMax <- TRUE

## Run fit request
fitRequest_none <- hzar.first.fitRequest.old.ML(model=clineModel_none,obsData=mydata,verbose=T)
fitRequest_right <- hzar.first.fitRequest.old.ML(model=clineModel_right,obsData=mydata,verbose=T)
fitRequest_left <- hzar.first.fitRequest.old.ML(model=clineModel_left,obsData=mydata,verbose=T)
fitRequest_both <- hzar.first.fitRequest.old.ML(model=clineModel_both,obsData=mydata,verbose=T)
fitRequest_mirror <- hzar.first.fitRequest.old.ML(model=clineModel_mirror,obsData=mydata,verbose=T)

## Fit models
myfitlist_none <- hzar.chain.doSeq(hzar.request=fitRequest_none,count=3,collapse=F)
myfitlist_right <- hzar.chain.doSeq(hzar.request=fitRequest_right,count=3,collapse=F)
myfitlist_left <- hzar.chain.doSeq(hzar.request=fitRequest_left,count=3,collapse=F)
myfitlist_both <- hzar.chain.doSeq(hzar.request=fitRequest_both,count=3,collapse=F)
myfitlist_mirror <- hzar.chain.doSeq(hzar.request=fitRequest_mirror,count=3,collapse=F)

## Can plot them to observe differences
hzar.plot.cline(myfitlist_none[[3]],xlab="Distance (km) from Cline Center",ylab="Proportion of Eastern Gray Fox Ancestry", main="hzar.plot.cline(tails.none)")
hzar.plot.cline(myfitlist_right[[3]],xlab="Distance (km) from Cline Center",ylab="Proportion of Eastern Gray Fox Ancestry", main="hzar.plot.cline(tails.right)")
hzar.plot.cline(myfitlist_left[[3]],xlab="Distance (km) from Cline Center",ylab="Proportion of Eastern Gray Fox Ancestry", main="hzar.plot.cline(tails.left)")
hzar.plot.cline(myfitlist_both[[3]],xlab="Distance (km) from Cline Center",ylab="Proportion of Eastern Gray Fox Ancestry", main="hzar.plot.cline(tails.both)")
hzar.plot.cline(myfitlist_mirror[[3]],xlab="Distance (km) from Cline Center",ylab="Proportion of Eastern Gray Fox Ancestry", main="hzar.plot.cline(tails.mirror)")

## Need to group data to run AIC analyses
groupedData_none <- hzar.dataGroup.add(myfitlist_none)
groupedData_right <- hzar.dataGroup.add(myfitlist_right)
groupedData_left <- hzar.dataGroup.add(myfitlist_left)
groupedData_both <- hzar.dataGroup.add(myfitlist_both)
groupedData_mirror <- hzar.dataGroup.add(myfitlist_mirror)

## Generate sample size corrected AIC (AICc)
AICc_none <- hzar.AICc.hzar.dataGroup(groupedData_none)
AICc_right <- hzar.AICc.hzar.dataGroup(groupedData_right)
AICc_left <- hzar.AICc.hzar.dataGroup(groupedData_left)
AICc_both <- hzar.AICc.hzar.dataGroup(groupedData_both)
AICc_mirror <- hzar.AICc.hzar.dataGroup(groupedData_mirror)

## Fit to data group to enable next steps
hzar.fit2DataGroup(myfitlist_none[[3]]) -> fit3_none
hzar.fit2DataGroup(myfitlist_right[[3]]) -> fit3_right
hzar.fit2DataGroup(myfitlist_left[[3]]) -> fit3_left
hzar.fit2DataGroup(myfitlist_both[[3]]) -> fit3_both
hzar.fit2DataGroup(myfitlist_mirror[[3]]) -> fit3_mirror

## Get cline location and width for top model
### Top model = no tails ###
hzar.get.ML.cline(myfitlist_none[[3]])$param.all$width
hzar.get.ML.cline(myfitlist_none[[3]])$param.all$center

## Extract center & +/- 2 LL range for top model
hzar.getLLCutParam(fit3_none,"width",2)
hzar.getLLCutParam(fit3_none,"center",2)

## PDF graph version of top model
pdf("GrayFox_hzarcline_Combined_EBK0.5_tails.none.pdf",6,6)
hzar.plot.fzCline(fit3_none,xlab="Distance (km) from Cline Center",ylab="Proportion of Eastern Gray Fox Ancestry",main="GF.hzar.plot.cline(tails = none)")
#hzar.plot.cline(myfitlist_mtDNA[[3]],xlab="Smoothed coastline distance in km (AK to CA)",ylab="American mtDNA",main="hzar.plot.cline()")
dev.off()

##############################
#####AHMM Modified Spline#####
##############################

## Read in Data
Nuclear <- read.csv("GrayFox_HZAR_Combined.csv", header=T, stringsAsFactors=F) 

## first create vector with distance based intervals
intervals <- seq(from = -2180, to = 1340, by =20)

## Add a new column to the matrix that assigns a location bin to the data
Nuclear$D.int <- findInterval(Nuclear$km_new, intervals, rightmost.closed = TRUE) 

## Quick plot to check data
plot(Nuclear$D.int, Nuclear$E.anc)

#####################################################################
### Make a distance file and summary matrix unique to this dataset###
#####################################################################

## Calculate the mean and Standard Error of each bin
unique <- unique(Nuclear$D.int)

## Make a list of locations
LocIDs <- Nuclear$D.int

p <- Nuclear$E.anc

## Read in distances of "sampling locations"
d <- read.csv("GrayFox_distances_ForHZAR_CombinedData_AHMMmodel.csv", header=T, stringsAsFactors=F)

## Output vector of distances
d$distance -> dists

## Get population means
sapply(unique(LocIDs),function(x){mean(p[which(LocIDs %in% x)])},USE.NAMES=F) %>% round(4) -> means
means

## Get population variances
sapply(unique(LocIDs),function(x){var(p[which(LocIDs %in% x)])},USE.NAMES=F) -> vars

## Get population counts
sapply(unique(LocIDs),function(x){length(p[which(LocIDs %in% x)])},USE.NAMES=F) -> counts

## Create data frame for input into Hzar
data.frame(LocIDs=unique(LocIDs),dists,means,vars,counts) -> m

## quick plot
plot(m$dists, m$means)

################
##Model Clines##
################

mydata <- hzar.doMolecularData1DPops(
  distance=m$dists,
  pObs=m$means,
  nEff=m$counts,
  siteID=m$LocIDs)

clineModel_none <- hzar.makeCline1DFreq(data=mydata,tails="none")
clineModel_none <- hzar.model.addCenterRange(clineModel_none, -500,500)
clineModel_none <- hzar.model.addMaxWidth(meta.model=clineModel_none,maxValue=1000)

clineModel_right <- hzar.makeCline1DFreq(data=mydata,tails="right")
clineModel_right<- hzar.model.addCenterRange(clineModel_right, -500,500)
clineModel_right <- hzar.model.addMaxWidth(meta.model=clineModel_right,maxValue=1000)

clineModel_left <- hzar.makeCline1DFreq(data=mydata,tails="left")
clineModel_left <- hzar.model.addCenterRange(clineModel_left, -500,500)
clineModel_left <- hzar.model.addMaxWidth(meta.model=clineModel_left,maxValue=1000)

clineModel_both <- hzar.makeCline1DFreq(data=mydata,tails="both")
clineModel_both <- hzar.model.addCenterRange(clineModel_both, -500,500)
clineModel_both <- hzar.model.addMaxWidth(meta.model=clineModel_both,maxValue=1000)

clineModel_mirror <- hzar.makeCline1DFreq(data=mydata,tails="mirror")
clineModel_mirror <- hzar.model.addCenterRange(clineModel_both, -500,500)
clineModel_mirror <- hzar.model.addMaxWidth(meta.model=clineModel_mirror,maxValue=1000)


## Set initial pMin and pMax

0 -> hzar.meta.init(clineModel_none)$pMin
1 -> hzar.meta.init(clineModel_none)$pMax

0 -> hzar.meta.init(clineModel_right)$pMin
1 -> hzar.meta.init(clineModel_right)$pMax

0 -> hzar.meta.init(clineModel_left)$pMin
1 -> hzar.meta.init(clineModel_left)$pMax

0 -> hzar.meta.init(clineModel_both)$pMin
1 -> hzar.meta.init(clineModel_both)$pMax

0 -> hzar.meta.init(clineModel_mirror)$pMin
1 -> hzar.meta.init(clineModel_mirror)$pMax

hzar.meta.fix(clineModel_none)$pMin <- TRUE
hzar.meta.fix(clineModel_none)$pMax <- TRUE

hzar.meta.fix(clineModel_right)$pMin <- TRUE
hzar.meta.fix(clineModel_right)$pMax <- TRUE

hzar.meta.fix(clineModel_left)$pMin <- TRUE
hzar.meta.fix(clineModel_left)$pMax <- TRUE

hzar.meta.fix(clineModel_both)$pMin <- TRUE
hzar.meta.fix(clineModel_both)$pMax <- TRUE

hzar.meta.fix(clineModel_mirror)$pMin <- TRUE
hzar.meta.fix(clineModel_mirror)$pMax <- TRUE

## Run fit request
fitRequest_none <- hzar.first.fitRequest.old.ML(model=clineModel_none,obsData=mydata,verbose=T)
fitRequest_right <- hzar.first.fitRequest.old.ML(model=clineModel_right,obsData=mydata,verbose=T)
fitRequest_left <- hzar.first.fitRequest.old.ML(model=clineModel_left,obsData=mydata,verbose=T)
fitRequest_both <- hzar.first.fitRequest.old.ML(model=clineModel_both,obsData=mydata,verbose=T)
fitRequest_mirror <- hzar.first.fitRequest.old.ML(model=clineModel_mirror,obsData=mydata,verbose=T)

## Fit models
myfitlist_none <- hzar.chain.doSeq(hzar.request=fitRequest_none,count=3,collapse=F)
myfitlist_right <- hzar.chain.doSeq(hzar.request=fitRequest_right,count=3,collapse=F)
myfitlist_left <- hzar.chain.doSeq(hzar.request=fitRequest_left,count=3,collapse=F)
myfitlist_both <- hzar.chain.doSeq(hzar.request=fitRequest_both,count=3,collapse=F)
myfitlist_mirror <- hzar.chain.doSeq(hzar.request=fitRequest_mirror,count=3,collapse=F)

## Can plot them to observe differences
hzar.plot.cline(myfitlist_none[[3]],xlab="Distance (km) from Cline Center",ylab="Proportion of Eastern Gray Fox Ancestry", main="hzar.plot.cline(tails.none)")
hzar.plot.cline(myfitlist_right[[3]],xlab="Distance (km) from Cline Center",ylab="Proportion of Eastern Gray Fox Ancestry", main="hzar.plot.cline(tails.right)")
hzar.plot.cline(myfitlist_left[[3]],xlab="Distance (km) from Cline Center",ylab="Proportion of Eastern Gray Fox Ancestry", main="hzar.plot.cline(tails.left)")
hzar.plot.cline(myfitlist_both[[3]],xlab="Distance (km) from Cline Center",ylab="Proportion of Eastern Gray Fox Ancestry", main="hzar.plot.cline(tails.both)")
hzar.plot.cline(myfitlist_mirror[[3]],xlab="Distance (km) from Cline Center",ylab="Proportion of Eastern Gray Fox Ancestry", main="hzar.plot.cline(tails.mirror)")

## Need to group data to run AIC analyses
groupedData_none <- hzar.dataGroup.add(myfitlist_none)
groupedData_right <- hzar.dataGroup.add(myfitlist_right)
groupedData_left <- hzar.dataGroup.add(myfitlist_left)
groupedData_both <- hzar.dataGroup.add(myfitlist_both)
groupedData_mirror <- hzar.dataGroup.add(myfitlist_mirror)

## Generate sample size corrected AIC (AICc)
AICc_none <- hzar.AICc.hzar.dataGroup(groupedData_none)
AICc_right <- hzar.AICc.hzar.dataGroup(groupedData_right)
AICc_left <- hzar.AICc.hzar.dataGroup(groupedData_left)
AICc_both <- hzar.AICc.hzar.dataGroup(groupedData_both)
AICc_mirror <- hzar.AICc.hzar.dataGroup(groupedData_mirror)

## Fit to data group to enable next steps
hzar.fit2DataGroup(myfitlist_none[[3]]) -> fit3_none
hzar.fit2DataGroup(myfitlist_right[[3]]) -> fit3_right
hzar.fit2DataGroup(myfitlist_left[[3]]) -> fit3_left
hzar.fit2DataGroup(myfitlist_both[[3]]) -> fit3_both
hzar.fit2DataGroup(myfitlist_mirror[[3]]) -> fit3_mirror

### Get cline location and width for top model
### Top model = no tails ###
hzar.get.ML.cline(myfitlist_none[[3]])$param.all$width
hzar.get.ML.cline(myfitlist_none[[3]])$param.all$center

## Extract center & +/- 2 LL range for top model
hzar.getLLCutParam(fit3_none,"width",2)
hzar.getLLCutParam(fit3_none,"center",2)

## PDF graph version of top model
pdf("GrayFox_hzarcline_Combined_AHMMshiftedcenter_tails.none.pdf",6,6)
hzar.plot.fzCline(fit3_none,xlab="Distance (km) from Cline Center",ylab="Proportion of Eastern Gray Fox Ancestry",main="GF.hzar.plot.cline(tails = none)")
#hzar.plot.cline(myfitlist_mtDNA[[3]],xlab="Smoothed coastline distance in km (AK to CA)",ylab="American mtDNA",main="hzar.plot.cline()")
dev.off()

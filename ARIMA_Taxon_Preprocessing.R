# Isaiah E. Smith
# For Analytical Paleobiology Class Research Project
# Summer 2021

# This script is followed by ARIMA_Taxon_Analysis.R

rm(list=ls())


library(divDyn)
library(chronosphere)
library(forecast)


# get data from PBDB, and specify version
dat <- fetch("pbdb", ver="20210530")

#backup<-dat
#dat<-backup

#This is a pre-processed version of the PBDB that was downloaded with this API call:
attributes(dat)$chronosphere$API

#################### Data processing #################### 

## Initial filtering

# The downloaded dataset needs to be subsetted to represent the taxon of interest:

### Cambrian Fauna
dat <- dat[dat$class=="Trilobita",]
#dat<- dat[dat$order=="Lingulida",]
#dat <- dat[dat$phylum=="Hyolitha",]
#dat <- dat[dat$class=="Monoplacophora",]
#dat <- dat[dat$class=="Eocrinoidea",]

### Paleozoic Fauna
#dat <- dat[dat$class=="Stenolaemata",]
#dat <- dat[dat$class=="Ostracoda",]
#dat <- dat[dat$class=="Cephalopoda",]
#dat <- dat[dat$class=="Anthozoa",]
#dat <- dat[dat$class=="Crinoidea",]
#dat <- dat[dat$order=="Orthida" | dat$order=="Paleotremata" | dat$order=="Terebratulida" | dat$order=="Pentamerida" 
#           | dat$order=="Rhynchonellida" | dat$order=="Strophomenida" | dat$order=="Spiriferida" ,]
#dat <- dat[dat$class=="Asteroidea" | dat$class=="Ophiuroidea" | dat$order=="Somasteroidea",]



# The analyses will be conducted at the level of genus. 
# The entries that have no valid genus names in the PBDB can be omitted at this stage.

# omit non-genus entries
dat<- dat[!dat$genus=="", ]


# Stratigraphic resolution
data(stages)
data(tens)
data(keys)


# Assigns listed interval names to number values
stgMin <- categorize(dat[ ,"early_interval"], keys$stgInt)
stgMax <- categorize(dat[ ,"late_interval"], keys$stgInt)


# convert to numeric
stgMin <- as.numeric(stgMin)
stgMax <- as.numeric(stgMax)


# empty container
dat$stg <- rep(NA, nrow(dat))

# select entries, where
stgCondition <- c(
  # the early and late interval fields indicate the same stg
  which(stgMax==stgMin),
  # or the late_interval field is empty
  which(stgMax==-1))

# in these entries, use the stg indicated by the early_interval
dat$stg[stgCondition] <- stgMin[stgCondition]	


# Sampling assessment

table(dat$stg)

sum(table(dat$stg))

# proportion of the total data.
sum(table(dat$stg))/nrow(dat)


# As we cannot use unresolved occurrences in any way, we omit them at this stage.  

# omit unbinned
dats <- dat[!is.na(dat$stg),]



# "This function [binstat] will return the basic sampling summaries of a dataset" (from function documentation)

bs <- binstat(dats, tax="genus", bin="stg", 
              coll="collection_no", ref="reference_no", duplicates=FALSE)

bs$occs

# binstat() function fills in missing bins (species absent) from first bin to first taxon occurrence with NAs, 
# but does not continue filling with NA's after species goes extinct. We must correct for this by creating those
# absent bins and filling them with zeros until the present (until each stage is accounted for).

# create rows and fill with NA's
bs[(length(bs$occs)+1):length(stages$mid),] <-rep(NA, ncol(bs))

bs$occs

#Replace NA's with 0 in occs column
bs$occs[is.na(bs$occs)] <- 0

bs$occs


## Plotting

# Back 250 MYA
tsplot(stages, boxes="sys", boxes.col="systemCol", 
       shading="series", xlim=c(250, 50), ylim=c(0,2000))

# Back 500 MYA
tsplot(stages, boxes="sys", boxes.col="systemCol", 
       shading="series", xlim=c(500, 50), ylim=c(0,6000), ylab="Number occurrences")

lines(stages$mid, bs$occs)

# save to function
tp <- function(...)	tsplot(stages, boxes="sys", boxes.col="systemCol", 
                           shading="series", xlim=1:95, ...)


tp(ylim=c(0,3000), ylab="Number of collections")	
lines(stages$mid, bs$colls)





# Richness through time
dd <- divDyn(dats, tax="genus", bin="stg")

# plot range-through diversity
tp(ylim=c(0,500), ylab="richness")	
lines(stages$mid[1:51], dd$divRT, lwd=2)



# Test for correlation

# For now, I will only look for correlation in bins where there are occurrences present
bs_occs_pres <- bs$occs[bs$occs>0]

bs_occs_pres

dd_divRT_pres <- dd$divRT[!is.na(dd$divRT)]

dd_divRT_pres

cor.test(dd_divRT_pres, bs_occs_pres, method="spearman", na.action(na.omit))

############## Subsampling #################
# There is a strong correlation with per-bin diversity and number of occurrences, 
# so we will subsample to account for that

########### Shareholder Quorum Subsampling (SQS): 

my_q <- 0.5

sqsDD <- subsample(dats, tax="genus", bin="stg", q=my_q, type="sqs",
                   duplicates=FALSE, coll="collection_no", iter=300)

sqsDD[(max(dats$stg)+1):95,]<-NA

sqsDD[is.na(sqsDD)] <- 0


# plot subsampled data:

# range-through diversity
tp(ylim=c(0, 200), ylab="Subsampled richness (SQS, 0.5)")	
lines(stages$mid, sqsDD$divRT, lwd=2)

# corrected-sampled-in-bin diversity
tp(ylim=c(0,200), ylab="Subsampled richness (SQS, 0.5)")	
lines(stages$mid, sqsDD$divCSIB, lwd=2)




########### Classical Rarefaction :

bs$colls

# find minimum value of samples to use for subsampling quota 

my_cr_q <- min(bs$colls, na.rm = T)

cr <- subsample(dats, tax="genus", bin="stg", 
                  q=my_cr_q, coll="collection_no", duplicates=FALSE, iter=10) 


# fill out missing stages with NA's
cr[(max(dats$stg)+1):95,]<-NA

cr

x11()

tp(ylim=c(0,100), ylab="Subsampled richness")	

lines(stages$mid, cr$divRT, lwd=2)







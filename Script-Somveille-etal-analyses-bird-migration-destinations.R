
#### Script to run the analyses of migration destinations presented in Somveille et al. (2018) Where the wild birds go: explaining the differences in migration destinations across terrestrial bird species. Ecography


## Load required libraries

library(rgdal)
library(geosphere)
library(fields)
library(MASS)
library(raster)
library(move)
library(vioplot)
library(igraph)
library(emdist)


setwd("~/Data")


################################################################
################################################################

#######      Prepare species and environmental data    #########

################################################################
################################################################



## Load presence-absence data for migratory bird species

PresAbs_BR_NH <- read.table("PresAbs_BR_NH.txt", header=T)  # Breeding range of species that breed during the northern summer
PresAbs_NB_NH <- read.table("PresAbs_NB_NH.txt", header=T)  # Non-breeding range of species that breed during the northern summer
PresAbs_BR_SH <- read.table("PresAbs_BR_SH.txt", header=T)  # Breeding range of species that breed during the northern winter
PresAbs_NB_SH <- read.table("PresAbs_NB_SH.txt", header=T)  # Non-breeding range of species that breed during the northern winter


## Remove species that do not have at least 50% of either their breeding or non-breeding range on land

sppinfo <- read.csv("List_species_both&coastal.csv", sep=";", header=T)
PresAbs_BR_NH <- PresAbs_BR_NH[,-match(c(as.character(sppinfo$CoastalBothBoth),as.character(sppinfo$CoastalNHEH),as.character(sppinfo$CoastalNHWH),as.character(sppinfo$CoastalBothEH)), colnames(PresAbs_BR_NH), nomatch=0)]
PresAbs_NB_NH <- PresAbs_NB_NH[,-match(c(as.character(sppinfo$CoastalBothBoth),as.character(sppinfo$CoastalNHEH),as.character(sppinfo$CoastalNHWH),as.character(sppinfo$CoastalBothEH)), colnames(PresAbs_NB_NH), nomatch=0)]
PresAbs_BR_SH <- PresAbs_BR_SH[,-match(c(as.character(sppinfo$CoastalSHEH),as.character(sppinfo$CoastalSHWH)), colnames(PresAbs_BR_SH), nomatch=0)]
PresAbs_NB_SH <- PresAbs_NB_SH[,-match(c(as.character(sppinfo$CoastalSHEH),as.character(sppinfo$CoastalSHWH)), colnames(PresAbs_NB_SH), nomatch=0)]


## Load the global hexagon grid

hexgrid <- readOGR("Hex_grid", "isea3h7_hexgrid", verbose=FALSE)
hexgridWH <- hexgrid[which(hexgrid@data[,4] <= -30),]
hexgridEH <- hexgrid[which(hexgrid@data[,4] > -30),]


## Load environmental data

envData <- read.csv("Env_data.csv")
envData2 <- envData[-which((envData$Tmean_NW=="#DIV/0!" | envData$Tmean_NS=="#DIV/0!") | is.na(envData$habitatCoverage) == T),]  # Remove hexagons without environmental data 
TempMeanNS <- as.numeric(as.vector(envData2$Tmean_NS))		# Mean temperature during the northern summer
TempMeanNW <- as.numeric(as.vector(envData2$Tmean_NW))		# Mean temperature during the northern winter
TempMinNS <- envData2$Tmin_NS								# Minimum temperature during the northern summer
TempMinNW <- envData2$Tmin_NW								# Minimum temperature during the northern winter
habitat <- envData2$habitatCoverage							# Habitat coverage
NDVI_NS <- ((envData2$NDVI_may + envData2$NDVI_jun + envData2$NDVI_jul + envData2$NDVI_aug)/4) * 100		# NDVI during the northern summer
NDVI_NS[NDVI_NS<1] <- 1
NDVI_NW <- ((envData2$NDVI_nov + envData2$NDVI_dec + envData2$NDVI_jan + envData2$NDVI_feb)/4) * 100		# NDVI during the northern winter
NDVI_NW[NDVI_NW<1] <- 1


## Get hexagons ID and update hexagon grid object

hexidWH <- PresAbs_BR_NH[match(envData2[which(envData2$LONGITUDE<=-30),1], PresAbs_BR_NH[,1]),1]
hexidEH <- PresAbs_BR_NH[match(envData2[which(envData2$LONGITUDE>-30),1], PresAbs_BR_NH[,1]),1]
hexgrid2 <- rbind(hexgridWH[match(hexidWH, hexgridWH@data[,1]),], hexgridEH[match(hexidEH, hexgridEH@data[,1]),])


## Geographic coordinates of hexagons

lonlat <- cbind(envData2$LONGITUDE, envData2$LATITUDE)
colnames(lonlat) <- c("LONGITUDE", "LATITUDE")
west_Hem <- lonlat[which(lonlat[,1]<=-30),]  # subset for Western Hemisphere
east_Hem <- lonlat[which(lonlat[,1]>-30),]   # subset for Eastern Hemisphere


## z-tranformation of climate data

sdTempMean <- sd(c(TempMeanNS,TempMeanNW))
mnTempMean <- mean(c(TempMeanNS,TempMeanNW))
sdTempMin <- sd(c(TempMinNS,TempMinNW))
mnTempMin <- mean(c(TempMinNS,TempMinNW))
TempMeanNS <- (TempMeanNS - mnTempMean) / sdTempMean
TempMeanNW <- (TempMeanNW - mnTempMean) / sdTempMean
TempMinNS <- (TempMinNS - mnTempMin) / sdTempMin
TempMinNW <- (TempMinNW - mnTempMin) / sdTempMin


## Split environmental data into western hemisphere (WH) and eastern hemisphere (EH)

TempMeanNS_WH <- TempMeanNS[which(envData2$LONGITUDE<=-30)]
TempMeanNW_WH <- TempMeanNW[which(envData2$LONGITUDE<=-30)]
TempMinNS_WH <- TempMinNS[which(envData2$LONGITUDE<=-30)]
TempMinNW_WH <- TempMinNW[which(envData2$LONGITUDE<=-30)]
NDVI_NS_WH <- NDVI_NS[which(envData2$LONGITUDE<=-30)]
NDVI_NW_WH <- NDVI_NW[which(envData2$LONGITUDE<=-30)]
habitat_WH <- habitat[which(envData2$LONGITUDE<=-30)]
TempMeanNS_EH <- TempMeanNS[which(envData2$LONGITUDE>-30)]
TempMeanNW_EH <- TempMeanNW[which(envData2$LONGITUDE>-30)]
TempMinNS_EH <- TempMinNS[which(envData2$LONGITUDE>-30)]
TempMinNW_EH <- TempMinNW[which(envData2$LONGITUDE>-30)]
NDVI_NS_EH <- NDVI_NS[which(envData2$LONGITUDE>-30)]
NDVI_NW_EH <- NDVI_NW[which(envData2$LONGITUDE>-30)]
habitat_EH <- habitat[which(envData2$LONGITUDE>-30)]


## Resource gain and resource scarcity

resourceGain_NS_WH <- NDVI_NS_WH - NDVI_NW_WH
resourceGain_NW_WH <- NDVI_NW_WH - NDVI_NS_WH
resourceScarcity_NS_WH <- -resourceGain_NS_WH
resourceScarcity_NW_WH <- -resourceGain_NW_WH
resourceGain_NS_EH <- NDVI_NS_EH - NDVI_NW_EH
resourceGain_NW_EH <- NDVI_NW_EH - NDVI_NS_EH
resourceScarcity_NS_EH <- -resourceGain_NS_EH
resourceScarcity_NW_EH <- -resourceGain_NW_EH


## Remove presence-absence for hexagons without environmental data 

PresAbs_BR_NH <- PresAbs_BR_NH[-which((envData$Tmean_NW=="#DIV/0!" & envData$Tmean_NS=="#DIV/0!") | is.na(envData$habitatCoverage) == T),]
PresAbs_NB_NH <- PresAbs_NB_NH[-which((envData$Tmean_NW=="#DIV/0!" & envData$Tmean_NS=="#DIV/0!") | is.na(envData$habitatCoverage) == T),]
PresAbs_BR_SH <- PresAbs_BR_SH[-which((envData$Tmean_NW=="#DIV/0!" & envData$Tmean_NS=="#DIV/0!") | is.na(envData$habitatCoverage) == T),]
PresAbs_NB_SH <- PresAbs_NB_SH[-which((envData$Tmean_NW=="#DIV/0!" & envData$Tmean_NS=="#DIV/0!") | is.na(envData$habitatCoverage) == T),]

# Western Hemisphere
PresAbs_BR_NH_WH <- PresAbs_BR_NH[match(envData2[which(envData2$LONGITUDE<=-30),1], PresAbs_BR_NH[,1]),]
PresAbs_NB_NH_WH <- PresAbs_NB_NH[match(envData2[which(envData2$LONGITUDE<=-30),1], PresAbs_NB_NH[,1]),]
PresAbs_BR_SH_WH <- PresAbs_BR_SH[match(envData2[which(envData2$LONGITUDE<=-30),1], PresAbs_BR_SH[,1]),]
PresAbs_NB_SH_WH <- PresAbs_NB_SH[match(envData2[which(envData2$LONGITUDE<=-30),1], PresAbs_NB_SH[,1]),]

# Eastern Hemisphere
PresAbs_BR_NH_EH <- PresAbs_BR_NH[match(envData2[which(envData2$LONGITUDE>-30),1], PresAbs_BR_NH[,1]),]
PresAbs_NB_NH_EH <- PresAbs_NB_NH[match(envData2[which(envData2$LONGITUDE>-30),1], PresAbs_NB_NH[,1]),]
PresAbs_BR_SH_EH <- PresAbs_BR_SH[match(envData2[which(envData2$LONGITUDE>-30),1], PresAbs_BR_SH[,1]),]
PresAbs_NB_SH_EH <- PresAbs_NB_SH[match(envData2[which(envData2$LONGITUDE>-30),1], PresAbs_NB_SH[,1]),]



## Compute seasonal range size (i.e. number of hexagons in which the species is present at a given season)

rangeSize_BR_NH_WH <- apply(PresAbs_BR_NH_WH[,-1], 2, sum)
rangeSize_BR_SH_WH <- apply(PresAbs_BR_SH_WH[,-1], 2, sum)
rangeSize_NB_NH_WH <- apply(PresAbs_NB_NH_WH[,-1], 2, sum)
rangeSize_NB_SH_WH <- apply(PresAbs_NB_SH_WH[,-1], 2, sum)
rangeSize_BR_NH_EH <- apply(PresAbs_BR_NH_EH[,-1], 2, sum)
rangeSize_BR_SH_EH <- apply(PresAbs_BR_SH_EH[,-1], 2, sum)
rangeSize_NB_NH_EH <- apply(PresAbs_NB_NH_EH[,-1], 2, sum)
rangeSize_NB_SH_EH <- apply(PresAbs_NB_SH_EH[,-1], 2, sum)
rangeSizes_WH <- c(rangeSize_BR_NH_WH, rangeSize_BR_SH_WH, rangeSize_NB_NH_WH, rangeSize_NB_SH_WH)
rangeSizes_EH <- c(rangeSize_BR_NH_EH, rangeSize_BR_SH_EH, rangeSize_NB_NH_EH, rangeSize_NB_SH_EH)


## Compute range overlap (i.e. number of hexagons in which the species is present during both breeding and non-breeding seasons divided by the number of hexagons in which the species is present during at least breeding or non-breeding)

# Western Hemisphere
perm_NH_WH <- PresAbs_BR_NH_WH * PresAbs_NB_NH_WH
perm_NH_WH <- apply(perm_NH_WH[,-1], 2, sum)
rangeOverlap_NH_WH <- perm_NH_WH / (rangeSize_BR_NH_WH + rangeSize_NB_NH_WH - perm_NH_WH)
rangeOverlap_NH_WH  <- rangeOverlap_NH_WH[-which(rangeOverlap_NH_WH == "NaN")]
perm_SH_WH <- PresAbs_BR_SH_WH * PresAbs_NB_SH_WH
perm_SH_WH <- apply(perm_SH_WH[,-1], 2, sum)
rangeOverlap_SH_WH <- perm_SH_WH / (rangeSize_BR_SH_WH + rangeSize_NB_SH_WH - perm_SH_WH)
rangeOverlap_SH_WH  <- rangeOverlap_SH_WH[-which(rangeOverlap_SH_WH == "NaN")]

# Eastern Hemisphere
perm_NH_EH <- PresAbs_BR_NH_EH * PresAbs_NB_NH_EH
perm_NH_EH <- apply(perm_NH_EH[,-1], 2, sum)
rangeOverlap_NH_EH <- perm_NH_EH / (rangeSize_BR_NH_EH + rangeSize_NB_NH_EH - perm_NH_EH)
rangeOverlap_NH_EH  <- rangeOverlap_NH_EH[-which(rangeOverlap_NH_EH == "NaN")]
perm_SH_EH <- PresAbs_BR_SH_EH * PresAbs_NB_SH_EH
perm_SH_EH <- apply(perm_SH_EH[,-1], 2, sum)
rangeOverlap_SH_EH <- perm_SH_EH / (rangeSize_BR_SH_EH + rangeSize_NB_SH_EH - perm_SH_EH)
rangeOverlap_SH_EH  <- rangeOverlap_SH_EH[-which(rangeOverlap_SH_EH == "NaN")]


## Migratory species with less than 20% overlap between their breeding and non-breeding ranges

selectedSpeciesNH_EH <- which(rangeOverlap_NH_EH < 0.2)
selectedSpeciesSH_EH <- which(rangeOverlap_SH_EH < 0.2)
selectedSpeciesNH_WH <- which(rangeOverlap_NH_WH < 0.2)
selectedSpeciesSH_WH <- which(rangeOverlap_SH_WH < 0.2)


## Correct some mistakes in coding seasons

#Anthus hoeschi
selectedSpeciesSH_EH <-  c(selectedSpeciesSH_EH, selectedSpeciesNH_EH[30])
selectedSpeciesNH_EH <- selectedSpeciesNH_EH[-30]

#Gallinago.nigripennis
selectedSpeciesNH_EH <- selectedSpeciesNH_EH[-162]

#Pitta.moluccensis and Sylvia.cantillans 
selectedSpeciesNH_EH <-  c(selectedSpeciesNH_EH, selectedSpeciesSH_EH[c(13,14)])
selectedSpeciesSH_EH <- selectedSpeciesSH_EH[-c(13,14)]


## Keep the presence-absence data only of selected species

PresAbs_BR_NH_WH2 <- PresAbs_BR_NH_WH[,match(names(selectedSpeciesNH_WH), colnames(PresAbs_BR_NH_WH))]
PresAbs_NB_NH_WH2 <- PresAbs_NB_NH_WH[,match(names(selectedSpeciesNH_WH), colnames(PresAbs_NB_NH_WH))]
PresAbs_BR_NH_EH2 <- PresAbs_BR_NH_EH[,match(names(selectedSpeciesNH_EH)[1:386], colnames(PresAbs_BR_NH_EH))]
PresAbs_NB_NH_EH2 <- PresAbs_NB_NH_EH[,match(names(selectedSpeciesNH_EH)[1:386], colnames(PresAbs_NB_NH_EH))]
PresAbs_BR_NH_EH2 <- cbind(PresAbs_BR_NH_EH2, PresAbs_BR_SH_EH[,match(names(selectedSpeciesNH_EH)[387:388], colnames(PresAbs_BR_SH_EH))])
PresAbs_NB_NH_EH2 <- cbind(PresAbs_NB_NH_EH2, PresAbs_NB_SH_EH[,match(names(selectedSpeciesNH_EH)[387:388], colnames(PresAbs_NB_SH_EH))])
PresAbs_BR_SH_WH2 <- PresAbs_BR_SH_WH[,match(names(selectedSpeciesSH_WH), colnames(PresAbs_BR_SH_WH))]
PresAbs_NB_SH_WH2 <- PresAbs_NB_SH_WH[,match(names(selectedSpeciesSH_WH), colnames(PresAbs_NB_SH_WH))]
PresAbs_BR_SH_EH2 <- PresAbs_BR_SH_EH[,match(names(selectedSpeciesSH_EH)[1:12], colnames(PresAbs_BR_SH_EH))]
PresAbs_NB_SH_EH2 <- PresAbs_NB_SH_EH[,match(names(selectedSpeciesSH_EH)[1:12], colnames(PresAbs_NB_SH_EH))]
PresAbs_BR_SH_EH2 <- cbind(PresAbs_BR_SH_EH2, PresAbs_BR_NH_EH[,match(names(selectedSpeciesSH_EH)[13], colnames(PresAbs_BR_NH_EH))])
PresAbs_NB_SH_EH2 <- cbind(PresAbs_NB_SH_EH2, PresAbs_NB_NH_EH[,match(names(selectedSpeciesSH_EH)[13], colnames(PresAbs_NB_NH_EH))])


## Remove species with 0 presences (or only 1 because it creates problems for the kernel estimation) during at least one season

toremove_NH_WH <- which(apply(PresAbs_BR_NH_WH2, 2, sum)<=1 | apply(PresAbs_NB_NH_WH2, 2, sum)<=1)
if(length(toremove_NH_WH) > 0){
	PresAbs_BR_NH_WH2 <- PresAbs_BR_NH_WH2[,-toremove_NH_WH]
	PresAbs_NB_NH_WH2 <- PresAbs_NB_NH_WH2[,-toremove_NH_WH]
}
toremove_SH_WH <- which(apply(PresAbs_BR_SH_WH2, 2, sum)<=1 | apply(PresAbs_NB_SH_WH2, 2, sum)<=1)
if(length(toremove_SH_WH) > 0){
	PresAbs_BR_SH_WH2 <- PresAbs_BR_SH_WH2[,-toremove_SH_WH]
	PresAbs_NB_SH_WH2 <- PresAbs_NB_SH_WH2[,-toremove_SH_WH]
}
toremove_NH_EH <- which(apply(PresAbs_BR_NH_EH2, 2, sum)<=1 | apply(PresAbs_NB_NH_EH2, 2, sum)<=1)
if(length(toremove_NH_EH) > 0){
	PresAbs_BR_NH_EH2 <- PresAbs_BR_NH_EH2[,-toremove_NH_EH]
	PresAbs_NB_NH_EH2 <- PresAbs_NB_NH_EH2[,-toremove_NH_EH]
}
toremove_SH_EH <- which(apply(PresAbs_BR_SH_EH2, 2, sum)<=1 | apply(PresAbs_NB_SH_EH2, 2, sum)<=1)
if(length(toremove_SH_EH) > 0){
	PresAbs_BR_SH_EH2 <- PresAbs_BR_SH_EH2[,-toremove_SH_EH]
	PresAbs_NB_SH_EH2 <- PresAbs_NB_SH_EH2[,-toremove_SH_EH]
}





################################################################
################################################################

######      Investigate trade-offs between variables    ########

################################################################
################################################################



###  Thermal distance between breeding and non-breeding grounds  ###


## Convert temperature values (mean and min) associated with presences into a density raster

nicheDensityRaster <- function(seasonalNiche){
	niche.kernel <- kde2d(seasonalNiche[,1], seasonalNiche[,2], n=20, h=1, lims=c(-3.4,1.6, -3.3,1.4)) 
	niche.kernel$z = niche.kernel$z/sum(niche.kernel$z)
	niche.raster <- raster(niche.kernel)
	threshold99=0; i=0
	while(threshold99 <= 0.99){
		i=i+1
		threshold99 = threshold99 + sort(as.vector(niche.raster), decreasing=T)[i]
	}
	niche.raster[which(as.vector(niche.raster) < sort(as.vector(niche.raster), decreasing=T)[i])] = 0
	niche.raster = niche.raster / sum(as.vector(niche.raster))	
	return(niche.raster)
}

breeding.thermal.nichesNH_WH <- apply(PresAbs_BR_NH_WH2, 2, function(x) nicheDensityRaster(cbind(TempMeanNS_WH[which(x==1)], TempMinNS_WH[which(x==1)])))
nonbreeding.thermal.nichesNH_WH <- apply(PresAbs_NB_NH_WH2, 2, function(x) nicheDensityRaster(cbind(TempMeanNW_WH[which(x==1)], TempMinNW_WH[which(x==1)])))
breeding.thermal.nichesSH_WH <- apply(PresAbs_BR_SH_WH2, 2, function(x) nicheDensityRaster(cbind(TempMeanNW_WH[which(x==1)], TempMinNW_WH[which(x==1)])))
nonbreeding.thermal.nichesSH_WH <- apply(PresAbs_NB_SH_WH2, 2, function(x) nicheDensityRaster(cbind(TempMeanNS_WH[which(x==1)], TempMinNS_WH[which(x==1)])))
breeding.thermal.nichesNH_EH <- apply(PresAbs_BR_NH_EH2, 2, function(x) nicheDensityRaster(cbind(TempMeanNS_EH[which(x==1)], TempMinNS_EH[which(x==1)])))
nonbreeding.thermal.nichesNH_EH <- apply(PresAbs_NB_NH_EH2, 2, function(x) nicheDensityRaster(cbind(TempMeanNW_EH[which(x==1)], TempMinNW_EH[which(x==1)])))
breeding.thermal.nichesSH_EH <- apply(PresAbs_BR_SH_EH2, 2, function(x) nicheDensityRaster(cbind(TempMeanNW_EH[which(x==1)], TempMinNW_EH[which(x==1)])))
nonbreeding.thermal.nichesSH_EH <- apply(PresAbs_NB_SH_EH2, 2, function(x) nicheDensityRaster(cbind(TempMeanNS_EH[which(x==1)], TempMinNS_EH[which(x==1)])))


## Compute thermal distance (using the earth mover's distance) for every migratory species

thermal.distancesNH_WH <- mapply(function(X,Y){ emd2d(as.matrix(X),as.matrix(Y)) }, X=breeding.thermal.nichesNH_WH, Y=nonbreeding.thermal.nichesNH_WH)
thermal.distancesSH_WH <- mapply(function(X,Y){ emd2d(as.matrix(X),as.matrix(Y)) }, X=breeding.thermal.nichesSH_WH, Y=nonbreeding.thermal.nichesSH_WH)
thermal.distancesNH_EH <- mapply(function(X,Y){ emd2d(as.matrix(X),as.matrix(Y)) }, X=breeding.thermal.nichesNH_EH, Y=nonbreeding.thermal.nichesNH_EH)
thermal.distancesSH_EH <- mapply(function(X,Y){ emd2d(as.matrix(X),as.matrix(Y)) }, X=breeding.thermal.nichesSH_EH, Y=nonbreeding.thermal.nichesSH_EH)
thermal.distances.obs = c(thermal.distancesNH_WH, thermal.distancesNH_EH, thermal.distancesSH_WH, thermal.distancesSH_EH)



###  Habitat distance between breeding and non-breeding grounds  ###

breeding.habitatNH_WH <- apply(PresAbs_BR_NH_WH2, 2, function(x) mean(habitat[which(x==1)]))
nonbreeding.habitatNH_WH <- apply(PresAbs_NB_NH_WH2, 2, function(x) mean(habitat[which(x==1)]))
breeding.habitatSH_WH <- apply(PresAbs_BR_SH_WH2, 2, function(x) mean(habitat[which(x==1)]))
nonbreeding.habitatSH_WH <- apply(PresAbs_NB_SH_WH2, 2, function(x) mean(habitat[which(x==1)]))
breeding.habitatNH_EH <- apply(PresAbs_BR_NH_EH2, 2, function(x) mean(habitat[which(x==1)]))
nonbreeding.habitatNH_EH <- apply(PresAbs_NB_NH_EH2, 2, function(x) mean(habitat[which(x==1)]))
breeding.habitatSH_EH <- apply(PresAbs_BR_SH_EH2, 2, function(x) mean(habitat[which(x==1)]))
nonbreeding.habitatSH_EH <- apply(PresAbs_NB_SH_EH2, 2, function(x) mean(habitat[which(x==1)]))

habitat.distancesNH_WH <- abs(breeding.habitatNH_WH - nonbreeding.habitatNH_WH)
habitat.distancesSH_WH <- abs(breeding.habitatSH_WH - nonbreeding.habitatSH_WH)
habitat.distancesNH_EH <- abs(breeding.habitatNH_EH - nonbreeding.habitatNH_EH)
habitat.distancesSH_EH <- abs(breeding.habitatSH_EH - nonbreeding.habitatSH_EH)
habitat.distances.obs = c(habitat.distancesNH_WH, habitat.distancesNH_EH, habitat.distancesSH_WH, habitat.distancesSH_EH)



###  Geographical distance between breeding and non-breeding grounds  ###

centroids.breeding.groundsNH_WH <- t(apply(PresAbs_BR_NH_WH2, 2, function(x) apply(west_Hem[which(x==1),], 2, mean)))
centroids.nonbreeding.groundsNH_WH <- t(apply(PresAbs_NB_NH_WH2, 2, function(x) apply(west_Hem[which(x==1),], 2, mean)))
centroids.breeding.groundsSH_WH <- t(apply(PresAbs_BR_SH_WH2, 2, function(x) apply(west_Hem[which(x==1),], 2, mean)))
centroids.nonbreeding.groundsSH_WH <- t(apply(PresAbs_NB_SH_WH2, 2, function(x) apply(west_Hem[which(x==1),], 2, mean)))
centroids.breeding.groundsNH_EH <- t(apply(PresAbs_BR_NH_EH2, 2, function(x) apply(east_Hem[which(x==1),], 2, mean)))
centroids.nonbreeding.groundsNH_EH <- t(apply(PresAbs_NB_NH_EH2, 2, function(x) apply(east_Hem[which(x==1),], 2, mean)))
centroids.breeding.groundsSH_EH <- t(apply(PresAbs_BR_SH_EH2, 2, function(x) apply(east_Hem[which(x==1),], 2, mean)))
centroids.nonbreeding.groundsSH_EH <- t(apply(PresAbs_NB_SH_EH2, 2, function(x) apply(east_Hem[which(x==1),], 2, mean)))

geo.distancesNH_WH <- rdist.earth.vec(centroids.breeding.groundsNH_WH, centroids.nonbreeding.groundsNH_WH, miles=F)
geo.distancesSH_WH <- rdist.earth.vec(centroids.breeding.groundsSH_WH, centroids.nonbreeding.groundsSH_WH, miles=F)
geo.distancesNH_EH <- rdist.earth.vec(centroids.breeding.groundsNH_EH, centroids.nonbreeding.groundsNH_EH, miles=F)
geo.distancesSH_EH <- rdist.earth.vec(centroids.breeding.groundsSH_EH, centroids.nonbreeding.groundsSH_EH, miles=F)
geo.distances.obs <- c(geo.distancesNH_WH, geo.distancesNH_EH, geo.distancesSH_WH, geo.distancesSH_EH)



###  Resource scarcity  ###

breeding.resource.scarcityNH_WH <- apply(PresAbs_BR_NH_WH2, 2, function(x) mean(resourceScarcity_NS_WH[which(x==1)]))
nonbreeding.resource.scarcityNH_WH <- apply(PresAbs_NB_NH_WH2, 2, function(x) mean(resourceScarcity_NW_WH[which(x==1)]))
resource.scarcityNH_WH <- breeding.resource.scarcityNH_WH + nonbreeding.resource.scarcityNH_WH
breeding.resource.scarcitySH_WH <- apply(PresAbs_BR_SH_WH2, 2, function(x) mean(resourceScarcity_NW_WH[which(x==1)]))
nonbreeding.resource.scarcitySH_WH <- apply(PresAbs_NB_SH_WH2, 2, function(x) mean(resourceScarcity_NS_WH[which(x==1)]))
resource.scarcitySH_WH <- breeding.resource.scarcitySH_WH + nonbreeding.resource.scarcitySH_WH
breeding.resource.scarcityNH_EH <- apply(PresAbs_BR_NH_EH2, 2, function(x) mean(resourceScarcity_NS_EH[which(x==1)]))
nonbreeding.resource.scarcityNH_EH <- apply(PresAbs_NB_NH_EH2, 2, function(x) mean(resourceScarcity_NW_EH[which(x==1)]))
resource.scarcityNH_EH <- breeding.resource.scarcityNH_EH + nonbreeding.resource.scarcityNH_EH
breeding.resource.scarcitySH_EH <- apply(PresAbs_BR_SH_EH2, 2, function(x) mean(resourceScarcity_NW_EH[which(x==1)]))
nonbreeding.resource.scarcitySH_EH <- apply(PresAbs_NB_SH_EH2, 2, function(x) mean(resourceScarcity_NS_EH[which(x==1)]))
resource.scarcitySH_EH <- breeding.resource.scarcitySH_EH + nonbreeding.resource.scarcitySH_EH
resource.scarcity.obs <- c(resource.scarcityNH_WH, resource.scarcityNH_EH, resource.scarcitySH_WH, resource.scarcitySH_EH)



####  Figure 1: Relationships between the geographical distance separating the breeding and non-breeding ranges, the year-round availability of surplis primary productivity and temperature difference between the breeding and non-breeding seasons (thermal distance) across migratory species  #### 

par(mfrow=c(1,1), mar=c(2.9,3,1.5,3), mgp=c(1.5,0.5,0))
hist(sqrt(geo.distances.obs), xlim=c(0,120), xlab="", ylab="", main="", xaxt="n", axes=F, col="light grey", border="grey")
axis(side=4)
axis(side=1)
par(new=T)
plot(sqrt(geo.distances.obs), resource.scarcity.obs, xlim=c(0,120), xlab="Geographical distance (square-root)", ylab="Resource scarcity", main="", xaxt="n", axes=F, pch=20, cex=0.9, cex.lab=1, col="yellow")
axis(side=2)
points(sqrt(geo.distances.obs)[which(thermal.distances.obs < quantile(thermal.distances.obs)[4])], resource.scarcity.obs[which(thermal.distances.obs < quantile(thermal.distances.obs)[4])], col="orange", cex= 0.9, cex.lab=1.1, pch=20)
points(sqrt(geo.distances.obs)[which(thermal.distances.obs < quantile(thermal.distances.obs)[3])], resource.scarcity.obs[which(thermal.distances.obs < quantile(thermal.distances.obs)[3])], col="red", cex= 0.9, cex.lab=1.1, pch=20)
points(sqrt(geo.distances.obs)[which(thermal.distances.obs < quantile(thermal.distances.obs)[2])], resource.scarcity.obs[which(thermal.distances.obs < quantile(thermal.distances.obs)[2])], col="brown4", cex= 0.9, cex.lab=1.1, pch=20)
mtext("Number of species", side=4, line=1.4, cex=1,las=0)
legend("topright", inset=.005, bg="white", box.col="white", title="Thermal distance", c(paste("<", round(quantile(thermal.distances.obs)[2], 2), sep=" "), paste(round(quantile(thermal.distances.obs)[2], 2), round(quantile(thermal.distances.obs)[3], 2), sep="–"), paste(round(quantile(thermal.distances.obs)[3], 2), round(quantile(thermal.distances.obs)[4], 2), sep="–"), paste(">", round(quantile(thermal.distances.obs)[4], 2), sep=" ")), fill=c("brown4", "red", "orange", "yellow"), cex=0.9)



####  Fig A3: Relationships between the four cost-benefit factors predicted to affect bird migration  ### 

par(mfrow=c(4,2), mar=c(2.9,3,1.5,3), mgp=c(1.5,0.5,0))

# Relationship between geographical distance and thermal distance
plot(sqrt(geo.distances.obs), thermal.distances.obs, xlim=c(0,120), xlab="Geographical distance (square-root)", ylab="Thermal distance", main="", xaxt="n", axes=F, pch=20, cex=0.6, cex.lab=1.1)
axis(side=2)
axis(side=1)
mod = lm(thermal.distances.obs ~ sqrt(geo.distances.obs) + I(sqrt(geo.distances.obs)^2))
lines(sort(sqrt(geo.distances.obs)), fitted(mod)[order(sqrt(geo.distances.obs))], type="l", col="red")
mtext("A", cex=1.3, side=3, line=0, at=-20)
mtext(bquote(R^2 == .(round(summary(mod)$r.squared,2))), cex=1, side=3, line=-1, at=20)

# Relationship between resource scarcity and thermal distance
plot(resource.scarcity.obs, thermal.distances.obs, xlab="Resource scarcity", ylab="Thermal distance", main="", xaxt="n", axes=F, pch=20, cex=0.6, cex.lab=1.1)
axis(side=2)
axis(side=1)
mod = lm(thermal.distances.obs ~ resource.scarcity.obs + I(resource.scarcity.obs^2))
lines(sort(resource.scarcity.obs), fitted(mod)[order(resource.scarcity.obs)], type="l", col="red")
mtext("B", cex=1.3, side=3, line=0, at=-80)
mtext(bquote(R^2 == .(round(summary(mod)$r.squared,2))), cex=1, side=3, line=-1, at=5)

# Relationship between geographical distance and habitat distance
plot(sqrt(geo.distances.obs), habitat.distances.obs, xlim=c(0,120), xlab="Geographical distance (square-root)", ylab="Habitat distance", main="", xaxt="n", axes=F, pch=20, cex=0.6, cex.lab=1.1)
axis(side=2)
axis(side=1)
mod = lm(habitat.distances.obs ~ sqrt(geo.distances.obs) + I(sqrt(geo.distances.obs)^2))
lines(sort(sqrt(geo.distances.obs)), fitted(mod)[order(sqrt(geo.distances.obs))], type="l", col="red")
mtext("C", cex=1.3, side=3, line=0, at=-20)
mtext(bquote(R^2 == .(round(summary(mod)$r.squared,2))), cex=1, side=3, line=-1, at=20)

# Relationship between resource scarcity and habitat distance
plot(resource.scarcity.obs, habitat.distances.obs, xlab="Resource scarcity", ylab="Habitat distance", main="", xaxt="n", axes=F, pch=20, cex=0.6, cex.lab=1.1)
axis(side=2)
axis(side=1)
mod = lm(habitat.distances.obs ~ resource.scarcity.obs + I(resource.scarcity.obs^2))
lines(sort(resource.scarcity.obs), fitted(mod)[order(resource.scarcity.obs)], type="l", col="red")
mtext("D", cex=1.3, side=3, line=0, at=-80)
mtext(bquote(R^2 == .(round(summary(mod)$r.squared,2))), cex=1, side=3, line=-1, at=-45)

# Relationship between geographical distance and resource scarcity
plot(sqrt(geo.distances.obs), resource.scarcity.obs, xlim=c(0,120), xlab="Geographical distance (square-root)", ylab="Resource scarcity", main="", xaxt="n", axes=F, pch=20, cex=0.6, cex.lab=1.1)
axis(side=2)
axis(side=1)
mod = lm(resource.scarcity.obs ~ sqrt(geo.distances.obs) + I(sqrt(geo.distances.obs)^2))
lines(sort(sqrt(geo.distances.obs)), fitted(mod)[order(sqrt(geo.distances.obs))], type="l", col="red")
mtext("E", cex=1.3, side=3, line=0, at=-20)
mtext(bquote(R^2 == .(round(summary(mod)$r.squared,2))), cex=1, side=1, line=-1.5, at=20)

# Relationship between habitat distance and thermal distance
plot(habitat.distances.obs, thermal.distances.obs, xlab="Habitat distance", ylab="Thermal distance", main="", xaxt="n", axes=F, pch=20, cex=0.6, cex.lab=1.1)
axis(side=2)
axis(side=1)
mod = lm(thermal.distances.obs ~ habitat.distances.obs + I(habitat.distances.obs^2))
lines(sort(habitat.distances.obs), fitted(mod)[order(habitat.distances.obs)], type="l", col="red")
mtext("F", cex=1.3, side=3, line=0, at=-0.1)
mtext(bquote(R^2 == .(round(summary(mod)$r.squared,2))), cex=1, side=3, line=-1.5, at=0.45)

# Relationship between three factors (geographical distance, resource scarcity and thermal distance)
hist(sqrt(geo.distances.obs), xlim=c(0,120), xlab="", ylab="", main="", xaxt="n", axes=F, col="light grey", border="grey")
axis(side=4)
axis(side=1)
par(new=T, mar=c(2.9,3,1.5,3), mgp=c(1.5,0.5,0))
plot(sqrt(geo.distances.obs), resource.scarcity.obs, xlim=c(0,120), xlab="Geographical distance (square-root)", ylab="Resource scarcity", main="", xaxt="n", axes=F, pch=20, cex=0.6, cex.lab=1.1, col="yellow")
axis(side=2)
points(sqrt(geo.distances.obs)[which(thermal.distances.obs < quantile(thermal.distances.obs)[4])], resource.scarcity.obs[which(thermal.distances.obs < quantile(thermal.distances.obs)[4])], col="orange", cex=0.6, cex.lab=1.1, pch=20)
points(sqrt(geo.distances.obs)[which(thermal.distances.obs < quantile(thermal.distances.obs)[3])], resource.scarcity.obs[which(thermal.distances.obs < quantile(thermal.distances.obs)[3])], col="red", cex=0.6, cex.lab=1.1, pch=20)
points(sqrt(geo.distances.obs)[which(thermal.distances.obs < quantile(thermal.distances.obs)[2])], resource.scarcity.obs[which(thermal.distances.obs < quantile(thermal.distances.obs)[2])], col="brown4", cex=0.6, cex.lab=1.1, pch=20)
mtext("G", cex=1.3, side=3, line=0, at=-20)
mtext("Number of species", side=4, line=1.4, cex=0.7,las=0)
legend("topright", inset=.005, bg="white", box.col="white", title="Thermal distance", c(paste("<", round(quantile(thermal.distances.obs)[2], 2), sep=" "), paste(round(quantile(thermal.distances.obs)[2], 2), round(quantile(thermal.distances.obs)[3], 2), sep="–"), paste(round(quantile(thermal.distances.obs)[3], 2), round(quantile(thermal.distances.obs)[4], 2), sep="–"), paste(">", round(quantile(thermal.distances.obs)[4], 2), sep=" ")), fill=c("brown4", "red", "orange", "yellow"), cex=0.7)

# Relationship between three factors (geographical distance, resource scarcity and habitat distance)
hist(sqrt(geo.distances.obs), xlim=c(0,120), xlab="", ylab="", main="", xaxt="n", axes=F, col="light grey", border="grey")
axis(side=4)
axis(side=1)
par(new=T, mar=c(2.9,3,1.5,3), mgp=c(1.5,0.5,0))
plot(sqrt(geo.distances.obs), resource.scarcity.obs, xlim=c(0,120), xlab="Geographical distance (square-root)", ylab="Resource scarcity", main="", xaxt="n", axes=F, pch=20, cex=0.6, cex.lab=1.1, col="yellow")
axis(side=2)
points(sqrt(geo.distances.obs)[which(habitat.distances.obs < quantile(habitat.distances.obs)[4])], resource.scarcity.obs[which(habitat.distances.obs < quantile(habitat.distances.obs)[4])], col="orange", cex=0.6, cex.lab=1.1, pch=20)
points(sqrt(geo.distances.obs)[which(habitat.distances.obs < quantile(habitat.distances.obs)[3])], resource.scarcity.obs[which(habitat.distances.obs < quantile(habitat.distances.obs)[3])], col="red", cex=0.6, cex.lab=1.1, pch=20)
points(sqrt(geo.distances.obs)[which(habitat.distances.obs < quantile(habitat.distances.obs)[2])], resource.scarcity.obs[which(habitat.distances.obs < quantile(habitat.distances.obs)[2])], col="brown4", cex=0.6, cex.lab=1.1, pch=20)
mtext("H", cex=1.3, side=3, line=0, at=-20)
mtext("Number of species", side=4, line=1.4, cex=0.7,las=0)
legend("topright", inset=.005, bg="white", box.col="white", title="Habitat distance", c(paste("<", round(quantile(habitat.distances.obs)[2], 2), sep=" "), paste(round(quantile(habitat.distances.obs)[2], 2), round(quantile(habitat.distances.obs)[3], 2), sep="–"), paste(round(quantile(habitat.distances.obs)[3], 2), round(quantile(habitat.distances.obs)[4], 2), sep="–"), paste(">", round(quantile(habitat.distances.obs)[4], 2), sep=" ")), fill=c("brown4", "red", "orange", "yellow"), cex=0.7)







################################################################
################################################################

#######      Simulate realistic geographical ranges     ########

################################################################
################################################################



## Compute the great circle distance between each pair of hexagons on the grid

# Western Hemisphere
pairwise.distance.WH <- rdist.earth(west_Hem, miles=F)
diag(pairwise.distance.WH) <- 0
pairwise.distance.WH.01 <- apply(pairwise.distance.WH, 2, function(x) ifelse(x<250 & x>0, 1, 0))
neig.list_WH <- apply(pairwise.distance.WH.01, 2, function(x) which(x==1))

# Eastern Hemisphere
pairwise.distance.EH <- rdist.earth(east_Hem, miles=F)
diag(pairwise.distance.EH) <- 0
pairwise.distance.EH.01 <- apply(pairwise.distance.EH, 2, function(x) ifelse(x<250 & x>0, 1, 0))
neig.list_EH <- apply(pairwise.distance.EH.01, 2, function(x) which(x==1))


## Compute the bearing of the rhumb line between every pairs of hexagons and only keep the neighbouring hexagons

pairwise.bearingsWH = t(apply(west_Hem, 1, function(x) bearingRhumb(x, west_Hem)))
pairwise.bearingsEH = t(apply(east_Hem, 1, function(x) bearingRhumb(x, east_Hem)))

pairwise.WH.01.bearings <- pairwise.distance.WH.01 * pairwise.bearingsWH
pairwise.WH.01.bearings[pairwise.WH.01.bearings==0] <- NA
pairwise.EH.01.bearings <- pairwise.distance.EH.01 * pairwise.bearingsEH
pairwise.EH.01.bearings[pairwise.EH.01.bearings==0] <- NA


## Probablity density function based on the bearing values

proba_bearing <- function(x, blockBearings){
	if(is.na(x) == T){prob <- 0}else{
		if(x>=0 & x<30){prob <- length(which(blockBearings>=0 & blockBearings<30))/length(blockBearings)}
		if(x>=30 & x<60){prob <- length(which(blockBearings>=30 & blockBearings<60))/length(blockBearings)}
		if(x>=60 & x<90){prob <- length(which(blockBearings>=60 & blockBearings<90))/length(blockBearings)}
		if(x>=90 & x<120){prob <- length(which(blockBearings>=90 & blockBearings<120))/length(blockBearings)}
		if(x>=120 & x<150){prob <- length(which(blockBearings>=120 & blockBearings<150))/length(blockBearings)}
		if(x>=150 & x<180){prob <- length(which(blockBearings>=150 & blockBearings<180))/length(blockBearings)}
		if(x>=180 & x<210){prob <- length(which(blockBearings>=180 & blockBearings<210))/length(blockBearings)}
		if(x>=210 & x<240){prob <- length(which(blockBearings>=210 & blockBearings<240))/length(blockBearings)}
		if(x>=240 & x<270){prob <- length(which(blockBearings>=240 & blockBearings<270))/length(blockBearings)}
		if(x>=270 & x<300){prob <- length(which(blockBearings>=270 & blockBearings<300))/length(blockBearings)}
		if(x>=300 & x<330){prob <- length(which(blockBearings>=300 & blockBearings<330))/length(blockBearings)}
		if(x>=330 & x<360){prob <- length(which(blockBearings>=330 & blockBearings<360))/length(blockBearings)}
	}
	return(prob)
}


## Function to simulate n realistic geographical ranges	

rangeSimulations <- function(presencesAbsences, latlong, neigList, pairwiseDists, pairwise01Bearings, n){
	#Range features
	occupied <- which(presencesAbsences == 1)
	blocks <- clusters(graph_from_adjacency_matrix(pairwise01Bearings[occupied, occupied]))
	blocks.number <- blocks$no
	blocks.id <- (1:blocks$no)[order(blocks$csize, decreasing=T)]
	blocks.sizes <- blocks$csize[order(blocks$csize, decreasing=T)]
	blocks.centroids <- list(); blocks.bearings <- list()
	for(i in 1:blocks.number){
		if(blocks.sizes[i] > 1){
			blocks.centroids[[i]] <- apply(latlong[occupied[which(blocks$membership == blocks.id[i])],], 2, mean)
		}else{
			blocks.centroids[[i]] <- latlong[occupied[which(blocks$membership == blocks.id[i])],]
		}
		blocks.bearings[[i]] <- bearingRhumb(blocks.centroids[[i]], latlong[occupied[which(blocks$membership == blocks.id[i])],])
	}
	blocks.dists <- rdist.earth(matrix(unlist(blocks.centroids), ncol=2, byrow=T), miles=F)[1,]	
	#Simulate range
	ranges.sim <- list()
	for(j in 1:n){	
		central.hexagon <- sample(1:length(pairwise01Bearings[1,]),1)
		range.simulated <- central.hexagon		
		for(k in 1:blocks.number){
			hex_possible <- order(abs(pairwiseDists[central.hexagon,] - blocks.dists[k]))
			hex_possible <- hex_possible[which(match(hex_possible, range.simulated, nomatch=0)==0)][1:20]
			block.simulated <- sample(hex_possible, 1)
			while(length(block.simulated) < blocks.sizes[k]){
				neigs <- unique(unlist(neigList[block.simulated]))
				neigs <- neigs[which(match(neigs, c(block.simulated,range.simulated), nomatch=0)==0)]			
				if(length(neigs)>1){
					neigs.bearings <- pairwise01Bearings[block.simulated, neigs]
					if(length(block.simulated)>1){
						neigs.proba <- apply(neigs.bearings, c(1,2), function(x) proba_bearing(x,blocks.bearings[[k]]))
						neigs.proba <- apply(neigs.proba, 2, sum)	
					}else{
					neigs.proba <- sapply(neigs.bearings, function(x) proba_bearing(x,blocks.bearings[[k]])) 
					}
					if(sum(neigs.proba) == 0){ neigs.proba <- rep(1,length(neigs.proba)) }
					neigs.proba <- neigs.proba / sum(neigs.proba)
					if(length(which(neigs.proba > 0)) >= ceiling(length(neigs)/4)){
						block.simulated <- c(block.simulated, sample(neigs, ceiling(length(neigs)/4), prob=neigs.proba))
					}else{
						block.simulated <- c(block.simulated, neigs[which(neigs.proba>0)])
					}
				}
				if(length(neigs)==1){ block.simulated <- c(block.simulated, neigs) }
				if(length(neigs)==0){ block.simulated <- c(block.simulated, (1:length(neigList))[-c(block.simulated,range.simulated)][order(pairwiseDists[block.simulated[1],-c(block.simulated,range.simulated)])[1]]) }
			}
			block.simulated <- block.simulated[1:blocks.sizes[k]]
			range.simulated <- unique(c(range.simulated, block.simulated))	
		}
		ranges.sim[[j]] <- range.simulated	
	}
	return(ranges.sim)
}
	
	
## Simulate 100 realistic geographical ranges for each observed seasonal range

simulated.ranges_BR_NH_WH <- apply(PresAbs_BR_NH_WH2, 2, function(x) rangeSimulations(x, west_Hem, neig.list_WH, pairwise.distance.WH, pairwise.WH.01.bearings, 100))
simulated.ranges_NB_NH_WH <- apply(PresAbs_NB_NH_WH2, 2, function(x) rangeSimulations(x, west_Hem, neig.list_WH, pairwise.distance.WH, pairwise.WH.01.bearings, 100))
simulated.ranges_BR_SH_WH <- apply(PresAbs_BR_SH_WH2, 2, function(x) rangeSimulations(x, west_Hem, neig.list_WH, pairwise.distance.WH, pairwise.WH.01.bearings, 100))
simulated.ranges_NB_SH_WH <- apply(PresAbs_NB_SH_WH2, 2, function(x) rangeSimulations(x, west_Hem, neig.list_WH, pairwise.distance.WH, pairwise.WH.01.bearings, 100))
simulated.ranges_BR_NH_EH <- apply(PresAbs_BR_NH_EH2, 2, function(x) rangeSimulations(x, east_Hem, neig.list_EH, pairwise.distance.EH, pairwise.EH.01.bearings, 100))
simulated.ranges_NB_NH_EH <- apply(PresAbs_NB_NH_EH2, 2, function(x) rangeSimulations(x, east_Hem, neig.list_EH, pairwise.distance.EH, pairwise.EH.01.bearings, 100))
simulated.ranges_BR_SH_EH <- apply(PresAbs_BR_SH_EH2, 2, function(x) rangeSimulations(x, east_Hem, neig.list_EH, pairwise.distance.EH, pairwise.EH.01.bearings, 100))
simulated.ranges_NB_SH_EH <- apply(PresAbs_NB_SH_EH2, 2, function(x) rangeSimulations(x, east_Hem, neig.list_EH, pairwise.distance.EH, pairwise.EH.01.bearings, 100))







################################################################
################################################################

######      Compare to alternative migration options    ########

################################################################
################################################################



#load("simuranges.RData")

## Compute thermal distances for alternative migration options

thermal.distances.simulatedNH_WH <- list()
for(i in 1:length(nonbreeding.thermal.nichesNH_WH)){
	breeding.thermal.niches.simulated_NH_WH <- lapply(simulated.ranges_BR_NH_WH[[i]], function(x) nicheDensityRaster(cbind(TempMeanNS_WH[x], TempMinNS_WH[x])))
	thermal.distances.simulatedBR_NH_WH <- mapply(function(X){ emd(X,nonbreeding.thermal.nichesNH_WH[[i]]) }, X=breeding.thermal.niches.simulated_NH_WH)
	nonbreeding.thermal.niches.simulated_NH_WH <- lapply(simulated.ranges_NB_NH_WH[[i]], function(x) nicheDensityRaster(cbind(TempMeanNW_WH[x], TempMinNW_WH[x])))
	thermal.distances.simulatedNB_NH_WH <- mapply(function(Y){ emd(breeding.thermal.nichesNH_WH[[i]],Y) }, Y=nonbreeding.thermal.niches.simulated_NH_WH)
	thermal.distances.simulatedNH_WH[[i]] <- c(unlist(thermal.distances.simulatedBR_NH_WH), unlist(thermal.distances.simulatedNB_NH_WH))
}
thermal.distances.simulatedSH_WH <- list()
for(i in 1:length(nonbreeding.thermal.nichesSH_WH)){
	breeding.thermal.niches.simulated_SH_WH <- lapply(simulated.ranges_BR_SH_WH[[i]], function(x) nicheDensityRaster(cbind(TempMeanNW_WH[x], TempMinNW_WH[x])))
	thermal.distances.simulatedBR_SH_WH <- mapply(function(X){ emd(X,nonbreeding.thermal.nichesSH_WH[[i]]) }, X=breeding.thermal.niches.simulated_SH_WH)
	nonbreeding.thermal.niches.simulated_SH_WH <- lapply(simulated.ranges_NB_SH_WH[[i]], function(x) nicheDensityRaster(cbind(TempMeanNS_WH[x], TempMinNS_WH[x])))
	thermal.distances.simulatedNB_SH_WH <- mapply(function(Y){ emd(breeding.thermal.nichesSH_WH[[i]],Y) }, Y=nonbreeding.thermal.niches.simulated_NH_WH)
	thermal.distances.simulatedSH_WH[[i]] <- c(unlist(thermal.distances.simulatedBR_SH_WH), unlist(thermal.distances.simulatedNB_SH_WH))
}
thermal.distances.simulatedNH_EH <- list()
for(i in 1:length(nonbreeding.thermal.nichesNH_EH)){
	breeding.thermal.niches.simulated_NH_EH <- lapply(simulated.ranges_BR_NH_EH[[i]], function(x) nicheDensityRaster(cbind(TempMeanNS_EH[x], TempMinNS_EH[x])))
	thermal.distances.simulatedBR_NH_EH <- mapply(function(X){ emd(X,nonbreeding.thermal.nichesNH_EH[[i]]) }, X=breeding.thermal.niches.simulated_NH_EH)
	nonbreeding.thermal.niches.simulated_NH_EH <- lapply(simulated.ranges_NB_NH_EH[[i]], function(x) nicheDensityRaster(cbind(TempMeanNW_EH[x], TempMinNW_EH[x])))
	thermal.distances.simulatedNB_NH_EH <- mapply(function(Y){ emd(breeding.thermal.nichesNH_EH[[i]],Y) }, Y=nonbreeding.thermal.niches.simulated_NH_EH)
	thermal.distances.simulatedNH_EH[[i]] <- c(unlist(thermal.distances.simulatedBR_NH_EH), unlist(thermal.distances.simulatedNB_NH_EH))
}
thermal.distances.simulatedSH_EH <- list()
for(i in 1:length(nonbreeding.thermal.nichesSH_EH)){
	breeding.thermal.niches.simulated_SH_EH <- lapply(simulated.ranges_BR_SH_EH[[i]], function(x) nicheDensityRaster(cbind(TempMeanNW_EH[x], TempMinNW_EH[x])))
	thermal.distances.simulatedBR_SH_EH <- mapply(function(X){ emd(X,nonbreeding.thermal.nichesSH_EH[[i]]) }, X=breeding.thermal.niches.simulated_SH_EH)
	nonbreeding.thermal.niches.simulated_SH_EH <- lapply(simulated.ranges_NB_SH_EH[[i]], function(x) nicheDensityRaster(cbind(TempMeanNS_EH[x], TempMinNS_EH[x])))
	thermal.distances.simulatedNB_SH_EH <- mapply(function(Y){ emd(breeding.thermal.nichesSH_EH[[i]],Y) }, Y=nonbreeding.thermal.niches.simulated_NH_EH)
	thermal.distances.simulatedSH_EH[[i]] <- c(unlist(thermal.distances.simulatedBR_SH_EH), unlist(thermal.distances.simulatedNB_SH_EH))
}
thermal.distances.simulated <- c(thermal.distances.simulatedNH_WH, thermal.distances.simulatedNH_EH, thermal.distances.simulatedSH_WH, thermal.distances.simulatedSH_EH)



## Compute habitat distances for alternative migration options

habitat.distances.simulatedNH_WH <- list()
for(i in 1:length(breeding.habitatNH_WH)){
	breeding.habitat.simulated_NH_WH <- lapply(simulated.ranges_BR_NH_WH[[i]], function(x) mean(habitat[x]))
	habitat.distances.simulatedBR_NH_WH <- abs(unlist(breeding.habitat.simulated_NH_WH) - nonbreeding.habitatNH_WH[i])
	nonbreeding.habitat.simulated_NH_WH <- lapply(simulated.ranges_NB_NH_WH[[i]], function(x) mean(habitat[x]))
	habitat.distances.simulatedNB_NH_WH <- abs(breeding.habitatNH_WH[i] - unlist(nonbreeding.habitat.simulated_NH_WH))
	habitat.distances.simulatedNH_WH[[i]] <- c(unlist(habitat.distances.simulatedBR_NH_WH), unlist(habitat.distances.simulatedNB_NH_WH))
}
habitat.distances.simulatedSH_WH <- list()
for(i in 1:length(breeding.habitatSH_WH)){
	breeding.habitat.simulated_SH_WH <- lapply(simulated.ranges_BR_SH_WH[[i]], function(x) mean(habitat[x]))
	habitat.distances.simulatedBR_SH_WH <- abs(unlist(breeding.habitat.simulated_SH_WH) - nonbreeding.habitatNH_WH[i])
	nonbreeding.habitat.simulated_SH_WH <- lapply(simulated.ranges_NB_SH_WH[[i]], function(x) mean(habitat[x]))
	habitat.distances.simulatedNB_SH_WH <- abs(breeding.habitatSH_WH[i] - unlist(nonbreeding.habitat.simulated_SH_WH))
	habitat.distances.simulatedSH_WH[[i]] <- c(unlist(habitat.distances.simulatedBR_SH_WH), unlist(habitat.distances.simulatedNB_SH_WH))
}
habitat.distances.simulatedNH_EH <- list()
for(i in 1:length(breeding.habitatNH_EH)){
	breeding.habitat.simulated_NH_EH <- lapply(simulated.ranges_BR_NH_EH[[i]], function(x) mean(habitat[x]))
	habitat.distances.simulatedBR_NH_EH <- abs(unlist(breeding.habitat.simulated_NH_EH) - nonbreeding.habitatNH_EH[i])
	nonbreeding.habitat.simulated_NH_EH <- lapply(simulated.ranges_NB_NH_EH[[i]], function(x) mean(habitat[x]))
	habitat.distances.simulatedNB_NH_EH <- abs(breeding.habitatNH_EH[i] - unlist(nonbreeding.habitat.simulated_NH_EH))
	habitat.distances.simulatedNH_EH[[i]] <- c(unlist(habitat.distances.simulatedBR_NH_EH), unlist(habitat.distances.simulatedNB_NH_EH))
}
habitat.distances.simulatedSH_EH <- list()
for(i in 1:length(breeding.habitatSH_EH)){
	breeding.habitat.simulated_SH_EH <- lapply(simulated.ranges_BR_SH_EH[[i]], function(x) mean(habitat[x]))
	habitat.distances.simulatedBR_SH_EH <- abs(unlist(breeding.habitat.simulated_SH_EH) - nonbreeding.habitatNH_EH[i])
	nonbreeding.habitat.simulated_SH_EH <- lapply(simulated.ranges_NB_SH_EH[[i]], function(x) mean(habitat[x]))
	habitat.distances.simulatedNB_SH_EH <- abs(breeding.habitatSH_EH[i] - unlist(nonbreeding.habitat.simulated_SH_EH))
	habitat.distances.simulatedSH_EH[[i]] <- c(unlist(habitat.distances.simulatedBR_SH_EH), unlist(habitat.distances.simulatedNB_SH_EH))
}
habitat.distances.simulated <- c(habitat.distances.simulatedNH_WH, habitat.distances.simulatedNH_EH, habitat.distances.simulatedSH_WH, habitat.distances.simulatedSH_EH)



## Compute geographical (migration) distance for alternative migration options

centroids.simulatedBR_NH_WH <- lapply(simulated.ranges_BR_NH_WH, function(x) lapply(x, function(x2) apply(west_Hem[x2,], 2, mean)))
centroids.simulatedNB_NH_WH <- lapply(simulated.ranges_NB_NH_WH, function(x) lapply(x, function(x2) apply(west_Hem[x2,], 2, mean)))
geo.distances.simulatedNH_WH <- list()
for(i in 1:length(centroids.simulatedBR_NH_WH)){
	geo.distances.simulatedBR_NH_WH <- rdist.earth(matrix(unlist(centroids.simulatedBR_NH_WH[[i]]), nrow=100, byrow=T), t(as.matrix(centroids.nonbreeding.groundsNH_WH[i,])), miles=F)
	geo.distances.simulatedNB_NH_WH <- rdist.earth(matrix(unlist(centroids.simulatedNB_NH_WH[[i]]), nrow=100, byrow=T), t(as.matrix(centroids.breeding.groundsNH_WH[i,])), miles=F)
	geo.distances.simulatedNH_WH[[i]] <- c(unlist(geo.distances.simulatedBR_NH_WH), unlist(geo.distances.simulatedNB_NH_WH))
}
centroids.simulatedBR_SH_WH <- lapply(simulated.ranges_BR_SH_WH, function(x) lapply(x, function(x2) apply(west_Hem[x2,], 2, mean)))
centroids.simulatedNB_SH_WH <- lapply(simulated.ranges_NB_SH_WH, function(x) lapply(x, function(x2) apply(west_Hem[x2,], 2, mean)))
geo.distances.simulatedSH_WH <- list()
for(i in 1:length(centroids.simulatedBR_SH_WH)){
	geo.distances.simulatedBR_SH_WH <- rdist.earth(matrix(unlist(centroids.simulatedBR_SH_WH[[i]]), nrow=100, byrow=T), t(as.matrix(centroids.nonbreeding.groundsSH_WH[i,])), miles=F)
	geo.distances.simulatedNB_SH_WH <- rdist.earth(matrix(unlist(centroids.simulatedNB_SH_WH[[i]]), nrow=100, byrow=T), t(as.matrix(centroids.breeding.groundsSH_WH[i,])), miles=F)
	geo.distances.simulatedSH_WH[[i]] <- c(unlist(geo.distances.simulatedBR_SH_WH), unlist(geo.distances.simulatedNB_SH_WH))
}
centroids.simulatedBR_NH_EH <- lapply(simulated.ranges_BR_NH_EH, function(x) lapply(x, function(x2) apply(east_Hem[x2,], 2, mean)))
centroids.simulatedNB_NH_EH <- lapply(simulated.ranges_NB_NH_EH, function(x) lapply(x, function(x2) apply(east_Hem[x2,], 2, mean)))
geo.distances.simulatedNH_EH <- list()
for(i in 1:length(centroids.simulatedBR_NH_EH)){
	geo.distances.simulatedBR_NH_EH <- rdist.earth(matrix(unlist(centroids.simulatedBR_NH_EH[[i]]), nrow=100, byrow=T), t(as.matrix(centroids.nonbreeding.groundsNH_EH[i,])), miles=F)
	geo.distances.simulatedNB_NH_EH <- rdist.earth(matrix(unlist(centroids.simulatedNB_NH_EH[[i]]), nrow=100, byrow=T), t(as.matrix(centroids.breeding.groundsNH_EH[i,])), miles=F)
	geo.distances.simulatedNH_EH[[i]] <- c(unlist(geo.distances.simulatedBR_NH_EH), unlist(geo.distances.simulatedNB_NH_EH))
}
centroids.simulatedBR_SH_EH <- lapply(simulated.ranges_BR_SH_EH, function(x) lapply(x, function(x2) apply(east_Hem[x2,], 2, mean)))
centroids.simulatedNB_SH_EH <- lapply(simulated.ranges_NB_SH_EH, function(x) lapply(x, function(x2) apply(east_Hem[x2,], 2, mean)))
geo.distances.simulatedSH_EH <- list()
for(i in 1:length(centroids.simulatedBR_SH_EH)){
	geo.distances.simulatedBR_SH_EH <- rdist.earth(matrix(unlist(centroids.simulatedBR_SH_EH[[i]]), nrow=100, byrow=T), t(as.matrix(centroids.nonbreeding.groundsSH_EH[i,])), miles=F)
	geo.distances.simulatedNB_SH_EH <- rdist.earth(matrix(unlist(centroids.simulatedNB_SH_EH[[i]]), nrow=100, byrow=T), t(as.matrix(centroids.breeding.groundsSH_EH[i,])), miles=F)
	geo.distances.simulatedSH_EH[[i]] <- c(unlist(geo.distances.simulatedBR_SH_EH), unlist(geo.distances.simulatedNB_SH_EH))
}
geo.distances.simulated <- c(geo.distances.simulatedNH_WH, geo.distances.simulatedNH_EH, geo.distances.simulatedSH_WH, geo.distances.simulatedSH_EH)



## Compute resource scarcity for alternative migration options

resource.scarcity.simulatedBR_NH_WH <- lapply(simulated.ranges_BR_NH_WH, function(x) lapply(x, function(x2) mean(resourceScarcity_NS_WH[x2])))
resource.scarcity.simulatedNB_NH_WH <- lapply(simulated.ranges_NB_NH_WH, function(x) lapply(x, function(x2) mean(resourceScarcity_NW_WH[x2])))
resource.scarcity.simulatedNH_WH <- list()
for(i in 1:length(resource.scarcity.simulatedBR_NH_WH)){
	resource.scarcity.simuBR_NH_WH <- unlist(resource.scarcity.simulatedBR_NH_WH[[i]]) + nonbreeding.resource.scarcityNH_WH[i]
	resource.scarcity.simuNB_NH_WH <- unlist(resource.scarcity.simulatedNB_NH_WH[[i]]) + breeding.resource.scarcityNH_WH[i]
	resource.scarcity.simulatedNH_WH[[i]] <- c(resource.scarcity.simuBR_NH_WH, resource.scarcity.simuNB_NH_WH)
}
resource.scarcity.simulatedBR_SH_WH <- lapply(simulated.ranges_BR_SH_WH, function(x) lapply(x, function(x2) mean(resourceScarcity_NW_WH[x2])))
resource.scarcity.simulatedNB_SH_WH <- lapply(simulated.ranges_NB_SH_WH, function(x) lapply(x, function(x2) mean(resourceScarcity_NS_WH[x2])))
resource.scarcity.simulatedSH_WH <- list()
for(i in 1:length(resource.scarcity.simulatedBR_SH_WH)){
	resource.scarcity.simuBR_SH_WH <- unlist(resource.scarcity.simulatedBR_SH_WH[[i]]) + nonbreeding.resource.scarcitySH_WH[i]
	resource.scarcity.simuNB_SH_WH <- unlist(resource.scarcity.simulatedNB_SH_WH[[i]]) + breeding.resource.scarcitySH_WH[i]
	resource.scarcity.simulatedSH_WH[[i]] <- c(resource.scarcity.simuBR_SH_WH, resource.scarcity.simuNB_SH_WH)
}
resource.scarcity.simulatedBR_NH_EH <- lapply(simulated.ranges_BR_NH_EH, function(x) lapply(x, function(x2) mean(resourceScarcity_NS_EH[x2])))
resource.scarcity.simulatedNB_NH_EH <- lapply(simulated.ranges_NB_NH_EH, function(x) lapply(x, function(x2) mean(resourceScarcity_NW_EH[x2])))
resource.scarcity.simulatedNH_EH <- list()
for(i in 1:length(resource.scarcity.simulatedBR_NH_EH)){
	resource.scarcity.simuBR_NH_EH <- unlist(resource.scarcity.simulatedBR_NH_EH[[i]]) + nonbreeding.resource.scarcityNH_EH[i]
	resource.scarcity.simuNB_NH_EH <- unlist(resource.scarcity.simulatedNB_NH_EH[[i]]) + breeding.resource.scarcityNH_EH[i]
	resource.scarcity.simulatedNH_EH[[i]] <- c(resource.scarcity.simuBR_NH_EH, resource.scarcity.simuNB_NH_EH)
}
resource.scarcity.simulatedBR_SH_EH <- lapply(simulated.ranges_BR_SH_EH, function(x) lapply(x, function(x2) mean(resourceScarcity_NW_EH[x2])))
resource.scarcity.simulatedNB_SH_EH <- lapply(simulated.ranges_NB_SH_EH, function(x) lapply(x, function(x2) mean(resourceScarcity_NS_EH[x2])))
resource.scarcity.simulatedSH_EH <- list()
for(i in 1:length(resource.scarcity.simulatedBR_SH_EH)){
	resource.scarcity.simuBR_SH_EH <- unlist(resource.scarcity.simulatedBR_SH_EH[[i]]) + nonbreeding.resource.scarcitySH_EH[i]
	resource.scarcity.simuNB_SH_EH <- unlist(resource.scarcity.simulatedNB_SH_EH[[i]]) + breeding.resource.scarcitySH_EH[i]
	resource.scarcity.simulatedSH_EH[[i]] <- c(resource.scarcity.simuBR_SH_EH, resource.scarcity.simuNB_SH_EH)
}
resource.scarcity.simulated <- c(resource.scarcity.simulatedNH_WH, resource.scarcity.simulatedNH_EH, resource.scarcity.simulatedSH_WH, resource.scarcity.simulatedSH_EH)



## Compute the rank of each observed migration among 200 migration alternatives

# Considering thermal distance only
ranks_thermal <- vector()
for(i in 1:length(thermal.distances.simulated)){
	ranks_thermal[i] <- length(which(thermal.distances.simulated[[i]] < thermal.distances.obs[i])) / length(thermal.distances.simulated[[i]])	
}

# Considering habitat distance only
ranks_habitat <- vector()
for(i in 1:length(habitat.distances.simulated)){
	ranks_habitat[i] <- length(which(habitat.distances.simulated[[i]] < habitat.distances.obs[i])) / length(habitat.distances.simulated[[i]])	
}

# Considering geographical distance only
ranks_geo <- vector()
for(i in 1:length(geo.distances.simulated)){
	ranks_geo[i] <- length(which(geo.distances.simulated[[i]] < geo.distances.obs[i])) / length(geo.distances.simulated[[i]])	
}

# Considering resource scarcity only
ranks_resource <- vector()
for(i in 1:length(resource.scarcity.simulated)){
	ranks_resource[i] <- length(which(resource.scarcity.simulated[[i]] < resource.scarcity.obs[i])) / length(resource.scarcity.simulated[[i]])	
}

# Considering geographical distance and thermal distance
ranks_geo_thermal <- vector()
for(i in 1:length(geo.distances.simulated)){
	geo_thermal <- scale(c(geo.distances.obs[i], geo.distances.simulated[[i]])) + scale(c(thermal.distances.obs[i], thermal.distances.simulated[[i]]))
	ranks_geo_thermal[i] <- length(which(geo_thermal[-1] < geo_thermal[1])) / length(geo_thermal[-1])	
}

# Considering geographical distance and resource scarcity
ranks_geo_resource <- vector()
for(i in 1:length(geo.distances.simulated)){
	geo_resource <- scale(c(geo.distances.obs[i], geo.distances.simulated[[i]])) + scale(c(resource.scarcity.obs[i], resource.scarcity.simulated[[i]]))
	ranks_geo_resource[i] <- length(which(geo_resource[-1] < geo_resource[1])) / length(geo_resource[-1])	
}

# Considering thermal distance and resource scarcity
ranks_thermal_resource <- vector()
for(i in 1:length(thermal.distances.simulated)){
	thermal_resource <- scale(c(thermal.distances.obs[i], thermal.distances.simulated[[i]])) + scale(c(resource.scarcity.obs[i], resource.scarcity.simulated[[i]]))
	ranks_thermal_resource[i] <- length(which(thermal_resource[-1] < thermal_resource[1])) / length(thermal_resource[-1])	
}

# Considering geographical distance and habitat distance
ranks_geo_habitat <- vector()
for(i in 1:length(geo.distances.simulated)){
	geo_habitat <- scale(c(geo.distances.obs[i], geo.distances.simulated[[i]])) + scale(c(habitat.distances.obs[i], habitat.distances.simulated[[i]]))
	ranks_geo_habitat[i] <- length(which(geo_habitat[-1] < geo_habitat[1])) / length(geo_habitat[-1])	
}

# Considering thermal distance and habitat distance
ranks_thermal_habitat <- vector()
for(i in 1:length(thermal.distances.simulated)){
	thermal_habitat <- scale(c(thermal.distances.obs[i], thermal.distances.simulated[[i]])) + scale(c(habitat.distances.obs[i], habitat.distances.simulated[[i]]))
	ranks_thermal_habitat[i] <- length(which(thermal_habitat[-1] < thermal_habitat[1])) / length(thermal_habitat[-1])	
}

# Considering habitat distance and resource scarcity
ranks_habitat_resource <- vector()
for(i in 1:length(habitat.distances.simulated)){
	habitat_resource <- scale(c(habitat.distances.obs[i], habitat.distances.simulated[[i]])) + scale(c(resource.scarcity.obs[i], resource.scarcity.simulated[[i]]))
	ranks_habitat_resource[i] <- length(which(habitat_resource[-1] < habitat_resource[1])) / length(habitat_resource[-1])	
}

# Considering geographical distance, thermal distance and resource scarcity
ranks_geo_thermal_resource <- vector()
for(i in 1:length(geo.distances.simulated)){
	geo_thermal_resource <- scale(c(geo.distances.obs[i], geo.distances.simulated[[i]])) + scale(c(thermal.distances.obs[i], thermal.distances.simulated[[i]])) + scale(c(resource.scarcity.obs[i], resource.scarcity.simulated[[i]]))
	ranks_geo_thermal_resource[i] <- length(which(geo_thermal_resource[-1] < geo_thermal_resource[1])) / length(geo_thermal_resource[-1])	
}

# Considering geographical distance, habitat distance and resource scarcity
ranks_geo_habitat_resource <- vector()
for(i in 1:length(geo.distances.simulated)){
	geo_habitat_resource <- scale(c(geo.distances.obs[i], geo.distances.simulated[[i]])) + scale(c(habitat.distances.obs[i], habitat.distances.simulated[[i]])) + scale(c(resource.scarcity.obs[i], resource.scarcity.simulated[[i]]))
	ranks_geo_habitat_resource[i] <- length(which(geo_habitat_resource[-1] < geo_habitat_resource[1])) / length(geo_habitat_resource[-1])	
}

# Considering geographical distance, thermal distance and habitat distance
ranks_geo_thermal_habitat <- vector()
for(i in 1:length(geo.distances.simulated)){
	geo_thermal_habitat <- scale(c(geo.distances.obs[i], geo.distances.simulated[[i]])) + scale(c(thermal.distances.obs[i], thermal.distances.simulated[[i]])) + scale(c(habitat.distances.obs[i], habitat.distances.simulated[[i]]))
	ranks_geo_thermal_habitat[i] <- length(which(geo_thermal_habitat[-1] < geo_thermal_habitat[1])) / length(geo_thermal_habitat[-1])	
}

# Considering thermal distance, habitat distance and resource scarcity
ranks_thermal_habitat_resource <- vector()
for(i in 1:length(habitat.distances.simulated)){
	habitat_thermal_resource <- scale(c(habitat.distances.obs[i], habitat.distances.simulated[[i]])) + scale(c(thermal.distances.obs[i], thermal.distances.simulated[[i]])) + scale(c(resource.scarcity.obs[i], resource.scarcity.simulated[[i]]))
	ranks_thermal_habitat_resource[i] <- length(which(habitat_thermal_resource[-1] < habitat_thermal_resource[1])) / length(habitat_thermal_resource[-1])	
}

# Considering geographical distance, thermal distance, habitat distance and resource scarcity
ranks_geo_thermal_habitat_resource <- vector()
for(i in 1:length(geo.distances.simulated)){
	geo_thermal_habitat_resource <- scale(c(geo.distances.obs[i], geo.distances.simulated[[i]])) + scale(c(thermal.distances.obs[i], thermal.distances.simulated[[i]])) + scale(c(habitat.distances.obs[i], habitat.distances.simulated[[i]])) + scale(c(resource.scarcity.obs[i], resource.scarcity.simulated[[i]]))
	ranks_geo_thermal_habitat_resource[i] <- length(which(geo_thermal_habitat_resource[-1] < geo_thermal_habitat_resource[1])) / length(geo_thermal_habitat_resource[-1])	
}




####  Figure 2 - Frequency distributions of the scaled ranks of each species' observed migration among alternative simulated options  ####

par(mfrow=c(5,3), mar=c(3.5,2.5,1,0.5), mgp=c(0.5,0.5,0))
hist(ranks_geo, main="", ylab="", xlab="", xlim=c(0,1), cex.lab=1, cex.axis=0.8, col="grey", border="white", breaks=10)
mtext("Number of species", cex=0.65, side=2, line=1.4, at=80)
mtext("Scaled rank for geographical distance", cex=0.65, side=1, line=1.7, at=0.5)
mtext("A", cex=1.2, side=3, line=-0.5, at=0.5)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-4, at=0.8)
hist(ranks_thermal, main="", ylab="", xlab="", xlim=c(0,1), cex.lab=1, cex.axis=0.8, col="grey", border="white", breaks=10)
mtext("Number of species", cex=0.65, side=2, line=1.4, at=70)
mtext("Scaled rank for thermal distance", cex=0.65, side=1, line=1.7, at=0.5)
mtext("B", cex=1.2, side=3, line=-0.5, at=0.5)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-4, at=0.8)
hist(ranks_habitat, main="", ylab="", xlab="", xlim=c(0,1), cex.lab=1, cex.axis=0.8, col="grey", border="white", breaks=10)
mtext("Number of species", cex=0.65, side=2, line=1.4, at=40)
mtext("Scaled rank for habitat distance", cex=0.65, side=1, line=1.7, at=0.5)
mtext("C", cex=1.2, side=3, line=-0.5, at=0.5)
abline(a=65.2, b=0)
mtext("P = 0.0002", cex=1, side=3, line=-2, at=0.2)
hist(ranks_resource, main="", ylab="", xlab="", xlim=c(0,1), cex.lab=1, cex.axis=0.8, col="grey", border="white", breaks=10)
mtext("Number of species", cex=0.65, side=2, line=1.4, at=80)
mtext("Scaled rank for resource scarcity", cex=0.65, side=1, line=1.7, at=0.5)
mtext("D", cex=1.2, side=3, line=-0.5, at=0.5)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-4, at=0.8)
hist(ranks_geo_thermal, main="", ylab="", xlab="", xlim=c(0,1), cex.lab=1, cex.axis=0.8, col="grey", border="white", breaks=10)
mtext("Number of species", cex=0.65, side=2, line=1.4, at=130)
mtext("Scaled rank for geographical distance\n+ thermal distance", cex=0.65, side=1, line=2.2, at=0.5)
mtext("E", cex=1.2, side=3, line=-0.5, at=0.5)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-5, at=0.8)
hist(ranks_geo_habitat, main="", ylab="", xlab="", xlim=c(0,1), cex.lab=1, cex.axis=0.8, col="grey", border="white", breaks=10)
mtext("Number of species", cex=0.65, side=2, line=1.4, at=60)
mtext("Scaled rank for geographical distance\n+ habitat distance", cex=0.65, side=1, line=2.2, at=0.5)
mtext("F", cex=1.2, side=3, line=-0.5, at=0.5)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-4, at=0.8)
hist(ranks_geo_resource, main="", ylab="", xlab="", xlim=c(0,1), cex.lab=1, cex.axis=0.8, col="grey", border="white", breaks=10)
mtext("Number of species", cex=0.65, side=2, line=1.4, at=110)
mtext("Scaled rank for geographical distance\n+ resource scarcity", cex=0.65, side=1, line=2.2, at=0.5)
mtext("G", cex=1.2, side=3, line=-0.5, at=0.5)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-5, at=0.8)
hist(ranks_thermal_habitat, main="", ylab="", xlab="", xlim=c(0,1), cex.lab=1, cex.axis=0.8, col="grey", border="white", breaks=10)
mtext("Number of species", cex=0.65, side=2, line=1.4, at=60)
mtext("Scaled rank for thermal distance\n+ habitat distance", cex=0.65, side=1, line=2.2, at=0.5)
mtext("H", cex=1.2, side=3, line=-0.5, at=0.5)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-4, at=0.8)
hist(ranks_thermal_resource, main="", ylab="", xlab="", xlim=c(0,1), cex.lab=1, cex.axis=0.8, col="grey", border="white", breaks=10)
mtext("Number of species", cex=0.65, side=2, line=1.4, at=80)
mtext("Scaled rank for thermal distance\n+ resource scarcity", cex=0.65, side=1, line=2.2, at=0.5)
mtext("I", cex=1.2, side=3, line=-0.5, at=0.5)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-5, at=0.8)
hist(ranks_habitat_resource, main="", ylab="", xlab="", xlim=c(0,1), cex.lab=1, cex.axis=0.8, col="grey", border="white", breaks=10)
mtext("Number of species", cex=0.65, side=2, line=1.4, at=60)
mtext("Scaled rank for habitat distance\n+ resource scarcity", cex=0.65, side=1, line=2.2, at=0.5)
mtext("J", cex=1.2, side=3, line=-0.5, at=0.5)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-4, at=0.8)
hist(ranks_geo_thermal_habitat, main="", ylab="", xlab="", xlim=c(0,1), cex.lab=1, cex.axis=0.8, col="grey", border="white", breaks=10)
mtext("Number of species", cex=0.65, side=2, line=1.4, at=100)
mtext("Scaled rank for geographical distance\n+ thermal distance + habitat distance", cex=0.65, side=1, line=2.2, at=0.5)
mtext("K", cex=1.2, side=3, line=-0.5, at=0.5)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-5, at=0.8)
hist(ranks_geo_thermal_resource, main="", ylab="", xlab="", xlim=c(0,1), cex.lab=1, cex.axis=0.8, col="grey", border="white", breaks=10)
mtext("Number of species", cex=0.65, side=2, line=1.4, at=175)
mtext("Scaled rank for geographical distance\n+ thermal distance + resource scarcity", cex=0.65, side=1, line=2.2, at=0.5)
mtext("L", cex=1.2, side=3, line=-0.5, at=0.5)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-5, at=0.8)
hist(ranks_geo_habitat_resource, main="", ylab="", xlab="", xlim=c(0,1), cex.lab=1, cex.axis=0.8, col="grey", border="white", breaks=10)
mtext("Number of species", cex=0.65, side=2, line=1.4, at=100)
mtext("Scaled rank for geographical distance\n+ habitat distance + resource scarcity", cex=0.65, side=1, line=2.2, at=0.5)
mtext("M", cex=1.2, side=3, line=-0.5, at=0.5)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-5, at=0.8)
hist(ranks_thermal_habitat_resource, main="", ylab="", xlab="", xlim=c(0,1), cex.lab=1, cex.axis=0.8, col="grey", border="white", breaks=10)
mtext("Number of species", cex=0.65, side=2, line=1.4, at=70)
mtext("Scaled rank for thermal distance\n+ habitat distance + resource scarcity", cex=0.65, side=1, line=2.2, at=0.5)
mtext("N", cex=1.2, side=3, line=-0.5, at=0.5)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-4, at=0.8)
hist(ranks_geo_thermal_habitat_resource, main="", ylab="", xlab="", xlim=c(0,1), cex.lab=1, cex.axis=0.8, col="grey", border="white", breaks=10)
mtext("Number of species", cex=0.65, side=2, line=1.4, at=110)
mtext("Scaled rank for geographical distance + thermal\ndistance + habitat distance + resource scarcity", cex=0.65, side=1, line=2.2, at=0.5)
mtext("O", cex=1.2, side=3, line=-0.5, at=0.5)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-5, at=0.8)

#Kolmogorov-Smirnov tests of whether the distribution of ranks is left-skewed wompared to the uniform distribution
ks.test(ranks_thermal, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_habitat, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_geo, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_resource, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_geo_thermal, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_geo_habitat, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_geo_resource, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_thermal_habitat, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_thermal_resource, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_geo_thermal_habitat, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_geo_thermal_resource, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_thermal_habitat_resource, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_geo_thermal_habitat_resource, runif(100000,0,1), alternative="greater")$p.value

ks.test(ranks_geo_thermal, ranks_geo, alternative="greater")

ks.test(ranks_geo_thermal_habitat_resource, ranks_geo_thermal_resource, alternative="greater")$p.value
ks.test(ranks_geo_thermal_habitat_resource, ranks_geo_thermal_resource, alternative="greater")$p.value






####  Figure 3: Global patterns in the apparent performance of species' migrations in relation to the alternative options available to them  ###


ranks_geo_thermal_resourceNH_WH <- vector()
for(i in 1:length(geo.distances.simulatedNH_WH)){
	geo_thermal_resource <- scale(c(geo.distancesNH_WH[i], geo.distances.simulatedNH_WH[[i]])) + scale(c(thermal.distancesNH_WH[i], thermal.distances.simulatedNH_WH[[i]])) + scale(c(resource.scarcityNH_WH[i], resource.scarcity.simulatedNH_WH[[i]]))
	ranks_geo_thermal_resourceNH_WH[i] <- length(which(geo_thermal_resource[-1] < geo_thermal_resource[1])) / length(geo_thermal_resource[-1])	
}
ranks_geo_thermal_resourceSH_WH <- vector()
for(i in 1:length(geo.distances.simulatedSH_WH)){
	geo_thermal_resource <- scale(c(geo.distancesSH_WH[i], geo.distances.simulatedSH_WH[[i]])) + scale(c(thermal.distancesSH_WH[i], thermal.distances.simulatedSH_WH[[i]])) + scale(c(resource.scarcitySH_WH[i], resource.scarcity.simulatedSH_WH[[i]]))
	ranks_geo_thermal_resourceSH_WH[i] <- length(which(geo_thermal_resource[-1] < geo_thermal_resource[1])) / length(geo_thermal_resource[-1])	
}
ranks_geo_thermal_resourceWH <- c(ranks_geo_thermal_resourceNH_WH, ranks_geo_thermal_resourceSH_WH)
ranks_geo_thermal_resourceNH_EH <- vector()
for(i in 1:length(geo.distances.simulatedNH_EH)){
	geo_thermal_resource <- scale(c(geo.distancesNH_EH[i], geo.distances.simulatedNH_EH[[i]])) + scale(c(thermal.distancesNH_EH[i], thermal.distances.simulatedNH_EH[[i]])) + scale(c(resource.scarcityNH_EH[i], resource.scarcity.simulatedNH_EH[[i]]))
	ranks_geo_thermal_resourceNH_EH[i] <- length(which(geo_thermal_resource[-1] < geo_thermal_resource[1])) / length(geo_thermal_resource[-1])	
}
ranks_geo_thermal_resourceSH_EH <- vector()
for(i in 1:length(geo.distances.simulatedSH_EH)){
	geo_thermal_resource <- scale(c(geo.distancesSH_EH[i], geo.distances.simulatedSH_EH[[i]])) + scale(c(thermal.distancesSH_EH[i], thermal.distances.simulatedSH_EH[[i]])) + scale(c(resource.scarcitySH_EH[i], resource.scarcity.simulatedSH_EH[[i]]))
	ranks_geo_thermal_resourceSH_EH[i] <- length(which(geo_thermal_resource[-1] < geo_thermal_resource[1])) / length(geo_thermal_resource[-1])	
}
ranks_geo_thermal_resourceEH <- c(ranks_geo_thermal_resourceNH_EH, ranks_geo_thermal_resourceSH_EH)
pres.list_BR_NH_WH <- apply(PresAbs_BR_NH_WH2, 2, function(x) hexidWH[which(x==1)])
pres.list_NB_NH_WH <- apply(PresAbs_NB_NH_WH2, 2, function(x) hexidWH[which(x==1)])
pres.list_BR_SH_WH <- apply(PresAbs_BR_SH_WH2, 2, function(x) hexidWH[which(x==1)])
pres.list_NB_SH_WH <- apply(PresAbs_NB_SH_WH2, 2, function(x) hexidWH[which(x==1)])
pres.list_BR_NH_EH <- apply(PresAbs_BR_NH_EH2, 2, function(x) hexidEH[which(x==1)])
pres.list_NB_NH_EH <- apply(PresAbs_NB_NH_EH2, 2, function(x) hexidEH[which(x==1)])
pres.list_BR_SH_EH <- apply(PresAbs_BR_SH_EH2, 2, function(x) hexidEH[which(x==1)])
pres.list_NB_SH_EH <- apply(PresAbs_NB_SH_EH2, 2, function(x) hexidEH[which(x==1)])
hex.weight_BR_WH <- vector()
for(i in 1:length(hexidWH)){
	hex.weight_BR_WH[i] <- sum(ranks_geo_thermal_resourceWH[c(unlist(lapply(pres.list_BR_NH_WH, function(x) is.element(hexidWH[i],x))), unlist(lapply(pres.list_BR_SH_WH, function(x) is.element(hexidWH[i],x))))]) / length(ranks_geo_thermal_resourceWH[c(unlist(lapply(pres.list_BR_NH_WH, function(x) is.element(hexidWH[i],x))), unlist(lapply(pres.list_BR_SH_WH, function(x) is.element(hexidWH[i],x))))])
}
hex.weight_NB_WH <- vector()
for(i in 1:length(hexidWH)){
	hex.weight_NB_WH[i] <- sum(ranks_geo_thermal_resourceWH[c(unlist(lapply(pres.list_NB_NH_WH, function(x) is.element(hexidWH[i],x))), unlist(lapply(pres.list_NB_SH_WH, function(x) is.element(hexidWH[i],x))))]) / length(ranks_geo_thermal_resourceWH[c(unlist(lapply(pres.list_NB_NH_WH, function(x) is.element(hexidWH[i],x))), unlist(lapply(pres.list_NB_SH_WH, function(x) is.element(hexidWH[i],x))))])
}
hex.weight_BR_EH <- vector()
for(i in 1:length(hexidEH)){
	hex.weight_BR_EH[i] <- sum(ranks_geo_thermal_resourceEH[c(unlist(lapply(pres.list_BR_NH_EH, function(x) is.element(hexidEH[i],x))), unlist(lapply(pres.list_BR_SH_EH, function(x) is.element(hexidEH[i],x))))]) / length(ranks_geo_thermal_resourceEH[c(unlist(lapply(pres.list_BR_NH_EH, function(x) is.element(hexidEH[i],x))), unlist(lapply(pres.list_BR_SH_EH, function(x) is.element(hexidEH[i],x))))])
}
hex.weight_NB_EH <- vector()
for(i in 1:length(hexidEH)){
	hex.weight_NB_EH[i] <- sum(ranks_geo_thermal_resourceEH[c(unlist(lapply(pres.list_NB_NH_EH, function(x) is.element(hexidEH[i],x))), unlist(lapply(pres.list_NB_SH_EH, function(x) is.element(hexidEH[i],x))))]) / length(ranks_geo_thermal_resourceEH[c(unlist(lapply(pres.list_NB_NH_EH, function(x) is.element(hexidEH[i],x))), unlist(lapply(pres.list_NB_SH_EH, function(x) is.element(hexidEH[i],x))))])
}
hex.weightBR <- c(hex.weight_BR_WH, hex.weight_BR_EH)
hex.weightNB <- c(hex.weight_NB_WH, hex.weight_NB_EH)

#Plot weighted richness
jpeg("Figure3.jpg", width=1200, height=800, quality=100)
par(mfrow=c(2,2), mar=c(0.1,0.1,1.4,0.1), mgp=c(1.5,0.5,0))
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(5)[as.numeric(cut(hex.weightBR, breaks=c(-0.1,0.1,0.2,0.3,0.4,1.1)))]
datcol[which(is.na(datcol)==T)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey")
mtext("A", cex=1.9, side=3, line=-2, at=-180)
mtext("Breeding season", cex=1.5, side=3, line=0.15, at=0)
legend("bottomleft", inset=.04, bg="grey", box.col="grey", title="Average rank\nof migrants", c("> 0.4","0.3–0.4", "0.2–0.3", "0.1–0.2", "0–0.1", "No species"), fill=c(rev(rbPal(5)),"white"), cex=1.5)
datcol <- rbPal(5)[as.numeric(cut(hex.weightNB, breaks=c(-0.1,0.1,0.2,0.3,0.4,1.1)))]
datcol[which(is.na(datcol)==T)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey")
mtext("B", cex=1.9, side=3, line=-2, at=-180)
mtext("Non-breeding season", cex=1.5, side=3, line=0.15, at=0)

#Plot observed richness in migrants
PresAbs_BR_WH <- cbind(PresAbs_BR_NH_WH2, PresAbs_BR_SH_WH2)
PresAbs_NB_WH <- cbind(PresAbs_NB_NH_WH2, PresAbs_NB_SH_WH2)
PresAbs_BR_EH <- cbind(PresAbs_BR_NH_EH2, PresAbs_BR_SH_EH2)
PresAbs_NB_EH <- cbind(PresAbs_NB_NH_EH2, PresAbs_NB_SH_EH2)
mbr <- c(apply(PresAbs_BR_WH, 1, sum), apply(PresAbs_BR_EH, 1, sum))
mnb <- c(apply(PresAbs_NB_WH, 1, sum), apply(PresAbs_NB_EH, 1, sum))
rbPal <- colorRampPalette(c("yellow3", "dark green"))
datcol <- rbPal(5)[as.numeric(cut(mbr, breaks=c(0.9,25,50,75,100,150)))]
datcol[which(is.na(datcol)==T)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey")
mtext("C", cex=1.9, side=3, line=-2, at=-180)
legend("bottomleft", inset=.04, bg="grey", box.col="grey", title="Richness\nin migrants", c("> 100","75–100", "50–75", "25–50", "1–25", "0"), fill=c(rev(rbPal(5)),"white"), cex=1.5)
datcol <- rbPal(5)[as.numeric(cut(mnb, breaks=c(0.9,25,50,75,100,150)))]
datcol[which(is.na(datcol)==T)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey")
mtext("D", cex=1.9, side=3, line=-2, at=-180)
dev.off()




####  Figure 4: Connections between the breeding and non-breeding ranges of species according to the apparent performance of their migration in relation to the alternatives available to them   ###


centroids.breeding.grounds <- rbind(centroids.breeding.groundsNH_WH, centroids.breeding.groundsNH_EH, centroids.breeding.groundsSH_WH, centroids.breeding.groundsSH_EH)
centroids.nonbreeding.grounds <- rbind(centroids.nonbreeding.groundsNH_WH, centroids.nonbreeding.groundsNH_EH, centroids.nonbreeding.groundsSH_WH, centroids.nonbreeding.groundsSH_EH)
#Split the species into 3 groups: good rank, intermerdiate rank and bad rank
centroids.breeding.grounds_1 <- centroids.breeding.grounds[which(ranks_geo_thermal_resource < 0.1),]
centroids.nonbreeding.grounds_1 <- centroids.nonbreeding.grounds[which(ranks_geo_thermal_resource < 0.1),]
centroids.breeding.grounds_2 <- centroids.breeding.grounds[which(ranks_geo_thermal_resource >= 0.1 & ranks_geo_thermal_resource <= 0.5),]
centroids.nonbreeding.grounds_2 <- centroids.nonbreeding.grounds[which(ranks_geo_thermal_resource >= 0.1 & ranks_geo_thermal_resource <= 0.5),]
centroids.breeding.grounds_3 <- centroids.breeding.grounds[which(ranks_geo_thermal_resource > 0.5),]
centroids.nonbreeding.grounds_3 <- centroids.nonbreeding.grounds[which(ranks_geo_thermal_resource > 0.5),]
#Plot
jpeg("Figure4.jpg", width=600, height=1000, quality=100)
par(mfrow=c(3,1), mar=c(0.1,0.1,0.1,0.1), mgp=c(1.5,0.5,0))
plot(hexgrid2, col= "dark grey", border = "dark grey", bg="light grey")
points(centroids.breeding.grounds_1, pch=20, col="red", cex=0.3)
points(centroids.nonbreeding.grounds_1, pch=20, col="blue", cex=0.3)
for(i in 1:length(centroids.breeding.grounds_1[,1])){
	inter <- gcIntermediate(centroids.breeding.grounds_1[i,], centroids.nonbreeding.grounds_1[i,], n=50, addStartEnd=T)
	lines(inter, lwd=1, col="yellow")
}
mtext("A", cex=1.9, side=3, line=-2.2, at=-180)
mtext("Scaled rank < 0.1", cex=1.5, side=1, line=-2.2, at=25)
plot(hexgrid2, col= "dark grey", border = "dark grey", bg="light grey")
points(centroids.breeding.grounds_2, pch=20, col="red", cex=0.3)
points(centroids.nonbreeding.grounds_2, pch=20, col="blue", cex=0.3)
for(i in 1:length(centroids.breeding.grounds_2[,1])){
	inter <- gcIntermediate(centroids.breeding.grounds_2[i,], centroids.nonbreeding.grounds_2[i,], n=50, addStartEnd=T)
	lines(inter, lwd=1, col="orange")
}
mtext("B", cex=1.9, side=3, line=-2.2, at=-180)
mtext("0.1 <= Scaled rank <= 0.5", cex=1.5, side=1, line=-2.2, at=25)
plot(hexgrid2, col= "dark grey", border = "dark grey", bg="light grey")
points(centroids.breeding.grounds_3, pch=20, col="red", cex=0.3)
points(centroids.nonbreeding.grounds_3, pch=20, col="blue", cex=0.3)
for(i in 1:length(centroids.breeding.grounds_3[,1])){
	inter <- gcIntermediate(centroids.breeding.grounds_3[i,], centroids.nonbreeding.grounds_3[i,], n=50, addStartEnd=T)
	lines(inter, lwd=1, col="brown4")
}
mtext("C", cex=1.9, side=3, line=-2.2, at=-180)
mtext("Scaled rank > 0.5", cex=1.5, side=1, line=-2.2, at=25)
dev.off()






### Figure A2 - illustrating how the ranks were computed using Schrenck's Bittern (Ixobrychus.eurhythmus) ###

jpeg("FigureA2.jpg", width=800, height=600, quality=100)
par(mfrow=c(2,2), mar=c(2.5,2.5,0.2,1.5), mgp=c(1.5,0.5,0))
plot(hexgridEH, col="grey", border = "grey")
plot(hexgridEH[match(PresAbs_BR_NH[which(PresAbs_BR_NH[,which(colnames(PresAbs_NB_NH) == "Ixobrychus.eurhythmus")] == 1),1], hexgridEH@data[,1]),], col="red3", border = "red3", add=T)
plot(hexgridEH[match(hexidEH[simulated.ranges_BR_NH_EH[[182]][[73]]], hexgridEH@data[,1]),], col="black", border = "black", add=T)
plot(hexgridEH[match(hexidEH[simulated.ranges_BR_NH_EH[[182]][[57]]], hexgridEH@data[,1]),], col="black", border = "black", add=T)
mtext("A", cex=1.8, side=3, line=-1.5, at=-50)
plot(scale(c(geo.distances.obs[406], geo.distances.simulated[[406]])), scale(c(thermal.distances.obs[406], thermal.distances.simulated[[406]])), pch=20, axes=F, xlab="Geographical distance", ylab = "Thermal distance", cex.lab=1.2, add=F)
axis(1,at=c(-2,-1,0,1,2), labels=c(-2,-1,0,1,2))
axis(2,at=c(-2,-1,0,1,2,3), labels=c(-2,-1,0,1,2,3))
abline(a = scale(c(thermal.distances.obs[406], thermal.distances.simulated[[406]]))[1], b = 0, col="orange")
abline(a = (scale(c(geo.distances.obs[406], geo.distances.simulated[[406]])) + scale(c(thermal.distances.obs[406], thermal.distances.simulated[[406]])))[1], b = -1, col="blue")
points(scale(c(geo.distances.obs[406], geo.distances.simulated[[406]]))[1], scale(c(thermal.distances.obs[406], thermal.distances.simulated[[406]]))[1], pch=20, col="red3", cex=2.5)
mtext("C", cex=1.8, side=3, line=-1.5, at=-2.35)
plot(hexgridEH, col="grey", border = "grey")
plot(hexgridEH[match(PresAbs_NB_NH[which(PresAbs_NB_NH[,which(colnames(PresAbs_NB_NH) == "Ixobrychus.eurhythmus")] == 1),1], hexgridEH@data[,1]),], col="red3", border = "red3", add=T)
plot(hexgridEH[match(hexidEH[simulated.ranges_NB_NH_EH[[182]][[77]]], hexgridEH@data[,1]),], col="black", border = "black", add=T)
plot(hexgridEH[match(hexidEH[simulated.ranges_NB_NH_EH[[182]][[66]]], hexgridEH@data[,1]),], col="black", border = "black", add=T)
mtext("B", cex=1.8, side=3, line=-1.5, at=-50)
geo_thermal_resource_schrenck <- scale(c(geo.distances.obs[406], geo.distances.simulated[[406]])) + scale(c(thermal.distances.obs[406], thermal.distances.simulated[[406]])) + scale(c(resource.scarcity.obs[406], resource.scarcity.simulated[[406]]))
plot(geo_thermal_resource_schrenck, rnorm(length(geo_thermal_resource_schrenck), 1, 0.1), pch=20, axes=F, xlab="Thermal distance + geographical distance + resource scarcity", ylab = "", ylim=c(0,2), cex.lab=1.2)
axis(1,at=c(-3,-2,-1,0,1,2,3), labels=c(-3,-2,-1,0,1,2,3))
abline(v=geo_thermal_resource_schrenck[1], col="green")
points(geo_thermal_resource_schrenck[1], 1, pch=20, ylim=c(0,2), col="red3", cex=2.5)
mtext("D", cex=1.8, side=3, line=-1.5, at=-3.7)
dev.off()





###  Export excel table for appendix of publication

species.names <- c(names(PresAbs_BR_NH_WH2), names(PresAbs_BR_NH_EH2), names(PresAbs_BR_SH_WH2), names(PresAbs_BR_SH_EH2))
longitudinal.hemisphere <- c(rep("WH", length(PresAbs_BR_NH_WH2[1,])), rep("EH", length(PresAbs_BR_NH_EH2[1,])), rep("WH", length(PresAbs_BR_SH_WH2[1,])), rep("EH", length(PresAbs_BR_SH_EH2[1,])))
xlsTable <- cbind(species.names, longitudinal.hemisphere, geo.distances.obs, thermal.distances.obs, habitat.distances.obs, resource.scarcity.obs, ranks_geo, ranks_thermal, ranks_habitat, ranks_resource, ranks_geo_thermal, ranks_geo_habitat, ranks_geo_resource, ranks_thermal_habitat, ranks_thermal_resource, ranks_habitat_resource, ranks_geo_thermal_habitat, ranks_geo_thermal_resource, ranks_geo_habitat_resource, ranks_thermal_habitat_resource, ranks_geo_thermal_habitat_resource)
colnames(xlsTable) <- c("species name", "Longitudinal hemisphere", "Geographical distance", "Thermal distance", "Habitat distance", "Resource scarcity", "Scaled rank geo distance", "Scaled rank thermal distance", "Scaled rank habitat distance", "Scaled rank resource scarcity", "Scaled rank geo distance + thermal distance", "Scaled rank geo distance + habitat distance", "Scaled rank geo distance + resource scarcity", "Scaled rank thermal distance + habitat distance", "Scaled rank thermal distance + resource scarcity", "Scaled rank habitat distance + resource scarcity", "Scaled rank geo distance + thermal distance + habitat distance", "Scaled rank geo distance + thermal distance + resource scarcity", "Scaled rank geo distance + habitat distance + resource scarcity", "Scaled rank thermal distance + habitat distance + resource scarcity", "Scaled rank geo distance + thermal distance + habitat distance + resource scarcity")
xlsTable <- xlsTable[order(species.names),]
write.csv(xlsTable, "Appendix1.csv")






####  Analysis of simulated migration destinations that are performing better than the observed one  ####


## Variables separated between western and eastern hemispheres

geo.distances.obs_WH <- c(geo.distancesNH_WH, geo.distancesSH_WH) 
geo.distances.obs_EH <- c(geo.distancesNH_EH, geo.distancesSH_EH)
thermal.distances.obs_WH <- c(thermal.distancesNH_WH, thermal.distancesSH_WH) 
thermal.distances.obs_EH <- c(thermal.distancesNH_EH, thermal.distancesSH_EH) 
resource.scarcity.obs_WH <- c(resource.scarcityNH_WH, resource.scarcitySH_WH) 
resource.scarcity.obs_EH <- c(resource.scarcityNH_EH, resource.scarcitySH_EH) 
geo.distances.simulated_WH <- c(geo.distances.simulatedNH_WH, geo.distances.simulatedSH_WH)
geo.distances.simulated_EH <- c(geo.distances.simulatedNH_EH, geo.distances.simulatedSH_EH)
thermal.distances.simulated_WH <- c(thermal.distances.simulatedNH_WH, thermal.distances.simulatedSH_WH)
thermal.distances.simulated_EH <- c(thermal.distances.simulatedNH_EH, thermal.distances.simulatedSH_EH)
resource.scarcity.simulated_WH <- c(resource.scarcity.simulatedNH_WH, resource.scarcity.simulatedSH_WH)
resource.scarcity.simulated_EH <- c(resource.scarcity.simulatedNH_EH, resource.scarcity.simulatedSH_EH)


## Select simulated migration destinations that are performing better than the observed one

geo_best_simu_brsim_WH <- list()
geo_best_simu_nbsim_WH <- list()
for(i in 1:length(geo.distances.simulated_WH)){
	geo_thermal_resource <- scale(c(geo.distances.obs_WH[i], geo.distances.simulated_WH[[i]])) + scale(c(thermal.distances.obs_WH[i], thermal.distances.simulated_WH[[i]])) + scale(c(resource.scarcity.obs_WH[i], resource.scarcity.simulated_WH[[i]]))
	geo_best_simu_brsim_WH[[i]] <- which(geo_thermal_resource[2:101] < geo_thermal_resource[1])
	geo_best_simu_nbsim_WH[[i]] <- which(geo_thermal_resource[102:201] < geo_thermal_resource[1])
}
geo_best_simu_brsim_EH <- list()
geo_best_simu_nbsim_EH <- list()
for(i in 1:length(geo.distances.simulated_EH)){
	geo_thermal_resource <- scale(c(geo.distances.obs_EH[i], geo.distances.simulated_EH[[i]])) + scale(c(thermal.distances.obs_EH[i], thermal.distances.simulated_EH[[i]])) + scale(c(resource.scarcity.obs_EH[i], resource.scarcity.simulated_EH[[i]]))
	geo_best_simu_brsim_EH[[i]] <- which(geo_thermal_resource[2:101] < geo_thermal_resource[1])
	geo_best_simu_nbsim_EH[[i]] <- which(geo_thermal_resource[102:201] < geo_thermal_resource[1])
}

geo_best_simu_brsim_dist_WH <- list()
geo_best_simu_nbsim_dist_WH <- list()
for(i in 1:length(geo.distances.simulated_WH)){
	geo_best_simu_brsim_dist_WH[[i]] <- geo.distances.simulated_WH[[i]][1:100][geo_best_simu_brsim_WH[[i]]]
	geo_best_simu_nbsim_dist_WH[[i]] <- geo.distances.simulated_WH[[i]][101:200][geo_best_simu_nbsim_WH[[i]]]
}
geo_best_simu_brsim_dist_EH <- list()
geo_best_simu_nbsim_dist_EH <- list()
for(i in 1:length(geo.distances.simulated_EH)){
	geo_best_simu_brsim_dist_EH[[i]] <- geo.distances.simulated_EH[[i]][1:100][geo_best_simu_brsim_EH[[i]]]
	geo_best_simu_nbsim_dist_EH[[i]] <- geo.distances.simulated_EH[[i]][101:200][geo_best_simu_nbsim_EH[[i]]]
}

geo_best_simu_brsim_dist_diff_WH <- list()
geo_best_simu_nbsim_dist_diff_WH <- list()
for(i in 1:length(geo_best_simu_brsim_dist_WH)){
	geo_best_simu_brsim_dist_diff_WH[[i]] <- geo_best_simu_brsim_dist_WH[[i]] - geo.distances.obs_WH[i]
	geo_best_simu_nbsim_dist_diff_WH[[i]] <- geo_best_simu_nbsim_dist_WH[[i]] - geo.distances.obs_WH[i]
}
geo_best_simu_brsim_dist_diff_EH <- list()
geo_best_simu_nbsim_dist_diff_EH <- list()
for(i in 1:length(geo_best_simu_brsim_dist_EH)){
	geo_best_simu_brsim_dist_diff_EH[[i]] <- geo_best_simu_brsim_dist_EH[[i]] - geo.distances.obs_EH[i]
	geo_best_simu_nbsim_dist_diff_EH[[i]] <- geo_best_simu_nbsim_dist_EH[[i]] - geo.distances.obs_EH[i]
}

geo_best_simu_brsim_dist_diff_mean_WH <- unlist(lapply(geo_best_simu_brsim_dist_diff_WH, mean))
geo_best_simu_nbsim_dist_diff_mean_WH <- unlist(lapply(geo_best_simu_nbsim_dist_diff_WH, mean))
geo_best_simu_brsim_dist_diff_mean_EH <- unlist(lapply(geo_best_simu_brsim_dist_diff_EH, mean))
geo_best_simu_nbsim_dist_diff_mean_EH <- unlist(lapply(geo_best_simu_nbsim_dist_diff_EH, mean))

best_simu_brsim_number_WH <- unlist(lapply(geo_best_simu_brsim_dist_diff_WH, length))
best_simu_nbsim_number_WH <- unlist(lapply(geo_best_simu_nbsim_dist_diff_WH, length))
best_simu_brsim_number_EH <- unlist(lapply(geo_best_simu_brsim_dist_diff_EH, length))
best_simu_nbsim_number_EH <- unlist(lapply(geo_best_simu_nbsim_dist_diff_EH, length))




## Migration longitude (computed as the average between the longitude of the centroid of the breeding range and  the longitude of the centroid of the non-breeding range)

# Observed migrations
lon_obs_NH_WH <- vector()
for(i in 1:length(centroids.nonbreeding.groundsNH_WH[,1])){
	lon_obs_NH_WH[i] <- centroids.breeding.groundsNH_WH[i,1] - centroids.nonbreeding.groundsNH_WH[i,1]
}
lon_obs_SH_WH <- vector()
for(i in 1:length(centroids.nonbreeding.groundsSH_WH[,1])){
	lon_obs_SH_WH[i] <- centroids.breeding.groundsSH_WH[i,1] - centroids.nonbreeding.groundsSH_WH[i,1]
}
lon_obs_NH_EH <- vector()
for(i in 1:length(centroids.nonbreeding.groundsNH_EH[,1])){
	lon_obs_NH_EH[i] <- centroids.breeding.groundsNH_EH[i,1] - centroids.nonbreeding.groundsNH_EH[i,1]
}
lon_obs_SH_EH <- vector()
for(i in 1:length(centroids.nonbreeding.groundsSH_EH[,1])){
	lon_obs_SH_EH[i] <- centroids.breeding.groundsSH_EH[i,1] - centroids.nonbreeding.groundsSH_EH[i,1]
}
lon_obs_WH <- c(lon_obs_NH_WH, lon_obs_SH_WH)
lon_obs_EH <- c(lon_obs_NH_EH, lon_obs_SH_EH)

# Simulated migrations
lon_sim_NH_WH_nbsim <- list()
lon_sim_NH_WH_brsim <- list()
for(i in 1:length(centroids.breeding.groundsNH_WH[,1])){	
	nbsim <- vector()
	brsim <- vector()
	for(j in 1:100){
		nbsim[j] <- centroids.simulatedNB_NH_WH[[i]][[j]][1] - centroids.nonbreeding.groundsNH_WH[i,1]
	}
	lon_sim_NH_WH_nbsim[[i]] <- nbsim
	for(j in 1:100){
		brsim[j] <- centroids.simulatedBR_NH_WH[[i]][[j]][1] - centroids.breeding.groundsNH_WH[i,1]
	}
	lon_sim_NH_WH_brsim[[i]] <- brsim
}
lon_sim_SH_WH_nbsim <- list()
lon_sim_SH_WH_brsim <- list()
for(i in 1:length(centroids.breeding.groundsSH_WH[,1])){	
	nbsim <- vector()
	brsim <- vector()
	for(j in 1:100){
		nbsim[j] <- centroids.simulatedNB_SH_WH[[i]][[j]][1] - centroids.nonbreeding.groundsSH_WH[i,1]
	}
	lon_sim_SH_WH_nbsim[[i]] <- nbsim
	for(j in 1:100){
		brsim[j] <- centroids.simulatedBR_SH_WH[[i]][[j]][1] - centroids.breeding.groundsSH_WH[i,1]
	}
	lon_sim_SH_WH_brsim[[i]] <- brsim
}
lon_sim_NH_EH_nbsim <- list()
lon_sim_NH_EH_brsim <- list()
for(i in 1:length(centroids.breeding.groundsNH_EH[,1])){	
	nbsim <- vector()
	brsim <- vector()
	for(j in 1:100){
		nbsim[j] <- centroids.simulatedNB_NH_EH[[i]][[j]][1] - centroids.nonbreeding.groundsNH_EH[i,1]
	}
	lon_sim_NH_EH_nbsim[[i]] <- nbsim
	for(j in 1:100){
		brsim[j] <- centroids.simulatedBR_NH_EH[[i]][[j]][1] - centroids.breeding.groundsNH_EH[i,1]
	}
	lon_sim_NH_EH_brsim[[i]] <- brsim
}
lon_sim_SH_EH_nbsim <- list()
lon_sim_SH_EH_brsim <- list()
for(i in 1:length(centroids.breeding.groundsSH_EH[,1])){	
	nbsim <- vector()
	brsim <- vector()
	for(j in 1:100){
		nbsim[j] <- centroids.simulatedNB_SH_EH[[i]][[j]][1] - centroids.nonbreeding.groundsSH_EH[i,1]
	}
	lon_sim_SH_EH_nbsim[[i]] <- nbsim
	for(j in 1:100){
		brsim[j] <- centroids.simulatedBR_SH_EH[[i]][[j]][1] - centroids.breeding.groundsSH_EH[i,1]
	}
	lon_sim_SH_EH_brsim[[i]] <- brsim
}
lon_sim_nbsim_WH <- c(lon_sim_NH_WH_nbsim, lon_sim_SH_WH_nbsim)
lon_sim_nbsim_EH <- c(lon_sim_NH_EH_nbsim, lon_sim_SH_EH_nbsim)
lon_sim_brsim_WH <- c(lon_sim_NH_WH_brsim, lon_sim_SH_WH_brsim)
lon_sim_brsim_EH <- c(lon_sim_NH_EH_brsim, lon_sim_SH_EH_brsim)

lon_nbsim_best_simu_sel_WH <- list()
lon_brsim_best_simu_sel_WH <- list()
for(i in 1:length(geo_best_simu_nbsim_WH)){
	lon_nbsim_best_simu_sel_WH[[i]] <- lon_sim_nbsim_WH[[i]][geo_best_simu_nbsim_WH[[i]]]
	lon_brsim_best_simu_sel_WH[[i]] <- lon_sim_brsim_WH[[i]][geo_best_simu_brsim_WH[[i]]]
}
lon_nbsim_best_simu_sel_EH <- list()
lon_brsim_best_simu_sel_EH <- list()
for(i in 1:length(geo_best_simu_nbsim_EH)){
	lon_nbsim_best_simu_sel_EH[[i]] <- lon_sim_nbsim_EH[[i]][geo_best_simu_nbsim_EH[[i]]]
	lon_brsim_best_simu_sel_EH[[i]] <- lon_sim_brsim_EH[[i]][geo_best_simu_brsim_EH[[i]]]
}
lon_nbsim_best_simu_mean_WH <- unlist(lapply(lon_nbsim_best_simu_sel_WH, mean))
lon_brsim_best_simu_mean_WH <- unlist(lapply(lon_brsim_best_simu_sel_WH, mean))
lon_nbsim_best_simu_mean_EH <- unlist(lapply(lon_nbsim_best_simu_sel_EH, mean))
lon_brsim_best_simu_mean_EH <- unlist(lapply(lon_brsim_best_simu_sel_EH, mean))



## Migration latitude (computed as the average between the longitude of the centroid of the breeding range and  the longitude of the centroid of the non-breeding range)

# Observed migrations
lat_obs_NH_WH <- vector()
for(i in 1:length(centroids.nonbreeding.groundsNH_WH[,1])){
	lat_obs_NH_WH[i] <- abs(centroids.breeding.groundsNH_WH[i,2] - centroids.nonbreeding.groundsNH_WH[i,2])
}
lat_obs_SH_WH <- vector()
for(i in 1:length(centroids.nonbreeding.groundsSH_WH[,1])){
	lat_obs_SH_WH[i] <- abs(centroids.breeding.groundsSH_WH[i,2] - centroids.nonbreeding.groundsSH_WH[i,2])
}
lat_obs_NH_EH <- vector()
for(i in 1:length(centroids.nonbreeding.groundsNH_EH[,1])){
	lat_obs_NH_EH[i] <- abs(centroids.breeding.groundsNH_EH[i,2] - centroids.nonbreeding.groundsNH_EH[i,2])
}
lat_obs_SH_EH <- vector()
for(i in 1:length(centroids.nonbreeding.groundsSH_EH[,1])){
	lat_obs_SH_EH[i] <- abs(centroids.breeding.groundsSH_EH[i,2] - centroids.nonbreeding.groundsSH_EH[i,2])
}
lat_obs_WH <- c(lat_obs_NH_WH, lat_obs_SH_WH)
lat_obs_EH <- c(lat_obs_NH_EH, lat_obs_SH_EH)

# Simulated migrations
lat_sim_NH_WH_nbsim <- list()
lat_sim_NH_WH_brsim <- list()
for(i in 1:length(centroids.breeding.groundsNH_WH[,1])){	
	nbsim <- vector()
	brsim <- vector()
	for(j in 1:100){
		nbsim[j] <- centroids.simulatedNB_NH_WH[[i]][[j]][2] - centroids.nonbreeding.groundsNH_WH[i,2]
	}
	lat_sim_NH_WH_nbsim[[i]] <- nbsim
	for(j in 1:100){
		brsim[j] <- centroids.simulatedBR_NH_WH[[i]][[j]][2] - centroids.breeding.groundsNH_WH[i,2]
	}
	lat_sim_NH_WH_brsim[[i]] <- brsim
}
lat_sim_SH_WH_nbsim <- list()
lat_sim_SH_WH_brsim <- list()
for(i in 1:length(centroids.breeding.groundsSH_WH[,1])){	
	nbsim <- vector()
	brsim <- vector()
	for(j in 1:100){
		nbsim[j] <- centroids.simulatedNB_SH_WH[[i]][[j]][2] - centroids.nonbreeding.groundsSH_WH[i,2]
	}
	lat_sim_SH_WH_nbsim[[i]] <- nbsim
	for(j in 1:100){
		brsim[j] <- centroids.simulatedBR_SH_WH[[i]][[j]][2] - centroids.breeding.groundsSH_WH[i,2]
	}
	lat_sim_SH_WH_brsim[[i]] <- brsim
}
lat_sim_NH_EH_nbsim <- list()
lat_sim_NH_EH_brsim <- list()
for(i in 1:length(centroids.breeding.groundsNH_EH[,1])){	
	nbsim <- vector()
	brsim <- vector()
	for(j in 1:100){
		nbsim[j] <- centroids.simulatedNB_NH_EH[[i]][[j]][2] - centroids.nonbreeding.groundsNH_EH[i,2]
	}
	lat_sim_NH_EH_nbsim[[i]] <- nbsim
	for(j in 1:100){
		brsim[j] <- centroids.simulatedBR_NH_EH[[i]][[j]][2] - centroids.breeding.groundsNH_EH[i,2]
	}
	lat_sim_NH_EH_brsim[[i]] <- brsim
}
lat_sim_SH_EH_nbsim <- list()
lat_sim_SH_EH_brsim <- list()
for(i in 1:length(centroids.breeding.groundsSH_EH[,1])){	
	nbsim <- vector()
	brsim <- vector()
	for(j in 1:100){
		nbsim[j] <- centroids.simulatedNB_SH_EH[[i]][[j]][2] - centroids.nonbreeding.groundsSH_EH[i,2]
	}
	lat_sim_SH_EH_nbsim[[i]] <- nbsim
	for(j in 1:100){
		brsim[j] <- centroids.simulatedBR_SH_EH[[i]][[j]][2] - centroids.breeding.groundsSH_EH[i,2]
	}
	lat_sim_SH_EH_brsim[[i]] <- brsim
}
lat_sim_nbsim_WH <- c(lat_sim_NH_WH_nbsim, lat_sim_SH_WH_nbsim)
lat_sim_nbsim_EH <- c(lat_sim_NH_EH_nbsim, lat_sim_SH_EH_nbsim)
lat_sim_brsim_WH <- c(lat_sim_NH_WH_brsim, lat_sim_SH_WH_brsim)
lat_sim_brsim_EH <- c(lat_sim_NH_EH_brsim, lat_sim_SH_EH_brsim)

lat_nbsim_best_simu_sel_WH <- list()
lat_brsim_best_simu_sel_WH <- list()
for(i in 1:length(geo_best_simu_nbsim_WH)){
	lat_nbsim_best_simu_sel_WH[[i]] <- lat_sim_nbsim_WH[[i]][geo_best_simu_nbsim_WH[[i]]]
	lat_brsim_best_simu_sel_WH[[i]] <- lat_sim_brsim_WH[[i]][geo_best_simu_brsim_WH[[i]]]
}
lat_nbsim_best_simu_sel_EH <- list()
lat_brsim_best_simu_sel_EH <- list()
for(i in 1:length(geo_best_simu_nbsim_EH)){
	lat_nbsim_best_simu_sel_EH[[i]] <- lat_sim_nbsim_EH[[i]][geo_best_simu_nbsim_EH[[i]]]
	lat_brsim_best_simu_sel_EH[[i]] <- lat_sim_brsim_EH[[i]][geo_best_simu_brsim_EH[[i]]]
}
lat_nbsim_best_simu_mean_WH <- unlist(lapply(lat_nbsim_best_simu_sel_WH, mean))
lat_brsim_best_simu_mean_WH <- unlist(lapply(lat_brsim_best_simu_sel_WH, mean))
lat_nbsim_best_simu_mean_EH <- unlist(lapply(lat_nbsim_best_simu_sel_EH, mean))
lat_brsim_best_simu_mean_EH <- unlist(lapply(lat_brsim_best_simu_sel_EH, mean))




# Are better performing simulated migratory destinations mostly driven by simulated breeding range or simulated non-breeding range?

plot(best_simu_brsim_number_WH, best_simu_nbsim_number_WH)
points(1:100,1:100,type="l")
plot(best_simu_brsim_number_EH, best_simu_nbsim_number_EH)
points(1:100,1:100,type="l")


# Comparing migration distance and bearing of better simulated migration distance with the observed ones

plot(geo_best_simu_brsim_dist_diff_mean_WH, bearing_brsim_best_simu_diff_mean_WH, pch=20)
plot(geo_best_simu_nbsim_dist_diff_mean_WH, bearing_nbsim_best_simu_diff_mean_WH, pch=20)
plot(geo_best_simu_brsim_dist_diff_mean_EH, bearing_brsim_best_simu_diff_mean_EH, pch=20)
plot(geo_best_simu_nbsim_dist_diff_mean_EH, bearing_nbsim_best_simu_diff_mean_EH, pch=20)





####   Fig A4: Contrast between the locations of the observed migratory destinations of each species, and of the simulated alterative migration options that performe better than the observed


lon_brsim_best_simu_mean_EH[which(lon_brsim_best_simu_mean_EH=="NaN")] <- 0
lat_brsim_best_simu_mean_EH[which(lat_brsim_best_simu_mean_EH=="NaN")] <- 0
lon_nbsim_best_simu_mean_EH[which(lon_nbsim_best_simu_mean_EH=="NaN")] <- 0
lat_nbsim_best_simu_mean_EH[which(lat_nbsim_best_simu_mean_EH=="NaN")] <- 0
lon_brsim_best_simu_mean_WH[which(lon_brsim_best_simu_mean_WH=="NaN")] <- 0
lat_brsim_best_simu_mean_WH[which(lat_brsim_best_simu_mean_WH=="NaN")] <- 0
lon_nbsim_best_simu_mean_WH[which(lon_nbsim_best_simu_mean_WH=="NaN")] <- 0
lat_nbsim_best_simu_mean_WH[which(lat_nbsim_best_simu_mean_WH=="NaN")] <- 0


distas_brsim_WH <- vector()
for(i in 1:length(lon_brsim_best_simu_mean_WH)){
	distas_brsim_WH[i] <- dist(rbind(c(0,0), c(lon_brsim_best_simu_mean_WH[i], lat_brsim_best_simu_mean_WH[i])))
}
distas_nbsim_WH <- vector()
for(i in 1:length(lon_nbsim_best_simu_mean_WH)){
	distas_nbsim_WH[i] <- dist(rbind(c(0,0), c(lon_nbsim_best_simu_mean_WH[i], lat_nbsim_best_simu_mean_WH[i])))
}
distas_brsim_EH <- vector()
for(i in 1:length(lon_brsim_best_simu_mean_EH)){
	distas_brsim_EH[i] <- dist(rbind(c(0,0), c(lon_brsim_best_simu_mean_EH[i], lat_brsim_best_simu_mean_EH[i])))
}
distas_nbsim_EH <- vector()
for(i in 1:length(lon_nbsim_best_simu_mean_EH)){
	distas_nbsim_EH[i] <- dist(rbind(c(0,0), c(lon_nbsim_best_simu_mean_EH[i], lat_nbsim_best_simu_mean_EH[i])))
}

par(mfrow=c(1,2), mar=c(3.5,3.5,3.5,0.1), mgp=c(1.5,0.5,0), bg="grey")

plot(c(lon_brsim_best_simu_mean_WH, lon_brsim_best_simu_mean_EH), c(lat_brsim_best_simu_mean_WH, lat_brsim_best_simu_mean_EH), pch=20, col="yellow", cex=0.7, xlab="Longitude difference", ylab="Latitude difference", main="Simulated breeding grounds", xlim=c(-110,110), ylim=c(-80,80))
points(c(lon_brsim_best_simu_mean_WH[which(distas_brsim_WH >= 15)], lon_brsim_best_simu_mean_EH[which(distas_brsim_EH >= 15)]), c(lat_brsim_best_simu_mean_WH[which(distas_brsim_WH >= 15)], lat_brsim_best_simu_mean_EH[which(distas_brsim_EH >= 15)]), pch=20, col="orange", cex=0.7)
points(c(lon_brsim_best_simu_mean_WH[which(distas_brsim_WH >= 30)], lon_brsim_best_simu_mean_EH[which(distas_brsim_EH >= 30)]), c(lat_brsim_best_simu_mean_WH[which(distas_brsim_WH >= 30)], lat_brsim_best_simu_mean_EH[which(distas_brsim_EH >= 30)]), pch=20, col="red", cex=0.7)
points(c(lon_brsim_best_simu_mean_WH[which(distas_brsim_WH >= 45)], lon_brsim_best_simu_mean_EH[which(distas_brsim_EH >= 45)]), c(lat_brsim_best_simu_mean_WH[which(distas_brsim_WH >= 45)], lat_brsim_best_simu_mean_EH[which(distas_brsim_EH >= 45)]), pch=20, col="brown4", cex=0.7)
abline(h=0)
abline(v=0)

plot(c(lon_nbsim_best_simu_mean_WH, lon_nbsim_best_simu_mean_EH), c(lat_nbsim_best_simu_mean_WH, lat_nbsim_best_simu_mean_EH), pch=20, col="yellow", cex=0.7, xlab="Longitude difference", ylab="Latitude difference", main="Simulated non-breeding grounds", xlim=c(-110,110), ylim=c(-80,80))
points(c(lon_nbsim_best_simu_mean_WH[which(distas_nbsim_WH >= 15)], lon_nbsim_best_simu_mean_EH[which(distas_nbsim_EH >= 15)]), c(lat_nbsim_best_simu_mean_WH[which(distas_nbsim_WH >= 15)], lat_nbsim_best_simu_mean_EH[which(distas_nbsim_EH >= 15)]), pch=20, col="orange", cex=0.7)
points(c(lon_nbsim_best_simu_mean_WH[which(distas_nbsim_WH >= 30)], lon_nbsim_best_simu_mean_EH[which(distas_nbsim_EH >= 30)]), c(lat_nbsim_best_simu_mean_WH[which(distas_nbsim_WH >= 30)], lat_nbsim_best_simu_mean_EH[which(distas_nbsim_EH >= 30)]), pch=20, col="red", cex=0.7)
points(c(lon_nbsim_best_simu_mean_WH[which(distas_nbsim_WH >= 45)], lon_nbsim_best_simu_mean_EH[which(distas_nbsim_EH >= 45)]), c(lat_nbsim_best_simu_mean_WH[which(distas_nbsim_WH >= 45)], lat_nbsim_best_simu_mean_EH[which(distas_nbsim_EH >= 45)]), pch=20, col="brown4", cex=0.7)
abline(h=0)
abline(v=0)




centroids.breeding.grounds_WH <- rbind(centroids.breeding.groundsNH_WH, centroids.breeding.groundsSH_WH)
centroids.nonbreeding.grounds_WH <- rbind(centroids.nonbreeding.groundsNH_WH, centroids.nonbreeding.groundsSH_WH)
centroids.breeding.grounds_EH <- rbind(centroids.breeding.groundsNH_EH, centroids.breeding.groundsSH_EH)
centroids.nonbreeding.grounds_EH <- rbind(centroids.nonbreeding.groundsNH_EH, centroids.nonbreeding.groundsSH_EH)

centroids.simu_nbsim_WH <- c(centroids.simulatedNB_NH_WH, centroids.simulatedNB_SH_WH)
centroids.simu_brsim_WH <- c(centroids.simulatedBR_NH_WH, centroids.simulatedBR_SH_WH)
centroids.simu_nbsim_EH <- c(centroids.simulatedNB_NH_EH, centroids.simulatedNB_SH_EH)
centroids.simu_brsim_EH <- c(centroids.simulatedBR_NH_EH, centroids.simulatedBR_SH_EH)

centroids.simu_nbsim_sel_WH <- list()
centroids.simu_brsim_sel_WH <- list()
for(i in 1:length(geo_best_simu_nbsim_WH)){
	centroids.simu_nbsim_sel_WH[[i]] <- centroids.simu_nbsim_WH[[i]][geo_best_simu_nbsim_WH[[i]]]
	centroids.simu_brsim_sel_WH[[i]] <- centroids.simu_brsim_WH[[i]][geo_best_simu_brsim_WH[[i]]]
}
centroids.simu_nbsim_sel_EH <- list()
centroids.simu_brsim_sel_EH <- list()
for(i in 1:length(geo_best_simu_nbsim_EH)){
	centroids.simu_nbsim_sel_EH[[i]] <- centroids.simu_nbsim_EH[[i]][geo_best_simu_nbsim_EH[[i]]]
	centroids.simu_brsim_sel_EH[[i]] <- centroids.simu_brsim_EH[[i]][geo_best_simu_brsim_EH[[i]]]
}



####   Fig A5: Mapping species migrations according to the performance of the breeding range compared to simulated alternative migrations

par(mfrow=c(2,2), mar=c(0.1,0.1,0.1,0.1), mgp=c(1.5,0.5,0))

#plot(NULL, xlim=c(-180,180), ylim=c(-90,90))
plot(hexgrid2, col= "dark grey", border = "dark grey", bg="light grey")
for(k in 1:length(centroids.breeding.grounds_EH[,1])){
	dista = dist(rbind(c(0,0), c(lon_brsim_best_simu_mean_EH[k], lat_brsim_best_simu_mean_EH[k])))
	if(lon_brsim_best_simu_mean_EH[k] < 0 & lat_brsim_best_simu_mean_EH[k] > 0 & dista < 15){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_EH[k]/30, col="yellow")
	}
	if(lon_brsim_best_simu_mean_EH[k] < 0 & lat_brsim_best_simu_mean_EH[k] > 0 & dista >= 15 & dista < 30){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_EH[k]/30, col="orange")
	}
	if(lon_brsim_best_simu_mean_EH[k] < 0 & lat_brsim_best_simu_mean_EH[k] > 0 & dista >= 30 & dista < 45){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_EH[k]/30, col="red")
	}
	if(lon_brsim_best_simu_mean_EH[k] < 0 & lat_brsim_best_simu_mean_EH[k] > 0 & dista >= 45){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_EH[k]/30, col="brown4")
	}
}
for(k in 1:length(centroids.breeding.grounds_WH[,1])){
	dista = dist(rbind(c(0,0), c(lon_brsim_best_simu_mean_WH[k], lat_brsim_best_simu_mean_WH[k])))
	if(lon_brsim_best_simu_mean_WH[k] < 0 & lat_brsim_best_simu_mean_WH[k] > 0 & dista < 15){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_WH[k]/30, col="yellow")
	}
	if(lon_brsim_best_simu_mean_WH[k] < 0 & lat_brsim_best_simu_mean_WH[k] > 0 & dista >= 15 & dista < 30){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_WH[k]/30, col="orange")
	}
	if(lon_brsim_best_simu_mean_WH[k] < 0 & lat_brsim_best_simu_mean_WH[k] > 0 & dista >= 30 & dista < 45){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_WH[k]/30, col="red")
	}
	if(lon_brsim_best_simu_mean_WH[k] < 0 & lat_brsim_best_simu_mean_WH[k] > 0 & dista >= 45){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_WH[k]/30, col="brown4")
	}
}
plot(hexgrid2, col= "dark grey", border = "dark grey", bg="light grey")
for(k in 1:length(centroids.breeding.grounds_EH[,1])){
	dista = dist(rbind(c(0,0), c(lon_brsim_best_simu_mean_EH[k], lat_brsim_best_simu_mean_EH[k])))
	if(lon_brsim_best_simu_mean_EH[k] > 0 & lat_brsim_best_simu_mean_EH[k] > 0 & dista < 15){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_EH[k]/30, col="yellow")
	}
	if(lon_brsim_best_simu_mean_EH[k] > 0 & lat_brsim_best_simu_mean_EH[k] > 0 & dista >= 15 & dista < 30){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_EH[k]/30, col="orange")
	}
	if(lon_brsim_best_simu_mean_EH[k] > 0 & lat_brsim_best_simu_mean_EH[k] > 0 & dista >= 30 & dista < 45){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_EH[k]/30, col="red")
	}
	if(lon_brsim_best_simu_mean_EH[k] > 0 & lat_brsim_best_simu_mean_EH[k] > 0 & dista >= 45){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_EH[k]/30, col="brown4")
	}
}
for(k in 1:length(centroids.breeding.grounds_WH[,1])){
	dista = dist(rbind(c(0,0), c(lon_brsim_best_simu_mean_WH[k], lat_brsim_best_simu_mean_WH[k])))
	if(lon_brsim_best_simu_mean_WH[k] > 0 & lat_brsim_best_simu_mean_WH[k] > 0 & dista < 15){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_WH[k]/30, col="yellow")
	}
	if(lon_brsim_best_simu_mean_WH[k] > 0 & lat_brsim_best_simu_mean_WH[k] > 0 & dista >= 15 & dista < 30){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_WH[k]/30, col="orange")
	}
	if(lon_brsim_best_simu_mean_WH[k] > 0 & lat_brsim_best_simu_mean_WH[k] > 0 & dista >= 30 & dista < 45){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_WH[k]/30, col="red")
	}
	if(lon_brsim_best_simu_mean_WH[k] > 0 & lat_brsim_best_simu_mean_WH[k] > 0 & dista >= 45){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_WH[k]/30, col="brown4")
	}
}
plot(hexgrid2, col= "dark grey", border = "dark grey", bg="light grey")
for(k in 1:length(centroids.breeding.grounds_EH[,1])){
	dista = dist(rbind(c(0,0), c(lon_brsim_best_simu_mean_EH[k], lat_brsim_best_simu_mean_EH[k])))
	if(lon_brsim_best_simu_mean_EH[k] < 0 & lat_brsim_best_simu_mean_EH[k] < 0 & dista < 15){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_EH[k]/30, col="yellow")
	}
	if(lon_brsim_best_simu_mean_EH[k] < 0 & lat_brsim_best_simu_mean_EH[k] < 0 & dista >= 15 & dista < 30){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_EH[k]/30, col="orange")
	}
	if(lon_brsim_best_simu_mean_EH[k] < 0 & lat_brsim_best_simu_mean_EH[k] < 0 & dista >= 30 & dista < 45){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_EH[k]/30, col="red")
	}
	if(lon_brsim_best_simu_mean_EH[k] < 0 & lat_brsim_best_simu_mean_EH[k] < 0 & dista >= 45){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_EH[k]/30, col="brown4")
	}
}
for(k in 1:length(centroids.breeding.grounds_WH[,1])){
	dista = dist(rbind(c(0,0), c(lon_brsim_best_simu_mean_WH[k], lat_brsim_best_simu_mean_WH[k])))
	if(lon_brsim_best_simu_mean_WH[k] < 0 & lat_brsim_best_simu_mean_WH[k] < 0 & dista < 15){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_WH[k]/30, col="yellow")
	}
	if(lon_brsim_best_simu_mean_WH[k] < 0 & lat_brsim_best_simu_mean_WH[k] < 0 & dista >= 15 & dista < 30){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_WH[k]/30, col="orange")
	}
	if(lon_brsim_best_simu_mean_WH[k] < 0 & lat_brsim_best_simu_mean_WH[k] < 0 & dista >= 30 & dista < 45){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_WH[k]/30, col="red")
	}
	if(lon_brsim_best_simu_mean_WH[k] < 0 & lat_brsim_best_simu_mean_WH[k] < 0 & dista >= 45){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_WH[k]/30, col="brown4")
	}
}
plot(hexgrid2, col= "dark grey", border = "dark grey", bg="light grey")
for(k in 1:length(centroids.breeding.grounds_EH[,1])){
	dista = dist(rbind(c(0,0), c(lon_brsim_best_simu_mean_EH[k], lat_brsim_best_simu_mean_EH[k])))
	if(lon_brsim_best_simu_mean_EH[k] > 0 & lat_brsim_best_simu_mean_EH[k] < 0 & dista < 15){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_EH[k]/30, col="yellow")
	}
	if(lon_brsim_best_simu_mean_EH[k] > 0 & lat_brsim_best_simu_mean_EH[k] < 0 & dista >= 15 & dista < 30){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_EH[k]/30, col="orange")
	}
	if(lon_brsim_best_simu_mean_EH[k] > 0 & lat_brsim_best_simu_mean_EH[k] < 0 & dista >= 30 & dista < 45){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_EH[k]/30, col="red")
	}
	if(lon_brsim_best_simu_mean_EH[k] > 0 & lat_brsim_best_simu_mean_EH[k] < 0 & dista >= 45){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_EH[k]/30, col="brown4")
	}
}
for(k in 1:length(centroids.breeding.grounds_WH[,1])){
	dista = dist(rbind(c(0,0), c(lon_brsim_best_simu_mean_WH[k], lat_brsim_best_simu_mean_WH[k])))
	if(lon_brsim_best_simu_mean_WH[k] > 0 & lat_brsim_best_simu_mean_WH[k] < 0 & dista < 15){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_WH[k]/30, col="yellow")
	}
	if(lon_brsim_best_simu_mean_WH[k] > 0 & lat_brsim_best_simu_mean_WH[k] < 0 & dista >= 15 & dista < 30){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_WH[k]/30, col="orange")
	}
	if(lon_brsim_best_simu_mean_WH[k] > 0 & lat_brsim_best_simu_mean_WH[k] < 0 & dista >= 30 & dista < 45){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_WH[k]/30, col="red")
	}
	if(lon_brsim_best_simu_mean_WH[k] > 0 & lat_brsim_best_simu_mean_WH[k] < 0 & dista >= 45){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_brsim_number_WH[k]/30, col="brown4")
	}
}


####   Fig A5: Mapping species migrations according to the performance of the non-breeding range compared to simulated alternative migrations

par(mfrow=c(2,2), mar=c(0.1,0.1,0.1,0.1), mgp=c(1.5,0.5,0))

plot(hexgrid2, col= "dark grey", border = "dark grey", bg="light grey")
for(k in 1:length(centroids.breeding.grounds_EH[,1])){
	dista = dist(rbind(c(0,0), c(lon_nbsim_best_simu_mean_EH[k], lat_nbsim_best_simu_mean_EH[k])))
	if(lon_nbsim_best_simu_mean_EH[k] < 0 & lat_nbsim_best_simu_mean_EH[k] > 0 & dista < 15){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_EH[k]/30, col="yellow")
	}
	if(lon_nbsim_best_simu_mean_EH[k] < 0 & lat_nbsim_best_simu_mean_EH[k] > 0 & dista >= 15 & dista < 30){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_EH[k]/30, col="orange")
	}
	if(lon_nbsim_best_simu_mean_EH[k] < 0 & lat_nbsim_best_simu_mean_EH[k] > 0 & dista >= 30 & dista < 45){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_EH[k]/30, col="red")
	}
	if(lon_nbsim_best_simu_mean_EH[k] < 0 & lat_nbsim_best_simu_mean_EH[k] > 0 & dista >= 45){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_EH[k]/30, col="brown4")
	}
}
for(k in 1:length(centroids.breeding.grounds_WH[,1])){
	dista = dist(rbind(c(0,0), c(lon_nbsim_best_simu_mean_WH[k], lat_nbsim_best_simu_mean_WH[k])))
	if(lon_nbsim_best_simu_mean_WH[k] < 0 & lat_nbsim_best_simu_mean_WH[k] > 0 & dista < 15){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_WH[k]/30, col="yellow")
	}
	if(lon_nbsim_best_simu_mean_WH[k] < 0 & lat_nbsim_best_simu_mean_WH[k] > 0 & dista >= 15 & dista < 30){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_WH[k]/30, col="orange")
	}
	if(lon_nbsim_best_simu_mean_WH[k] < 0 & lat_nbsim_best_simu_mean_WH[k] > 0 & dista >= 30 & dista < 45){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_WH[k]/30, col="red")
	}
	if(lon_nbsim_best_simu_mean_WH[k] < 0 & lat_nbsim_best_simu_mean_WH[k] > 0 & dista >= 45){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_WH[k]/30, col="brown4")
	}
}
plot(hexgrid2, col= "dark grey", border = "dark grey", bg="light grey")
for(k in 1:length(centroids.breeding.grounds_EH[,1])){
	dista = dist(rbind(c(0,0), c(lon_nbsim_best_simu_mean_EH[k], lat_nbsim_best_simu_mean_EH[k])))
	if(lon_nbsim_best_simu_mean_EH[k] > 0 & lat_nbsim_best_simu_mean_EH[k] > 0 & dista < 15){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_EH[k]/30, col="yellow")
	}
	if(lon_nbsim_best_simu_mean_EH[k] > 0 & lat_nbsim_best_simu_mean_EH[k] > 0 & dista >= 15 & dista < 30){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_EH[k]/30, col="orange")
	}
	if(lon_nbsim_best_simu_mean_EH[k] > 0 & lat_nbsim_best_simu_mean_EH[k] > 0 & dista >= 30 & dista < 45){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_EH[k]/30, col="red")
	}
	if(lon_nbsim_best_simu_mean_EH[k] > 0 & lat_nbsim_best_simu_mean_EH[k] > 0 & dista >= 45){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_EH[k]/30, col="brown4")
	}
}
for(k in 1:length(centroids.breeding.grounds_WH[,1])){
	dista = dist(rbind(c(0,0), c(lon_nbsim_best_simu_mean_WH[k], lat_nbsim_best_simu_mean_WH[k])))
	if(lon_nbsim_best_simu_mean_WH[k] > 0 & lat_nbsim_best_simu_mean_WH[k] > 0 & dista < 15){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_WH[k]/30, col="yellow")
	}
	if(lon_nbsim_best_simu_mean_WH[k] > 0 & lat_nbsim_best_simu_mean_WH[k] > 0 & dista >= 15 & dista < 30){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_WH[k]/30, col="orange")
	}
	if(lon_nbsim_best_simu_mean_WH[k] > 0 & lat_nbsim_best_simu_mean_WH[k] > 0 & dista >= 30 & dista < 45){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_WH[k]/30, col="red")
	}
	if(lon_nbsim_best_simu_mean_WH[k] > 0 & lat_nbsim_best_simu_mean_WH[k] > 0 & dista >= 45){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_WH[k]/30, col="brown4")
	}
}
plot(hexgrid2, col= "dark grey", border = "dark grey", bg="light grey")
for(k in 1:length(centroids.breeding.grounds_EH[,1])){
	dista = dist(rbind(c(0,0), c(lon_nbsim_best_simu_mean_EH[k], lat_nbsim_best_simu_mean_EH[k])))
	if(lon_nbsim_best_simu_mean_EH[k] < 0 & lat_nbsim_best_simu_mean_EH[k] < 0 & dista < 15){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_EH[k]/30, col="yellow")
	}
	if(lon_nbsim_best_simu_mean_EH[k] < 0 & lat_nbsim_best_simu_mean_EH[k] < 0 & dista >= 15 & dista < 30){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_EH[k]/30, col="orange")
	}
	if(lon_nbsim_best_simu_mean_EH[k] < 0 & lat_nbsim_best_simu_mean_EH[k] < 0 & dista >= 30 & dista < 45){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_EH[k]/30, col="red")
	}
	if(lon_nbsim_best_simu_mean_EH[k] < 0 & lat_nbsim_best_simu_mean_EH[k] < 0 & dista >= 45){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_EH[k]/30, col="brown4")
	}
}
for(k in 1:length(centroids.breeding.grounds_WH[,1])){
	dista = dist(rbind(c(0,0), c(lon_nbsim_best_simu_mean_WH[k], lat_nbsim_best_simu_mean_WH[k])))
	if(lon_nbsim_best_simu_mean_WH[k] < 0 & lat_nbsim_best_simu_mean_WH[k] < 0 & dista < 15){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_WH[k]/30, col="yellow")
	}
	if(lon_nbsim_best_simu_mean_WH[k] < 0 & lat_nbsim_best_simu_mean_WH[k] < 0 & dista >= 15 & dista < 30){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_WH[k]/30, col="orange")
	}
	if(lon_nbsim_best_simu_mean_WH[k] < 0 & lat_nbsim_best_simu_mean_WH[k] < 0 & dista >= 30 & dista < 45){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_WH[k]/30, col="red")
	}
	if(lon_nbsim_best_simu_mean_WH[k] < 0 & lat_nbsim_best_simu_mean_WH[k] < 0 & dista >= 45){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_WH[k]/30, col="brown4")
	}
}
plot(hexgrid2, col= "dark grey", border = "dark grey", bg="light grey")
for(k in 1:length(centroids.breeding.grounds_EH[,1])){
	dista = dist(rbind(c(0,0), c(lon_nbsim_best_simu_mean_EH[k], lat_nbsim_best_simu_mean_EH[k])))
	if(lon_nbsim_best_simu_mean_EH[k] > 0 & lat_nbsim_best_simu_mean_EH[k] < 0 & dista < 15){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_EH[k]/30, col="yellow")
	}
	if(lon_nbsim_best_simu_mean_EH[k] > 0 & lat_nbsim_best_simu_mean_EH[k] < 0 & dista >= 15 & dista < 30){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_EH[k]/30, col="orange")
	}
	if(lon_nbsim_best_simu_mean_EH[k] > 0 & lat_nbsim_best_simu_mean_EH[k] < 0 & dista >= 30 & dista < 45){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_EH[k]/30, col="red")
	}
	if(lon_nbsim_best_simu_mean_EH[k] > 0 & lat_nbsim_best_simu_mean_EH[k] < 0 & dista >= 45){
		inter <- gcIntermediate(centroids.breeding.grounds_EH[k,], centroids.nonbreeding.grounds_EH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_EH[k]/30, col="brown4")
	}
}
for(k in 1:length(centroids.breeding.grounds_WH[,1])){
	dista = dist(rbind(c(0,0), c(lon_nbsim_best_simu_mean_WH[k], lat_nbsim_best_simu_mean_WH[k])))
	if(lon_nbsim_best_simu_mean_WH[k] > 0 & lat_nbsim_best_simu_mean_WH[k] < 0 & dista < 15){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_WH[k]/30, col="yellow")
	}
	if(lon_nbsim_best_simu_mean_WH[k] > 0 & lat_nbsim_best_simu_mean_WH[k] < 0 & dista >= 15 & dista < 30){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_WH[k]/30, col="orange")
	}
	if(lon_nbsim_best_simu_mean_WH[k] > 0 & lat_nbsim_best_simu_mean_WH[k] < 0 & dista >= 30 & dista < 45){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_WH[k]/30, col="red")
	}
	if(lon_nbsim_best_simu_mean_WH[k] > 0 & lat_nbsim_best_simu_mean_WH[k] < 0 & dista >= 45){
		inter <- gcIntermediate(centroids.breeding.grounds_WH[k,], centroids.nonbreeding.grounds_WH[k,], n=50, addStartEnd=T)
		lines(inter, lwd= best_simu_nbsim_number_WH[k]/30, col="brown4")
	}
}




















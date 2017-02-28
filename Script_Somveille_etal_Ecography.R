library(maptools)

library(sp)
library(rgdal)
library(ape)
library(nlme)






library(geosphere)
library(fields)
library(MASS)
library(raster)
library(move)
library(vioplot)
library(igraph)

setwd("/Users/mariussomveille/Desktop/PhD/Chapter 3 – niche/Submission process/Ecography/Data")




################################################################
################################################################

#######      Prepare species and environmental data    #########

################################################################
################################################################



#Presence-absence data for migratory bird species
PresAbs_BR_NH <- read.table("PresAbs_BR_NH.txt", header=T)  # Breeding range of species breeding during the northern summer
PresAbs_NB_NH <- read.table("PresAbs_NB_NH.txt", header=T)  # Non-breeding range of species breeding during the northern summer
PresAbs_BR_SH <- read.table("PresAbs_BR_SH.txt", header=T)  # Breeding range of species breeding during the northern winter
PresAbs_NB_SH <- read.table("PresAbs_NB_SH.txt", header=T)  # Non-breeding range of species breeding during the northern winter

#Remove species that do not have at least 50% of either their breeding or non-breeding range on land
sppinfo <- read.csv("List_species_both&coastal.csv", sep=";", header=T)
PresAbs_BR_NH <- PresAbs_BR_NH[,-match(c(as.character(sppinfo$CoastalBothBoth),as.character(sppinfo$CoastalNHEH),as.character(sppinfo$CoastalNHWH),as.character(sppinfo$CoastalBothEH)), colnames(PresAbs_BR_NH), nomatch=0)]
PresAbs_NB_NH <- PresAbs_NB_NH[,-match(c(as.character(sppinfo$CoastalBothBoth),as.character(sppinfo$CoastalNHEH),as.character(sppinfo$CoastalNHWH),as.character(sppinfo$CoastalBothEH)), colnames(PresAbs_NB_NH), nomatch=0)]
PresAbs_BR_SH <- PresAbs_BR_SH[,-match(c(as.character(sppinfo$CoastalSHEH),as.character(sppinfo$CoastalSHWH)), colnames(PresAbs_BR_SH), nomatch=0)]
PresAbs_NB_SH <- PresAbs_NB_SH[,-match(c(as.character(sppinfo$CoastalSHEH),as.character(sppinfo$CoastalSHWH)), colnames(PresAbs_NB_SH), nomatch=0)]


#Environmental data
envData <- read.csv("Env_data.csv")
envData2 <- envData[-which(envData$Tmean_NW==0 & envData$Tmean_NS==0),]
PrecNS <- log(envData2$Prec_NS+1)
PrecNW <- log(envData2$Prec_NW+1)
TempNS <- envData2$Tmean_NS/10
TempNW <- envData2$Tmean_NW/10
NDVI_NS <- ((envData2$NDVI_may + envData2$NDVI_jun + envData2$NDVI_jul + envData2$NDVI_aug)/4) * 100
NDVI_NW <- ((envData2$NDVI_nov + envData2$NDVI_dec + envData2$NDVI_jan + envData2$NDVI_feb)/4) * 100
NDVI_NS[NDVI_NS<1] <- 1
NDVI_NW[NDVI_NW<1] <- 1


#z-tranform climate data to make axes comparable
sdTemp <- sd(c(TempNS,TempNW))
mnTemp <- mean(c(TempNS,TempNW))
sdPrec <- sd(c(PrecNS,PrecNW))
mnPrec <- mean(c(PrecNS,PrecNW))
TempNS <- (TempNS - mnTemp) / sdTemp
TempNW <- (TempNW - mnTemp) / sdTemp
PrecNS <- (PrecNS - mnPrec) / sdPrec
PrecNW <- (PrecNW - mnPrec) / sdPrec


#Split environmental data into western hemisphere (WH) and eastern hemisphere (EH)
TempNS_WH <- TempNS[which(envData2$LONGITUDE<=-30)]
TempNW_WH <- TempNW[which(envData2$LONGITUDE<=-30)]
PrecNS_WH <- PrecNS[which(envData2$LONGITUDE<=-30)]
PrecNW_WH <- PrecNW[which(envData2$LONGITUDE<=-30)]
NDVI_NS_WH <- NDVI_NS[which(envData2$LONGITUDE<=-30)]
NDVI_NW_WH <- NDVI_NW[which(envData2$LONGITUDE<=-30)]
resourceGain_NS_WH <- NDVI_NS_WH - NDVI_NW_WH
resourceGain_NW_WH <- NDVI_NW_WH - NDVI_NS_WH
resourceScarcity_NS_WH <- -resourceGain_NS_WH
resourceScarcity_NW_WH <- -resourceGain_NW_WH
TempNS_EH <- TempNS[which(envData2$LONGITUDE>-30)]
TempNW_EH <- TempNW[which(envData2$LONGITUDE>-30)]
PrecNS_EH <- PrecNS[which(envData2$LONGITUDE>-30)]
PrecNW_EH <- PrecNW[which(envData2$LONGITUDE>-30)]
NDVI_NS_EH <- NDVI_NS[which(envData2$LONGITUDE>-30)]
NDVI_NW_EH <- NDVI_NW[which(envData2$LONGITUDE>-30)]
resourceGain_NS_EH <- NDVI_NS_EH - NDVI_NW_EH
resourceGain_NW_EH <- NDVI_NW_EH - NDVI_NS_EH
resourceScarcity_NS_EH <- -resourceGain_NS_EH
resourceScarcity_NW_EH <- -resourceGain_NW_EH


#Geographic coordinates
lonlat <- cbind(envData2$LONGITUDE, envData2$LATITUDE)
colnames(lonlat) <- c("LONGITUDE", "LATITUDE")
west_Hem <- lonlat[which(lonlat[,1]<=-30),]  # subset for Western Hemisphere
east_Hem <- lonlat[which(lonlat[,1]>-30),]   # subset for Eastern Hemisphere


#Clean presence-absence data
PresAbs_BR_NH <- PresAbs_BR_NH[-which(envData$Tmean_NW==0 & envData$Tmean_NS==0),]
PresAbs_NB_NH <- PresAbs_NB_NH[-which(envData$Tmean_NW==0 & envData$Tmean_NS==0),]
PresAbs_BR_SH <- PresAbs_BR_SH[-which(envData$Tmean_NW==0 & envData$Tmean_NS==0),]
PresAbs_NB_SH <- PresAbs_NB_SH[-which(envData$Tmean_NW==0 & envData$Tmean_NS==0),]

PresAbs_BR_NH_EH <- PresAbs_BR_NH[match(envData2[which(envData2$LONGITUDE>-30),1], PresAbs_BR_NH[,1]),]
PresAbs_NB_NH_EH <- PresAbs_NB_NH[match(envData2[which(envData2$LONGITUDE>-30),1], PresAbs_NB_NH[,1]),]
PresAbs_BR_NH_WH <- PresAbs_BR_NH[match(envData2[which(envData2$LONGITUDE<=-30),1], PresAbs_BR_NH[,1]),]
PresAbs_NB_NH_WH <- PresAbs_NB_NH[match(envData2[which(envData2$LONGITUDE<=-30),1], PresAbs_NB_NH[,1]),]

PresAbs_BR_SH_EH <- PresAbs_BR_SH[match(envData2[which(envData2$LONGITUDE>-30),1], PresAbs_BR_SH[,1]),]
PresAbs_NB_SH_EH <- PresAbs_NB_SH[match(envData2[which(envData2$LONGITUDE>-30),1], PresAbs_NB_SH[,1]),]
PresAbs_BR_SH_WH <- PresAbs_BR_SH[match(envData2[which(envData2$LONGITUDE<=-30),1], PresAbs_BR_SH[,1]),]
PresAbs_NB_SH_WH <- PresAbs_NB_SH[match(envData2[which(envData2$LONGITUDE<=-30),1], PresAbs_NB_SH[,1]),]



#Computing range size
rangeSize_migra_BR_NH_WH <- apply(PresAbs_BR_NH_WH[,-1], 2, sum)
rangeSize_migra_BR_SH_WH <- apply(PresAbs_BR_SH_WH[,-1], 2, sum)
rangeSize_migra_NB_NH_WH <- apply(PresAbs_NB_NH_WH[,-1], 2, sum)
rangeSize_migra_NB_SH_WH <- apply(PresAbs_NB_SH_WH[,-1], 2, sum)
rangeSize_migra_BR_NH_EH <- apply(PresAbs_BR_NH_EH[,-1], 2, sum)
rangeSize_migra_BR_SH_EH <- apply(PresAbs_BR_SH_EH[,-1], 2, sum)
rangeSize_migra_NB_NH_EH <- apply(PresAbs_NB_NH_EH[,-1], 2, sum)
rangeSize_migra_NB_SH_EH <- apply(PresAbs_NB_SH_EH[,-1], 2, sum)
rangeSizes_WH <- c(rangeSize_migra_BR_NH_WH, rangeSize_migra_BR_SH_WH, rangeSize_migra_NB_NH_WH, rangeSize_migra_NB_SH_WH)
rangeSizes_EH <- c(rangeSize_migra_BR_NH_EH, rangeSize_migra_BR_SH_EH, rangeSize_migra_NB_NH_EH, rangeSize_migra_NB_SH_EH)


#Computing range overlap
perm_NH_WH <- PresAbs_BR_NH_WH * PresAbs_NB_NH_WH
perm_NH_WH <- apply(perm_NH_WH[,-1], 2, sum)
rangeOverlap_NH_WH <- perm_NH_WH / (rangeSize_migra_BR_NH_WH + rangeSize_migra_NB_NH_WH - perm_NH_WH)
rangeOverlap_NH_WH  <- rangeOverlap_NH_WH[-which(rangeOverlap_NH_WH == "NaN")]
perm_SH_WH <- PresAbs_BR_SH_WH * PresAbs_NB_SH_WH
perm_SH_WH <- apply(perm_SH_WH[,-1], 2, sum)
rangeOverlap_SH_WH <- perm_SH_WH / (rangeSize_migra_BR_SH_WH + rangeSize_migra_NB_SH_WH - perm_SH_WH)
rangeOverlap_SH_WH  <- rangeOverlap_SH_WH[-which(rangeOverlap_SH_WH == "NaN")]

perm_NH_EH <- PresAbs_BR_NH_EH * PresAbs_NB_NH_EH
perm_NH_EH <- apply(perm_NH_EH[,-1], 2, sum)
rangeOverlap_NH_EH <- perm_NH_EH / (rangeSize_migra_BR_NH_EH + rangeSize_migra_NB_NH_EH - perm_NH_EH)
rangeOverlap_NH_EH  <- rangeOverlap_NH_EH[-which(rangeOverlap_NH_EH == "NaN")]
perm_SH_EH <- PresAbs_BR_SH_EH * PresAbs_NB_SH_EH
perm_SH_EH <- apply(perm_SH_EH[,-1], 2, sum)
rangeOverlap_SH_EH <- perm_SH_EH / (rangeSize_migra_BR_SH_EH + rangeSize_migra_NB_SH_EH - perm_SH_EH)
rangeOverlap_SH_EH  <- rangeOverlap_SH_EH[-which(rangeOverlap_SH_EH == "NaN")]


#Select migratory species with less than 20% overlap between their breeding and non-breeding ranges
selectedSpeciesNH_EH <- which(rangeOverlap_NH_EH < 0.2)
selectedSpeciesSH_EH <- which(rangeOverlap_SH_EH < 0.2)
selectedSpeciesNH_WH <- which(rangeOverlap_NH_WH < 0.2)
selectedSpeciesSH_WH <- which(rangeOverlap_SH_WH < 0.2)

#correct mistakes in coding seasons
#Anthus hoeschi
selectedSpeciesSH_EH <-  c(selectedSpeciesSH_EH, selectedSpeciesNH_EH[29])
selectedSpeciesNH_EH <- selectedSpeciesNH_EH[-29]
#Gallinago.nigripennis
selectedSpeciesNH_EH <- selectedSpeciesNH_EH[-159]
#Pitta.moluccensis and Sylvia.cantillans 
selectedSpeciesNH_EH <-  c(selectedSpeciesNH_EH, selectedSpeciesSH_EH[c(12,13)])
selectedSpeciesSH_EH <- selectedSpeciesSH_EH[-c(12,13)]


#Keep the presence-absence data only of selected species 
PresAbs_BR_NH_WH2 <- PresAbs_BR_NH_WH[,match(names(selectedSpeciesNH_WH), colnames(PresAbs_BR_NH_WH))]
PresAbs_NB_NH_WH2 <- PresAbs_NB_NH_WH[,match(names(selectedSpeciesNH_WH), colnames(PresAbs_NB_NH_WH))]
PresAbs_BR_NH_EH2 <- PresAbs_BR_NH_EH[,match(names(selectedSpeciesNH_EH)[1:380], colnames(PresAbs_BR_NH_EH))]
PresAbs_NB_NH_EH2 <- PresAbs_NB_NH_EH[,match(names(selectedSpeciesNH_EH)[1:380], colnames(PresAbs_NB_NH_EH))]
PresAbs_BR_NH_EH2 <- cbind(PresAbs_BR_NH_EH2, PresAbs_BR_SH_EH[,match(names(selectedSpeciesNH_EH)[381:382], colnames(PresAbs_BR_SH_EH))])
PresAbs_NB_NH_EH2 <- cbind(PresAbs_NB_NH_EH2, PresAbs_NB_SH_EH[,match(names(selectedSpeciesNH_EH)[381:382], colnames(PresAbs_NB_SH_EH))])
PresAbs_BR_SH_WH2 <- PresAbs_BR_SH_WH[,match(names(selectedSpeciesSH_WH), colnames(PresAbs_BR_SH_WH))]
PresAbs_NB_SH_WH2 <- PresAbs_NB_SH_WH[,match(names(selectedSpeciesSH_WH), colnames(PresAbs_NB_SH_WH))]
PresAbs_BR_SH_EH2 <- PresAbs_BR_SH_EH[,match(names(selectedSpeciesSH_EH)[1:11], colnames(PresAbs_BR_SH_EH))]
PresAbs_NB_SH_EH2 <- PresAbs_NB_SH_EH[,match(names(selectedSpeciesSH_EH)[1:11], colnames(PresAbs_NB_SH_EH))]
PresAbs_BR_SH_EH2 <- cbind(PresAbs_BR_SH_EH2, PresAbs_BR_NH_EH[,match(names(selectedSpeciesSH_EH)[12], colnames(PresAbs_BR_NH_EH))])
PresAbs_NB_SH_EH2 <- cbind(PresAbs_NB_SH_EH2, PresAbs_NB_NH_EH[,match(names(selectedSpeciesSH_EH)[12], colnames(PresAbs_NB_NH_EH))])

#Remove species with 0 presences (or only 1 because it creates problems for the kernel estimation) for at least one season
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


#Species migrating in both longitudinal hemispheres
#selectedSpeciesNH_both <- selectedSpeciesNH_EH[which(is.element(names(selectedSpeciesNH_EH), names(selectedSpeciesNH_WH)))]
#aa = apply(PresAbs_BR_NH_WH[,match(names(selectedSpeciesNH_both), colnames(PresAbs_BR_NH_WH))], 2, sum)
#bb = apply(PresAbs_NB_NH_WH[,match(names(selectedSpeciesNH_both), colnames(PresAbs_NB_NH_WH))], 2, sum)
#cc = apply(PresAbs_BR_NH_EH[,match(names(selectedSpeciesNH_both), colnames(PresAbs_BR_NH_EH))], 2, sum)
#dd = apply(PresAbs_NB_NH_EH[,match(names(selectedSpeciesNH_both), colnames(PresAbs_NB_NH_EH))], 2, sum)
#which(aa>0 & bb>0 & cc>0 & dd>0)
#names(selectedSpeciesNH_both)





################################################################
################################################################

######      Investigate trade-offs between variables    ########

################################################################
################################################################


#Convert the temperature and precipitation values recorded for presences into a density raster
nicheDensityRaster <- function(seasonalNiche){
	niche.kernel <- kde2d(seasonalNiche[,1], seasonalNiche[,2], n=20, h=1, lims=c(-3.7,1.5, -2.7,2.5)) 
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

#Compute niche distances (using the earth mover's distance) for every migratory species
breeding.nichesNH_WH <- apply(PresAbs_BR_NH_WH2, 2, function(x) nicheDensityRaster(cbind(TempNS_WH[which(x==1)], PrecNS_WH[which(x==1)])))
nonbreeding.nichesNH_WH <- apply(PresAbs_NB_NH_WH2, 2, function(x) nicheDensityRaster(cbind(TempNW_WH[which(x==1)], PrecNW_WH[which(x==1)])))
breeding.nichesSH_WH <- apply(PresAbs_BR_SH_WH2, 2, function(x) nicheDensityRaster(cbind(TempNW_WH[which(x==1)], PrecNW_WH[which(x==1)])))
nonbreeding.nichesSH_WH <- apply(PresAbs_NB_SH_WH2, 2, function(x) nicheDensityRaster(cbind(TempNS_WH[which(x==1)], PrecNS_WH[which(x==1)])))
breeding.nichesNH_EH <- apply(PresAbs_BR_NH_EH2, 2, function(x) nicheDensityRaster(cbind(TempNS_EH[which(x==1)], PrecNS_EH[which(x==1)])))
nonbreeding.nichesNH_EH <- apply(PresAbs_NB_NH_EH2, 2, function(x) nicheDensityRaster(cbind(TempNW_EH[which(x==1)], PrecNW_EH[which(x==1)])))
breeding.nichesSH_EH <- apply(PresAbs_BR_SH_EH2, 2, function(x) nicheDensityRaster(cbind(TempNW_EH[which(x==1)], PrecNW_EH[which(x==1)])))
nonbreeding.nichesSH_EH <- apply(PresAbs_NB_SH_EH2, 2, function(x) nicheDensityRaster(cbind(TempNS_EH[which(x==1)], PrecNS_EH[which(x==1)])))

niche.distancesNH_WH <- mapply(function(X,Y){ emd(X,Y) }, X=breeding.nichesNH_WH, Y=nonbreeding.nichesNH_WH)
niche.distancesSH_WH <- mapply(function(X,Y){ emd(X,Y) }, X=breeding.nichesSH_WH, Y=nonbreeding.nichesSH_WH)
niche.distancesNH_EH <- mapply(function(X,Y){ emd(X,Y) }, X=breeding.nichesNH_EH, Y=nonbreeding.nichesNH_EH)
niche.distancesSH_EH <- mapply(function(X,Y){ emd(X,Y) }, X=breeding.nichesSH_EH, Y=nonbreeding.nichesSH_EH)
niche.distances.obs = c(niche.distancesNH_WH, niche.distancesNH_EH, niche.distancesSH_WH, niche.distancesSH_EH)


#Computing geographical (migration) distance for every migratory species
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


#Computing resource scarcity for every migratory species
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


### Figure 3 - trade-off ### 

par(mfrow=c(2,2), mar=c(2.5,3,1.5,3), mgp=c(1.5,0.5,0))
# Relationship between geographical distance and niche distance
plot(sqrt(geo.distances.obs), niche.distances.obs, xlim=c(0,120), ylim=c(0,3.1), xlab="Geographical distance (square-root)", ylab="Niche distance", main="", xaxt="n", axes=F, pch=20, cex=0.7)
axis(side=2)
axis(side=1)
mod = lm(niche.distances.obs ~ sqrt(geo.distances.obs) + I(sqrt(geo.distances.obs)^2))
lines(sort(sqrt(geo.distances.obs)), fitted(mod)[order(sqrt(geo.distances.obs))], type="l", col="red")
mtext("A", cex=1.3, side=3, line=0, at=-22)
mtext(bquote(R^2 == .(round(summary(mod)$r.squared,2))), cex=1, side=3, line=-1.6, at=30)
# Relationship between niche distance and resource scarcity
plot(niche.distances.obs, resource.scarcity.obs, xlim=c(0,3.1), ylim=c(-63,21), xlab="Niche distance", ylab="Resource scarcity", main="", xaxt="n", axes=F, pch=20, cex=0.7)
axis(side=2)
axis(side=1)
mod = lm(resource.scarcity.obs ~ niche.distances.obs + I(niche.distances.obs^2))
lines(sort(niche.distances.obs), fitted(mod)[order(niche.distances.obs)], type="l", col="red")
mtext("B", cex=1.3, side=3, line=0, at=-0.55)
mtext(bquote(R^2 == .(round(summary(mod)$r.squared,2))), cex=1, side=3, line=-6, at=2.6)
# Relationship between geographical distance and resource scarcity
plot(sqrt(geo.distances.obs), resource.scarcity.obs, xlim=c(0,120), ylim=c(-63,21), xlab="Geographical distance (square-root)", ylab="Resource scarcity", main="", xaxt="n", axes=F, pch=20, cex=0.7)
axis(side=2)
axis(side=1)
mod = lm(resource.scarcity.obs ~ sqrt(geo.distances.obs) + I(sqrt(geo.distances.obs)^2))
lines(sort(sqrt(geo.distances.obs)), fitted(mod)[order(sqrt(geo.distances.obs))], type="l", col="red")
mtext("C", cex=1.3, side=3, line=0, at=-22)
mtext(bquote(R^2 == .(round(summary(mod)$r.squared,2))), cex=1, side=3, line=-1.5, at=90)
# Relationship between the three factors
hist(sqrt(geo.distances.obs), xlim=c(0,120), xlab="", ylab="", main="", xaxt="n", axes=F, col="light grey", border="grey")
axis(side=4)
axis(side=1)
par(new=T, mar=c(2.5,3,1.5,3), mgp=c(1.5,0.5,0))
plot(sqrt(geo.distances.obs), resource.scarcity.obs, xlim=c(0,120), ylim=c(-63,21), xlab="Geographical distance (square-root)", ylab="Resource scarcity", main="", xaxt="n", axes=F, pch=20, cex=0.7, col="yellow")
axis(side=2)
points(sqrt(geo.distances.obs)[which(niche.distances.obs < quantile(niche.distances.obs)[4])], resource.scarcity.obs[which(niche.distances.obs < quantile(niche.distances.obs)[4])], col="orange", cex=0.7, pch=20)
points(sqrt(geo.distances.obs)[which(niche.distances.obs < quantile(niche.distances.obs)[3])], resource.scarcity.obs[which(niche.distances.obs < quantile(niche.distances.obs)[3])], col="red", cex=0.7, pch=20)
points(sqrt(geo.distances.obs)[which(niche.distances.obs < quantile(niche.distances.obs)[2])], resource.scarcity.obs[which(niche.distances.obs < quantile(niche.distances.obs)[2])], col="brown4", cex=0.7, pch=20)
mtext("D", cex=1.3, side=3, line=0, at=-22)
mtext("Number of species", side=4, line=1.4, cex=0.85,las=0)
legend("topright", inset=.005, bg="white", box.col="white", title="Niche distance", c(paste("<", round(quantile(niche.distances.obs)[2], 2), sep=" "), paste(round(quantile(niche.distances.obs)[2], 2), round(quantile(niche.distances.obs)[3], 2), sep="–"), paste(round(quantile(niche.distances.obs)[3], 2), round(quantile(niche.distances.obs)[4], 2), sep="–"), paste(">", round(quantile(niche.distances.obs)[4], 2), sep=" ")), fill=c("brown4", "red", "orange", "yellow"), cex=0.8)






################################################################
################################################################

#######      Comparing to if species did not migrate    ########

################################################################
################################################################


#Niche distance if migratory species were resident
breeding.niches.residentinNB.NH_WH <- apply(PresAbs_NB_NH_WH2, 2, function(x) nicheDensityRaster(cbind(TempNS_WH[which(x==1)], PrecNS_WH[which(x==1)])))
nonbreeding.niches.residentinBR.NH_WH <- apply(PresAbs_BR_NH_WH2, 2, function(x) nicheDensityRaster(cbind(TempNW_WH[which(x==1)], PrecNW_WH[which(x==1)])))
breeding.niches.residentinNB.SH_WH <- apply(PresAbs_NB_SH_WH2, 2, function(x) nicheDensityRaster(cbind(TempNW_WH[which(x==1)], PrecNW_WH[which(x==1)])))
nonbreeding.niches.residentinBR.SH_WH <- apply(PresAbs_BR_SH_WH2, 2, function(x) nicheDensityRaster(cbind(TempNS_WH[which(x==1)], PrecNS_WH[which(x==1)])))
breeding.niches.residentinNB.NH_EH <- apply(PresAbs_NB_NH_EH2, 2, function(x) nicheDensityRaster(cbind(TempNS_EH[which(x==1)], PrecNS_EH[which(x==1)])))
nonbreeding.niches.residentinBR.NH_EH <- apply(PresAbs_BR_NH_EH2, 2, function(x) nicheDensityRaster(cbind(TempNW_EH[which(x==1)], PrecNW_EH[which(x==1)])))
breeding.niches.residentinNB.SH_EH <- apply(PresAbs_NB_SH_EH2, 2, function(x) nicheDensityRaster(cbind(TempNW_EH[which(x==1)], PrecNW_EH[which(x==1)])))
nonbreeding.niches.residentinBR.SH_EH <- apply(PresAbs_BR_SH_EH2, 2, function(x) nicheDensityRaster(cbind(TempNS_EH[which(x==1)], PrecNS_EH[which(x==1)])))

niche.distances.residentinNB.NH_WH <- mapply(function(X,Y){ emd(X,Y) }, X=breeding.niches.residentinNB.NH_WH, Y=nonbreeding.nichesNH_WH)
niche.distances.residentinBR.NH_WH <- mapply(function(X,Y){ emd(X,Y) }, X=breeding.nichesNH_WH, Y=nonbreeding.niches.residentinBR.NH_WH)
niche.distances.residentinNB.SH_WH <- mapply(function(X,Y){ emd(X,Y) }, X=breeding.niches.residentinNB.SH_WH, Y=nonbreeding.nichesSH_WH)
niche.distances.residentinBR.SH_WH <- mapply(function(X,Y){ emd(X,Y) }, X=breeding.nichesSH_WH, Y=nonbreeding.niches.residentinBR.SH_WH)
niche.distances.residentinNB.NH_EH <- mapply(function(X,Y){ emd(X,Y) }, X=breeding.niches.residentinNB.NH_EH, Y=nonbreeding.nichesNH_EH)
niche.distances.residentinBR.NH_EH <- mapply(function(X,Y){ emd(X,Y) }, X=breeding.nichesNH_EH, Y=nonbreeding.niches.residentinBR.NH_EH)
niche.distances.residentinNB.SH_EH <- mapply(function(X,Y){ emd(X,Y) }, X=breeding.niches.residentinNB.SH_EH, Y=nonbreeding.nichesSH_EH)
niche.distances.residentinBR.SH_EH <- mapply(function(X,Y){ emd(X,Y) }, X=breeding.nichesSH_EH, Y=nonbreeding.niches.residentinBR.SH_EH)
niche.distances.residentinNB = c(niche.distances.residentinNB.NH_WH, niche.distances.residentinNB.NH_EH, niche.distances.residentinNB.SH_WH, niche.distances.residentinNB.SH_EH)
niche.distances.residentinBR = c(niche.distances.residentinBR.NH_WH, niche.distances.residentinBR.NH_EH, niche.distances.residentinBR.SH_WH, niche.distances.residentinBR.SH_EH)


### Figure 4 - comparing niche distance to if resident ###

par(mfrow=c(1,1), mar=c(3,4,1.5,3), mgp=c(1.5,0.5,0), las=1, bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=range(c((niche.distances.obs - niche.distances.residentinNB), (niche.distances.obs - niche.distances.residentinBR))), axes=F, ann=F)
vioplot((niche.distances.obs - niche.distances.residentinNB), (niche.distances.obs - niche.distances.residentinBR), col="grey", names=c("Versus staying resident\nin non-breeding grounds", "Versus staying resident\nin breeding grounds"), add=T)
abline(h=0, lty = 3)
axis(side=1, at=1:2, labels=c("Versus staying resident\nin non-breeding grounds", "Versus staying resident\nin breeding grounds"), line=1, lwd=0)
axis(side=2, at=c(-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5), labels=c(-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5))
mtext("Difference in niche distance", side=2, line=2.5, cex=1.2,las=0)






################################################################
################################################################

#######      Simulate realistic geographical ranges     ########

################################################################
################################################################


#Great circle distance between each pair of hexagons
pairwise.distance.WH <- rdist.earth(west_Hem, miles=F)
diag(pairwise.distance.WH) <- 0
pairwise.distance.WH.01 <- apply(pairwise.distance.WH, 2, function(x) ifelse(x<250 & x>0, 1, 0))
neig.list_WH <- apply(pairwise.distance.WH.01, 2, function(x) which(x==1))
pairwise.distance.EH <- rdist.earth(east_Hem, miles=F)
diag(pairwise.distance.EH) <- 0
pairwise.distance.EH.01 <- apply(pairwise.distance.EH, 2, function(x) ifelse(x<250 & x>0, 1, 0))
neig.list_EH <- apply(pairwise.distance.EH.01, 2, function(x) which(x==1))


#Compute bearing between every pairs of hexagons
pairwise.bearingsWH = t(apply(west_Hem, 1, function(x) bearingRhumb(x, west_Hem)))
pairwise.bearingsEH = t(apply(east_Hem, 1, function(x) bearingRhumb(x, east_Hem)))

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


pairwise.WH.01.bearings <- pairwise.distance.WH.01 * pairwise.bearingsWH
pairwise.WH.01.bearings[pairwise.WH.01.bearings==0] <- NA
pairwise.EH.01.bearings <- pairwise.distance.EH.01 * pairwise.bearingsEH
pairwise.EH.01.bearings[pairwise.EH.01.bearings==0] <- NA

#Function to simulate n realistic geographical ranges	
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
					if(length(which(neigs.proba > 0)) >= ceiling(length(neigs)/2)){
						block.simulated <- c(block.simulated, sample(neigs, ceiling(length(neigs)/2), prob=neigs.proba))
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
	
#Simulate 100 realistic geographical ranges for each observed seasonal range
simulated.ranges_BR_NH_WH <- apply(PresAbs_BR_NH_WH2, 2, function(x) rangeSimulations(x, west_Hem, neig.list_WH, pairwise.distance.WH, pairwise.WH.01.bearings, 100))
simulated.ranges_NB_NH_WH <- apply(PresAbs_NB_NH_WH2, 2, function(x) rangeSimulations(x, west_Hem, neig.list_WH, pairwise.distance.WH, pairwise.WH.01.bearings, 100))
simulated.ranges_BR_SH_WH <- apply(PresAbs_BR_SH_WH2, 2, function(x) rangeSimulations(x, west_Hem, neig.list_WH, pairwise.distance.WH, pairwise.WH.01.bearings, 100))
simulated.ranges_NB_SH_WH <- apply(PresAbs_NB_SH_WH2, 2, function(x) rangeSimulations(x, west_Hem, neig.list_WH, pairwise.distance.WH, pairwise.WH.01.bearings, 100))
simulated.ranges_BR_NH_EH <- apply(PresAbs_BR_NH_EH2, 2, function(x) rangeSimulations(x, east_Hem, neig.list_EH, pairwise.distance.EH, pairwise.EH.01.bearings, 100))
simulated.ranges_NB_NH_EH <- apply(PresAbs_NB_NH_EH2, 2, function(x) rangeSimulations(x, east_Hem, neig.list_EH, pairwise.distance.EH, pairwise.EH.01.bearings, 100))
simulated.ranges_BR_SH_EH <- apply(PresAbs_BR_SH_EH2, 2, function(x) rangeSimulations(x, east_Hem, neig.list_EH, pairwise.distance.EH, pairwise.EH.01.bearings, 100))
simulated.ranges_NB_SH_EH <- apply(PresAbs_NB_SH_EH2, 2, function(x) rangeSimulations(x, east_Hem, neig.list_EH, pairwise.distance.EH, pairwise.EH.01.bearings, 100))
#load("simulated_ranges.RData")





################################################################
################################################################

######      Compare to alternative migration options    ########

################################################################
################################################################



#Compute niche distances (using the earth mover's distance) for alternative migration options
niche.distances.simulatedBR_NH_WH <- list()
niche.distances.simulatedNB_NH_WH <- list()
for(i in 1:length(nonbreeding.nichesNH_WH)){
	breeding.niches.simulated_NH_WH <- lapply(simulated.ranges_BR_NH_WH[[i]], function(x) nicheDensityRaster(cbind(TempNS_WH[x], PrecNS_WH[x])))
	niche.distances.simulatedBR_NH_WH[[i]] <- mapply(function(X){ emd(X,nonbreeding.nichesNH_WH[[i]]) }, X=breeding.niches.simulated_NH_WH)
	nonbreeding.niches.simulated_NH_WH <- lapply(simulated.ranges_NB_NH_WH[[i]], function(x) nicheDensityRaster(cbind(TempNW_WH[x], PrecNW_WH[x])))
	niche.distances.simulatedNB_NH_WH[[i]] <- mapply(function(Y){ emd(breeding.nichesNH_WH[[i]],Y) }, Y=nonbreeding.niches.simulated_NH_WH)
}
niche.distances.simulatedBR_SH_WH <- list()
niche.distances.simulatedNB_SH_WH <- list()
for(i in 1:length(nonbreeding.nichesSH_WH)){
	breeding.niches.simulated_SH_WH <- lapply(simulated.ranges_BR_SH_WH[[i]], function(x) nicheDensityRaster(cbind(TempNW_WH[x], PrecNW_WH[x])))
	niche.distances.simulatedBR_SH_WH[[i]] <- mapply(function(X){ emd(X,nonbreeding.nichesSH_WH[[i]]) }, X=breeding.niches.simulated_SH_WH)
	nonbreeding.niches.simulated_SH_WH <- lapply(simulated.ranges_NB_SH_WH[[i]], function(x) nicheDensityRaster(cbind(TempNS_WH[x], PrecNS_WH[x])))
	niche.distances.simulatedNB_SH_WH[[i]] <- mapply(function(Y){ emd(breeding.nichesSH_WH[[i]],Y) }, Y=nonbreeding.niches.simulated_SH_WH)
}
niche.distances.simulatedBR_NH_EH <- list()
niche.distances.simulatedNB_NH_EH <- list()
for(i in 1:length(nonbreeding.nichesNH_EH)){
	breeding.niches.simulated_NH_EH <- lapply(simulated.ranges_BR_NH_EH[[i]], function(x) nicheDensityRaster(cbind(TempNS_EH[x], PrecNS_EH[x])))
	niche.distances.simulatedBR_NH_EH[[i]] <- mapply(function(X){ emd(X,nonbreeding.nichesNH_EH[[i]]) }, X=breeding.niches.simulated_NH_EH)
	nonbreeding.niches.simulated_NH_EH <- lapply(simulated.ranges_NB_NH_EH[[i]], function(x) nicheDensityRaster(cbind(TempNW_EH[x], PrecNW_EH[x])))
	niche.distances.simulatedNB_NH_EH[[i]] <- mapply(function(Y){ emd(breeding.nichesNH_EH[[i]],Y) }, Y=nonbreeding.niches.simulated_NH_EH)
}
niche.distances.simulatedBR_SH_EH <- list()
niche.distances.simulatedNB_SH_EH <- list()
for(i in 1:length(nonbreeding.nichesSH_EH)){
	breeding.niches.simulated_SH_EH <- lapply(simulated.ranges_BR_SH_EH[[i]], function(x) nicheDensityRaster(cbind(TempNW_EH[x], PrecNW_EH[x])))
	niche.distances.simulatedBR_SH_EH[[i]] <- mapply(function(X){ emd(X,nonbreeding.nichesSH_EH[[i]]) }, X=breeding.niches.simulated_SH_EH)
	nonbreeding.niches.simulated_SH_EH <- lapply(simulated.ranges_NB_SH_EH[[i]], function(x) nicheDensityRaster(cbind(TempNS_EH[x], PrecNS_EH[x])))
	niche.distances.simulatedNB_SH_EH[[i]] <- mapply(function(Y){ emd(breeding.nichesSH_EH[[i]],Y) }, Y=nonbreeding.niches.simulated_SH_EH)
}
#load("simulated_niche_distances.RData")
niche.distances.simulatedNH_WH <- mapply(function(X,Y){ c(X,Y) }, X=niche.distances.simulatedBR_NH_WH, Y=niche.distances.simulatedNB_NH_WH)
niche.distances.simulatedNH_WH <- split(niche.distances.simulatedNH_WH, rep(1:ncol(niche.distances.simulatedNH_WH),each=nrow(niche.distances.simulatedNH_WH)))
niche.distances.simulatedSH_WH <- mapply(function(X,Y){ c(X,Y) }, X=niche.distances.simulatedBR_SH_WH, Y=niche.distances.simulatedNB_SH_WH)
niche.distances.simulatedSH_WH <- split(niche.distances.simulatedSH_WH, rep(1:ncol(niche.distances.simulatedSH_WH),each=nrow(niche.distances.simulatedSH_WH)))
niche.distances.simulatedNH_EH <- mapply(function(X,Y){ c(X,Y) }, X=niche.distances.simulatedBR_NH_EH, Y=niche.distances.simulatedNB_NH_EH)
niche.distances.simulatedNH_EH <- split(niche.distances.simulatedNH_EH, rep(1:ncol(niche.distances.simulatedNH_EH),each=nrow(niche.distances.simulatedNH_EH)))
niche.distances.simulatedSH_EH <- mapply(function(X,Y){ c(X,Y) }, X=niche.distances.simulatedBR_SH_EH, Y=niche.distances.simulatedNB_SH_EH)
niche.distances.simulatedSH_EH <- split(niche.distances.simulatedSH_EH, rep(1:ncol(niche.distances.simulatedSH_EH),each=nrow(niche.distances.simulatedSH_EH)))
niche.distances.simulated <- c(niche.distances.simulatedNH_WH, niche.distances.simulatedNH_EH, niche.distances.simulatedSH_WH, niche.distances.simulatedSH_EH)


#Compute geographical (migration) distance for alternative migration options
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


#Compute resource scarcity for alternative migration options
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



#Compute ranks
ranks_niche <- vector()
for(i in 1:length(niche.distances.simulated)){
	ranks_niche[i] <- length(which(niche.distances.simulated[[i]] < niche.distances.obs[i])) / length(niche.distances.simulated[[i]])	
}
ranks_geo <- vector()
for(i in 1:length(geo.distances.simulated)){
	ranks_geo[i] <- length(which(geo.distances.simulated[[i]] < geo.distances.obs[i])) / length(geo.distances.simulated[[i]])	
}
ranks_resource <- vector()
for(i in 1:length(resource.scarcity.simulated)){
	ranks_resource[i] <- length(which(resource.scarcity.simulated[[i]] < resource.scarcity.obs[i])) / length(resource.scarcity.simulated[[i]])	
}
ranks_geo_niche <- vector()
for(i in 1:length(geo.distances.simulated)){
	geo_niche <- scale(c(geo.distances.obs[i], geo.distances.simulated[[i]])) + scale(c(niche.distances.obs[i], niche.distances.simulated[[i]]))
	ranks_geo_niche[i] <- length(which(geo_niche[-1] < geo_niche[1])) / length(geo_niche[-1])	
}
ranks_geo_resource <- vector()
for(i in 1:length(geo.distances.simulated)){
	geo_resource <- scale(c(geo.distances.obs[i], geo.distances.simulated[[i]])) + scale(c(resource.scarcity.obs[i], resource.scarcity.simulated[[i]]))
	ranks_geo_resource[i] <- length(which(geo_resource[-1] < geo_resource[1])) / length(geo_resource[-1])	
}
ranks_niche_resource <- vector()
for(i in 1:length(niche.distances.simulated)){
	niche_resource <- scale(c(niche.distances.obs[i], niche.distances.simulated[[i]])) + scale(c(resource.scarcity.obs[i], resource.scarcity.simulated[[i]]))
	ranks_niche_resource[i] <- length(which(niche_resource[-1] < niche_resource[1])) / length(niche_resource[-1])	
}
ranks_geo_niche_resource <- vector()
for(i in 1:length(geo.distances.simulated)){
	geo_niche_resource <- scale(c(geo.distances.obs[i], geo.distances.simulated[[i]])) + scale(c(niche.distances.obs[i], niche.distances.simulated[[i]])) + scale(c(resource.scarcity.obs[i], resource.scarcity.simulated[[i]]))
	ranks_geo_niche_resource[i] <- length(which(geo_niche_resource[-1] < geo_niche_resource[1])) / length(geo_niche_resource[-1])	
}




### Figure 2 - illustrating how we compute the ranks using Schrenck's Bittern (Ixobrychus.eurhythmus) ###

#Load shapefile of the global hexagon grid
hexgrid <- readOGR("Hex_grid", "isea3h7_analyses_clean", verbose=FALSE)
hexgrid <- hexgrid[,c(1,2,15,16)]
hexgridWH <- hexgrid[which(hexgrid@data[,4] <= -30),]
hexgridEH <- hexgrid[which(hexgrid@data[,4] > -30),]

#Hexagons ID in eastern hemisphere
hexidEH <- PresAbs_BR_NH[match(envData2[which(envData2$LONGITUDE>-30),1], PresAbs_BR_NH[,1]),1]

#Plot
par(mfrow=c(2,2), mar=c(2.5,2.5,0.2,1.5), mgp=c(1.5,0.5,0))
plot(hexgridEH, col="grey", border = "grey")
plot(hexgridEH[match(PresAbs_BR_NH[which(PresAbs_BR_NH[,which(colnames(PresAbs_NB_NH) == "Ixobrychus.eurhythmus")] == 1),1], hexgridEH@data[,1]),], col="red3", border = "red3", add=T)
plot(hexgridEH[match(hexidEH[simulated.ranges_BR_NH_EH[[181]][[12]]], hexgridEH@data[,1]),], col="black", border = "black", add=T)
plot(hexgridEH[match(hexidEH[simulated.ranges_BR_NH_EH[[181]][[20]]], hexgridEH@data[,1]),], col="black", border = "black", add=T)
mtext("A", cex=1.8, side=3, line=-1.5, at=-50)
plot(scale(c(geo.distances.obs[181], geo.distances.simulated[[181]])), scale(c(niche.distances.obs[181], niche.distances.simulated[[181]])), pch=20, axes=F, xlab="Geographic distance", ylab = "Niche distance", cex.lab=1.2, add=F)
axis(1,at=c(-2,-1,0,1,2), labels=c(-2,-1,0,1,2))
axis(2,at=c(-2,-1,0,1,2,3), labels=c(-2,-1,0,1,2,3))
abline(a = scale(c(niche.distances.obs[181], niche.distances.simulated[[181]]))[1], b = 0, col="orange")
abline(a = (scale(c(geo.distances.obs[181], geo.distances.simulated[[181]])) + scale(c(niche.distances.obs[181], niche.distances.simulated[[181]])))[1], b = -1, col="blue")
points(scale(c(geo.distances.obs[181], geo.distances.simulated[[181]]))[1], scale(c(niche.distances.obs[181], niche.distances.simulated[[181]]))[1], pch=20, col="red3", cex=2.5)
mtext("C", cex=1.8, side=3, line=-1.5, at=-2.35)
plot(hexgridEH, col="grey", border = "grey")
plot(hexgridEH[match(PresAbs_NB_NH[which(PresAbs_NB_NH[,which(colnames(PresAbs_NB_NH) == "Ixobrychus.eurhythmus")] == 1),1], hexgridEH@data[,1]),], col="red3", border = "red3", add=T)
plot(hexgridEH[match(hexidEH[simulated.ranges_NB_NH_EH[[181]][[1]]], hexgridEH@data[,1]),], col="black", border = "black", add=T)
plot(hexgridEH[match(hexidEH[simulated.ranges_NB_NH_EH[[181]][[7]]], hexgridEH@data[,1]),], col="black", border = "black", add=T)
mtext("B", cex=1.8, side=3, line=-1.5, at=-50)
geo_niche_resource_schrenck <- scale(c(geo.distances.obs[181], geo.distances.simulated[[181]])) + scale(c(niche.distances.obs[181], niche.distances.simulated[[181]])) + scale(c(resource.scarcity.obs[181], resource.scarcity.simulated[[181]]))
plot(geo_niche_resource_schrenck, rnorm(length(geo_niche_resource_schrenck), 1, 0.1), pch=20, axes=F, xlab="Niche distance + geographic distance + resource scarcity", ylab = "", ylim=c(0,2), cex.lab=1.2)
axis(1,at=c(-3,-2,-1,0,1,2,3), labels=c(-3,-2,-1,0,1,2,3))
abline(v=geo_niche_resource_schrenck[1], col="green")
points(geo_niche_resource_schrenck[1], 1, pch=20, ylim=c(0,2), col="red3", cex=2.5)
mtext("D", cex=1.8, side=3, line=-1.5, at=-3.5)



### Figure 5 - histograms of ranks ###

par(mfrow=c(3,3), mar=c(4.5,3.5,1.5,0.5), mgp=c(2,0.5,0))
hist(ranks_niche, main="", ylab="Number of species", xlab="Scaled rank for niche distance", xlim=c(0,1), cex.lab=1.3, col="grey", border="white")
mtext("A", cex=1.3, side=3, line=-0.25, at=-0.2)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-2.5, at=0.7)
hist(ranks_geo, main="", ylab="", xlab="Scaled rank for geographic distance", xlim=c(0,1), cex.lab=1.3, col="grey", border="white")
mtext("B", cex=1.3, side=3, line=-0.25, at=-0.2)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-2.5, at=0.7)
hist(ranks_resource, main="", ylab="", xlab="Scaled rank for resource scarcity", xlim=c(0,1), cex.lab=1.3, col="grey", border="white")
mtext("C", cex=1.3, side=3, line=-0.25, at=-0.2)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-2.5, at=0.7)
par(new=F, mar=c(4.5,3.5,1.5,0.5), mgp=c(3,0.5,0))
hist(ranks_geo_niche, main="", ylab="", xlab="Scaled rank for niche distance \n+ geographic distance", xlim=c(0,1), cex.lab=1.3, col="grey", border="white")
mtext("Number of species", cex=0.9, side=2, line=2, at=135)
mtext("D", cex=1.3, side=3, line=-0.25, at=-0.2)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-2.5, at=0.7)
hist(ranks_niche_resource, main="", ylab="", xlab="Scaled rank for niche distance \n+ resource scarcity", xlim=c(0,1), cex.lab=1.3, col="grey", border="white")
mtext("E", cex=1.3, side=3, line=-0.25, at=-0.2)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-2.5, at=0.7)
hist(ranks_geo_resource, main="", ylab="", xlab="Scaled rank for geographic distance \n+ resource scarcity", xlim=c(0,1), cex.lab=1.3, col="grey", border="white")
mtext("F", cex=1.3, side=3, line=-0.25, at=-0.2)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-2.5, at=0.7)
hist(ranks_geo_niche_resource, main="", ylab="", xlab="Scaled rank for niche distance \n+ geographic distance + resource scarcity", xlim=c(0,1), cex.lab=1.3, col="grey", border="white")
mtext("G", cex=1.3, side=3, line=-0.25, at=-0.2)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-2.5, at=0.7)
mtext("Number of species", cex=0.9, side=2, line=2, at=155)

#Kolmogorov-Smirnov tests of skewdness of ranks distributions
ks.test(ranks_niche, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_geo, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_resource, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_geo_niche, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_geo_resource, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_niche_resource, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_geo_niche_resources, runif(100000,0,1), alternative="greater")$p.value





### Figure 6 - Weighted Richness Maps ###



#### MAPS

PresAbs_BR_NH_WH_compet <- PresAbs_BR_NH_WH[,match(names(selectedSpeciesNH_WH), colnames(PresAbs_BR_NH_WH))]
PresAbs_NB_NH_WH_compet <- PresAbs_NB_NH_WH[,match(names(selectedSpeciesNH_WH), colnames(PresAbs_NB_NH_WH))]
PresAbs_BR_SH_WH_compet <- PresAbs_BR_SH_WH[,match(names(selectedSpeciesSH_WH), colnames(PresAbs_BR_SH_WH))]
PresAbs_NB_SH_WH_compet <- PresAbs_NB_SH_WH[,match(names(selectedSpeciesSH_WH), colnames(PresAbs_NB_SH_WH))]
PresAbs_BR_WH_compet <- cbind(PresAbs_BR_NH_WH_compet, PresAbs_BR_SH_WH_compet)
PresAbs_NB_WH_compet <- cbind(PresAbs_NB_NH_WH_compet, PresAbs_NB_SH_WH_compet)

PresAbs_BR_NH_EH_compet <- PresAbs_BR_NH_EH[,match(names(selectedSpeciesNH_EH), colnames(PresAbs_BR_NH_EH))[-which(is.na(match(names(selectedSpeciesNH_EH), colnames(PresAbs_BR_NH_EH)))==T)]]
PresAbs_NB_NH_EH_compet <- PresAbs_NB_NH_EH[,match(names(selectedSpeciesNH_EH), colnames(PresAbs_NB_NH_EH))[-which(is.na(match(names(selectedSpeciesNH_EH), colnames(PresAbs_NB_NH_EH)))==T)]]
for(j in 380:381){
PresAbs_BR_NH_EH_compet[,j] <- PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]
PresAbs_NB_NH_EH_compet[,j] <- PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]
}
PresAbs_BR_SH_EH_compet <- PresAbs_BR_SH_EH[,match(names(selectedSpeciesSH_EH), colnames(PresAbs_BR_SH_EH))[-which(is.na(match(names(selectedSpeciesSH_EH), colnames(PresAbs_BR_SH_EH)))==T)]]
PresAbs_NB_SH_EH_compet <- PresAbs_NB_SH_EH[,match(names(selectedSpeciesSH_EH), colnames(PresAbs_NB_SH_EH))[-which(is.na(match(names(selectedSpeciesSH_EH), colnames(PresAbs_BR_SH_EH)))==T)]]
PresAbs_BR_SH_EH_compet[,13] <- PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[13]))]
PresAbs_NB_SH_EH_compet[,13] <- PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[13]))]
PresAbs_BR_EH_compet <- cbind(PresAbs_BR_NH_EH_compet, PresAbs_BR_SH_EH_compet)
PresAbs_NB_EH_compet <- cbind(PresAbs_NB_NH_EH_compet, PresAbs_NB_SH_EH_compet)


#####   Figure 6 – Weighted Richness Maps

par(mfrow=c(4,2), mar=c(0.1,0.1,0.1,0.1), mgp=c(1.5,0.5,0))

## Weighted richness maps for Niche

weighted.richness.br_WH <- apply(PresAbs_BR_WH_compet, 1, function(x) sum(((c(rank_nicheDistObs_WH, rank_nicheDistObs_SHWH)-1)/200)*x)) / apply(PresAbs_BR_WH_compet, 1, sum)
weighted.richness.nb_WH <- apply(PresAbs_NB_WH_compet, 1, function(x) sum(((c(rank_nicheDistObs_WH, rank_nicheDistObs_SHWH)-1)/200)*x)) / apply(PresAbs_NB_WH_compet, 1, sum)
weighted.richness.br_WH[which(weighted.richness.br_WH == "NaN")] <- 0
weighted.richness.nb_WH[which(weighted.richness.nb_WH == "NaN")] <- 0
weighted.richness.br_EH <- apply(PresAbs_BR_EH_compet, 1, function(x) sum(((c(rank_nicheDistObs_EH, rank_nicheDistObs_SHEH)-1)/200)*x)) / apply(PresAbs_BR_EH_compet, 1, sum)
weighted.richness.nb_EH <- apply(PresAbs_NB_EH_compet, 1, function(x) sum(((c(rank_nicheDistObs_EH, rank_nicheDistObs_SHEH)-1)/200)*x)) / apply(PresAbs_NB_EH_compet, 1, sum)
weighted.richness.br_EH[which(weighted.richness.br_EH == "NaN")] <- 0
weighted.richness.nb_EH[which(weighted.richness.nb_EH == "NaN")] <- 0
weighted.richness.br <- c(weighted.richness.br_WH, weighted.richness.br_EH)
weighted.richness.nb <- c(weighted.richness.nb_WH, weighted.richness.nb_EH)

rbPal <- colorRampPalette(c("yellow", "red3"))

datcol <- rbPal(5)[as.numeric(cut(weighted.richness.br, breaks=c(-0.1,0.1,0.2,0.3,0.4,1.1)))]
plot(hexgrid[match(c(PresAbs_BR_NH_WH[,1], PresAbs_BR_NH_EH[,1]), hexgrid@data[,1]),], col= datcol, border = datcol, bg="grey")
mtext("A", cex=1.3, side=3, line=-2, at=-180)
legend("bottomleft", inset=.04, bg="grey", box.col="grey", title="Average rank\nof migrants", c("> 0.4","0.3–0.4", "0.2–0.3", "0.1–0.2", "0–0.1"), fill=rev(rbPal(5)), cex=0.8)

datcol <- rbPal(5)[as.numeric(cut(weighted.richness.nb, breaks=c(-0.1,0.1,0.2,0.3,0.4,1.1)))]
plot(hexgrid[match(c(PresAbs_NB_NH_WH[,1], PresAbs_NB_NH_EH[,1]), hexgrid@data[,1]),], col= datcol, border = datcol, bg="grey")
mtext("B", cex=1.3, side=3, line=-2, at=-180)




## Weighted richness maps for Resources + GeoDist

weighted.richness.br_WH <- apply(PresAbs_BR_WH_compet, 1, function(x) sum(((c(rank_resources_geoDist_Obs_WH, rank_resources_geoDist_Obs_SHWH)-1)/200)*x)) / apply(PresAbs_BR_WH_compet, 1, sum)
weighted.richness.nb_WH <- apply(PresAbs_NB_WH_compet, 1, function(x) sum(((c(rank_resources_geoDist_Obs_WH, rank_resources_geoDist_Obs_SHWH)-1)/200)*x)) / apply(PresAbs_NB_WH_compet, 1, sum)
weighted.richness.br_WH[which(weighted.richness.br_WH == "NaN")] <- 0
weighted.richness.nb_WH[which(weighted.richness.nb_WH == "NaN")] <- 0
weighted.richness.br_EH <- apply(PresAbs_BR_EH_compet, 1, function(x) sum(((c(rank_resources_geoDist_Obs_EH, rank_resources_geoDist_Obs_SHEH)-1)/200)*x)) / apply(PresAbs_BR_EH_compet, 1, sum)
weighted.richness.nb_EH <- apply(PresAbs_NB_EH_compet, 1, function(x) sum(((c(rank_resources_geoDist_Obs_EH, rank_resources_geoDist_Obs_SHEH)-1)/200)*x)) / apply(PresAbs_NB_EH_compet, 1, sum)
weighted.richness.br_EH[which(weighted.richness.br_EH == "NaN")] <- 0
weighted.richness.nb_EH[which(weighted.richness.nb_EH == "NaN")] <- 0
weighted.richness.br <- c(weighted.richness.br_WH, weighted.richness.br_EH)
weighted.richness.nb <- c(weighted.richness.nb_WH, weighted.richness.nb_EH)

datcol <- rbPal(5)[as.numeric(cut(weighted.richness.br, breaks=c(-0.1,0.1,0.2,0.3,0.4,1.1)))]
plot(hexgrid[match(c(PresAbs_BR_NH_WH[,1], PresAbs_BR_NH_EH[,1]), hexgrid@data[,1]),], col= datcol, border = datcol, bg="grey")
mtext("C", cex=1.3, side=3, line=-2, at=-180)

datcol <- rbPal(5)[as.numeric(cut(weighted.richness.nb, breaks=c(-0.1,0.1,0.2,0.3,0.4,1.1)))]
plot(hexgrid[match(c(PresAbs_NB_NH_WH[,1], PresAbs_NB_NH_EH[,1]), hexgrid@data[,1]),], col= datcol, border = datcol, bg="grey")
mtext("D", cex=1.3, side=3, line=-2, at=-180)




## Weighted richness maps for Niche + Resources + GeoDist

weighted.richness.br_WH <- apply(PresAbs_BR_WH_compet, 1, function(x) sum(((c(rank_nicheDist_resources_geoDist_Obs_WH, rank_nicheDist_resources_geoDist_Obs_SHWH)-1)/200)*x)) / apply(PresAbs_BR_WH_compet, 1, sum)
weighted.richness.nb_WH <- apply(PresAbs_NB_WH_compet, 1, function(x) sum(((c(rank_nicheDist_resources_geoDist_Obs_WH, rank_nicheDist_resources_geoDist_Obs_SHWH)-1)/200)*x)) / apply(PresAbs_NB_WH_compet, 1, sum)
weighted.richness.br_WH[which(weighted.richness.br_WH == "NaN")] <- 0
weighted.richness.nb_WH[which(weighted.richness.nb_WH == "NaN")] <- 0
weighted.richness.br_EH <- apply(PresAbs_BR_EH_compet, 1, function(x) sum(((c(rank_nicheDist_resources_geoDist_Obs_EH, rank_nicheDist_resources_geoDist_Obs_SHEH)-1)/200)*x)) / apply(PresAbs_BR_EH_compet, 1, sum)
weighted.richness.nb_EH <- apply(PresAbs_NB_EH_compet, 1, function(x) sum(((c(rank_nicheDist_resources_geoDist_Obs_EH, rank_nicheDist_resources_geoDist_Obs_SHEH)-1)/200)*x)) / apply(PresAbs_NB_EH_compet, 1, sum)
weighted.richness.br_EH[which(weighted.richness.br_EH == "NaN")] <- 0
weighted.richness.nb_EH[which(weighted.richness.nb_EH == "NaN")] <- 0
weighted.richness.br <- c(weighted.richness.br_WH, weighted.richness.br_EH)
weighted.richness.nb <- c(weighted.richness.nb_WH, weighted.richness.nb_EH)

datcol <- rbPal(5)[as.numeric(cut(weighted.richness.br, breaks=c(-0.1,0.1,0.2,0.3,0.4,1.1)))]
plot(hexgrid[match(c(PresAbs_BR_NH_WH[,1], PresAbs_BR_NH_EH[,1]), hexgrid@data[,1]),], col= datcol, border = datcol, bg="grey")
mtext("E", cex=1.3, side=3, line=-2, at=-180)

datcol <- rbPal(5)[as.numeric(cut(weighted.richness.nb, breaks=c(-0.1,0.1,0.2,0.3,0.4,1.1)))]
plot(hexgrid[match(c(PresAbs_NB_NH_WH[,1], PresAbs_NB_NH_EH[,1]), hexgrid@data[,1]),], col= datcol, border = datcol, bg="grey")
mtext("F", cex=1.3, side=3, line=-2, at=-180)




## richness in migrants

mbr_NHWH <- apply(PresAbs_BR_NH_WH_compet, 1, sum)
mnb_NHWH <- apply(PresAbs_NB_NH_WH_compet, 1, sum)
mbr_SHWH <- apply(PresAbs_BR_SH_WH_compet, 1, sum)
mnb_SHWH <- apply(PresAbs_NB_SH_WH_compet, 1, sum)
mbr_WH <- mbr_NHWH + mbr_SHWH
mnb_WH <- mnb_NHWH + mnb_SHWH
mbr_NHEH <- apply(PresAbs_BR_NH_EH_compet, 1, sum)
mnb_NHEH <- apply(PresAbs_NB_NH_EH_compet, 1, sum)
mbr_SHEH <- apply(PresAbs_BR_SH_EH_compet, 1, sum)
mnb_SHEH <- apply(PresAbs_NB_SH_EH_compet, 1, sum)
mbr_EH <- mbr_NHEH + mbr_SHEH
mnb_EH <- mnb_NHEH + mnb_SHEH

mbr <- c(mbr_WH, mbr_EH)
mnb <- c(mnb_WH, mnb_EH)


rbPal <- colorRampPalette(c("yellow3", "dark green"))

datcol <- rbPal(5)[as.numeric(cut(mbr, breaks=c(-0.1,25,50,75,100,150)))]
plot(hexgrid[match(c(PresAbs_BR_NH_WH[,1], PresAbs_BR_NH_EH[,1]), hexgrid@data[,1]),], col= datcol, border = datcol, bg="grey")
mtext("G", cex=1.3, side=3, line=-2, at=-180)
legend("bottomleft", inset=.04, bg="grey", box.col="grey", title="Richness\nin migrants", c("> 100","75–100", "50–75", "25–50", "0–25"), fill=rev(rbPal(5)), cex=0.8)

datcol <- rbPal(5)[as.numeric(cut(mnb, breaks=c(-0.1,25,50,75,100,150)))]
plot(hexgrid[match(c(PresAbs_NB_NH_WH[,1], PresAbs_NB_NH_EH[,1]), hexgrid@data[,1]),], col= datcol, border = datcol, bg="grey")
mtext("H", cex=1.3, side=3, line=-2, at=-180)






###  Figure 7 – Weighted arrows maps 


centroidsWH <- matrix(nrow=length(selectedSpeciesNH_WH), ncol=4)
centroidsEH <- matrix(nrow=length(selectedSpeciesNH_EH), ncol=4)
centroidsSHWH <- matrix(nrow=length(selectedSpeciesSH_WH), ncol=4)
centroidsSHEH <- matrix(nrow=length(selectedSpeciesSH_EH), ncol=4)
for(j in 1:length(selectedSpeciesNH_WH)){ 
centroidsWH[j,1] <- mean(west_Hem[which(PresAbs_BR_NH_WH[,match(names(selectedSpeciesNH_WH)[j], colnames(PresAbs_BR_NH_WH))]==1),1])
centroidsWH[j,2] <- mean(west_Hem[which(PresAbs_BR_NH_WH[,match(names(selectedSpeciesNH_WH)[j], colnames(PresAbs_BR_NH_WH))]==1),2])
centroidsWH[j,3] <- mean(west_Hem[which(PresAbs_NB_NH_WH[,match(names(selectedSpeciesNH_WH)[j], colnames(PresAbs_NB_NH_WH))]==1),1])
centroidsWH[j,4] <- mean(west_Hem[which(PresAbs_NB_NH_WH[,match(names(selectedSpeciesNH_WH)[j], colnames(PresAbs_NB_NH_WH))]==1),2])
}
for(j in 1:length(selectedSpeciesSH_WH)){ 
centroidsSHWH[j,1] <- mean(west_Hem[which(PresAbs_BR_SH_WH[,match(names(selectedSpeciesSH_WH)[j], colnames(PresAbs_BR_SH_WH))]==1),1])
centroidsSHWH[j,2] <- mean(west_Hem[which(PresAbs_BR_SH_WH[,match(names(selectedSpeciesSH_WH)[j], colnames(PresAbs_BR_SH_WH))]==1),2])
centroidsSHWH[j,3] <- mean(west_Hem[which(PresAbs_NB_SH_WH[,match(names(selectedSpeciesSH_WH)[j], colnames(PresAbs_NB_SH_WH))]==1),1])
centroidsSHWH[j,4] <- mean(west_Hem[which(PresAbs_NB_SH_WH[,match(names(selectedSpeciesSH_WH)[j], colnames(PresAbs_NB_SH_WH))]==1),2])
}
for(j in 1:length(selectedSpeciesNH_EH)){ 
centroidsEH[j,1] <- mean(east_Hem[which(PresAbs_BR_NH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_BR_NH_EH))]==1),1])
centroidsEH[j,2] <- mean(east_Hem[which(PresAbs_BR_NH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_BR_NH_EH))]==1),2])
centroidsEH[j,3] <- mean(east_Hem[which(PresAbs_NB_NH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_NB_NH_EH))]==1),1])
centroidsEH[j,4] <- mean(east_Hem[which(PresAbs_NB_NH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_NB_NH_EH))]==1),2])
}
for(j in 380:381){ 
centroidsEH[j,1] <- mean(east_Hem[which(PresAbs_BR_SH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_BR_SH_EH))]==1),1])
centroidsEH[j,2] <- mean(east_Hem[which(PresAbs_BR_SH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_BR_SH_EH))]==1),2])
centroidsEH[j,3] <- mean(east_Hem[which(PresAbs_NB_SH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_NB_SH_EH))]==1),1])
centroidsEH[j,4] <- mean(east_Hem[which(PresAbs_NB_SH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_NB_SH_EH))]==1),2])
}
for(j in 1:length(selectedSpeciesSH_EH)){ 
centroidsSHEH[j,1] <- mean(east_Hem[which(PresAbs_BR_SH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_BR_SH_EH))]==1),1])
centroidsSHEH[j,2] <- mean(east_Hem[which(PresAbs_BR_SH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_BR_SH_EH))]==1),2])
centroidsSHEH[j,3] <- mean(east_Hem[which(PresAbs_NB_SH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_NB_SH_EH))]==1),1])
centroidsSHEH[j,4] <- mean(east_Hem[which(PresAbs_NB_SH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_NB_SH_EH))]==1),2])
}
j=13
centroidsSHEH[j,1] <- mean(east_Hem[which(PresAbs_BR_NH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_BR_NH_EH))]==1),1])
centroidsSHEH[j,2] <- mean(east_Hem[which(PresAbs_BR_NH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_BR_NH_EH))]==1),2])
centroidsSHEH[j,3] <- mean(east_Hem[which(PresAbs_NB_NH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_NB_NH_EH))]==1),1])
centroidsSHEH[j,4] <- mean(east_Hem[which(PresAbs_NB_NH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_NB_NH_EH))]==1),2])

centroids <- rbind(centroidsWH, centroidsSHWH, centroidsEH, centroidsSHEH)

weights.NHWH <- (rank_nicheDistObs_WH-1)/200
weights.SHWH <- (rank_nicheDistObs_SHWH-1)/200
weights.NHEH <- (rank_nicheDistObs_EH-1)/200
weights.SHEH <- (rank_nicheDistObs_SHEH-1)/200
weights.niche <- c(weights.NHWH, weights.SHWH, weights.NHEH, weights.SHEH)
centroids.niche <- centroids[order(weights.niche),]
weights.niche <- weights.niche[order(weights.niche)]
weights.niche1 <- weights.niche[which(weights.niche < 0.1)]
centroids.niche1 <- centroids.niche[which(weights.niche < 0.1),]
weights.niche2 <- weights.niche[which(weights.niche >= 0.1 & weights.niche < 0.5)]
centroids.niche2 <- centroids.niche[which(weights.niche >= 0.1 & weights.niche < 0.5),]
weights.niche3 <- weights.niche[which(weights.niche >= 0.5)]
centroids.niche3 <- centroids.niche[which(weights.niche >= 0.5),]

weights.NHWH <- (rank_resources_geoDist_Obs_WH-1)/200
weights.SHWH <- (rank_resources_geoDist_Obs_SHWH-1)/200
weights.NHEH <- (rank_resources_geoDist_Obs_EH-1)/200
weights.SHEH <- (rank_resources_geoDist_Obs_SHEH-1)/200
weights.ResGeo <- c(weights.NHWH, weights.SHWH, weights.NHEH, weights.SHEH)
centroids.ResGeo <- centroids[order(weights.ResGeo),]
weights.ResGeo <- weights.ResGeo[order(weights.ResGeo)]
weights.ResGeo1 <- weights.ResGeo[which(weights.ResGeo < 0.1)]
centroids.ResGeo1 <- centroids.ResGeo[which(weights.ResGeo < 0.1),]
weights.ResGeo2 <- weights.ResGeo[which(weights.ResGeo >= 0.1 & weights.ResGeo < 0.5)]
centroids.ResGeo2 <- centroids.ResGeo[which(weights.ResGeo >= 0.1 & weights.ResGeo < 0.5),]
weights.ResGeo3 <- weights.ResGeo[which(weights.ResGeo >= 0.5)]
centroids.ResGeo3 <- centroids.ResGeo[which(weights.ResGeo >= 0.5),]

weights.NHWH <- (rank_nicheDist_resources_geoDist_Obs_WH-1)/200
weights.SHWH <- (rank_nicheDist_resources_geoDist_Obs_SHWH-1)/200
weights.NHEH <- (rank_nicheDist_resources_geoDist_Obs_EH-1)/200
weights.SHEH <- (rank_nicheDist_resources_geoDist_Obs_SHEH-1)/200
weights.nicheResGeo <- c(weights.NHWH, weights.SHWH, weights.NHEH, weights.SHEH)
centroids.nicheResGeo <- centroids[order(weights.nicheResGeo),]
weights.nicheResGeo <- weights.nicheResGeo[order(weights.nicheResGeo)]
weights.nicheResGeo1 <- weights.nicheResGeo[which(weights.nicheResGeo < 0.1)]
centroids.nicheResGeo1 <- centroids.nicheResGeo[which(weights.nicheResGeo < 0.1),]
weights.nicheResGeo2 <- weights.nicheResGeo[which(weights.nicheResGeo >= 0.1 & weights.nicheResGeo < 0.5)]
centroids.nicheResGeo2 <- centroids.nicheResGeo[which(weights.nicheResGeo >= 0.1 & weights.nicheResGeo < 0.5),]
weights.nicheResGeo3 <- weights.nicheResGeo[which(weights.nicheResGeo >= 0.5)]
centroids.nicheResGeo3 <- centroids.nicheResGeo[which(weights.nicheResGeo >= 0.5),]

library(maps)
library(geosphere)

#rbPal <- colorRampPalette(c("yellow", "brown4"))
#rbPal1 <- colorRampPalette(c("yellow", "orange"))
#rbPal2 <- colorRampPalette(c("orange", "brown4"))

par(mfrow=c(3,3), mar=c(0.1,0.1,0.1,0.1), mgp=c(1.5,0.5,0))

plot(hexgrid, col= "grey", border = "grey")
for(i in 1:length(weights.niche1)){
	inter <- gcIntermediate(c(centroids.niche1[i,1],centroids.niche1[i,2]), c(centroids.niche1[i,3],centroids.niche1[i,4]), n=50, addStartEnd=T)
	#colcol <- rbPal1(10)[as.numeric(cut(weights.niche1, breaks=10))]
	#lines(inter, lwd=(weights.niche1[i]*1)+1, col=colcol[i])
	lines(inter, lwd=1, col="yellow")
	points(centroids.niche1[i,1], centroids.niche1[i,2], pch=20, col="red", cex=0.3)
	points(centroids.niche1[i,3], centroids.niche1[i,4], pch=20, col="blue", cex=0.3)
}
mtext("A", cex=1.3, side=3, line=-2.2, at=-180)
#legend("bottomleft", inset=c(0, 0.1), bg="white", box.col="white", title="", c("Rank < 20", "Breeding centroid", "Non-breeding centroid"), lty=c(1,NA,NA), pch=c(NA,20,20), col=c("yellow","red", "blue"))

plot(hexgrid, col= "grey", border = "grey")
for(i in 1:length(weights.niche2)){
	inter <- gcIntermediate(c(centroids.niche2[i,1],centroids.niche2[i,2]), c(centroids.niche2[i,3],centroids.niche2[i,4]), n=50, addStartEnd=T)
	#colcol <- rbPal(10)[as.numeric(cut(weights.niche2, breaks=10))]
	#lines(inter, lwd=(weights.niche2[i]*1)+1, col=colcol[i])
	lines(inter, lwd=1, col="orange")
	points(centroids.niche2[i,1], centroids.niche2[i,2], pch=20, col="red", cex=0.3)
	points(centroids.niche2[i,3], centroids.niche2[i,4], pch=20, col="blue", cex=0.3)
}
mtext("B", cex=1.3, side=3, line=-2.2, at=-180)
#legend("bottomleft", inset=c(0, 0.1), bg="white", box.col="white", title="", c("20 <= Rank < 100 ", "Breeding centroid", "Non-breeding centroid"), lty=c(1,NA,NA), pch=c(NA,20,20), col=c("orange","red", "blue"))

plot(hexgrid, col= "grey", border = "grey")
for(i in 1:length(weights.niche3)){
	inter <- gcIntermediate(c(centroids.niche3[i,1],centroids.niche3[i,2]), c(centroids.niche3[i,3],centroids.niche3[i,4]), n=50, addStartEnd=T)
	#colcol <- rbPal(10)[as.numeric(cut(weights.niche3, breaks=10))]
	#lines(inter, lwd=(weights.niche3[i]*1)+1, col=colcol[i])
	lines(inter, lwd=1, col="brown4")
	points(centroids.niche3[i,1], centroids.niche3[i,2], pch=20, col="red", cex=0.3)
	points(centroids.niche3[i,3], centroids.niche3[i,4], pch=20, col="blue", cex=0.3)
}
mtext("C", cex=1.3, side=3, line=-2.2, at=-180)
#legend("bottomleft", inset=c(0, 0.1), bg="white", box.col="white", title="", c("Rank >= 100", "Breeding centroid", "Non-breeding centroid"), lty=c(1,NA,NA), pch=c(NA,20,20), col=c("brown4","red", "blue"))


plot(hexgrid, col= "grey", border = "grey")
for(i in 1:length(weights.ResGeo1)){
	inter <- gcIntermediate(c(centroids.ResGeo1[i,1],centroids.ResGeo1[i,2]), c(centroids.ResGeo1[i,3],centroids.ResGeo1[i,4]), n=50, addStartEnd=T)
	#colcol <- rbPal1(10)[as.numeric(cut(weights.ResGeo1, breaks=10))]
	#lines(inter, lwd=(weights.ResGeo1[i]*1)+1, col=colcol[i])
	lines(inter, lwd=1, col="yellow")
	points(centroids.ResGeo1[i,1], centroids.ResGeo1[i,2], pch=20, col="red", cex=0.3)
	points(centroids.ResGeo1[i,3], centroids.ResGeo1[i,4], pch=20, col="blue", cex=0.3)
}
mtext("D", cex=1.3, side=3, line=-2.2, at=-180)
#legend("bottomleft", inset=c(0, 0.1), bg="white", box.col="white", title="", c("Rank < 20", "Breeding centroid", "Non-breeding centroid"), lty=c(1,NA,NA), pch=c(NA,20,20), col=c("yellow","red", "blue"))

plot(hexgrid, col= "grey", border = "grey")
for(i in 1:length(weights.ResGeo2)){
	inter <- gcIntermediate(c(centroids.ResGeo2[i,1],centroids.ResGeo2[i,2]), c(centroids.ResGeo2[i,3],centroids.ResGeo2[i,4]), n=50, addStartEnd=T)
	#colcol <- rbPal(10)[as.numeric(cut(weights.ResGeo2, breaks=10))]
	#lines(inter, lwd=(weights.ResGeo2[i]*1)+1, col=colcol[i])
	lines(inter, lwd=1, col="orange")
	points(centroids.ResGeo2[i,1], centroids.ResGeo2[i,2], pch=20, col="red", cex=0.3)
	points(centroids.ResGeo2[i,3], centroids.ResGeo2[i,4], pch=20, col="blue", cex=0.3)
}
mtext("E", cex=1.3, side=3, line=-2.2, at=-180)
#legend("bottomleft", inset=c(0, 0.1), bg="white", box.col="white", title="", c("20 <= Rank < 100 ", "Breeding centroid", "Non-breeding centroid"), lty=c(1,NA,NA), pch=c(NA,20,20), col=c("orange","red", "blue"))

plot(hexgrid, col= "grey", border = "grey")
for(i in 1:length(weights.ResGeo3)){
	inter <- gcIntermediate(c(centroids.ResGeo3[i,1],centroids.ResGeo3[i,2]), c(centroids.ResGeo3[i,3],centroids.ResGeo3[i,4]), n=50, addStartEnd=T)
	#colcol <- rbPal(10)[as.numeric(cut(weights.ResGeo3, breaks=10))]
	#lines(inter, lwd=(weights.ResGeo3[i]*1)+1, col=colcol[i])
	lines(inter, lwd=1, col="brown4")
	points(centroids.ResGeo3[i,1], centroids.ResGeo3[i,2], pch=20, col="red", cex=0.3)
	points(centroids.ResGeo3[i,3], centroids.ResGeo3[i,4], pch=20, col="blue", cex=0.3)
}
mtext("F", cex=1.3, side=3, line=-2.2, at=-180)
#legend("bottomleft", inset=c(0, 0.1), bg="white", box.col="white", title="", c("Rank >= 100", "Breeding centroid", "Non-breeding cetroid"), lty=c(1,NA,NA), pch=c(NA,20,20), col=c("brown4","red", "blue"))


plot(hexgrid, col= "grey", border = "grey")
for(i in 1:length(weights.nicheResGeo1)){
	inter <- gcIntermediate(c(centroids.nicheResGeo1[i,1],centroids.nicheResGeo1[i,2]), c(centroids.nicheResGeo1[i,3],centroids.nicheResGeo1[i,4]), n=50, addStartEnd=T)
	#colcol <- rbPal1(10)[as.numeric(cut(weights.nicheResGeo1, breaks=10))]
	#lines(inter, lwd=(weights.nicheResGeo1[i]*1)+1, col=colcol[i])
	lines(inter, lwd=1, col="yellow")
	points(centroids.nicheResGeo1[i,1], centroids.nicheResGeo1[i,2], pch=20, col="red", cex=0.3)
	points(centroids.nicheResGeo1[i,3], centroids.nicheResGeo1[i,4], pch=20, col="blue", cex=0.3)
}
mtext("G", cex=1.3, side=3, line=-2.2, at=-180)
#legend("bottomleft", inset=c(0, 0.1), bg="white", box.col="white", title="", c("Rank < 20", "Breeding centroid", "Non-breeding centroid"), lty=c(1,NA,NA), pch=c(NA,20,20), col=c("yellow","red", "blue"))

plot(hexgrid, col= "grey", border = "grey")
for(i in 1:length(weights.nicheResGeo2)){
	inter <- gcIntermediate(c(centroids.nicheResGeo2[i,1],centroids.nicheResGeo2[i,2]), c(centroids.nicheResGeo2[i,3],centroids.nicheResGeo2[i,4]), n=50, addStartEnd=T)
	#colcol <- rbPal(10)[as.numeric(cut(weights.nicheResGeo2, breaks=10))]
	#lines(inter, lwd=(weights.nicheResGeo2[i]*1)+1, col=colcol[i])
	lines(inter, lwd=1, col="orange")
	points(centroids.nicheResGeo2[i,1], centroids.nicheResGeo2[i,2], pch=20, col="red", cex=0.3)
	points(centroids.nicheResGeo2[i,3], centroids.nicheResGeo2[i,4], pch=20, col="blue", cex=0.3)
}
mtext("H", cex=1.3, side=3, line=-2.2, at=-180)
#legend("bottomleft", inset=c(0, 0.1), bg="white", box.col="white", title="", c("20 <= Rank < 100 ", "Breeding centroid", "Non-breeding centroid"), lty=c(1,NA,NA), pch=c(NA,20,20), col=c("orange","red", "blue"))

plot(hexgrid, col= "grey", border = "grey")
for(i in 1:length(weights.nicheResGeo3)){
	inter <- gcIntermediate(c(centroids.nicheResGeo3[i,1],centroids.nicheResGeo3[i,2]), c(centroids.nicheResGeo3[i,3],centroids.nicheResGeo3[i,4]), n=50, addStartEnd=T)
	#colcol <- rbPal(10)[as.numeric(cut(weights.nicheResGeo3, breaks=10))]
	#lines(inter, lwd=(weights.nicheResGeo3[i]*1)+1, col=colcol[i])
	lines(inter, lwd=1, col="brown4")
	points(centroids.nicheResGeo3[i,1], centroids.nicheResGeo3[i,2], pch=20, col="red", cex=0.3)
	points(centroids.nicheResGeo3[i,3], centroids.nicheResGeo3[i,4], pch=20, col="blue", cex=0.3)
}
mtext("I", cex=1.3, side=3, line=-2.2, at=-180)
#legend("bottomleft", inset=c(0, 0.1), bg="white", box.col="white", title="", c("Rank >= 100", "Breeding centroid", "Non-breeding centroid"), lty=c(1,NA,NA), pch=c(NA,20,20), col=c("brown4","red", "blue"))



































## Run the models (separately for each combination of longitudinal hemisphere and breeding season)

rescale <- function(x) (x-min(x,na.rm=T))/(max(x,na.rm=T) - min(x,na.rm=T))

# Western Hemisphere #
#load("D:/shares/Marius_Somveille/niche analyses/ranges_smlWH.RData")
#range.sim_WH <- range.sim
#range.simNB_WH <- range.simNB

rank_nicheDistObs_WH <- vector() 
rank_nicheDistResNB_WH <- vector() 
rank_nicheDistResBR_WH <- vector() 
rank_geoDistObs_WH <- vector() 
rank_geoDistResNB_WH <- vector() 
rank_geoDistResBR_WH <- vector() 
rank_resourcesObs_WH <- vector() 
rank_resourcesResNB_WH <- vector() 
rank_resourcesResBR_WH <- vector() 
rank_nicheDist_geoDist_Obs_WH <- vector() 
rank_nicheDist_geoDist_ResNB_WH <- vector() 
rank_nicheDist_geoDist_ResBR_WH <- vector() 
rank_nicheDist_resources_Obs_WH <- vector() 
rank_nicheDist_resources_ResNB_WH <- vector() 
rank_nicheDist_resources_ResBR_WH <- vector() 
rank_resources_geoDist_Obs_WH <- vector() 
rank_resources_geoDist_ResNB_WH <- vector() 
rank_resources_geoDist_ResBR_WH <- vector() 
rank_nicheDist_resources_geoDist_Obs_WH <- vector() 
rank_nicheDist_resources_geoDist_ResNB_WH <- vector() 
rank_nicheDist_resources_geoDist_ResBR_WH <- vector() 
rank_nicheDist_random_WH <- vector() 
rank_geoDist_random_WH <- vector() 
rank_resources_random_WH <- vector() 
rank_nicheDist_geoDist_random_WH <- vector() 
rank_nicheDist_resources_random_WH <- vector() 
rank_resources_geoDist_random_WH <- vector() 
rank_nicheDist_resources_geoDist_random_WH <- vector() 
nicheDistObs_WH <- vector() 
geoDistObs_WH <- vector() 
resourcesObs_WH <- vector() 
nicheDistObs_rescaled_WH <- vector() 
nicheDistResNB_rescaled_WH <- vector() 
nicheDistResBR_rescaled_WH <- vector() 
geoDistObs_rescaled_WH <- vector() 
geoDistResNB_rescaled_WH <- vector() 
geoDistResBR_rescaled_WH <- vector() 
resourcesObs_rescaled_WH <- vector() 
resourcesResNB_rescaled_WH <- vector() 
resourcesResBR_rescaled_WH <- vector() 

for(j in 1:length(selectedSpeciesNH_WH)){

#resourcesScarcity_julyWH <- -((0.75*NDVI_ns_WH) - NDVI_nw_WH)
#resourcesScarcity_januaryWH <- -(NDVI_nw_WH - (0.75*NDVI_ns_WH))

## BR sim, NB obs

  nicheDistSim1 <- vector()
  geoDistSim1 <- vector()
  resourcesSim1 <- vector()
  #competSim1 <- vector()
  for(k in 1:99){ #length(range.sim_WH[[j]])){
    range.simu.BR <- range.sim_WH[[j]][[k]]

	nicheDistSim1[k] <- dist(rbind(c(mean(TempNW_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]), mean(PrecNW_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])), c(mean(TempNS_WH[range.simu.BR]), mean(PrecNS_WH[range.simu.BR]))))

	## breeding niche ##
	br <- cbind(TempNS_WH[range.simu.BR], PrecNS_WH[range.simu.BR])
	br_kernel <- kde2d(br[,1], br[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
	# Convert the kernel into a raster
	br_kernel_2_mat <- rep(br_kernel[[1]][1],50)
	for(i in 2:50){
		br_kernel_2_mat <- c(br_kernel_2_mat, rep(br_kernel[[1]][i],50))
	}
	br_kernel_2_mat2 <- rep(br_kernel[[2]],50)
	br_kernel_2_mat3 <- br_kernel[[3]][1,]
	for(i in 2:50){
		br_kernel_2_mat3 <- c(br_kernel_2_mat3, br_kernel[[3]][i,])
	}
	br_kernel_mat <- as.data.frame(cbind(br_kernel_2_mat, br_kernel_2_mat2, br_kernel_2_mat3))
	breeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(breeding.niche.raster) <- extent(c(-3, 3, -3,3))
	breeding.niche.raster <- rasterize(br_kernel_mat[,1:2], breeding.niche.raster, br_kernel_mat[,3])
	breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))
	# Keep only the top 99% of the kernel, set the rest to 0
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(breeding.niche.raster), decreasing=T)[i]
	}
	breeding.niche.raster[which(as.vector(breeding.niche.raster) < sort(as.vector(breeding.niche.raster), decreasing=T)[i])] = 0
	breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))

	## non-breeding niche ##
	nb <- cbind(TempNW_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)], PrecNW_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])
	nb_kernel <- kde2d(nb[,1], nb[,2], h=0.5, n=50, lims=c(-3, 3, -3,3))
	nb_kernel_2_mat <- rep(nb_kernel[[1]][1],50)
	for(i in 2:50){
		nb_kernel_2_mat <- c(nb_kernel_2_mat, rep(nb_kernel[[1]][i],50))
	}
	nb_kernel_2_mat2 <- rep(nb_kernel[[2]],50)
	nb_kernel_2_mat3 <- nb_kernel[[3]][1,]
	for(i in 2:50){
		nb_kernel_2_mat3 <- c(nb_kernel_2_mat3, nb_kernel[[3]][i,])
	}
	nb_kernel_mat <- as.data.frame(cbind(nb_kernel_2_mat, nb_kernel_2_mat2, nb_kernel_2_mat3))
	nonbreeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(nonbreeding.niche.raster) <- extent(c(-3, 3, -3,3))
	nonbreeding.niche.raster <- rasterize(nb_kernel_mat[,1:2], nonbreeding.niche.raster, nb_kernel_mat[,3])
	nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i]
	}
	nonbreeding.niche.raster[which(as.vector(nonbreeding.niche.raster) < sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i])] = 0
	nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))

	# Compute observed distance
	nicheDistSim1[k] <- emd(breeding.niche.raster, nonbreeding.niche.raster, threshold=2)


    geoDistSim1[k] <- rdist.earth(t(as.matrix(c(mean(west_Hem[range.simu.BR,1]), mean(west_Hem[range.simu.BR,2])))),t(as.matrix(c(mean(west_Hem[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1),1]), mean(west_Hem[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1),2])))), miles = F)

  res_JulySim <- resourcesScarcity_julyWH[range.simu.BR]
  resourcesSim1[k] <- mean(res_JulySim) + mean(resourcesScarcity_januaryWH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])
  
  #competSim1[k] <- mean(compet_july_WH[range.simu.BR]) + mean(compet_january_WH[which(PresAbs_NB_NH_WH_compet[,j] == 1)])
  
  }

## NB sim, BR obs

  nicheDistSim2 <- vector()
  geoDistSim2 <- vector()
  resourcesSim2 <- vector()
  #competSim2 <- vector()
  for(k in 1:99){ #length(range.simNB_WH[[j]])){
    #range.simu.BR <- match(range.simNB_WH[[j]][[k]], ID_WH)
    #range.simu.BR <- range.simu.BR[which(is.na(range.simu.BR)==F)]
    range.simu.NB <- range.simNB_WH[[j]][[k]]

    
	## breeding niche ##
	br <- cbind(TempNS_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)], PrecNS_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])
	br_kernel <- kde2d(br[,1], br[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
	# Convert the kernel into a raster
	br_kernel_2_mat <- rep(br_kernel[[1]][1],50)
	for(i in 2:50){
		br_kernel_2_mat <- c(br_kernel_2_mat, rep(br_kernel[[1]][i],50))
	}
	br_kernel_2_mat2 <- rep(br_kernel[[2]],50)
	br_kernel_2_mat3 <- br_kernel[[3]][1,]
	for(i in 2:50){
		br_kernel_2_mat3 <- c(br_kernel_2_mat3, br_kernel[[3]][i,])
	}
	br_kernel_mat <- as.data.frame(cbind(br_kernel_2_mat, br_kernel_2_mat2, br_kernel_2_mat3))
	breeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(breeding.niche.raster) <- extent(c(-3, 3, -3,3))
	breeding.niche.raster <- rasterize(br_kernel_mat[,1:2], breeding.niche.raster, br_kernel_mat[,3])
	breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))
	# Keep only the top 99% of the kernel, set the rest to 0
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(breeding.niche.raster), decreasing=T)[i]
	}
	breeding.niche.raster[which(as.vector(breeding.niche.raster) < sort(as.vector(breeding.niche.raster), decreasing=T)[i])] = 0
	breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))

	## non-breeding niche ##
	nb <- cbind(TempNW_WH[range.simu.NB], PrecNW_WH[range.simu.NB])
	nb_kernel <- kde2d(nb[,1], nb[,2], h=0.5, n=50, lims=c(-3, 3, -3,3))
	nb_kernel_2_mat <- rep(nb_kernel[[1]][1],50)
	for(i in 2:50){
		nb_kernel_2_mat <- c(nb_kernel_2_mat, rep(nb_kernel[[1]][i],50))
	}
	nb_kernel_2_mat2 <- rep(nb_kernel[[2]],50)
	nb_kernel_2_mat3 <- nb_kernel[[3]][1,]
	for(i in 2:50){
		nb_kernel_2_mat3 <- c(nb_kernel_2_mat3, nb_kernel[[3]][i,])
	}
	nb_kernel_mat <- as.data.frame(cbind(nb_kernel_2_mat, nb_kernel_2_mat2, nb_kernel_2_mat3))
	nonbreeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(nonbreeding.niche.raster) <- extent(c(-3, 3, -3,3))
	nonbreeding.niche.raster <- rasterize(nb_kernel_mat[,1:2], nonbreeding.niche.raster, nb_kernel_mat[,3])
	nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i]
	}
	nonbreeding.niche.raster[which(as.vector(nonbreeding.niche.raster) < sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i])] = 0
	nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))

	# Compute observed distance
	nicheDistSim2[k] <- emd(breeding.niche.raster, nonbreeding.niche.raster, threshold=2)


    geoDistSim2[k] <- rdist.earth(t(as.matrix(c(mean(west_Hem[range.simu.NB,1]), mean(west_Hem[range.simu.NB,2])))),t(as.matrix(c(mean(west_Hem[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1),1]), mean(west_Hem[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1),2])))), miles = F)

  res_JanuarySim <- resourcesScarcity_januaryWH[range.simu.NB]
  resourcesSim2[k] <- mean(res_JanuarySim) + mean(resourcesScarcity_julyWH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])
  
    #competSim2[k] <- mean(compet_july_WH[which(PresAbs_BR_NH_WH_compet[,j] == 1)]) + mean(compet_january_WH[range.simu.NB])
  
  }


## breeding niche ##
br <- cbind(TempNS_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)], PrecNS_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])
br_kernel <- kde2d(br[,1], br[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
# Convert the kernel into a raster
br_kernel_2_mat <- rep(br_kernel[[1]][1],50)
for(i in 2:50){
	br_kernel_2_mat <- c(br_kernel_2_mat, rep(br_kernel[[1]][i],50))
}	
br_kernel_2_mat2 <- rep(br_kernel[[2]],50)
br_kernel_2_mat3 <- br_kernel[[3]][1,]
for(i in 2:50){
	br_kernel_2_mat3 <- c(br_kernel_2_mat3, br_kernel[[3]][i,])
}
br_kernel_mat <- as.data.frame(cbind(br_kernel_2_mat, br_kernel_2_mat2, br_kernel_2_mat3))
breeding.niche.raster <- raster(ncol=50, nrow=50)
extent(breeding.niche.raster) <- extent(c(-3, 3, -3,3))
breeding.niche.raster <- rasterize(br_kernel_mat[,1:2], breeding.niche.raster, br_kernel_mat[,3])
breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))
# Keep only the top 99% of the kernel, set the rest to 0
thres = 0
i=0
while(thres <= 0.99){
	i = i+1
	thres = thres + sort(as.vector(breeding.niche.raster), decreasing=T)[i]
}
breeding.niche.raster[which(as.vector(breeding.niche.raster) < sort(as.vector(breeding.niche.raster), decreasing=T)[i])] = 0
breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))

## non-breeding niche ##
nb <- cbind(TempNW_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)], PrecNW_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])
nb_kernel <- kde2d(nb[,1], nb[,2], h=0.5, n=50, lims=c(-3, 3, -3,3))
nb_kernel_2_mat <- rep(nb_kernel[[1]][1],50)
for(i in 2:50){
	nb_kernel_2_mat <- c(nb_kernel_2_mat, rep(nb_kernel[[1]][i],50))
}
nb_kernel_2_mat2 <- rep(nb_kernel[[2]],50)
nb_kernel_2_mat3 <- nb_kernel[[3]][1,]
for(i in 2:50){
	nb_kernel_2_mat3 <- c(nb_kernel_2_mat3, nb_kernel[[3]][i,])
}
nb_kernel_mat <- as.data.frame(cbind(nb_kernel_2_mat, nb_kernel_2_mat2, nb_kernel_2_mat3))
nonbreeding.niche.raster <- raster(ncol=50, nrow=50)
extent(nonbreeding.niche.raster) <- extent(c(-3, 3, -3,3))
nonbreeding.niche.raster <- rasterize(nb_kernel_mat[,1:2], nonbreeding.niche.raster, nb_kernel_mat[,3])
nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))
thres = 0
i=0
while(thres <= 0.99){
	i = i+1
	thres = thres + sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i]
}
nonbreeding.niche.raster[which(as.vector(nonbreeding.niche.raster) < sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i])] = 0
nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))

# Compute observed distance
nicheDistObs <- emd(breeding.niche.raster, nonbreeding.niche.raster, threshold=2)
	
	
## stay resident on non-breeding ground ##
resnb <- cbind(TempNS_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)], PrecNS_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])
resnb_kernel <- kde2d(resnb[,1], resnb[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
# Convert the kernel into a raster
resnb_kernel_2_mat <- rep(resnb_kernel[[1]][1],50)
for(i in 2:50){
	resnb_kernel_2_mat <- c(resnb_kernel_2_mat, rep(resnb_kernel[[1]][i],50))
}
resnb_kernel_2_mat2 <- rep(resnb_kernel[[2]],50)
resnb_kernel_2_mat3 <- resnb_kernel[[3]][1,]
for(i in 2:50){
	resnb_kernel_2_mat3 <- c(resnb_kernel_2_mat3, resnb_kernel[[3]][i,])
}
resnb_kernel_mat <- as.data.frame(cbind(resnb_kernel_2_mat, resnb_kernel_2_mat2, resnb_kernel_2_mat3))
resident.nonbreeding.niche.raster <- raster(ncol=50, nrow=50)
extent(resident.nonbreeding.niche.raster) <- extent(c(-3, 3, -3,3))
resident.nonbreeding.niche.raster <- rasterize(resnb_kernel_mat[,1:2], resident.nonbreeding.niche.raster, resnb_kernel_mat[,3])
resident.nonbreeding.niche.raster = resident.nonbreeding.niche.raster / sum(as.vector(resident.nonbreeding.niche.raster))
# Keep only the top 99% of the kernel, set the rest to 0
thres = 0
i=0
while(thres <= 0.99){
	i = i+1
	thres = thres + sort(as.vector(resident.nonbreeding.niche.raster), decreasing=T)[i]
}
resident.nonbreeding.niche.raster[which(as.vector(resident.nonbreeding.niche.raster) < sort(as.vector(resident.nonbreeding.niche.raster), decreasing=T)[i])] = 0
resident.nonbreeding.niche.raster = resident.nonbreeding.niche.raster / sum(as.vector(resident.nonbreeding.niche.raster))
	
## stay resident on breeding ground ##
resbr <- cbind(TempNW_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)], PrecNW_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])
resbr_kernel <- kde2d(resbr[,1], resbr[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
# Convert the kernel into a raster
resbr_kernel_2_mat <- rep(resbr_kernel[[1]][1],50)
for(i in 2:50){
	resbr_kernel_2_mat <- c(resbr_kernel_2_mat, rep(resbr_kernel[[1]][i],50))
}
resbr_kernel_2_mat2 <- rep(resbr_kernel[[2]],50)
resbr_kernel_2_mat3 <- resbr_kernel[[3]][1,]
for(i in 2:50){
	resbr_kernel_2_mat3 <- c(resbr_kernel_2_mat3, resbr_kernel[[3]][i,])
}
resbr_kernel_mat <- as.data.frame(cbind(resbr_kernel_2_mat, resbr_kernel_2_mat2, resbr_kernel_2_mat3))
resident.breeding.niche.raster <- raster(ncol=50, nrow=50)
extent(resident.breeding.niche.raster) <- extent(c(-3, 3, -3,3))
resident.breeding.niche.raster <- rasterize(resbr_kernel_mat[,1:2], resident.breeding.niche.raster, resbr_kernel_mat[,3])
resident.breeding.niche.raster = resident.breeding.niche.raster / sum(as.vector(resident.breeding.niche.raster))
# Keep only the top 99% of the kernel, set the rest to 0
thres = 0
i=0
while(thres <= 0.99){
	i = i+1
	thres = thres + sort(as.vector(resident.breeding.niche.raster), decreasing=T)[i]
}
resident.breeding.niche.raster[which(as.vector(resident.breeding.niche.raster) < sort(as.vector(resident.breeding.niche.raster), decreasing=T)[i])] = 0
resident.breeding.niche.raster = resident.breeding.niche.raster / sum(as.vector(resident.breeding.niche.raster))

# Compute distances if resident
nicheDistResNB <- emd(resident.nonbreeding.niche.raster, nonbreeding.niche.raster, threshold=2)
nicheDistResBR <- emd(breeding.niche.raster, resident.breeding.niche.raster, threshold=2)

}


all_niche <- rescale(c(nicheDistObs, nicheDistResNB, nicheDistResBR, nicheDistSim1, nicheDistSim2))







nicheDistObs <- dist(rbind(c(mean(TempNW_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]), mean(PrecNW_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])), c(mean(TempNS_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]), mean(PrecNS_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]))))

nicheDistResNB <- dist(rbind(c(mean(TempNW_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]), mean(PrecNW_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])), c(mean(TempNS_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]), mean(PrecNS_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]))))

nicheDistResBR <- dist(rbind(c(mean(TempNW_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]), mean(PrecNW_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])), c(mean(TempNS_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]), mean(PrecNS_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]))))

all_niche <- rescale(c(nicheDistObs, nicheDistResNB, nicheDistResBR, nicheDistSim1, nicheDistSim2))


geoDistObs <- rdist.earth(t(as.matrix(c(mean(west_Hem[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1),1]), mean(west_Hem[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1),2])))),t(as.matrix(c(mean(west_Hem[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1),1]), mean(west_Hem[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1),2])))), miles = F)

all_geo <- rescale(c(geoDistObs, 0, 0, geoDistSim1, geoDistSim2))


resourcesObs <- mean(resourcesScarcity_julyWH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]) + mean(resourcesScarcity_januaryWH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])

resourcesResNB <- mean(resourcesScarcity_julyWH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]) + mean(resourcesScarcity_januaryWH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])

resourcesResBR <- mean(resourcesScarcity_julyWH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]) + mean(resourcesScarcity_januaryWH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])

all_resources <- rescale(c(resourcesObs, resourcesResNB, resourcesResBR, resourcesSim1, resourcesSim2))

nngg <- all_niche + all_geo
rrgg <- all_resources + all_geo
nnrr <- all_niche + all_resources
nnrrgg <- all_niche + all_resources + all_geo

nicheDistObs_WH[j] <- nicheDistObs
geoDistObs_WH[j] <- geoDistObs
resourcesObs_WH[j] <- resourcesObs

nicheDistObs_rescaled_WH[j] <- all_niche[1]
nicheDistResNB_rescaled_WH[j] <- all_niche[2]
nicheDistResBR_rescaled_WH[j] <- all_niche[3]
geoDistObs_rescaled_WH[j] <- all_geo[1]
geoDistResNB_rescaled_WH[j] <- 0
geoDistResBR_rescaled_WH[j] <- 0
resourcesObs_rescaled_WH[j] <- all_resources[1]
resourcesResNB_rescaled_WH[j] <- all_resources[2]
resourcesResBR_rescaled_WH[j] <- all_resources[3]

rank_nicheDistObs_WH[j] <- rank(all_niche)[1]
rank_nicheDistResNB_WH[j] <- rank(all_niche)[2]
rank_nicheDistResBR_WH[j] <- rank(all_niche)[3]
rank_geoDistObs_WH[j] <- rank(all_geo)[1]
rank_geoDistResNB_WH[j] <- rank(all_geo)[2]
rank_geoDistResBR_WH[j] <- rank(all_geo)[3]
rank_resourcesObs_WH[j] <- rank(all_resources)[1]
rank_resourcesResNB_WH[j] <- rank(all_resources)[2]
rank_resourcesResBR_WH[j] <- rank(all_resources)[3]
rank_nicheDist_geoDist_Obs_WH[j] <- rank(nngg)[1]
rank_nicheDist_geoDist_ResNB_WH[j] <- rank(nngg)[2]
rank_nicheDist_geoDist_ResBR_WH[j] <- rank(nngg)[3]
rank_nicheDist_resources_Obs_WH[j] <- rank(nnrr)[1]
rank_nicheDist_resources_ResNB_WH[j] <- rank(nnrr)[2]
rank_nicheDist_resources_ResBR_WH[j] <- rank(nnrr)[3]
rank_resources_geoDist_Obs_WH[j] <- rank(rrgg)[1]
rank_resources_geoDist_ResNB_WH[j] <- rank(rrgg)[2]
rank_resources_geoDist_ResBR_WH[j] <- rank(rrgg)[3]
rank_nicheDist_resources_geoDist_Obs_WH[j] <- rank(nnrrgg)[1]
rank_nicheDist_resources_geoDist_ResNB_WH[j] <- rank(nnrrgg)[2]
rank_nicheDist_resources_geoDist_ResBR_WH[j] <- rank(nnrrgg)[3]

rank_nicheDist_random_WH[j] <- sample(rank(all_niche)[4:201], 1)
rank_geoDist_random_WH[j] <- sample(rank(all_geo)[4:201], 1)
rank_resources_random_WH[j] <- sample(rank(all_resources)[4:201], 1)
rank_nicheDist_geoDist_random_WH[j] <- sample(rank(nngg)[4:201], 1)
rank_nicheDist_resources_random_WH[j] <- sample(rank(nnrr)[4:201], 1)
rank_resources_geoDist_random_WH[j] <- sample(rank(rrgg)[4:201], 1)
rank_nicheDist_resources_geoDist_random_WH[j] <- sample(rank(nnrrgg)[4:201], 1)

print(j)
  }
   
   


######  EH


load("D:/shares/Marius_Somveille/niche analyses/ranges_smlEH.RData")
range.sim_EH <- range.sim
range.simNB_EH <- range.simNB

rank_nicheDistObs_EH <- vector() 
rank_nicheDistResNB_EH <- vector() 
rank_nicheDistResBR_EH <- vector() 
rank_geoDistObs_EH <- vector() 
rank_geoDistResNB_EH <- vector() 
rank_geoDistResBR_EH <- vector() 
rank_resourcesObs_EH <- vector() 
rank_resourcesResNB_EH <- vector() 
rank_resourcesResBR_EH <- vector() 
rank_nicheDist_geoDist_Obs_EH <- vector() 
rank_nicheDist_geoDist_ResNB_EH <- vector() 
rank_nicheDist_geoDist_ResBR_EH <- vector() 
rank_nicheDist_resources_Obs_EH <- vector() 
rank_nicheDist_resources_ResNB_EH <- vector() 
rank_nicheDist_resources_ResBR_EH <- vector() 
rank_resources_geoDist_Obs_EH <- vector() 
rank_resources_geoDist_ResNB_EH <- vector() 
rank_resources_geoDist_ResBR_EH <- vector() 
rank_nicheDist_resources_geoDist_Obs_EH <- vector() 
rank_nicheDist_resources_geoDist_ResNB_EH <- vector() 
rank_nicheDist_resources_geoDist_ResBR_EH <- vector() 
rank_nicheDist_random_EH <- vector() 
rank_geoDist_random_EH <- vector() 
rank_resources_random_EH <- vector() 
rank_nicheDist_geoDist_random_EH <- vector() 
rank_nicheDist_resources_random_EH <- vector() 
rank_resources_geoDist_random_EH <- vector() 
rank_nicheDist_resources_geoDist_random_EH <- vector()
nicheDistObs_EH <- vector() 
geoDistObs_EH <- vector() 
resourcesObs_EH <- vector() 
nicheDistObs_rescaled_EH <- vector() 
nicheDistResNB_rescaled_EH <- vector() 
nicheDistResBR_rescaled_EH <- vector() 
geoDistObs_rescaled_EH <- vector() 
geoDistResNB_rescaled_EH <- vector() 
geoDistResBR_rescaled_EH <- vector() 
resourcesObs_rescaled_EH <- vector() 
resourcesResNB_rescaled_EH <- vector() 
resourcesResBR_rescaled_EH <- vector() 

for(j in 1:length(selectedSpeciesNH_EH)){

resourcesScarcity_julyEH <- -((0.75*NDVI_ns_EH) - NDVI_nw_EH)
resourcesScarcity_januaryEH <- -(NDVI_nw_EH - (0.75*NDVI_ns_EH))

## BR sim, NB obs

  nicheDistSim1 <- vector()
  geoDistSim1 <- vector()
  resourcesSim1 <- vector()
  #competSim1 <- vector()
  for(k in 1:99){ #length(range.sim_EH[[j]])){
    #range.simu.BR <- match(range.sim_EH[[j]][[k]], ID_EH)
    #range.simu.BR <- range.simu.BR[which(is.na(range.simu.BR)==F)]
    range.simu.BR <- range.sim_EH[[j]][[k]]


nicheDistSim1[k] <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])), c(mean(TempNS_EH[range.simu.BR]), mean(PrecNS_EH[range.simu.BR]))))

    geoDistSim1[k] <- rdist.earth(t(as.matrix(c(mean(east_Hem[range.simu.BR,1]), mean(east_Hem[range.simu.BR,2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1),1]), mean(east_Hem[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1),2])))), miles = F)

  res_JulySim <- resourcesScarcity_julyEH[range.simu.BR]
  resourcesSim1[k] <- mean(res_JulySim) + mean(resourcesScarcity_januaryEH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])
  
  #competSim1[k] <- mean(compet_july_EH[range.simu.BR]) + mean(compet_january_EH[which(PresAbs_NB_NH_EH_compet[,j] == 1)])
  }

## NB sim, BR obs

  nicheDistSim2 <- vector()
  geoDistSim2 <- vector()
  resourcesSim2 <- vector()
  #competSim2 <- vector()
  for(k in 1:99){ #length(range.simNB_EH[[j]])){
    #range.simu.BR <- match(range.simNB_EH[[j]][[k]], ID_EH)
    #range.simu.BR <- range.simu.BR[which(is.na(range.simu.BR)==F)]
    range.simu.NB <- range.simNB_EH[[j]][[k]]

    nicheDistSim2[k] <- dist(rbind(c(mean(TempNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])), c(mean(TempNW_EH[range.simu.NB]), mean(PrecNW_EH[range.simu.NB]))))
    
    geoDistSim2[k] <- rdist.earth(t(as.matrix(c(mean(east_Hem[range.simu.NB,1]), mean(east_Hem[range.simu.NB,2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1),1]), mean(east_Hem[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1),2])))), miles = F)

  res_JanuarySim <- resourcesScarcity_januaryEH[range.simu.NB]
  resourcesSim2[k] <- mean(res_JanuarySim) + mean(resourcesScarcity_julyEH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])
  
  #competSim2[k] <- mean(compet_july_EH[which(PresAbs_BR_NH_EH_compet[,j] == 1)]) + mean(compet_january_EH[range.simu.NB])
  }


nicheDistObs <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])), c(mean(TempNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]))))

nicheDistResNB <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])), c(mean(TempNS_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]))))

nicheDistResBR <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])), c(mean(TempNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]))))

all_niche <- rescale(c(nicheDistObs, nicheDistResNB, nicheDistResBR, nicheDistSim1, nicheDistSim2))


geoDistObs <- rdist.earth(t(as.matrix(c(mean(east_Hem[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1),1]), mean(east_Hem[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1),2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1),1]), mean(east_Hem[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1),2])))), miles = F)

all_geo <- rescale(c(geoDistObs, 0, 0, geoDistSim1, geoDistSim2))


resourcesObs <- mean(resourcesScarcity_julyEH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]) + mean(resourcesScarcity_januaryEH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])

resourcesResNB <- mean(resourcesScarcity_julyEH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]) + mean(resourcesScarcity_januaryEH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])

resourcesResBR <- mean(resourcesScarcity_julyEH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]) + mean(resourcesScarcity_januaryEH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])

all_resources <- rescale(c(resourcesObs, resourcesResNB, resourcesResBR, resourcesSim1, resourcesSim2))

nngg <- all_niche + all_geo
rrgg <- all_resources + all_geo
nnrr <- all_niche + all_resources
nnrrgg <- all_niche + all_resources + all_geo

nicheDistObs_EH[j] <- nicheDistObs
geoDistObs_EH[j] <- geoDistObs
resourcesObs_EH[j] <- resourcesObs

nicheDistObs_rescaled_EH[j] <- all_niche[1]
nicheDistResNB_rescaled_EH[j] <- all_niche[2]
nicheDistResBR_rescaled_EH[j] <- all_niche[3]
geoDistObs_rescaled_EH[j] <- all_geo[1]
geoDistResNB_rescaled_EH[j] <- 0
geoDistResBR_rescaled_EH[j] <- 0
resourcesObs_rescaled_EH[j] <- all_resources[1]
resourcesResNB_rescaled_EH[j] <- all_resources[2]
resourcesResBR_rescaled_EH[j] <- all_resources[3]

rank_nicheDistObs_EH[j] <- rank(all_niche)[1]
rank_nicheDistResNB_EH[j] <- rank(all_niche)[2]
rank_nicheDistResBR_EH[j] <- rank(all_niche)[3]
rank_geoDistObs_EH[j] <- rank(all_geo)[1]
rank_geoDistResNB_EH[j] <- rank(all_geo)[2]
rank_geoDistResBR_EH[j] <- rank(all_geo)[3]
rank_resourcesObs_EH[j] <- rank(all_resources)[1]
rank_resourcesResNB_EH[j] <- rank(all_resources)[2]
rank_resourcesResBR_EH[j] <- rank(all_resources)[3]
rank_nicheDist_geoDist_Obs_EH[j] <- rank(nngg)[1]
rank_nicheDist_geoDist_ResNB_EH[j] <- rank(nngg)[2]
rank_nicheDist_geoDist_ResBR_EH[j] <- rank(nngg)[3]
rank_nicheDist_resources_Obs_EH[j] <- rank(nnrr)[1]
rank_nicheDist_resources_ResNB_EH[j] <- rank(nnrr)[2]
rank_nicheDist_resources_ResBR_EH[j] <- rank(nnrr)[3]
rank_resources_geoDist_Obs_EH[j] <- rank(rrgg)[1]
rank_resources_geoDist_ResNB_EH[j] <- rank(rrgg)[2]
rank_resources_geoDist_ResBR_EH[j] <- rank(rrgg)[3]
rank_nicheDist_resources_geoDist_Obs_EH[j] <- rank(nnrrgg)[1]
rank_nicheDist_resources_geoDist_ResNB_EH[j] <- rank(nnrrgg)[2]
rank_nicheDist_resources_geoDist_ResBR_EH[j] <- rank(nnrrgg)[3]

rank_nicheDist_random_EH[j] <- sample(rank(all_niche)[4:201], 1)
rank_geoDist_random_EH[j] <- sample(rank(all_geo)[4:201], 1)
rank_resources_random_EH[j] <- sample(rank(all_resources)[4:201], 1)
rank_nicheDist_geoDist_random_EH[j] <- sample(rank(nngg)[4:201], 1)
rank_nicheDist_resources_random_EH[j] <- sample(rank(nnrr)[4:201], 1)
rank_resources_geoDist_random_EH[j] <- sample(rank(rrgg)[4:201], 1)
rank_nicheDist_resources_geoDist_random_EH[j] <- sample(rank(nnrrgg)[4:201], 1)

print(j)
  }
  
  
for(j in 380:381){

resourcesScarcity_julyEH <- -((0.75*NDVI_ns_EH) - NDVI_nw_EH)
resourcesScarcity_januaryEH <- -(NDVI_nw_EH - (0.75*NDVI_ns_EH))

## BR sim, NB obs

  nicheDistSim1 <- vector()
  geoDistSim1 <- vector()
  resourcesSim1 <- vector()
  #competSim1 <- vector()
  for(k in 1:99){ #length(range.sim_EH[[j]])){
    #range.simu.BR <- match(range.sim_EH[[j]][[k]], ID_EH)
    #range.simu.BR <- range.simu.BR[which(is.na(range.simu.BR)==F)]
    range.simu.BR <- range.sim_EH[[j]][[k]]


nicheDistSim1[k] <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])), c(mean(TempNS_EH[range.simu.BR]), mean(PrecNS_EH[range.simu.BR]))))

    geoDistSim1[k] <- rdist.earth(t(as.matrix(c(mean(east_Hem[range.simu.BR,1]), mean(east_Hem[range.simu.BR,2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1),1]), mean(east_Hem[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1),2])))), miles = F)

  res_JulySim <- resourcesScarcity_julyEH[range.simu.BR]
  resourcesSim1[k] <- mean(res_JulySim) + mean(resourcesScarcity_januaryEH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])
  
  #competSim1[k] <- mean(compet_july_EH[range.simu.BR]) + mean(compet_january_EH[which(PresAbs_NB_SH_EH_compet[,j] == 1)])
  }

## NB sim, BR obs

  nicheDistSim2 <- vector()
  geoDistSim2 <- vector()
  resourcesSim2 <- vector()
  #competSim2 <- vector()
  for(k in 1:99){ #length(range.simNB_EH[[j]])){
    #range.simu.BR <- match(range.simNB_EH[[j]][[k]], ID_EH)
    #range.simu.BR <- range.simu.BR[which(is.na(range.simu.BR)==F)]
    range.simu.NB <- range.simNB_EH[[j]][[k]]

    nicheDistSim2[k] <- dist(rbind(c(mean(TempNS_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])), c(mean(TempNW_EH[range.simu.NB]), mean(PrecNW_EH[range.simu.NB]))))
    
    geoDistSim2[k] <- rdist.earth(t(as.matrix(c(mean(east_Hem[range.simu.NB,1]), mean(east_Hem[range.simu.NB,2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1),1]), mean(east_Hem[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1),2])))), miles = F)

  res_JanuarySim <- resourcesScarcity_januaryEH[range.simu.NB]
  resourcesSim2[k] <- mean(res_JanuarySim) + mean(resourcesScarcity_julyEH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])
  
  #competSim2[k] <- mean(compet_july_EH[which(PresAbs_BR_SH_EH_compet[,j] == 1)]) + mean(compet_january_EH[range.simu.NB])
  }


nicheDistObs <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])), c(mean(TempNS_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]))))

nicheDistResNB <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])), c(mean(TempNS_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]))))

nicheDistResBR <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])), c(mean(TempNS_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]))))

all_niche <- rescale(c(nicheDistObs, nicheDistResNB, nicheDistResBR, nicheDistSim1, nicheDistSim2))


geoDistObs <- rdist.earth(t(as.matrix(c(mean(east_Hem[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1),1]), mean(east_Hem[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1),2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1),1]), mean(east_Hem[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1),2])))), miles = F)

all_geo <- rescale(c(geoDistObs, 0, 0, geoDistSim1, geoDistSim2))


resourcesObs <- mean(resourcesScarcity_julyEH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]) + mean(resourcesScarcity_januaryEH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])

resourcesResNB <- mean(resourcesScarcity_julyEH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]) + mean(resourcesScarcity_januaryEH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])

resourcesResBR <- mean(resourcesScarcity_julyEH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]) + mean(resourcesScarcity_januaryEH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])

all_resources <- rescale(c(resourcesObs, resourcesResNB, resourcesResBR, resourcesSim1, resourcesSim2))

nngg <- all_niche + all_geo
rrgg <- all_resources + all_geo
nnrr <- all_niche + all_resources
nnrrgg <- all_niche + all_resources + all_geo

nicheDistObs_EH[j] <- nicheDistObs
geoDistObs_EH[j] <- geoDistObs
resourcesObs_EH[j] <- resourcesObs

nicheDistObs_rescaled_EH[j] <- all_niche[1]
nicheDistResNB_rescaled_EH[j] <- all_niche[2]
nicheDistResBR_rescaled_EH[j] <- all_niche[3]
geoDistObs_rescaled_EH[j] <- all_geo[1]
geoDistResNB_rescaled_EH[j] <- 0
geoDistResBR_rescaled_EH[j] <- 0
resourcesObs_rescaled_EH[j] <- all_resources[1]
resourcesResNB_rescaled_EH[j] <- all_resources[2]
resourcesResBR_rescaled_EH[j] <- all_resources[3]

rank_nicheDistObs_EH[j] <- rank(all_niche)[1]
rank_nicheDistResNB_EH[j] <- rank(all_niche)[2]
rank_nicheDistResBR_EH[j] <- rank(all_niche)[3]
rank_geoDistObs_EH[j] <- rank(all_geo)[1]
rank_geoDistResNB_EH[j] <- rank(all_geo)[2]
rank_geoDistResBR_EH[j] <- rank(all_geo)[3]
rank_resourcesObs_EH[j] <- rank(all_resources)[1]
rank_resourcesResNB_EH[j] <- rank(all_resources)[2]
rank_resourcesResBR_EH[j] <- rank(all_resources)[3]
rank_nicheDist_geoDist_Obs_EH[j] <- rank(nngg)[1]
rank_nicheDist_geoDist_ResNB_EH[j] <- rank(nngg)[2]
rank_nicheDist_geoDist_ResBR_EH[j] <- rank(nngg)[3]
rank_nicheDist_resources_Obs_EH[j] <- rank(nnrr)[1]
rank_nicheDist_resources_ResNB_EH[j] <- rank(nnrr)[2]
rank_nicheDist_resources_ResBR_EH[j] <- rank(nnrr)[3]
rank_resources_geoDist_Obs_EH[j] <- rank(rrgg)[1]
rank_resources_geoDist_ResNB_EH[j] <- rank(rrgg)[2]
rank_resources_geoDist_ResBR_EH[j] <- rank(rrgg)[3]
rank_nicheDist_resources_geoDist_Obs_EH[j] <- rank(nnrrgg)[1]
rank_nicheDist_resources_geoDist_ResNB_EH[j] <- rank(nnrrgg)[2]
rank_nicheDist_resources_geoDist_ResBR_EH[j] <- rank(nnrrgg)[3]

rank_nicheDist_random_EH[j] <- sample(rank(all_niche)[4:201], 1)
rank_geoDist_random_EH[j] <- sample(rank(all_geo)[4:201], 1)
rank_resources_random_EH[j] <- sample(rank(all_resources)[4:201], 1)
rank_nicheDist_geoDist_random_EH[j] <- sample(rank(nngg)[4:201], 1)
rank_nicheDist_resources_random_EH[j] <- sample(rank(nnrr)[4:201], 1)
rank_resources_geoDist_random_EH[j] <- sample(rank(rrgg)[4:201], 1)
rank_nicheDist_resources_geoDist_random_EH[j] <- sample(rank(nnrrgg)[4:201], 1)

print(j)
  }
   
   

######  SH WH

load("D:/shares/Marius_Somveille/niche analyses/ranges_smlSH.RData")


rank_nicheDistObs_SHWH <- vector() 
rank_nicheDistResNB_SHWH <- vector() 
rank_nicheDistResBR_SHWH <- vector() 
rank_geoDistObs_SHWH <- vector() 
rank_geoDistResNB_SHWH <- vector() 
rank_geoDistResBR_SHWH <- vector() 
rank_resourcesObs_SHWH <- vector() 
rank_resourcesResNB_SHWH <- vector() 
rank_resourcesResBR_SHWH <- vector() 
rank_nicheDist_geoDist_Obs_SHWH <- vector() 
rank_nicheDist_geoDist_ResNB_SHWH <- vector() 
rank_nicheDist_geoDist_ResBR_SHWH <- vector() 
rank_nicheDist_resources_Obs_SHWH <- vector() 
rank_nicheDist_resources_ResNB_SHWH <- vector() 
rank_nicheDist_resources_ResBR_SHWH <- vector() 
rank_resources_geoDist_Obs_SHWH <- vector() 
rank_resources_geoDist_ResNB_SHWH <- vector() 
rank_resources_geoDist_ResBR_SHWH <- vector() 
rank_nicheDist_resources_geoDist_Obs_SHWH <- vector() 
rank_nicheDist_resources_geoDist_ResNB_SHWH <- vector() 
rank_nicheDist_resources_geoDist_ResBR_SHWH <- vector()
rank_nicheDist_random_SHWH <- vector() 
rank_geoDist_random_SHWH <- vector() 
rank_resources_random_SHWH <- vector() 
rank_nicheDist_geoDist_random_SHWH <- vector() 
rank_nicheDist_resources_random_SHWH <- vector() 
rank_resources_geoDist_random_SHWH <- vector() 
rank_nicheDist_resources_geoDist_random_SHWH <- vector() 
nicheDistObs_SHWH <- vector() 
geoDistObs_SHWH <- vector() 
resourcesObs_SHWH <- vector() 
nicheDistObs_rescaled_SHWH <- vector() 
nicheDistResNB_rescaled_SHWH <- vector() 
nicheDistResBR_rescaled_SHWH <- vector() 
geoDistObs_rescaled_SHWH <- vector() 
geoDistResNB_rescaled_SHWH <- vector() 
geoDistResBR_rescaled_SHWH <- vector() 
resourcesObs_rescaled_SHWH <- vector() 
resourcesResNB_rescaled_SHWH <- vector() 
resourcesResBR_rescaled_SHWH <- vector() 

for(j in 1:length(selectedSpeciesSH_WH)){

resourcesScarcity_januaryWH <- -((0.75*NDVI_nw_WH) - NDVI_ns_WH)
resourcesScarcity_julyWH <- -(NDVI_ns_WH - (0.75*NDVI_nw_WH))

## BR sim, NB obs

  nicheDistSim1 <- vector()
  geoDistSim1 <- vector()
  resourcesSim1 <- vector()
  for(k in 1:99){ #length(range.sim_WH[[j]])){
    range.simu.BR <- range.sim_SHWH[[j]][[k]]


nicheDistSim1[k] <- dist(rbind(c(mean(TempNS_WH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]), mean(PrecNS_WH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)])), c(mean(TempNW_WH[range.simu.BR]), mean(PrecNW_WH[range.simu.BR]))))

    geoDistSim1[k] <- rdist.earth(t(as.matrix(c(mean(west_Hem[range.simu.BR,1]), mean(west_Hem[range.simu.BR,2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1),1]), mean(east_Hem[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1),2])))), miles = F)

  res_JanuarySim <- resourcesScarcity_januaryWH[range.simu.BR]
  resourcesSim1[k] <- mean(res_JanuarySim) + mean(resourcesScarcity_julyWH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)])
  
  }

## NB sim, BR obs

  nicheDistSim2 <- vector()
  geoDistSim2 <- vector()
  resourcesSim2 <- vector()
  for(k in 1:99){ #length(range.simNB_WH[[j]])){
    #range.simu.BR <- match(range.simNB_WH[[j]][[k]], ID_WH)
    #range.simu.BR <- range.simu.BR[which(is.na(range.simu.BR)==F)]
    range.simu.NB <- range.simNB_WH[[j]][[k]]

    nicheDistSim2[k] <- dist(rbind(c(mean(TempNW_WH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]), mean(PrecNW_WH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)])), c(mean(TempNS_WH[range.simu.NB]), mean(PrecNS_WH[range.simu.NB]))))
    
    geoDistSim2[k] <- rdist.earth(t(as.matrix(c(mean(west_Hem[range.simu.NB,1]), mean(west_Hem[range.simu.NB,2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1),1]), mean(east_Hem[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1),2])))), miles = F)

  res_JulySim <- resourcesScarcity_julyWH[range.simu.NB]
  resourcesSim2[k] <- mean(res_JulySim) + mean(resourcesScarcity_januaryWH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)])
  }

nicheDistObs <- dist(rbind(c(mean(TempNS_WH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]), mean(PrecNS_WH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)])), c(mean(TempNW_WH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]), mean(PrecNW_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]))))

nicheDistResNB <- dist(rbind(c(mean(TempNW_WH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]), mean(PrecNW_WH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)])), c(mean(TempNS_WH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]), mean(PrecNS_WH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]))))

nicheDistResBR <- dist(rbind(c(mean(TempNW_WH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]), mean(PrecNW_WH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)])), c(mean(TempNS_WH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]), mean(PrecNS_WH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]))))

all_niche <- rescale(c(nicheDistObs, nicheDistResNB, nicheDistResBR, nicheDistSim1, nicheDistSim2))


geoDistObs <- rdist.earth(t(as.matrix(c(mean(west_Hem[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1),1]), mean(west_Hem[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1),2])))),t(as.matrix(c(mean(west_Hem[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1),1]), mean(west_Hem[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1),2])))), miles = F)

all_geo <- rescale(c(geoDistObs, 0, 0, geoDistSim1, geoDistSim2))


resourcesObs <- mean(resourcesScarcity_julyWH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]) + mean(resourcesScarcity_januaryWH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)])

resourcesResNB <- mean(resourcesScarcity_julyWH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]) + mean(resourcesScarcity_januaryWH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)])

resourcesResBR <- mean(resourcesScarcity_julyWH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]) + mean(resourcesScarcity_januaryWH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)])

all_resources <- rescale(c(resourcesObs, resourcesResNB, resourcesResBR, resourcesSim1, resourcesSim2))

nngg <- all_niche + all_geo
rrgg <- all_resources + all_geo
nnrr <- all_niche + all_resources
nnrrgg <- all_niche + all_resources + all_geo

nicheDistObs_SHWH[j] <- nicheDistObs
geoDistObs_SHWH[j] <- geoDistObs
resourcesObs_SHWH[j] <- resourcesObs

nicheDistObs_rescaled_SHWH[j] <- all_niche[1]
nicheDistResNB_rescaled_SHWH[j] <- all_niche[2]
nicheDistResBR_rescaled_SHWH[j] <- all_niche[3]
geoDistObs_rescaled_SHWH[j] <- all_geo[1]
geoDistResNB_rescaled_SHWH[j] <- 0
geoDistResBR_rescaled_SHWH[j] <- 0
resourcesObs_rescaled_SHWH[j] <- all_resources[1]
resourcesResNB_rescaled_SHWH[j] <- all_resources[2]
resourcesResBR_rescaled_SHWH[j] <- all_resources[3]

rank_nicheDistObs_SHWH[j] <- rank(all_niche)[1]
rank_nicheDistResNB_SHWH[j] <- rank(all_niche)[2]
rank_nicheDistResBR_SHWH[j] <- rank(all_niche)[3]
rank_geoDistObs_SHWH[j] <- rank(all_geo)[1]
rank_geoDistResNB_SHWH[j] <- rank(all_geo)[2]
rank_geoDistResBR_SHWH[j] <- rank(all_geo)[3]
rank_resourcesObs_SHWH[j] <- rank(all_resources)[1]
rank_resourcesResNB_SHWH[j] <- rank(all_resources)[2]
rank_resourcesResBR_SHWH[j] <- rank(all_resources)[3]
rank_nicheDist_geoDist_Obs_SHWH[j] <- rank(nngg)[1]
rank_nicheDist_geoDist_ResNB_SHWH[j] <- rank(nngg)[2]
rank_nicheDist_geoDist_ResBR_SHWH[j] <- rank(nngg)[3]
rank_nicheDist_resources_Obs_SHWH[j] <- rank(nnrr)[1]
rank_nicheDist_resources_ResNB_SHWH[j] <- rank(nnrr)[2]
rank_nicheDist_resources_ResBR_SHWH[j] <- rank(nnrr)[3]
rank_resources_geoDist_Obs_SHWH[j] <- rank(rrgg)[1]
rank_resources_geoDist_ResNB_SHWH[j] <- rank(rrgg)[2]
rank_resources_geoDist_ResBR_SHWH[j] <- rank(rrgg)[3]
rank_nicheDist_resources_geoDist_Obs_SHWH[j] <- rank(nnrrgg)[1]
rank_nicheDist_resources_geoDist_ResNB_SHWH[j] <- rank(nnrrgg)[2]
rank_nicheDist_resources_geoDist_ResBR_SHWH[j] <- rank(nnrrgg)[3]

rank_nicheDist_random_SHWH[j] <- sample(rank(all_niche)[4:201], 1)
rank_geoDist_random_SHWH[j] <- sample(rank(all_geo)[4:201], 1)
rank_resources_random_SHWH[j] <- sample(rank(all_resources)[4:201], 1)
rank_nicheDist_geoDist_random_SHWH[j] <- sample(rank(nngg)[4:201], 1)
rank_nicheDist_resources_random_SHWH[j] <- sample(rank(nnrr)[4:201], 1)
rank_resources_geoDist_random_SHWH[j] <- sample(rank(rrgg)[4:201], 1)
rank_nicheDist_resources_geoDist_random_SHWH[j] <- sample(rank(nnrrgg)[4:201], 1)

print(j)
  }
   
   


######  SH EH


rank_nicheDistObs_SHEH <- vector() 
rank_nicheDistResNB_SHEH <- vector() 
rank_nicheDistResBR_SHEH <- vector() 
rank_geoDistObs_SHEH <- vector() 
rank_geoDistResNB_SHEH <- vector() 
rank_geoDistResBR_SHEH <- vector() 
rank_resourcesObs_SHEH <- vector() 
rank_resourcesResNB_SHEH <- vector() 
rank_resourcesResBR_SHEH <- vector() 
rank_nicheDist_geoDist_Obs_SHEH <- vector() 
rank_nicheDist_geoDist_ResNB_SHEH <- vector() 
rank_nicheDist_geoDist_ResBR_SHEH <- vector() 
rank_nicheDist_resources_Obs_SHEH <- vector() 
rank_nicheDist_resources_ResNB_SHEH <- vector() 
rank_nicheDist_resources_ResBR_SHEH <- vector() 
rank_resources_geoDist_Obs_SHEH <- vector() 
rank_resources_geoDist_ResNB_SHEH <- vector() 
rank_resources_geoDist_ResBR_SHEH <- vector() 
rank_nicheDist_resources_geoDist_Obs_SHEH <- vector() 
rank_nicheDist_resources_geoDist_ResNB_SHEH <- vector() 
rank_nicheDist_resources_geoDist_ResBR_SHEH <- vector()
rank_nicheDist_random_SHEH <- vector() 
rank_geoDist_random_SHEH <- vector() 
rank_resources_random_SHEH <- vector() 
rank_nicheDist_geoDist_random_SHEH <- vector() 
rank_nicheDist_resources_random_SHEH <- vector() 
rank_resources_geoDist_random_SHEH <- vector() 
rank_nicheDist_resources_geoDist_random_SHEH <- vector()
nicheDistObs_SHEH <- vector()
geoDistObs_SHEH <- vector()
resourcesObs_SHEH <- vector()
nicheDistObs_rescaled_SHEH <- vector() 
nicheDistResNB_rescaled_SHEH <- vector() 
nicheDistResBR_rescaled_SHEH <- vector() 
geoDistObs_rescaled_SHEH <- vector() 
geoDistResNB_rescaled_SHEH <- vector() 
geoDistResBR_rescaled_SHEH <- vector() 
resourcesObs_rescaled_SHEH <- vector() 
resourcesResNB_rescaled_SHEH <- vector() 
resourcesResBR_rescaled_SHEH <- vector() 

for(j in 1:length(selectedSpeciesSH_EH)){

resourcesScarcity_januaryEH <- -((0.75*NDVI_nw_EH) - NDVI_ns_EH)
resourcesScarcity_julyEH <- -(NDVI_ns_EH - (0.75*NDVI_nw_EH))

## BR sim, NB obs

  nicheDistSim1 <- vector()
  geoDistSim1 <- vector()
  resourcesSim1 <- vector()
  for(k in 1:99){ #length(range.sim_EH[[j]])){
    #range.simu.BR <- match(range.sim_EH[[j]][[k]], ID_EH)
    #range.simu.BR <- range.simu.BR[which(is.na(range.simu.BR)==F)]
    range.simu.BR <- range.sim_SHEH[[j]][[k]]


	nicheDistSim1[k] <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])), c(mean(TempNS_EH[range.simu.BR]), mean(PrecNS_EH[range.simu.BR]))))

    geoDistSim1[k] <- rdist.earth(t(as.matrix(c(mean(east_Hem[range.simu.BR,1]), mean(east_Hem[range.simu.BR,2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1),1]), mean(east_Hem[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1),2])))), miles = F)

  res_JanuarySim <- resourcesScarcity_januaryEH[range.simu.BR]
  resourcesSim1[k] <- mean(res_JanuarySim) + mean(resourcesScarcity_julyEH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])
  }

## NB sim, BR obs

  nicheDistSim2 <- vector()
  geoDistSim2 <- vector()
  resourcesSim2 <- vector()
  for(k in 1:99){ #length(range.simNB_EH[[j]])){
    #range.simu.BR <- match(range.simNB_EH[[j]][[k]], ID_EH)
    #range.simu.BR <- range.simu.BR[which(is.na(range.simu.BR)==F)]
    range.simu.NB <- range.simNB_EH[[j]][[k]]

    nicheDistSim2[k] <- dist(rbind(c(mean(TempNS_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])), c(mean(TempNW_EH[range.simu.NB]), mean(PrecNW_EH[range.simu.NB]))))
    
    geoDistSim2[k] <- rdist.earth(t(as.matrix(c(mean(east_Hem[range.simu.NB,1]), mean(east_Hem[range.simu.NB,2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1),1]), mean(east_Hem[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1),2])))), miles = F)

  res_JulySim <- resourcesScarcity_julyEH[range.simu.NB]
  resourcesSim2[k] <- mean(res_JulySim) + mean(resourcesScarcity_januaryEH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])
  }

nicheDistObs <- dist(rbind(c(mean(TempNS_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])), c(mean(TempNW_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]))))

nicheDistResNB <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])), c(mean(TempNS_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]))))

nicheDistResBR <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])), c(mean(TempNS_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]))))

all_niche <- rescale(c(nicheDistObs, nicheDistResNB, nicheDistResBR, nicheDistSim1, nicheDistSim2))


geoDistObs <- rdist.earth(t(as.matrix(c(mean(east_Hem[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1),1]), mean(east_Hem[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1),2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1),1]), mean(east_Hem[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1),2])))), miles = F)

all_geo <- rescale(c(geoDistObs, 0, 0, geoDistSim1, geoDistSim2))


resourcesObs <- mean(resourcesScarcity_julyEH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]) + mean(resourcesScarcity_januaryEH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])

resourcesResNB <- mean(resourcesScarcity_julyEH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]) + mean(resourcesScarcity_januaryEH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])

resourcesResBR <- mean(resourcesScarcity_julyEH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]) + mean(resourcesScarcity_januaryEH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])

all_resources <- rescale(c(resourcesObs, resourcesResNB, resourcesResBR, resourcesSim1, resourcesSim2))

nngg <- all_niche + all_geo
rrgg <- all_resources + all_geo
nnrr <- all_niche + all_resources
nnrrgg <- all_niche + all_resources + all_geo

nicheDistObs_SHEH[j] <- nicheDistObs
geoDistObs_SHEH[j] <- geoDistObs
resourcesObs_SHEH[j] <- resourcesObs

nicheDistObs_rescaled_SHEH[j] <- all_niche[1]
nicheDistResNB_rescaled_SHEH[j] <- all_niche[2]
nicheDistResBR_rescaled_SHEH[j] <- all_niche[3]
geoDistObs_rescaled_SHEH[j] <- all_geo[1]
geoDistResNB_rescaled_SHEH[j] <- 0
geoDistResBR_rescaled_SHEH[j] <- 0
resourcesObs_rescaled_SHEH[j] <- all_resources[1]
resourcesResNB_rescaled_SHEH[j] <- all_resources[2]
resourcesResBR_rescaled_SHEH[j] <- all_resources[3]

rank_nicheDistObs_SHEH[j] <- rank(all_niche)[1]
rank_nicheDistResNB_SHEH[j] <- rank(all_niche)[2]
rank_nicheDistResBR_SHEH[j] <- rank(all_niche)[3]
rank_geoDistObs_SHEH[j] <- rank(all_geo)[1]
rank_geoDistResNB_SHEH[j] <- rank(all_geo)[2]
rank_geoDistResBR_SHEH[j] <- rank(all_geo)[3]
rank_resourcesObs_SHEH[j] <- rank(all_resources)[1]
rank_resourcesResNB_SHEH[j] <- rank(all_resources)[2]
rank_resourcesResBR_SHEH[j] <- rank(all_resources)[3]
rank_nicheDist_geoDist_Obs_SHEH[j] <- rank(nngg)[1]
rank_nicheDist_geoDist_ResNB_SHEH[j] <- rank(nngg)[2]
rank_nicheDist_geoDist_ResBR_SHEH[j] <- rank(nngg)[3]
rank_nicheDist_resources_Obs_SHEH[j] <- rank(nnrr)[1]
rank_nicheDist_resources_ResNB_SHEH[j] <- rank(nnrr)[2]
rank_nicheDist_resources_ResBR_SHEH[j] <- rank(nnrr)[3]
rank_resources_geoDist_Obs_SHEH[j] <- rank(rrgg)[1]
rank_resources_geoDist_ResNB_SHEH[j] <- rank(rrgg)[2]
rank_resources_geoDist_ResBR_SHEH[j] <- rank(rrgg)[3]
rank_nicheDist_resources_geoDist_Obs_SHEH[j] <- rank(nnrrgg)[1]
rank_nicheDist_resources_geoDist_ResNB_SHEH[j] <- rank(nnrrgg)[2]
rank_nicheDist_resources_geoDist_ResBR_SHEH[j] <- rank(nnrrgg)[3]

rank_nicheDist_random_SHEH[j] <- sample(rank(all_niche)[4:201], 1)
rank_geoDist_random_SHEH[j] <- sample(rank(all_geo)[4:201], 1)
rank_resources_random_SHEH[j] <- sample(rank(all_resources)[4:201], 1)
rank_nicheDist_geoDist_random_SHEH[j] <- sample(rank(nngg)[4:201], 1)
rank_nicheDist_resources_random_SHEH[j] <- sample(rank(nnrr)[4:201], 1)
rank_resources_geoDist_random_SHEH[j] <- sample(rank(rrgg)[4:201], 1)
rank_nicheDist_resources_geoDist_random_SHEH[j] <- sample(rank(nnrrgg)[4:201], 1)

print(j)
  }


j=13

resourcesScarcity_januaryEH <- -((0.75*NDVI_nw_EH) - NDVI_ns_EH)
resourcesScarcity_julyEH <- -(NDVI_ns_EH - (0.75*NDVI_nw_EH))

## BR sim, NB obs

  nicheDistSim1 <- vector()
  geoDistSim1 <- vector()
  resourcesSim1 <- vector()
  for(k in 1:99){ #length(range.sim_EH[[j]])){
    #range.simu.BR <- match(range.sim_EH[[j]][[k]], ID_EH)
    #range.simu.BR <- range.simu.BR[which(is.na(range.simu.BR)==F)]
    range.simu.BR <- range.sim_SHEH[[j]][[k]]


	nicheDistSim1[k] <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])), c(mean(TempNS_EH[range.simu.BR]), mean(PrecNS_EH[range.simu.BR]))))

    geoDistSim1[k] <- rdist.earth(t(as.matrix(c(mean(east_Hem[range.simu.BR,1]), mean(east_Hem[range.simu.BR,2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1),1]), mean(east_Hem[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1),2])))), miles = F)

  res_JanuarySim <- resourcesScarcity_januaryEH[range.simu.BR]
  resourcesSim1[k] <- mean(res_JanuarySim) + mean(resourcesScarcity_julyEH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])
  }

## NB sim, BR obs

  nicheDistSim2 <- vector()
  geoDistSim2 <- vector()
  resourcesSim2 <- vector()
  for(k in 1:99){ #length(range.simNB_EH[[j]])){
    #range.simu.BR <- match(range.simNB_EH[[j]][[k]], ID_EH)
    #range.simu.BR <- range.simu.BR[which(is.na(range.simu.BR)==F)]
    range.simu.NB <- range.simNB_EH[[j]][[k]]

    nicheDistSim2[k] <- dist(rbind(c(mean(TempNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])), c(mean(TempNW_EH[range.simu.NB]), mean(PrecNW_EH[range.simu.NB]))))
    
    geoDistSim2[k] <- rdist.earth(t(as.matrix(c(mean(east_Hem[range.simu.NB,1]), mean(east_Hem[range.simu.NB,2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1),1]), mean(east_Hem[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1),2])))), miles = F)

  res_JulySim <- resourcesScarcity_julyEH[range.simu.NB]
  resourcesSim2[k] <- mean(res_JulySim) + mean(resourcesScarcity_januaryEH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])
  }

nicheDistObs <- dist(rbind(c(mean(TempNS_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])), c(mean(TempNW_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]))))

nicheDistResNB <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])), c(mean(TempNS_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]))))

nicheDistResBR <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])), c(mean(TempNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]))))

all_niche <- rescale(c(nicheDistObs, nicheDistResNB, nicheDistResBR, nicheDistSim1, nicheDistSim2))


geoDistObs <- rdist.earth(t(as.matrix(c(mean(east_Hem[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1),1]), mean(east_Hem[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1),2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1),1]), mean(east_Hem[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1),2])))), miles = F)

all_geo <- rescale(c(geoDistObs, 0, 0, geoDistSim1, geoDistSim2))


resourcesObs <- mean(resourcesScarcity_julyEH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]) + mean(resourcesScarcity_januaryEH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])

resourcesResNB <- mean(resourcesScarcity_julyEH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]) + mean(resourcesScarcity_januaryEH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])

resourcesResBR <- mean(resourcesScarcity_julyEH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]) + mean(resourcesScarcity_januaryEH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])

all_resources <- rescale(c(resourcesObs, resourcesResNB, resourcesResBR, resourcesSim1, resourcesSim2))

nngg <- all_niche + all_geo
rrgg <- all_resources + all_geo
nnrr <- all_niche + all_resources
nnrrgg <- all_niche + all_resources + all_geo

nicheDistObs_SHEH[j] <- nicheDistObs
geoDistObs_SHEH[j] <- geoDistObs
resourcesObs_SHEH[j] <- resourcesObs

nicheDistObs_rescaled_SHEH[j] <- all_niche[1]
nicheDistResNB_rescaled_SHEH[j] <- all_niche[2]
nicheDistResBR_rescaled_SHEH[j] <- all_niche[3]
geoDistObs_rescaled_SHEH[j] <- all_geo[1]
geoDistResNB_rescaled_SHEH[j] <- 0
geoDistResBR_rescaled_SHEH[j] <- 0
resourcesObs_rescaled_SHEH[j] <- all_resources[1]
resourcesResNB_rescaled_SHEH[j] <- all_resources[2]
resourcesResBR_rescaled_SHEH[j] <- all_resources[3]

rank_nicheDistObs_SHEH[j] <- rank(all_niche)[1]
rank_nicheDistResNB_SHEH[j] <- rank(all_niche)[2]
rank_nicheDistResBR_SHEH[j] <- rank(all_niche)[3]
rank_geoDistObs_SHEH[j] <- rank(all_geo)[1]
rank_geoDistResNB_SHEH[j] <- rank(all_geo)[2]
rank_geoDistResBR_SHEH[j] <- rank(all_geo)[3]
rank_resourcesObs_SHEH[j] <- rank(all_resources)[1]
rank_resourcesResNB_SHEH[j] <- rank(all_resources)[2]
rank_resourcesResBR_SHEH[j] <- rank(all_resources)[3]
rank_nicheDist_geoDist_Obs_SHEH[j] <- rank(nngg)[1]
rank_nicheDist_geoDist_ResNB_SHEH[j] <- rank(nngg)[2]
rank_nicheDist_geoDist_ResBR_SHEH[j] <- rank(nngg)[3]
rank_nicheDist_resources_Obs_SHEH[j] <- rank(nnrr)[1]
rank_nicheDist_resources_ResNB_SHEH[j] <- rank(nnrr)[2]
rank_nicheDist_resources_ResBR_SHEH[j] <- rank(nnrr)[3]
rank_resources_geoDist_Obs_SHEH[j] <- rank(rrgg)[1]
rank_resources_geoDist_ResNB_SHEH[j] <- rank(rrgg)[2]
rank_resources_geoDist_ResBR_SHEH[j] <- rank(rrgg)[3]
rank_nicheDist_resources_geoDist_Obs_SHEH[j] <- rank(nnrrgg)[1]
rank_nicheDist_resources_geoDist_ResNB_SHEH[j] <- rank(nnrrgg)[2]
rank_nicheDist_resources_geoDist_ResBR_SHEH[j] <- rank(nnrrgg)[3]

rank_nicheDist_random_SHEH[j] <- sample(rank(all_niche)[4:201], 1)
rank_geoDist_random_SHEH[j] <- sample(rank(all_geo)[4:201], 1)
rank_resources_random_SHEH[j] <- sample(rank(all_resources)[4:201], 1)
rank_nicheDist_geoDist_random_SHEH[j] <- sample(rank(nngg)[4:201], 1)
rank_nicheDist_resources_random_SHEH[j] <- sample(rank(nnrr)[4:201], 1)
rank_resources_geoDist_random_SHEH[j] <- sample(rank(rrgg)[4:201], 1)
rank_nicheDist_resources_geoDist_random_SHEH[j] <- sample(rank(nnrrgg)[4:201], 1)




nicheDist <- c(nicheDistObs_EH, nicheDistObs_WH, nicheDistObs_SHEH, nicheDistObs_SHWH)
geoDist <- c(geoDistObs_EH, geoDistObs_WH, geoDistObs_SHEH, geoDistObs_SHWH)
resources <- c(resourcesObs_EH, resourcesObs_WH, resourcesObs_SHEH, resourcesObs_SHWH)

nicheDist_rescaled <- c(nicheDistObs_rescaled_EH, nicheDistObs_rescaled_WH, nicheDistObs_rescaled_SHEH, nicheDistObs_rescaled_SHWH)
geoDist_rescaled <- c(geoDistObs_rescaled_EH, geoDistObs_rescaled_WH, geoDistObs_rescaled_SHEH, geoDistObs_rescaled_SHWH)
resources_rescaled <- c(resourcesObs_rescaled_EH, resourcesObs_rescaled_WH, resourcesObs_rescaled_SHEH, resourcesObs_rescaled_SHWH)
nicheDistResNB_rescaled <- c(nicheDistResNB_rescaled_EH, nicheDistResNB_rescaled_WH, nicheDistResNB_rescaled_SHEH, nicheDistResNB_rescaled_SHWH)
geoDistResNB_rescaled <- c(geoDistResNB_rescaled_EH, geoDistResNB_rescaled_WH, geoDistResNB_rescaled_SHEH, geoDistResNB_rescaled_SHWH)
resourcesResNB_rescaled <- c(resourcesResNB_rescaled_EH, resourcesResNB_rescaled_WH, resourcesResNB_rescaled_SHEH, resourcesResNB_rescaled_SHWH)
nicheDistResBR_rescaled <- c(nicheDistResBR_rescaled_EH, nicheDistResBR_rescaled_WH, nicheDistResBR_rescaled_SHEH, nicheDistResBR_rescaled_SHWH)
geoDistResBR_rescaled <- c(geoDistResBR_rescaled_EH, geoDistResBR_rescaled_WH, geoDistResBR_rescaled_SHEH, geoDistResBR_rescaled_SHWH)
resourcesResBR_rescaled <- c(resourcesResBR_rescaled_EH, resourcesResBR_rescaled_WH, resourcesResBR_rescaled_SHEH, resourcesResBR_rescaled_SHWH)

rank_niche <- c(rank_nicheDistObs_EH, rank_nicheDistObs_WH, rank_nicheDistObs_SHEH, rank_nicheDistObs_SHWH)
rank_geo <- c(rank_geoDistObs_EH, rank_geoDistObs_WH, rank_geoDistObs_SHEH, rank_geoDistObs_SHWH)
rank_resources <- c(rank_resourcesObs_EH, rank_resourcesObs_WH, rank_resourcesObs_SHEH, rank_resourcesObs_SHWH)
rank_niche_geo <- c(rank_nicheDist_geoDist_Obs_EH, rank_nicheDist_geoDist_Obs_WH, rank_nicheDist_geoDist_Obs_SHEH, rank_nicheDist_geoDist_Obs_SHWH) 
rank_niche_resources <- c(rank_nicheDist_resources_Obs_EH, rank_nicheDist_resources_Obs_WH, rank_nicheDist_resources_Obs_SHEH, rank_nicheDist_resources_Obs_SHWH) 
rank_resources_geo <- c(rank_resources_geoDist_Obs_EH, rank_resources_geoDist_Obs_WH, rank_resources_geoDist_Obs_SHEH, rank_resources_geoDist_Obs_SHWH) 
rank_niche_geo_resources <- c(rank_nicheDist_resources_geoDist_Obs_EH, rank_nicheDist_resources_geoDist_Obs_WH, rank_nicheDist_resources_geoDist_Obs_SHEH, rank_nicheDist_resources_geoDist_Obs_SHWH) 

rank_niche <- (rank_niche - 1)/200
rank_geo <- (rank_geo - 1)/200
rank_resources <- (rank_resources - 1)/200
rank_niche_geo <- (rank_niche_geo - 1)/200
rank_niche_resources <- (rank_niche_resources - 1)/200
rank_resources_geo <- (rank_resources_geo - 1)/200
rank_niche_geo_resources <- (rank_niche_geo_resources - 1)/200


rank_niche_resNB <- c(rank_nicheDistResNB_EH, rank_nicheDistResNB_WH, rank_nicheDistResNB_SHEH, rank_nicheDistResNB_SHWH)
rank_geo_resNB <- c(rank_geoDistResNB_EH, rank_geoDistResNB_WH, rank_geoDistResNB_SHEH, rank_geoDistResNB_SHWH)
rank_resources_resNB <- c(rank_resourcesResNB_EH, rank_resourcesResNB_WH, rank_resourcesResNB_SHEH, rank_resourcesResNB_SHWH)
rank_niche_geo_resNB <- c(rank_nicheDist_geoDist_ResNB_EH, rank_nicheDist_geoDist_ResNB_WH, rank_nicheDist_geoDist_ResNB_SHEH, rank_nicheDist_geoDist_ResNB_SHWH) 
rank_niche_resources_resNB <- c(rank_nicheDist_resources_ResNB_EH, rank_nicheDist_resources_ResNB_WH, rank_nicheDist_resources_ResNB_SHEH, rank_nicheDist_resources_ResNB_SHWH) 
rank_resources_geo_resNB <- c(rank_resources_geoDist_ResNB_EH, rank_resources_geoDist_ResNB_WH, rank_resources_geoDist_ResNB_SHEH, rank_resources_geoDist_ResNB_SHWH) 
rank_niche_geo_resources_resNB <- c(rank_nicheDist_resources_geoDist_ResNB_EH, rank_nicheDist_resources_geoDist_ResNB_WH, rank_nicheDist_resources_geoDist_ResNB_SHEH, rank_nicheDist_resources_geoDist_ResNB_SHWH) 
rank_niche_resNB <- (rank_niche_resNB - 1)/200
rank_geo_resNB <- (rank_geo_resNB - 1)/200
rank_resources_resNB <- (rank_resources_resNB - 1)/200
rank_niche_geo_resNB <- (rank_niche_geo_resNB - 1)/200
rank_niche_resources_resNB <- (rank_niche_resources_resNB - 1)/200
rank_resources_geo_resNB <- (rank_resources_geo_resNB - 1)/200
rank_niche_geo_resources_resNB <- (rank_niche_geo_resources_resNB - 1)/200

rank_niche_resBR <- c(rank_nicheDistResBR_EH, rank_nicheDistResBR_WH, rank_nicheDistResBR_SHEH, rank_nicheDistResBR_SHWH)
rank_geo_resBR <- c(rank_geoDistResBR_EH, rank_geoDistResBR_WH, rank_geoDistResBR_SHEH, rank_geoDistResBR_SHWH)
rank_resources_resBR <- c(rank_resourcesResBR_EH, rank_resourcesResBR_WH, rank_resourcesResBR_SHEH, rank_resourcesResBR_SHWH)
rank_niche_geo_resBR <- c(rank_nicheDist_geoDist_ResBR_EH, rank_nicheDist_geoDist_ResBR_WH, rank_nicheDist_geoDist_ResBR_SHEH, rank_nicheDist_geoDist_ResBR_SHWH) 
rank_niche_resources_resBR <- c(rank_nicheDist_resources_ResBR_EH, rank_nicheDist_resources_ResBR_WH, rank_nicheDist_resources_ResBR_SHEH, rank_nicheDist_resources_ResBR_SHWH) 
rank_resources_geo_resBR <- c(rank_resources_geoDist_ResBR_EH, rank_resources_geoDist_ResBR_WH, rank_resources_geoDist_ResBR_SHEH, rank_resources_geoDist_ResBR_SHWH) 
rank_niche_geo_resources_resBR <- c(rank_nicheDist_resources_geoDist_ResBR_EH, rank_nicheDist_resources_geoDist_ResBR_WH, rank_nicheDist_resources_geoDist_ResBR_SHEH, rank_nicheDist_resources_geoDist_ResBR_SHWH) 
rank_niche_resBR <- (rank_niche_resBR - 1)/200
rank_geo_resBR <- (rank_geo_resBR - 1)/200
rank_resources_resBR <- (rank_resources_resBR - 1)/200
rank_niche_geo_resBR <- (rank_niche_geo_resBR - 1)/200
rank_niche_resources_resBR <- (rank_niche_resources_resBR - 1)/200
rank_resources_geo_resBR <- (rank_resources_geo_resBR - 1)/200
rank_niche_geo_resources_resBR <- (rank_niche_geo_resources_resBR - 1)/200

rank_niche_random <- c(rank_nicheDist_random_EH, rank_nicheDist_random_WH, rank_nicheDist_random_SHEH, rank_nicheDist_random_SHWH)
rank_geo_random <- c(rank_geoDist_random_EH, rank_geoDist_random_WH, rank_geoDist_random_SHEH, rank_geoDist_random_SHWH)
rank_resources_random <- c(rank_resources_random_EH, rank_resources_random_WH, rank_resources_random_SHEH, rank_resources_random_SHWH)
rank_niche_geo_random <- c(rank_nicheDist_geoDist_random_EH, rank_nicheDist_geoDist_random_WH, rank_nicheDist_geoDist_random_SHEH, rank_nicheDist_geoDist_random_SHWH) 
rank_niche_resources_random <- c(rank_nicheDist_resources_random_EH, rank_nicheDist_resources_random_WH, rank_nicheDist_resources_random_SHEH, rank_nicheDist_resources_random_SHWH) 
rank_resources_geo_random <- c(rank_resources_geoDist_random_EH, rank_resources_geoDist_random_WH, rank_resources_geoDist_random_SHEH, rank_resources_geoDist_random_SHWH) 
rank_niche_geo_resources_random <- c(rank_nicheDist_resources_geoDist_random_EH, rank_nicheDist_resources_geoDist_random_WH, rank_nicheDist_resources_geoDist_random_SHEH, rank_nicheDist_resources_geoDist_random_SHWH) 



#Export excel table for appendix of publication

species.names <- names(c(selectedSpeciesNH_EH, selectedSpeciesNH_WH, selectedSpeciesSH_EH, selectedSpeciesSH_WH))
hem <- c(rep("EH", length(selectedSpeciesNH_EH)), rep("WH", length(selectedSpeciesNH_WH)), rep("EH", length(selectedSpeciesSH_EH)), rep("WH", length(selectedSpeciesSH_WH)))

xlsTable <- cbind(species.names, hem, geoDist, nicheDist, resources, geoDist_rescaled, nicheDist_rescaled, resources_rescaled, rank_geo, rank_niche, rank_resources, rank_niche_geo, rank_niche_resources, rank_resources_geo, rank_niche_geo_resources)
colnames(xlsTable) <- c("species name", "Longitudinal hemisphere", "Geographical distance", "Niche distance", "Resource scarcity", "Geographical distance rescaled", "Niche distance rescaled", "Resource scarcity rescaled", "Rank geo distance", "Rank niche distance", "Rank resource scarcity", "Rank niche distance + geo distance", "Rank niche distance + resource scarcity", "Rank resource scarcity + geo distance", "Rank niche distance + resource scarcity + geo distance")

xlsTable <- xlsTable[order(species.names),]

write.csv(xlsTable, "/Users/mariussomveille/Desktop/PhD/Chapter 3 – niche/Submission process/Proceedings b/Revisions/Appendix1.csv")


































datasett <- read.csv("env_data_6months.csv")
datasett2 <- datasett[-which(datasett$Temp_NS==0 & datasett$Prec_NS==0 & datasett$Temp_NW==0 & datasett$Prec_NW==0),]


#hexgrid <- readOGR("/Users/mariussomveille/Desktop/Hex grid", "isea3h7_analyses_clean", verbose=FALSE)
#hexgrid <- hexgrid[,c(1,2,15,16)]
#hexgridWH <- hexgrid[which(hexgrid@data[,4] <= -30),]
#hexgridEH <- hexgrid[which(hexgrid@data[,4] > -30),]


hexid <- dataset2[,1]
















#### MIGRATION DISTANCE







geoDistObsWH_ranges <- vector()
geoDistObsEH_ranges <- vector()
geoDistObsSHWH_ranges <- vector()
geoDistObsSHEH_ranges <- vector()
for(j in 1:length(selectedSpeciesNH_WH)){ 
geoDistObsWH_ranges[j] <- as.vector(rdist.earth(t(as.matrix(c(mean(west_Hem[which(PresAbs_NB_NH_WH[,match(names(selectedSpeciesNH_WH)[j], colnames(PresAbs_NB_NH_WH))]==1),1]), mean(west_Hem[which(PresAbs_NB_NH_WH[,match(names(selectedSpeciesNH_WH)[j], colnames(PresAbs_NB_NH_WH))]==1),2])))),t(as.matrix(c(mean(west_Hem[which(PresAbs_BR_NH_WH[,match(names(selectedSpeciesNH_WH)[j], colnames(PresAbs_BR_NH_WH))]==1),1]), mean(west_Hem[which(PresAbs_BR_NH_WH[,match(names(selectedSpeciesNH_WH)[j], colnames(PresAbs_BR_NH_WH))]==1),2])))), miles = F)) 
}
for(j in 1:length(selectedSpeciesNH_EH)){ 
geoDistObsEH_ranges[j] <- as.vector(rdist.earth(t(as.matrix(c(mean(east_Hem[which(PresAbs_NB_NH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_NB_NH_EH))]==1),1]), mean(east_Hem[which(PresAbs_NB_NH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_NB_NH_EH))]==1),2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_BR_NH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_BR_NH_EH))]==1),1]), mean(east_Hem[which(PresAbs_BR_NH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_BR_NH_EH))]==1),2])))), miles = F)) 
}
j=387
geoDistObsEH_ranges[j] <- as.vector(rdist.earth(t(as.matrix(c(mean(east_Hem[which(PresAbs_NB_SH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_NB_SH_EH))]==1),1]), mean(east_Hem[which(PresAbs_NB_SH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_NB_SH_EH))]==1),2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_BR_SH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_BR_SH_EH))]==1),1]), mean(east_Hem[which(PresAbs_BR_SH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_BR_SH_EH))]==1),2])))), miles = F))
j=388
geoDistObsEH_ranges[j] <- as.vector(rdist.earth(t(as.matrix(c(mean(east_Hem[which(PresAbs_NB_SH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_NB_SH_EH))]==1),1]), mean(east_Hem[which(PresAbs_NB_SH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_NB_SH_EH))]==1),2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_BR_SH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_BR_SH_EH))]==1),1]), mean(east_Hem[which(PresAbs_BR_SH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_BR_SH_EH))]==1),2])))), miles = F))
for(j in 1:length(selectedSpeciesSH_WH)){ 
geoDistObsSHWH_ranges[j] <- as.vector(rdist.earth(t(as.matrix(c(mean(west_Hem[which(PresAbs_NB_SH_WH[,match(names(selectedSpeciesSH_WH)[j], colnames(PresAbs_NB_SH_WH))]==1),1]), mean(west_Hem[which(PresAbs_NB_SH_WH[,match(names(selectedSpeciesSH_WH)[j], colnames(PresAbs_NB_SH_WH))]==1),2])))),t(as.matrix(c(mean(west_Hem[which(PresAbs_BR_SH_WH[,match(names(selectedSpeciesSH_WH)[j], colnames(PresAbs_BR_SH_WH))]==1),1]), mean(west_Hem[which(PresAbs_BR_SH_WH[,match(names(selectedSpeciesSH_WH)[j], colnames(PresAbs_BR_SH_WH))]==1),2])))), miles = F)) 
}
for(j in 1:length(selectedSpeciesSH_EH)){ 
geoDistObsSHEH_ranges[j] <- as.vector(rdist.earth(t(as.matrix(c(mean(east_Hem[which(PresAbs_NB_SH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_NB_SH_EH))]==1),1]), mean(east_Hem[which(PresAbs_NB_SH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_NB_SH_EH))]==1),2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_BR_SH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_BR_SH_EH))]==1),1]), mean(east_Hem[which(PresAbs_BR_SH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_BR_SH_EH))]==1),2])))), miles = F)) 
}
j=13
geoDistObsSHEH_ranges[j] <- as.vector(rdist.earth(t(as.matrix(c(mean(east_Hem[which(PresAbs_NB_NH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_NB_NH_EH))]==1),1]), mean(east_Hem[which(PresAbs_NB_NH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_NB_NH_EH))]==1),2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_BR_NH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_BR_NH_EH))]==1),1]), mean(east_Hem[which(PresAbs_BR_NH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_BR_NH_EH))]==1),2])))), miles = F))



selectedSpeciesNH_WH <- selectedSpeciesNH_WH[which(is.na(geoDistObsWH_ranges) == F)]
selectedSpeciesNH_EH <- selectedSpeciesNH_EH[which(is.na(geoDistObsEH_ranges) == F)]
selectedSpeciesSH_WH <- selectedSpeciesSH_WH[which(is.na(geoDistObsSHWH_ranges) == F)]
selectedSpeciesSH_EH <- selectedSpeciesSH_EH[which(is.na(geoDistObsSHEH_ranges) == F)]

geoDistObsWH_ranges <- geoDistObsWH_ranges[which(is.na(geoDistObsWH_ranges) == F)]
geoDistObsEH_ranges <- geoDistObsEH_ranges[which(is.na(geoDistObsEH_ranges) == F)]
geoDistObsSHWH_ranges <- geoDistObsSHWH_ranges[which(is.na(geoDistObsSHWH_ranges) == F)]
geoDistObsSHEH_ranges <- geoDistObsSHEH_ranges[which(is.na(geoDistObsSHEH_ranges) == F)]

migrationDistanceObs <- sqrt(c(geoDistObsWH_ranges, geoDistObsEH_ranges, geoDistObsSHWH_ranges, geoDistObsSHEH_ranges))







#### NICHE DISTANCE & RESOURCE SCARCITY (CHANGE TO GAIN FOR PLOT)


nicheDistObsWH_ranges <- vector()
resourcesObsWH_ranges <- vector()
nicheDistObsWH_resNB_ranges <- vector()
nicheDistObsWH_resBR_ranges <- vector()

for(j in 1:length(selectedSpeciesNH_WH)){

  nicheDistObsWH_ranges[j] <- dist(rbind(c(mean(TempNW_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]), mean(PrecNW_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])), c(mean(TempNS_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]), mean(PrecNS_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]))))
  nicheDistObsWH_resNB_ranges[j] <- dist(rbind(c(mean(TempNW_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]), mean(PrecNW_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])), c(mean(TempNS_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]), mean(PrecNS_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]))))
  nicheDistObsWH_resBR_ranges[j] <- dist(rbind(c(mean(TempNW_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]), mean(PrecNW_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])), c(mean(TempNS_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]), mean(PrecNS_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]))))

res_JulyObs <- resourcesScarcity_julyWH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]
res_JanuaryObs <- resourcesScarcity_januaryWH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)]
  
resourcesObsWH_ranges[j] <- mean(res_JulyObs) + mean(res_JanuaryObs)

}




nicheDistObsEH_ranges <- vector()
resourcesObsEH_ranges <- vector()
nicheDistObsEH_resNB_ranges <- vector()
nicheDistObsEH_resBR_ranges <- vector()

for(j in 1:length(selectedSpeciesNH_EH)){

## observed migratory + if residents

  nicheDistObsEH_ranges[j] <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])), c(mean(TempNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]))))
  nicheDistObsEH_resNB_ranges[j] <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])), c(mean(TempNS_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]))))
  nicheDistObsEH_resBR_ranges[j] <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])), c(mean(TempNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]))))

res_JulyObs <- resourcesScarcity_julyEH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]
res_JanuaryObs <- resourcesScarcity_januaryEH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]
  
resourcesObsEH_ranges[j] <- mean(res_JulyObs) + mean(res_JanuaryObs)

}


for(j in 380:381){

## observed migratory + if residents

  nicheDistObsEH_ranges[j] <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])), c(mean(TempNS_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]))))
  nicheDistObsEH_resNB_ranges[j] <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])), c(mean(TempNS_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]))))
  nicheDistObsEH_resBR_ranges[j] <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])), c(mean(TempNS_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]))))

res_JulyObs <- resourcesScarcity_julyEH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]
res_JanuaryObs <- resourcesScarcity_januaryEH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)]
  
resourcesObsEH_ranges[j] <- mean(res_JulyObs) + mean(res_JanuaryObs)

}



nicheDistObsWHSH_ranges <- vector()
resourcesObsWHSH_ranges <- vector()
nicheDistObsWHSH_resNB_ranges <- vector()
nicheDistObsWHSH_resBR_ranges <- vector()

for(j in 1:length(selectedSpeciesSH_WH)){

  nicheDistObsWHSH_ranges[j] <- dist(rbind(c(mean(TempNS_WH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]), mean(PrecNS_WH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)])), c(mean(TempNW_WH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]), mean(PrecNW_WH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]))))
  nicheDistObsWHSH_resNB_ranges[j] <- dist(rbind(c(mean(TempNW_WH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]), mean(PrecNW_WH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)])), c(mean(TempNS_WH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]), mean(PrecNS_WH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]))))
  nicheDistObsWHSH_resBR_ranges[j] <- dist(rbind(c(mean(TempNW_WH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]), mean(PrecNW_WH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)])), c(mean(TempNS_WH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]), mean(PrecNS_WH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]))))

res_JulyObs <- resourcesScarcity_julyWH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]
res_JanuaryObs <- resourcesScarcity_januaryWH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)]
  
resourcesObsWHSH_ranges[j] <- mean(res_JulyObs) + mean(res_JanuaryObs)

}



nicheDistObsEHSH_ranges <- vector()
resourcesObsEHSH_ranges <- vector()
nicheDistObsEHSH_resNB_ranges <- vector()
nicheDistObsEHSH_resBR_ranges <- vector()

for(j in 1:length(selectedSpeciesSH_EH)){

  nicheDistObsEHSH_ranges[j] <- dist(rbind(c(mean(TempNS_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])), c(mean(TempNW_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]))))
  nicheDistObsEHSH_resNB_ranges[j] <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])), c(mean(TempNS_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]))))
  nicheDistObsEHSH_resBR_ranges[j] <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])), c(mean(TempNS_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]))))

res_JulyObs <- resourcesScarcity_julyEH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]
res_JanuaryObs <- resourcesScarcity_januaryEH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]
  
resourcesObsEHSH_ranges[j] <- mean(res_JulyObs) + mean(res_JanuaryObs)

}

j = 13

  nicheDistObsEHSH_ranges[j] <- dist(rbind(c(mean(TempNS_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])), c(mean(TempNW_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]))))
  nicheDistObsEHSH_resNB_ranges[j] <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])), c(mean(TempNS_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]))))
  nicheDistObsEHSH_resBR_ranges[j] <- dist(rbind(c(mean(TempNW_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNW_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])), c(mean(TempNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]), mean(PrecNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]))))

res_JulyObs <- resourcesScarcity_julyEH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]
res_JanuaryObs <- resourcesScarcity_januaryEH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)]
  
resourcesObsEHSH_ranges[j] <- mean(res_JulyObs) + mean(res_JanuaryObs)



nicheDistObs <- c(nicheDistObsWH_ranges, nicheDistObsEH_ranges, nicheDistObsWHSH_ranges, nicheDistObsEHSH_ranges)
nicheDistObs_resNB <- c(nicheDistObsWH_resNB_ranges, nicheDistObsEH_resNB_ranges, nicheDistObsWHSH_resNB_ranges, nicheDistObsEHSH_resNB_ranges)
nicheDistObs_resBR <- c(nicheDistObsWH_resBR_ranges, nicheDistObsEH_resBR_ranges, nicheDistObsWHSH_resBR_ranges, nicheDistObsEHSH_resBR_ranges)
resourcesObs <- c(resourcesObsWH_ranges, resourcesObsEH_ranges, resourcesObsWHSH_ranges, resourcesObsEHSH_ranges)






### COMPETITION BETWEEN MIGRANTS

PresAbs_BR_NH_WH_compet <- PresAbs_BR_NH_WH[,match(names(selectedSpeciesNH_WH), colnames(PresAbs_BR_NH_WH))]
PresAbs_NB_NH_WH_compet <- PresAbs_NB_NH_WH[,match(names(selectedSpeciesNH_WH), colnames(PresAbs_NB_NH_WH))]
PresAbs_BR_SH_WH_compet <- PresAbs_BR_SH_WH[,match(names(selectedSpeciesSH_WH), colnames(PresAbs_BR_SH_WH))]
PresAbs_NB_SH_WH_compet <- PresAbs_NB_SH_WH[,match(names(selectedSpeciesSH_WH), colnames(PresAbs_NB_SH_WH))]

PresAbs_BR_NH_EH_compet <- PresAbs_BR_NH_EH[,match(names(selectedSpeciesNH_EH), colnames(PresAbs_BR_NH_EH))[-which(is.na(match(names(selectedSpeciesNH_EH), colnames(PresAbs_BR_NH_EH)))==T)]]
PresAbs_NB_NH_EH_compet <- PresAbs_NB_NH_EH[,match(names(selectedSpeciesNH_EH), colnames(PresAbs_NB_NH_EH))[-which(is.na(match(names(selectedSpeciesNH_EH), colnames(PresAbs_NB_NH_EH)))==T)]]
for(j in 380:381){
PresAbs_BR_NH_EH_compet[,j] <- PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]
PresAbs_NB_NH_EH_compet[,j] <- PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]
}
PresAbs_BR_SH_EH_compet <- PresAbs_BR_SH_EH[,match(names(selectedSpeciesSH_EH), colnames(PresAbs_BR_SH_EH))[-which(is.na(match(names(selectedSpeciesSH_EH), colnames(PresAbs_BR_SH_EH)))==T)]]
PresAbs_NB_SH_EH_compet <- PresAbs_NB_SH_EH[,match(names(selectedSpeciesSH_EH), colnames(PresAbs_NB_SH_EH))[-which(is.na(match(names(selectedSpeciesSH_EH), colnames(PresAbs_BR_SH_EH)))==T)]]
PresAbs_BR_SH_EH_compet[,13] <- PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[13]))]
PresAbs_NB_SH_EH_compet[,13] <- PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[13]))]


compet.with.migrants_NH_WH <- vector()
for(i in 1:dim(PresAbs_BR_NH_WH_compet)[2]){
	compet.with.migrants_NH_WH[i] <- (sum(apply(PresAbs_BR_NH_WH_compet[which(PresAbs_BR_NH_WH_compet[,i] == 1),], 1, sum) - 1) / length(which(PresAbs_BR_NH_WH_compet[,i] == 1))) + (sum(apply(PresAbs_NB_NH_WH_compet[which(PresAbs_NB_NH_WH_compet[,i] == 1),], 1, sum) - 1) / length(which(PresAbs_NB_NH_WH_compet[,i] == 1)))
}
compet.with.migrants_SH_WH <- vector()
for(i in 1:dim(PresAbs_BR_SH_WH_compet)[2]){
	compet.with.migrants_SH_WH[i] <- (sum(apply(PresAbs_BR_SH_WH_compet[which(PresAbs_BR_SH_WH_compet[,i] == 1),], 1, sum) - 1) / length(which(PresAbs_BR_SH_WH_compet[,i] == 1))) + (sum(apply(PresAbs_NB_SH_WH_compet[which(PresAbs_NB_SH_WH_compet[,i] == 1),], 1, sum) - 1) / length(which(PresAbs_NB_SH_WH_compet[,i] == 1)))
}
compet.with.migrants_NH_EH <- vector()
for(i in 1:dim(PresAbs_BR_NH_EH_compet)[2]){
	compet.with.migrants_NH_EH[i] <- (sum(apply(PresAbs_BR_NH_EH_compet[which(PresAbs_BR_NH_EH_compet[,i] == 1),], 1, sum) - 1) / length(which(PresAbs_BR_NH_EH_compet[,i] == 1))) + (sum(apply(PresAbs_NB_NH_EH_compet[which(PresAbs_NB_NH_EH_compet[,i] == 1),], 1, sum) - 1) / length(which(PresAbs_NB_NH_EH_compet[,i] == 1)))
}
compet.with.migrants_SH_EH <- vector()
for(i in 1:dim(PresAbs_BR_SH_EH_compet)[2]){
	compet.with.migrants_SH_EH[i] <- (sum(apply(PresAbs_BR_SH_EH_compet[which(PresAbs_BR_SH_EH_compet[,i] == 1),], 1, sum) - 1) / length(which(PresAbs_BR_SH_EH_compet[,i] == 1))) + (sum(apply(PresAbs_NB_SH_EH_compet[which(PresAbs_NB_SH_EH_compet[,i] == 1),], 1, sum) - 1) / length(which(PresAbs_NB_SH_EH_compet[,i] == 1)))
}

compet.with.migrants <- c(compet.with.migrants_NH_WH, compet.with.migrants_NH_EH, compet.with.migrants_SH_WH, compet.with.migrants_SH_EH)




## New measure of competition using species ranges

rescale <- function(x) (x-min(x,na.rm=T))/(max(x,na.rm=T) - min(x,na.rm=T)) * 100


competMigr_julyWH <- apply(PresAbs_BR_NH_WH_compet, 1, sum) + 1
compet_july_WH <- vector()
for(i in 1:length(resourcesGain_julyWH)){
	compet_july_WH[i] <- ifelse(competMigr_julyWH[i]>0, resourcesGain_julyWH[i] / competMigr_julyWH[i], 0)
}
competMigr_januaryWH <- apply(PresAbs_NB_NH_WH_compet, 1, sum) + 1
compet_january_WH <- vector()
for(i in 1:length(resourcesGain_januaryWH)){
	compet_january_WH[i] <- ifelse(competMigr_januaryWH[i]>0, resourcesGain_januaryWH[i] / competMigr_januaryWH[i], 0)
}

competMigr_julyEH <- apply(PresAbs_BR_NH_EH_compet, 1, sum) + 1
compet_july_EH <- vector()
for(i in 1:length(resourcesGain_julyEH)){
	compet_july_EH[i] <- ifelse(competMigr_julyEH[i]>0, resourcesGain_julyEH[i] / competMigr_julyEH[i], 0)
}
competMigr_januaryEH <- apply(PresAbs_NB_NH_EH_compet, 1, sum) + 1 
compet_january_EH <- vector()
for(i in 1:length(resourcesGain_januaryEH)){
	compet_january_EH[i] <- ifelse(competMigr_januaryEH[i]>0, resourcesGain_januaryEH[i] / competMigr_januaryEH[i], 0)
}




compet_NH_WH <- vector()
for(i in 1:dim(PresAbs_BR_NH_WH_compet)[2]){
	compet_NH_WH[i] <- mean(compet_july_WH[which(PresAbs_BR_NH_WH_compet[,i] == 1)]) + mean(compet_january_WH[which(PresAbs_NB_NH_WH_compet[,i] == 1)])
}

compet_NH_EH <- vector()
for(i in 1:dim(PresAbs_BR_NH_EH_compet)[2]){
	compet_NH_EH[i] <- mean(compet_july_EH[which(PresAbs_BR_NH_EH_compet[,i] == 1)]) + mean(compet_january_EH[which(PresAbs_NB_NH_EH_compet[,i] == 1)])
}


rescale(c(geoDistObsWH_ranges, geoDistObsEH_ranges))

plot(sqrt(c(geoDistObsWH_ranges, geoDistObsEH_ranges)), c(compet_NH_WH, compet_NH_EH), pch=20)


compet_NH_WH <- vector()
for(i in 1:dim(PresAbs_BR_NH_WH_compet)[2]){
	compet_NH_WH[i] <- (sum(resourcesGain_julyWH[which(PresAbs_BR_NH_WH_compet[,i] == 1)]) / sum(apply(PresAbs_BR_NH_WH_compet[which(PresAbs_BR_NH_WH_compet[,i] == 1),], 1, sum) - 1)) + (sum(resourcesGain_januaryWH[which(PresAbs_NB_NH_WH_compet[,i] == 1)]) / sum(apply(PresAbs_NB_NH_WH_compet[which(PresAbs_NB_NH_WH_compet[,i] == 1),], 1, sum) - 1))
}



compet_NH_WH_july <- vector()
compet_NH_WH_january <- vector()
for(i in 1:length(selectedSpeciesNH_WH)){
	compet_NH_WH_july[i] <- sum(NDVI_ns_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[i]))]==1)]) / sum(apply(PresAbs_BR_NH_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[i]))] == 1),], 1, sum) - 1)
	compet_NH_WH_january[i] <- sum(NDVI_nw_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[i]))]==1)]) / sum(apply(PresAbs_NB_NH_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[i]))] == 1),], 1, sum) - 1)
}
compet_NH_WH <- compet_NH_WH_january + compet_NH_WH_july
	
	 / length(which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[i]))] == 1))) 
	
	+ (sum(apply(PresAbs_NB_NH_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[i]))] == 1),], 1, sum) - 1) / length(which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[i]))] == 1)))
}


# Figure 3 color-coded by competition between migrants instead of niche tracking
par(mar=c(2.5,2.5,0,2.5), mgp=c(1.5,0.5,0))
hist(migrationDistanceObs, xlim=c(0,120), xlab="", ylab="", main="", xaxt="n", axes=F, col="light grey", border="grey")
axis(side=4)
axis(side=1, at=c(0,20,40,60,80,100,120), labels = as.character(c(0,20,40,60,80,100,120)^2))
par(new=T, mar=c(2.5,2.5,0,2.5), mgp=c(1.5,0.5,0))
plot(migrationDistanceObs, resourcesObs, pch=20, ylab="Resource gain", xlab="Geographic distance (Km)", xlim=c(0,120), col="yellow", cex=1.2, axes=F)
points(migrationDistanceObs[which(compet.with.migrants > quantile(compet.with.migrants)[2])], resourcesObs[which(compet.with.migrants > quantile(compet.with.migrants)[2])], pch=20, xlim=c(0,120), col="orange", cex=1.2)
points(migrationDistanceObs[which(compet.with.migrants > quantile(compet.with.migrants)[3])], resourcesObs[which(compet.with.migrants > quantile(compet.with.migrants)[3])], pch=20, xlim=c(0,120), col="red", cex=1.2)
points(migrationDistanceObs[which(compet.with.migrants > quantile(compet.with.migrants)[4])], resourcesObs[which(compet.with.migrants > quantile(compet.with.migrants)[4])], pch=20, xlim=c(0,120), col="brown4", cex=1.2)
axis(side=2)
mtext("Number of species", side=4, line=1.5, cex.lab=1,las=0)
legend("topleft", inset=.05, title="Competition\nwith migrants", c(paste(">",round(quantile(compet.with.migrants)[4])), paste(round(quantile(compet.with.migrants)[3]),"–",round(quantile(compet.with.migrants)[4]), sep=""), paste(round(quantile(compet.with.migrants)[2]),"–",round(quantile(compet.with.migrants)[3]), sep=""), paste("<",round(quantile(compet.with.migrants)[2]))), fill=c("brown4", "red", "orange", "yellow"), bty="n")

summary(lm(resourcesObs ~ migrationDistanceObs*compet.with.migrants))





### EFFECT OF RANGE SIZE?

range.sizes_BR_NH_WH <- apply(PresAbs_BR_NH_WH_compet, 2, sum)
range.sizes_NB_NH_WH <- apply(PresAbs_NB_NH_WH_compet, 2, sum)
range.sizes_BR_SH_WH <- apply(PresAbs_BR_SH_WH_compet, 2, sum)
range.sizes_NB_SH_WH <- apply(PresAbs_NB_SH_WH_compet, 2, sum)
range.sizes_BR_NH_EH <- apply(PresAbs_BR_NH_EH_compet, 2, sum)
range.sizes_NB_NH_EH <- apply(PresAbs_NB_NH_EH_compet, 2, sum)
range.sizes_BR_SH_EH <- apply(PresAbs_BR_SH_EH_compet, 2, sum)
range.sizes_NB_SH_EH <- apply(PresAbs_NB_SH_EH_compet, 2, sum)

range.sizes_BR <- c(range.sizes_BR_NH_WH, range.sizes_BR_NH_EH, range.sizes_BR_SH_WH, range.sizes_BR_SH_EH)
range.sizes_NB <- c(range.sizes_NB_NH_WH, range.sizes_NB_NH_EH, range.sizes_NB_SH_WH, range.sizes_NB_SH_EH)






###  Earth mover distance in niche space ###

# create a raster grid across niche space (extent in temp: -4, 2; and in prec: -3,3)

# compute density kernels for breeding and non-breeding niches and fill the cells of grid with probability values

nicheDistObs_NHWH <- vector()
nicheDistResNB_NHWH <- vector()
nicheDistResBR_NHWH <- vector()

for(j in 1:length(selectedSpeciesNH_WH)){
	
	## breeding niche ##
	br <- cbind(TempNS_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)], PrecNS_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])
	br_kernel <- kde2d(br[,1], br[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
	# Convert the kernel into a raster
	br_kernel_2_mat <- rep(br_kernel[[1]][1],50)
	for(i in 2:50){
		br_kernel_2_mat <- c(br_kernel_2_mat, rep(br_kernel[[1]][i],50))
	}
	br_kernel_2_mat2 <- rep(br_kernel[[2]],50)
	br_kernel_2_mat3 <- br_kernel[[3]][1,]
	for(i in 2:50){
		br_kernel_2_mat3 <- c(br_kernel_2_mat3, br_kernel[[3]][i,])
	}
	br_kernel_mat <- as.data.frame(cbind(br_kernel_2_mat, br_kernel_2_mat2, br_kernel_2_mat3))
	breeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(breeding.niche.raster) <- extent(c(-3, 3, -3,3))
	breeding.niche.raster <- rasterize(br_kernel_mat[,1:2], breeding.niche.raster, br_kernel_mat[,3])
	breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))
	# Keep only the top 99% of the kernel, set the rest to 0
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(breeding.niche.raster), decreasing=T)[i]
	}
	breeding.niche.raster[which(as.vector(breeding.niche.raster) < sort(as.vector(breeding.niche.raster), decreasing=T)[i])] = 0
	breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))

	## non-breeding niche ##
	nb <- cbind(TempNW_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)], PrecNW_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])
	nb_kernel <- kde2d(nb[,1], nb[,2], h=0.5, n=50, lims=c(-3, 3, -3,3))
	nb_kernel_2_mat <- rep(nb_kernel[[1]][1],50)
	for(i in 2:50){
		nb_kernel_2_mat <- c(nb_kernel_2_mat, rep(nb_kernel[[1]][i],50))
	}
	nb_kernel_2_mat2 <- rep(nb_kernel[[2]],50)
	nb_kernel_2_mat3 <- nb_kernel[[3]][1,]
	for(i in 2:50){
		nb_kernel_2_mat3 <- c(nb_kernel_2_mat3, nb_kernel[[3]][i,])
	}
	nb_kernel_mat <- as.data.frame(cbind(nb_kernel_2_mat, nb_kernel_2_mat2, nb_kernel_2_mat3))
	nonbreeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(nonbreeding.niche.raster) <- extent(c(-3, 3, -3,3))
	nonbreeding.niche.raster <- rasterize(nb_kernel_mat[,1:2], nonbreeding.niche.raster, nb_kernel_mat[,3])
	nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i]
	}
	nonbreeding.niche.raster[which(as.vector(nonbreeding.niche.raster) < sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i])] = 0
	nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))

	# Compute observed distance
	nicheDistObs_NHWH[j] <- emd(breeding.niche.raster, nonbreeding.niche.raster, threshold=2)
	
	
	## stay resident on non-breeding ground ##
	resnb <- cbind(TempNS_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)], PrecNS_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])
	resnb_kernel <- kde2d(resnb[,1], resnb[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
	# Convert the kernel into a raster
	resnb_kernel_2_mat <- rep(resnb_kernel[[1]][1],50)
	for(i in 2:50){
		resnb_kernel_2_mat <- c(resnb_kernel_2_mat, rep(resnb_kernel[[1]][i],50))
	}
	resnb_kernel_2_mat2 <- rep(resnb_kernel[[2]],50)
	resnb_kernel_2_mat3 <- resnb_kernel[[3]][1,]
	for(i in 2:50){
		resnb_kernel_2_mat3 <- c(resnb_kernel_2_mat3, resnb_kernel[[3]][i,])
	}
	resnb_kernel_mat <- as.data.frame(cbind(resnb_kernel_2_mat, resnb_kernel_2_mat2, resnb_kernel_2_mat3))
	resident.nonbreeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(resident.nonbreeding.niche.raster) <- extent(c(-3, 3, -3,3))
	resident.nonbreeding.niche.raster <- rasterize(resnb_kernel_mat[,1:2], resident.nonbreeding.niche.raster, resnb_kernel_mat[,3])
	resident.nonbreeding.niche.raster = resident.nonbreeding.niche.raster / sum(as.vector(resident.nonbreeding.niche.raster))
	# Keep only the top 99% of the kernel, set the rest to 0
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(resident.nonbreeding.niche.raster), decreasing=T)[i]
	}
	resident.nonbreeding.niche.raster[which(as.vector(resident.nonbreeding.niche.raster) < sort(as.vector(resident.nonbreeding.niche.raster), decreasing=T)[i])] = 0
	resident.nonbreeding.niche.raster = resident.nonbreeding.niche.raster / sum(as.vector(resident.nonbreeding.niche.raster))
	
	
	## stay resident on breeding ground ##
	resbr <- cbind(TempNW_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)], PrecNW_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[j]))]==1)])
	resbr_kernel <- kde2d(resbr[,1], resbr[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
	# Convert the kernel into a raster
	resbr_kernel_2_mat <- rep(resbr_kernel[[1]][1],50)
	for(i in 2:50){
		resbr_kernel_2_mat <- c(resbr_kernel_2_mat, rep(resbr_kernel[[1]][i],50))
	}
	resbr_kernel_2_mat2 <- rep(resbr_kernel[[2]],50)
	resbr_kernel_2_mat3 <- resbr_kernel[[3]][1,]
	for(i in 2:50){
		resbr_kernel_2_mat3 <- c(resbr_kernel_2_mat3, resbr_kernel[[3]][i,])
	}
	resbr_kernel_mat <- as.data.frame(cbind(resbr_kernel_2_mat, resbr_kernel_2_mat2, resbr_kernel_2_mat3))
	resident.breeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(resident.breeding.niche.raster) <- extent(c(-3, 3, -3,3))
	resident.breeding.niche.raster <- rasterize(resbr_kernel_mat[,1:2], resident.breeding.niche.raster, resbr_kernel_mat[,3])
	resident.breeding.niche.raster = resident.breeding.niche.raster / sum(as.vector(resident.breeding.niche.raster))
	# Keep only the top 99% of the kernel, set the rest to 0
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(resident.breeding.niche.raster), decreasing=T)[i]
	}
	resident.breeding.niche.raster[which(as.vector(resident.breeding.niche.raster) < sort(as.vector(resident.breeding.niche.raster), decreasing=T)[i])] = 0
	resident.breeding.niche.raster = resident.breeding.niche.raster / sum(as.vector(resident.breeding.niche.raster))
	
	# Compute distances if resident
	nicheDistResNB_NHWH[j] <- emd(resident.nonbreeding.niche.raster, nonbreeding.niche.raster, threshold=2)
	nicheDistResBR_NHWH[j] <- emd(breeding.niche.raster, resident.breeding.niche.raster, threshold=2)

}



nicheDistObs_SHWH <- vector()
nicheDistResNB_SHWH <- vector()
nicheDistResBR_SHWH <- vector()

for(j in 1:length(selectedSpeciesSH_WH)){
	
	## breeding niche ##
	br <- cbind(TempNW_WH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)], PrecNW_WH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)])
	br_kernel <- kde2d(br[,1], br[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
	# Convert the kernel into a raster
	br_kernel_2_mat <- rep(br_kernel[[1]][1],50)
	for(i in 2:50){
		br_kernel_2_mat <- c(br_kernel_2_mat, rep(br_kernel[[1]][i],50))
	}
	br_kernel_2_mat2 <- rep(br_kernel[[2]],50)
	br_kernel_2_mat3 <- br_kernel[[3]][1,]
	for(i in 2:50){
		br_kernel_2_mat3 <- c(br_kernel_2_mat3, br_kernel[[3]][i,])
	}
	br_kernel_mat <- as.data.frame(cbind(br_kernel_2_mat, br_kernel_2_mat2, br_kernel_2_mat3))
	breeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(breeding.niche.raster) <- extent(c(-3, 3, -3,3))
	breeding.niche.raster <- rasterize(br_kernel_mat[,1:2], breeding.niche.raster, br_kernel_mat[,3])
	breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))
	# Keep only the top 99% of the kernel, set the rest to 0
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(breeding.niche.raster), decreasing=T)[i]
	}
	breeding.niche.raster[which(as.vector(breeding.niche.raster) < sort(as.vector(breeding.niche.raster), decreasing=T)[i])] = 0
	breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))

	## non-breeding niche ##
	nb <- cbind(TempNS_WH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)], PrecNS_WH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)])
	nb_kernel <- kde2d(nb[,1], nb[,2], h=0.5, n=50, lims=c(-3, 3, -3,3))
	nb_kernel_2_mat <- rep(nb_kernel[[1]][1],50)
	for(i in 2:50){
		nb_kernel_2_mat <- c(nb_kernel_2_mat, rep(nb_kernel[[1]][i],50))
	}
	nb_kernel_2_mat2 <- rep(nb_kernel[[2]],50)
	nb_kernel_2_mat3 <- nb_kernel[[3]][1,]
	for(i in 2:50){
		nb_kernel_2_mat3 <- c(nb_kernel_2_mat3, nb_kernel[[3]][i,])
	}
	nb_kernel_mat <- as.data.frame(cbind(nb_kernel_2_mat, nb_kernel_2_mat2, nb_kernel_2_mat3))
	nonbreeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(nonbreeding.niche.raster) <- extent(c(-3, 3, -3,3))
	nonbreeding.niche.raster <- rasterize(nb_kernel_mat[,1:2], nonbreeding.niche.raster, nb_kernel_mat[,3])
	nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i]
	}
	nonbreeding.niche.raster[which(as.vector(nonbreeding.niche.raster) < sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i])] = 0
	nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))

	# Compute observed distance
	nicheDistObs_SHWH[j] <- emd(breeding.niche.raster, nonbreeding.niche.raster, threshold=2)
	
	
	## stay resident on non-breeding ground ##
	resnb <- cbind(TempNW_WH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)], PrecNW_WH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)])
	resnb_kernel <- kde2d(resnb[,1], resnb[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
	# Convert the kernel into a raster
	resnb_kernel_2_mat <- rep(resnb_kernel[[1]][1],50)
	for(i in 2:50){
		resnb_kernel_2_mat <- c(resnb_kernel_2_mat, rep(resnb_kernel[[1]][i],50))
	}
	resnb_kernel_2_mat2 <- rep(resnb_kernel[[2]],50)
	resnb_kernel_2_mat3 <- resnb_kernel[[3]][1,]
	for(i in 2:50){
		resnb_kernel_2_mat3 <- c(resnb_kernel_2_mat3, resnb_kernel[[3]][i,])
	}
	resnb_kernel_mat <- as.data.frame(cbind(resnb_kernel_2_mat, resnb_kernel_2_mat2, resnb_kernel_2_mat3))
	resident.nonbreeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(resident.nonbreeding.niche.raster) <- extent(c(-3, 3, -3,3))
	resident.nonbreeding.niche.raster <- rasterize(resnb_kernel_mat[,1:2], resident.nonbreeding.niche.raster, resnb_kernel_mat[,3])
	resident.nonbreeding.niche.raster = resident.nonbreeding.niche.raster / sum(as.vector(resident.nonbreeding.niche.raster))
	# Keep only the top 99% of the kernel, set the rest to 0
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(resident.nonbreeding.niche.raster), decreasing=T)[i]
	}
	resident.nonbreeding.niche.raster[which(as.vector(resident.nonbreeding.niche.raster) < sort(as.vector(resident.nonbreeding.niche.raster), decreasing=T)[i])] = 0
	resident.nonbreeding.niche.raster = resident.nonbreeding.niche.raster / sum(as.vector(resident.nonbreeding.niche.raster))
	
	
	## stay resident on breeding ground ##
	resbr <- cbind(TempNS_WH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)], PrecNS_WH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[j]))]==1)])
	resbr_kernel <- kde2d(resbr[,1], resbr[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
	# Convert the kernel into a raster
	resbr_kernel_2_mat <- rep(resbr_kernel[[1]][1],50)
	for(i in 2:50){
		resbr_kernel_2_mat <- c(resbr_kernel_2_mat, rep(resbr_kernel[[1]][i],50))
	}
	resbr_kernel_2_mat2 <- rep(resbr_kernel[[2]],50)
	resbr_kernel_2_mat3 <- resbr_kernel[[3]][1,]
	for(i in 2:50){
		resbr_kernel_2_mat3 <- c(resbr_kernel_2_mat3, resbr_kernel[[3]][i,])
	}
	resbr_kernel_mat <- as.data.frame(cbind(resbr_kernel_2_mat, resbr_kernel_2_mat2, resbr_kernel_2_mat3))
	resident.breeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(resident.breeding.niche.raster) <- extent(c(-3, 3, -3,3))
	resident.breeding.niche.raster <- rasterize(resbr_kernel_mat[,1:2], resident.breeding.niche.raster, resbr_kernel_mat[,3])
	resident.breeding.niche.raster = resident.breeding.niche.raster / sum(as.vector(resident.breeding.niche.raster))
	# Keep only the top 99% of the kernel, set the rest to 0
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(resident.breeding.niche.raster), decreasing=T)[i]
	}
	resident.breeding.niche.raster[which(as.vector(resident.breeding.niche.raster) < sort(as.vector(resident.breeding.niche.raster), decreasing=T)[i])] = 0
	resident.breeding.niche.raster = resident.breeding.niche.raster / sum(as.vector(resident.breeding.niche.raster))
	
	# Compute distances if resident
	nicheDistResNB_SHWH[j] <- emd(resident.nonbreeding.niche.raster, nonbreeding.niche.raster, threshold=2)
	nicheDistResBR_SHWH[j] <- emd(breeding.niche.raster, resident.breeding.niche.raster, threshold=2)

}





nicheDistObs_NHEH <- vector()
nicheDistResNB_NHEH <- vector()
nicheDistResBR_NHEH <- vector()

for(j in 1:length(selectedSpeciesNH_EH)){
	
	## breeding niche ##
	br <- cbind(TempNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)], PrecNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])
	br_kernel <- kde2d(br[,1], br[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
	# Convert the kernel into a raster
	br_kernel_2_mat <- rep(br_kernel[[1]][1],50)
	for(i in 2:50){
		br_kernel_2_mat <- c(br_kernel_2_mat, rep(br_kernel[[1]][i],50))
	}
	br_kernel_2_mat2 <- rep(br_kernel[[2]],50)
	br_kernel_2_mat3 <- br_kernel[[3]][1,]
	for(i in 2:50){
		br_kernel_2_mat3 <- c(br_kernel_2_mat3, br_kernel[[3]][i,])
	}
	br_kernel_mat <- as.data.frame(cbind(br_kernel_2_mat, br_kernel_2_mat2, br_kernel_2_mat3))
	breeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(breeding.niche.raster) <- extent(c(-3, 3, -3,3))
	breeding.niche.raster <- rasterize(br_kernel_mat[,1:2], breeding.niche.raster, br_kernel_mat[,3])
	breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))
	# Keep only the top 99% of the kernel, set the rest to 0
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(breeding.niche.raster), decreasing=T)[i]
	}
	breeding.niche.raster[which(as.vector(breeding.niche.raster) < sort(as.vector(breeding.niche.raster), decreasing=T)[i])] = 0
	breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))

	## non-breeding niche ##
	nb <- cbind(TempNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)], PrecNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])
	nb_kernel <- kde2d(nb[,1], nb[,2], h=0.5, n=50, lims=c(-3, 3, -3,3))
	nb_kernel_2_mat <- rep(nb_kernel[[1]][1],50)
	for(i in 2:50){
		nb_kernel_2_mat <- c(nb_kernel_2_mat, rep(nb_kernel[[1]][i],50))
	}
	nb_kernel_2_mat2 <- rep(nb_kernel[[2]],50)
	nb_kernel_2_mat3 <- nb_kernel[[3]][1,]
	for(i in 2:50){
		nb_kernel_2_mat3 <- c(nb_kernel_2_mat3, nb_kernel[[3]][i,])
	}
	nb_kernel_mat <- as.data.frame(cbind(nb_kernel_2_mat, nb_kernel_2_mat2, nb_kernel_2_mat3))
	nonbreeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(nonbreeding.niche.raster) <- extent(c(-3, 3, -3,3))
	nonbreeding.niche.raster <- rasterize(nb_kernel_mat[,1:2], nonbreeding.niche.raster, nb_kernel_mat[,3])
	nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i]
	}
	nonbreeding.niche.raster[which(as.vector(nonbreeding.niche.raster) < sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i])] = 0
	nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))

	# Compute observed distance
	nicheDistObs_NHEH[j] <- emd(breeding.niche.raster, nonbreeding.niche.raster, threshold=2)
	
	
	## stay resident on non-breeding ground ##
	resnb <- cbind(TempNS_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)], PrecNS_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])
	resnb_kernel <- kde2d(resnb[,1], resnb[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
	# Convert the kernel into a raster
	resnb_kernel_2_mat <- rep(resnb_kernel[[1]][1],50)
	for(i in 2:50){
		resnb_kernel_2_mat <- c(resnb_kernel_2_mat, rep(resnb_kernel[[1]][i],50))
	}
	resnb_kernel_2_mat2 <- rep(resnb_kernel[[2]],50)
	resnb_kernel_2_mat3 <- resnb_kernel[[3]][1,]
	for(i in 2:50){
		resnb_kernel_2_mat3 <- c(resnb_kernel_2_mat3, resnb_kernel[[3]][i,])
	}
	resnb_kernel_mat <- as.data.frame(cbind(resnb_kernel_2_mat, resnb_kernel_2_mat2, resnb_kernel_2_mat3))
	resident.nonbreeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(resident.nonbreeding.niche.raster) <- extent(c(-3, 3, -3,3))
	resident.nonbreeding.niche.raster <- rasterize(resnb_kernel_mat[,1:2], resident.nonbreeding.niche.raster, resnb_kernel_mat[,3])
	resident.nonbreeding.niche.raster = resident.nonbreeding.niche.raster / sum(as.vector(resident.nonbreeding.niche.raster))
	# Keep only the top 99% of the kernel, set the rest to 0
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(resident.nonbreeding.niche.raster), decreasing=T)[i]
	}
	resident.nonbreeding.niche.raster[which(as.vector(resident.nonbreeding.niche.raster) < sort(as.vector(resident.nonbreeding.niche.raster), decreasing=T)[i])] = 0
	resident.nonbreeding.niche.raster = resident.nonbreeding.niche.raster / sum(as.vector(resident.nonbreeding.niche.raster))
	
	
	## stay resident on breeding ground ##
	resbr <- cbind(TempNW_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)], PrecNW_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])
	resbr_kernel <- kde2d(resbr[,1], resbr[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
	# Convert the kernel into a raster
	resbr_kernel_2_mat <- rep(resbr_kernel[[1]][1],50)
	for(i in 2:50){
		resbr_kernel_2_mat <- c(resbr_kernel_2_mat, rep(resbr_kernel[[1]][i],50))
	}
	resbr_kernel_2_mat2 <- rep(resbr_kernel[[2]],50)
	resbr_kernel_2_mat3 <- resbr_kernel[[3]][1,]
	for(i in 2:50){
		resbr_kernel_2_mat3 <- c(resbr_kernel_2_mat3, resbr_kernel[[3]][i,])
	}
	resbr_kernel_mat <- as.data.frame(cbind(resbr_kernel_2_mat, resbr_kernel_2_mat2, resbr_kernel_2_mat3))
	resident.breeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(resident.breeding.niche.raster) <- extent(c(-3, 3, -3,3))
	resident.breeding.niche.raster <- rasterize(resbr_kernel_mat[,1:2], resident.breeding.niche.raster, resbr_kernel_mat[,3])
	resident.breeding.niche.raster = resident.breeding.niche.raster / sum(as.vector(resident.breeding.niche.raster))
	# Keep only the top 99% of the kernel, set the rest to 0
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(resident.breeding.niche.raster), decreasing=T)[i]
	}
	resident.breeding.niche.raster[which(as.vector(resident.breeding.niche.raster) < sort(as.vector(resident.breeding.niche.raster), decreasing=T)[i])] = 0
	resident.breeding.niche.raster = resident.breeding.niche.raster / sum(as.vector(resident.breeding.niche.raster))
	
	# Compute distances if resident
	nicheDistResNB_NHEH[j] <- emd(resident.nonbreeding.niche.raster, nonbreeding.niche.raster, threshold=2)
	nicheDistResBR_NHEH[j] <- emd(breeding.niche.raster, resident.breeding.niche.raster, threshold=2)

}

for(j in 380:381){
	
	## breeding niche ##
	br <- cbind(TempNS_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)], PrecNS_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])
	br_kernel <- kde2d(br[,1], br[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
	# Convert the kernel into a raster
	br_kernel_2_mat <- rep(br_kernel[[1]][1],50)
	for(i in 2:50){
		br_kernel_2_mat <- c(br_kernel_2_mat, rep(br_kernel[[1]][i],50))
	}
	br_kernel_2_mat2 <- rep(br_kernel[[2]],50)
	br_kernel_2_mat3 <- br_kernel[[3]][1,]
	for(i in 2:50){
		br_kernel_2_mat3 <- c(br_kernel_2_mat3, br_kernel[[3]][i,])
	}
	br_kernel_mat <- as.data.frame(cbind(br_kernel_2_mat, br_kernel_2_mat2, br_kernel_2_mat3))
	breeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(breeding.niche.raster) <- extent(c(-3, 3, -3,3))
	breeding.niche.raster <- rasterize(br_kernel_mat[,1:2], breeding.niche.raster, br_kernel_mat[,3])
	breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))
	# Keep only the top 99% of the kernel, set the rest to 0
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(breeding.niche.raster), decreasing=T)[i]
	}
	breeding.niche.raster[which(as.vector(breeding.niche.raster) < sort(as.vector(breeding.niche.raster), decreasing=T)[i])] = 0
	breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))

	## non-breeding niche ##
	nb <- cbind(TempNW_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)], PrecNW_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])
	nb_kernel <- kde2d(nb[,1], nb[,2], h=0.5, n=50, lims=c(-3, 3, -3,3))
	nb_kernel_2_mat <- rep(nb_kernel[[1]][1],50)
	for(i in 2:50){
		nb_kernel_2_mat <- c(nb_kernel_2_mat, rep(nb_kernel[[1]][i],50))
	}
	nb_kernel_2_mat2 <- rep(nb_kernel[[2]],50)
	nb_kernel_2_mat3 <- nb_kernel[[3]][1,]
	for(i in 2:50){
		nb_kernel_2_mat3 <- c(nb_kernel_2_mat3, nb_kernel[[3]][i,])
	}
	nb_kernel_mat <- as.data.frame(cbind(nb_kernel_2_mat, nb_kernel_2_mat2, nb_kernel_2_mat3))
	nonbreeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(nonbreeding.niche.raster) <- extent(c(-3, 3, -3,3))
	nonbreeding.niche.raster <- rasterize(nb_kernel_mat[,1:2], nonbreeding.niche.raster, nb_kernel_mat[,3])
	nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i]
	}
	nonbreeding.niche.raster[which(as.vector(nonbreeding.niche.raster) < sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i])] = 0
	nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))

	# Compute observed distance
	nicheDistObs_NHEH[j] <- emd(breeding.niche.raster, nonbreeding.niche.raster, threshold=2)
	
	
	## stay resident on non-breeding ground ##
	resnb <- cbind(TempNS_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)], PrecNS_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])
	resnb_kernel <- kde2d(resnb[,1], resnb[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
	# Convert the kernel into a raster
	resnb_kernel_2_mat <- rep(resnb_kernel[[1]][1],50)
	for(i in 2:50){
		resnb_kernel_2_mat <- c(resnb_kernel_2_mat, rep(resnb_kernel[[1]][i],50))
	}
	resnb_kernel_2_mat2 <- rep(resnb_kernel[[2]],50)
	resnb_kernel_2_mat3 <- resnb_kernel[[3]][1,]
	for(i in 2:50){
		resnb_kernel_2_mat3 <- c(resnb_kernel_2_mat3, resnb_kernel[[3]][i,])
	}
	resnb_kernel_mat <- as.data.frame(cbind(resnb_kernel_2_mat, resnb_kernel_2_mat2, resnb_kernel_2_mat3))
	resident.nonbreeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(resident.nonbreeding.niche.raster) <- extent(c(-3, 3, -3,3))
	resident.nonbreeding.niche.raster <- rasterize(resnb_kernel_mat[,1:2], resident.nonbreeding.niche.raster, resnb_kernel_mat[,3])
	resident.nonbreeding.niche.raster = resident.nonbreeding.niche.raster / sum(as.vector(resident.nonbreeding.niche.raster))
	# Keep only the top 99% of the kernel, set the rest to 0
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(resident.nonbreeding.niche.raster), decreasing=T)[i]
	}
	resident.nonbreeding.niche.raster[which(as.vector(resident.nonbreeding.niche.raster) < sort(as.vector(resident.nonbreeding.niche.raster), decreasing=T)[i])] = 0
	resident.nonbreeding.niche.raster = resident.nonbreeding.niche.raster / sum(as.vector(resident.nonbreeding.niche.raster))
	
	
	## stay resident on breeding ground ##
	resbr <- cbind(TempNW_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)], PrecNW_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]==1)])
	resbr_kernel <- kde2d(resbr[,1], resbr[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
	# Convert the kernel into a raster
	resbr_kernel_2_mat <- rep(resbr_kernel[[1]][1],50)
	for(i in 2:50){
		resbr_kernel_2_mat <- c(resbr_kernel_2_mat, rep(resbr_kernel[[1]][i],50))
	}
	resbr_kernel_2_mat2 <- rep(resbr_kernel[[2]],50)
	resbr_kernel_2_mat3 <- resbr_kernel[[3]][1,]
	for(i in 2:50){
		resbr_kernel_2_mat3 <- c(resbr_kernel_2_mat3, resbr_kernel[[3]][i,])
	}
	resbr_kernel_mat <- as.data.frame(cbind(resbr_kernel_2_mat, resbr_kernel_2_mat2, resbr_kernel_2_mat3))
	resident.breeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(resident.breeding.niche.raster) <- extent(c(-3, 3, -3,3))
	resident.breeding.niche.raster <- rasterize(resbr_kernel_mat[,1:2], resident.breeding.niche.raster, resbr_kernel_mat[,3])
	resident.breeding.niche.raster = resident.breeding.niche.raster / sum(as.vector(resident.breeding.niche.raster))
	# Keep only the top 99% of the kernel, set the rest to 0
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(resident.breeding.niche.raster), decreasing=T)[i]
	}
	resident.breeding.niche.raster[which(as.vector(resident.breeding.niche.raster) < sort(as.vector(resident.breeding.niche.raster), decreasing=T)[i])] = 0
	resident.breeding.niche.raster = resident.breeding.niche.raster / sum(as.vector(resident.breeding.niche.raster))
	
	# Compute distances if resident
	nicheDistResNB_NHEH[j] <- emd(resident.nonbreeding.niche.raster, nonbreeding.niche.raster, threshold=2)
	nicheDistResBR_NHEH[j] <- emd(breeding.niche.raster, resident.breeding.niche.raster, threshold=2)

}



nicheDistObs_SHEH <- vector()
nicheDistResNB_SHEH <- vector()
nicheDistResBR_SHEH <- vector()

for(j in 1:length(selectedSpeciesSH_EH)){
	
	## breeding niche ##
	br <- cbind(TempNW_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)], PrecNW_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])
	br_kernel <- kde2d(br[,1], br[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
	# Convert the kernel into a raster
	br_kernel_2_mat <- rep(br_kernel[[1]][1],50)
	for(i in 2:50){
		br_kernel_2_mat <- c(br_kernel_2_mat, rep(br_kernel[[1]][i],50))
	}
	br_kernel_2_mat2 <- rep(br_kernel[[2]],50)
	br_kernel_2_mat3 <- br_kernel[[3]][1,]
	for(i in 2:50){
		br_kernel_2_mat3 <- c(br_kernel_2_mat3, br_kernel[[3]][i,])
	}
	br_kernel_mat <- as.data.frame(cbind(br_kernel_2_mat, br_kernel_2_mat2, br_kernel_2_mat3))
	breeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(breeding.niche.raster) <- extent(c(-3, 3, -3,3))
	breeding.niche.raster <- rasterize(br_kernel_mat[,1:2], breeding.niche.raster, br_kernel_mat[,3])
	breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))
	# Keep only the top 99% of the kernel, set the rest to 0
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(breeding.niche.raster), decreasing=T)[i]
	}
	breeding.niche.raster[which(as.vector(breeding.niche.raster) < sort(as.vector(breeding.niche.raster), decreasing=T)[i])] = 0
	breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))

	## non-breeding niche ##
	nb <- cbind(TempNS_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)], PrecNS_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])
	nb_kernel <- kde2d(nb[,1], nb[,2], h=0.5, n=50, lims=c(-3, 3, -3,3))
	nb_kernel_2_mat <- rep(nb_kernel[[1]][1],50)
	for(i in 2:50){
		nb_kernel_2_mat <- c(nb_kernel_2_mat, rep(nb_kernel[[1]][i],50))
	}
	nb_kernel_2_mat2 <- rep(nb_kernel[[2]],50)
	nb_kernel_2_mat3 <- nb_kernel[[3]][1,]
	for(i in 2:50){
		nb_kernel_2_mat3 <- c(nb_kernel_2_mat3, nb_kernel[[3]][i,])
	}
	nb_kernel_mat <- as.data.frame(cbind(nb_kernel_2_mat, nb_kernel_2_mat2, nb_kernel_2_mat3))
	nonbreeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(nonbreeding.niche.raster) <- extent(c(-3, 3, -3,3))
	nonbreeding.niche.raster <- rasterize(nb_kernel_mat[,1:2], nonbreeding.niche.raster, nb_kernel_mat[,3])
	nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i]
	}
	nonbreeding.niche.raster[which(as.vector(nonbreeding.niche.raster) < sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i])] = 0
	nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))

	# Compute observed distance
	nicheDistObs_SHEH[j] <- emd(breeding.niche.raster, nonbreeding.niche.raster, threshold=2)
	
	
	## stay resident on non-breeding ground ##
	resnb <- cbind(TempNW_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)], PrecNW_EH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])
	resnb_kernel <- kde2d(resnb[,1], resnb[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
	# Convert the kernel into a raster
	resnb_kernel_2_mat <- rep(resnb_kernel[[1]][1],50)
	for(i in 2:50){
		resnb_kernel_2_mat <- c(resnb_kernel_2_mat, rep(resnb_kernel[[1]][i],50))
	}
	resnb_kernel_2_mat2 <- rep(resnb_kernel[[2]],50)
	resnb_kernel_2_mat3 <- resnb_kernel[[3]][1,]
	for(i in 2:50){
		resnb_kernel_2_mat3 <- c(resnb_kernel_2_mat3, resnb_kernel[[3]][i,])
	}
	resnb_kernel_mat <- as.data.frame(cbind(resnb_kernel_2_mat, resnb_kernel_2_mat2, resnb_kernel_2_mat3))
	resident.nonbreeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(resident.nonbreeding.niche.raster) <- extent(c(-3, 3, -3,3))
	resident.nonbreeding.niche.raster <- rasterize(resnb_kernel_mat[,1:2], resident.nonbreeding.niche.raster, resnb_kernel_mat[,3])
	resident.nonbreeding.niche.raster = resident.nonbreeding.niche.raster / sum(as.vector(resident.nonbreeding.niche.raster))
	# Keep only the top 99% of the kernel, set the rest to 0
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(resident.nonbreeding.niche.raster), decreasing=T)[i]
	}
	resident.nonbreeding.niche.raster[which(as.vector(resident.nonbreeding.niche.raster) < sort(as.vector(resident.nonbreeding.niche.raster), decreasing=T)[i])] = 0
	resident.nonbreeding.niche.raster = resident.nonbreeding.niche.raster / sum(as.vector(resident.nonbreeding.niche.raster))
	
	
	## stay resident on breeding ground ##
	resbr <- cbind(TempNS_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)], PrecNS_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])
	resbr_kernel <- kde2d(resbr[,1], resbr[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
	# Convert the kernel into a raster
	resbr_kernel_2_mat <- rep(resbr_kernel[[1]][1],50)
	for(i in 2:50){
		resbr_kernel_2_mat <- c(resbr_kernel_2_mat, rep(resbr_kernel[[1]][i],50))
	}
	resbr_kernel_2_mat2 <- rep(resbr_kernel[[2]],50)
	resbr_kernel_2_mat3 <- resbr_kernel[[3]][1,]
	for(i in 2:50){
		resbr_kernel_2_mat3 <- c(resbr_kernel_2_mat3, resbr_kernel[[3]][i,])
	}
	resbr_kernel_mat <- as.data.frame(cbind(resbr_kernel_2_mat, resbr_kernel_2_mat2, resbr_kernel_2_mat3))
	resident.breeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(resident.breeding.niche.raster) <- extent(c(-3, 3, -3,3))
	resident.breeding.niche.raster <- rasterize(resbr_kernel_mat[,1:2], resident.breeding.niche.raster, resbr_kernel_mat[,3])
	resident.breeding.niche.raster = resident.breeding.niche.raster / sum(as.vector(resident.breeding.niche.raster))
	# Keep only the top 99% of the kernel, set the rest to 0
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(resident.breeding.niche.raster), decreasing=T)[i]
	}
	resident.breeding.niche.raster[which(as.vector(resident.breeding.niche.raster) < sort(as.vector(resident.breeding.niche.raster), decreasing=T)[i])] = 0
	resident.breeding.niche.raster = resident.breeding.niche.raster / sum(as.vector(resident.breeding.niche.raster))
	
	# Compute distances if resident
	nicheDistResNB_SHEH[j] <- emd(resident.nonbreeding.niche.raster, nonbreeding.niche.raster, threshold=2)
	nicheDistResBR_SHEH[j] <- emd(breeding.niche.raster, resident.breeding.niche.raster, threshold=2)

}


j = 13

## breeding niche ##
	br <- cbind(TempNW_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)], PrecNW_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])
	br_kernel <- kde2d(br[,1], br[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
	# Convert the kernel into a raster
	br_kernel_2_mat <- rep(br_kernel[[1]][1],50)
	for(i in 2:50){
		br_kernel_2_mat <- c(br_kernel_2_mat, rep(br_kernel[[1]][i],50))
	}
	br_kernel_2_mat2 <- rep(br_kernel[[2]],50)
	br_kernel_2_mat3 <- br_kernel[[3]][1,]
	for(i in 2:50){
		br_kernel_2_mat3 <- c(br_kernel_2_mat3, br_kernel[[3]][i,])
	}
	br_kernel_mat <- as.data.frame(cbind(br_kernel_2_mat, br_kernel_2_mat2, br_kernel_2_mat3))
	breeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(breeding.niche.raster) <- extent(c(-3, 3, -3,3))
	breeding.niche.raster <- rasterize(br_kernel_mat[,1:2], breeding.niche.raster, br_kernel_mat[,3])
	breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))
	# Keep only the top 99% of the kernel, set the rest to 0
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(breeding.niche.raster), decreasing=T)[i]
	}
	breeding.niche.raster[which(as.vector(breeding.niche.raster) < sort(as.vector(breeding.niche.raster), decreasing=T)[i])] = 0
	breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))

	## non-breeding niche ##
	nb <- cbind(TempNS_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)], PrecNS_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])
	nb_kernel <- kde2d(nb[,1], nb[,2], h=0.5, n=50, lims=c(-3, 3, -3,3))
	nb_kernel_2_mat <- rep(nb_kernel[[1]][1],50)
	for(i in 2:50){
		nb_kernel_2_mat <- c(nb_kernel_2_mat, rep(nb_kernel[[1]][i],50))
	}
	nb_kernel_2_mat2 <- rep(nb_kernel[[2]],50)
	nb_kernel_2_mat3 <- nb_kernel[[3]][1,]
	for(i in 2:50){
		nb_kernel_2_mat3 <- c(nb_kernel_2_mat3, nb_kernel[[3]][i,])
	}
	nb_kernel_mat <- as.data.frame(cbind(nb_kernel_2_mat, nb_kernel_2_mat2, nb_kernel_2_mat3))
	nonbreeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(nonbreeding.niche.raster) <- extent(c(-3, 3, -3,3))
	nonbreeding.niche.raster <- rasterize(nb_kernel_mat[,1:2], nonbreeding.niche.raster, nb_kernel_mat[,3])
	nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i]
	}
	nonbreeding.niche.raster[which(as.vector(nonbreeding.niche.raster) < sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i])] = 0
	nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))

	# Compute observed distance
	nicheDistObs_SHEH[j] <- emd(breeding.niche.raster, nonbreeding.niche.raster, threshold=2)
	
	
	## stay resident on non-breeding ground ##
	resnb <- cbind(TempNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)], PrecNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])
	resnb_kernel <- kde2d(resnb[,1], resnb[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
	# Convert the kernel into a raster
	resnb_kernel_2_mat <- rep(resnb_kernel[[1]][1],50)
	for(i in 2:50){
		resnb_kernel_2_mat <- c(resnb_kernel_2_mat, rep(resnb_kernel[[1]][i],50))
	}
	resnb_kernel_2_mat2 <- rep(resnb_kernel[[2]],50)
	resnb_kernel_2_mat3 <- resnb_kernel[[3]][1,]
	for(i in 2:50){
		resnb_kernel_2_mat3 <- c(resnb_kernel_2_mat3, resnb_kernel[[3]][i,])
	}
	resnb_kernel_mat <- as.data.frame(cbind(resnb_kernel_2_mat, resnb_kernel_2_mat2, resnb_kernel_2_mat3))
	resident.nonbreeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(resident.nonbreeding.niche.raster) <- extent(c(-3, 3, -3,3))
	resident.nonbreeding.niche.raster <- rasterize(resnb_kernel_mat[,1:2], resident.nonbreeding.niche.raster, resnb_kernel_mat[,3])
	resident.nonbreeding.niche.raster = resident.nonbreeding.niche.raster / sum(as.vector(resident.nonbreeding.niche.raster))
	# Keep only the top 99% of the kernel, set the rest to 0
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(resident.nonbreeding.niche.raster), decreasing=T)[i]
	}
	resident.nonbreeding.niche.raster[which(as.vector(resident.nonbreeding.niche.raster) < sort(as.vector(resident.nonbreeding.niche.raster), decreasing=T)[i])] = 0
	resident.nonbreeding.niche.raster = resident.nonbreeding.niche.raster / sum(as.vector(resident.nonbreeding.niche.raster))
	
	
	## stay resident on breeding ground ##
	resbr <- cbind(TempNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)], PrecNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[j]))]==1)])
	resbr_kernel <- kde2d(resbr[,1], resbr[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
	# Convert the kernel into a raster
	resbr_kernel_2_mat <- rep(resbr_kernel[[1]][1],50)
	for(i in 2:50){
		resbr_kernel_2_mat <- c(resbr_kernel_2_mat, rep(resbr_kernel[[1]][i],50))
	}
	resbr_kernel_2_mat2 <- rep(resbr_kernel[[2]],50)
	resbr_kernel_2_mat3 <- resbr_kernel[[3]][1,]
	for(i in 2:50){
		resbr_kernel_2_mat3 <- c(resbr_kernel_2_mat3, resbr_kernel[[3]][i,])
	}
	resbr_kernel_mat <- as.data.frame(cbind(resbr_kernel_2_mat, resbr_kernel_2_mat2, resbr_kernel_2_mat3))
	resident.breeding.niche.raster <- raster(ncol=50, nrow=50)
	extent(resident.breeding.niche.raster) <- extent(c(-3, 3, -3,3))
	resident.breeding.niche.raster <- rasterize(resbr_kernel_mat[,1:2], resident.breeding.niche.raster, resbr_kernel_mat[,3])
	resident.breeding.niche.raster = resident.breeding.niche.raster / sum(as.vector(resident.breeding.niche.raster))
	# Keep only the top 99% of the kernel, set the rest to 0
	thres = 0
	i=0
	while(thres <= 0.99){
		i = i+1
		thres = thres + sort(as.vector(resident.breeding.niche.raster), decreasing=T)[i]
	}
	resident.breeding.niche.raster[which(as.vector(resident.breeding.niche.raster) < sort(as.vector(resident.breeding.niche.raster), decreasing=T)[i])] = 0
	resident.breeding.niche.raster = resident.breeding.niche.raster / sum(as.vector(resident.breeding.niche.raster))
	
	# Compute distances if resident
	nicheDistResNB_SHEH[j] <- emd(resident.nonbreeding.niche.raster, nonbreeding.niche.raster, threshold=2)
	nicheDistResBR_SHEH[j] <- emd(breeding.niche.raster, resident.breeding.niche.raster, threshold=2)







nicheDistObs = c(nicheDistObs_NHEH, nicheDistObs_NHWH, nicheDistObs_SHEH, nicheDistObs_SHWH)
nicheDistResNB = c(nicheDistResNB_NHEH, nicheDistResNB_NHWH, nicheDistResNB_SHEH, nicheDistResNB_SHWH)
nicheDistResBR = c(nicheDistResBR_NHEH, nicheDistResBR_NHWH, nicheDistResBR_SHEH, nicheDistResBR_SHWH)











#########################
######   FIGURES   ######
#########################



nicheDist_diff_resNB <- nicheDistResNB_rescaled - nicheDist_rescaled
nicheDist_diff_resBR <- nicheDistResBR_rescaled - nicheDist_rescaled
nicheRank_diff_resNB <- rank_niche_resNB - rank_niche
nicheRank_diff_resBR <- rank_niche_resBR - rank_niche
geoRank_diff_resNB <- rank_geo_resNB - rank_geo
resourcesRank_diff_resNB <- rank_resources_resNB - rank_resources



###  Figure 2

## re-run model for j=182 (Bobolink)

hexidEH <- PresAbs_BR_NH_EH[,1]

par(mfrow=c(2,2), mar=c(2.5,2.5,0.2,1.5), mgp=c(1.5,0.5,0))

plot(hexgridEH, col="grey", border = "grey")
plot(hexgridEH[match(PresAbs_BR_NH[which(PresAbs_BR_NH[,which(colnames(PresAbs_NB_NH) == "Ixobrychus.eurhythmus")] == 1),1], hexgridEH@data[,1]),], col="red3", border = "red3", add=T)
plot(hexgridEH[match(hexidEH[range.sim_EH[[182]][[3]]], hexgridEH@data[,1]),], col="black", border = "black", add=T)
plot(hexgridEH[match(hexidEH[range.sim_EH[[182]][[17]]], hexgridEH@data[,1]),], col="black", border = "black", add=T)

mtext("A", cex=1.8, side=3, line=-1.5, at=-50)

plot(ggg, nnn, pch=20, axes=F, xlab="Geographic distance", ylab = "Niche distance", cex.lab=1.3, add=F)
 axis(1)
 axis(2)
 abline(a = nobs, b = 0, col="orange")
 abline(a = (nobs + gobs), b = -1, col="blue")
 points(gobs, nobs, pch=20, col="red3", cex=2.5)

mtext("C", cex=1.8, side=3, line=-1.5, at=-0.2)

plot(hexgridEH, col="grey", border = "grey")
plot(hexgridEH[match(PresAbs_NB_NH[which(PresAbs_NB_NH[,which(colnames(PresAbs_NB_NH) == "Ixobrychus.eurhythmus")] == 1),1], hexgridEH@data[,1]),], col="red3", border = "red3", add=T)
plot(hexgridEH[match(hexidEH[range.simNB_EH[[182]][[1]]], hexgridEH@data[,1]),], col="black", border = "black", add=T)
plot(hexgridEH[match(hexidEH[range.simNB_EH[[182]][[4]]], hexgridEH@data[,1]),], col="black", border = "black", add=T)

mtext("B", cex=1.8, side=3, line=-1.5, at=-50)

plot(nnnrrrggg, rnorm(length(nnnrrrggg), 1, 0.1), pch=20, axes=F, xlab="Niche distance + geographic distance + resource scarcity", ylab = "", ylim=c(0,2), xlim = c(0.5, 2.5), cex.lab=1.3)
axis(1)
abline(v=(nobs + robs + gobs), col="green")
points(nobs + robs + gobs, 1, pch=20, ylim=c(0,2), col="red3", cex=2.5)

mtext("D", cex=1.8, side=3, line=-1.5, at=0.1)





####  Figure 3 - single factor analyses

par(mfrow=c(2,2), mar=c(2.5,3,1,2), mgp=c(1.5,0.5,0))

# Niche tracking - histograms
hist(nicheRank_diff_resNB, xlim=c(-1,1), xlab="Difference in niche distance with if resident in NB", ylab="frequency", main="", xaxt="n", axes=F, breaks=10)
abline(v=0, col="red")
axis(side=2)
axis(side=1, at=c(-1,-0.5,0,0.5,1), labels = as.character(c(-1,-0.5,0,0.5,1)))
mtext("A", cex=1.5, side=3, line=-0.5, at=-1.35)

hist(nicheRank_diff_resBR, xlim=c(-1,1), xlab="Difference in niche distance with if resident in BR", ylab="frequency", main="", xaxt="n", axes=F, breaks=10)
abline(v=0, col="red")
axis(side=2)
axis(side=1, at=c(-1,-0.5,0,0.5,1), labels = as.character(c(-1,-0.5,0,0.5,1)))
mtext("B", cex=1.5, side=3, line=-0.5, at=-1.35)

# geo distance - histograms
hist(geoRank_diff_resNB, xlim=c(-1,1), xlab="Difference in geographical distance with if resident", ylab="frequency", main="", xaxt="n", axes=F, breaks=5)
abline(v=0, col="red")
axis(side=2)
axis(side=1, at=c(-1,-0.5,0,0.5,1), labels = as.character(c(-1,-0.5,0,0.5,1)))
mtext("C", cex=1.5, side=3, line=-0.5, at=-1.35)

hist(resourcesRank_diff_resNB, xlim=c(-1,1), xlab="Difference in resource scarcity with if resident", ylab="frequency", main="", xaxt="n", axes=F, breaks=10)
abline(v=0, col="red")
axis(side=2)
axis(side=1, at=c(-1,-0.5,0,0.5,1), labels = as.character(c(-1,-0.5,0,0.5,1)))
mtext("D", cex=1.5, side=3, line=-0.5, at=-1.35)




####  Figure 4 - trade-off

par(mfrow=c(2,2), mar=c(2.5,3,1.5,3), mgp=c(1.5,0.5,0))

# Relationship between geographical distance and niche distance
plot(geoDist_rescaled, nicheDist_rescaled, ylim=c(0,1), xlim=c(0,1), xlab="Scaled  geographical distance", ylab="Scaled niche distance", main="", xaxt="n", axes=F, pch=20, cex=0.7)
axis(side=2)
axis(side=1, at=c(0,0.2,0.4,0.6,0.8,1), labels = as.character(c(0,0.2,0.4,0.6,0.8,1)))
mod = lm(nicheDist_rescaled ~ geoDist_rescaled + I(geoDist_rescaled ^2))
lines(sort(geoDist_rescaled), fitted(mod)[order(geoDist_rescaled)], type="l", col="red")
mtext("A", cex=1.5, side=3, line=0, at=-0.15)
mtext(bquote(R^2 == .(round(summary(mod)$r.squared,2))), cex=1, side=3, line=-1.3, at=0.6)

# Relationship between niche distance and resource scarcity
plot(nicheDist_rescaled, resources_rescaled, ylim=c(0,1), xlim=c(0,1), xlab="Scaled niche distance", ylab="Scaled resource scarcity", main="", xaxt="n", axes=F, pch=20, cex=0.7)
axis(side=2)
axis(side=1, at=c(0,0.2,0.4,0.6,0.8,1), labels = as.character(c(0,0.2,0.4,0.6,0.8,1)))
mod = lm(resources_rescaled ~ nicheDist_rescaled + I(nicheDist_rescaled ^2))
lines(sort(nicheDist_rescaled), fitted(mod)[order(nicheDist_rescaled)], type="l", col="red")
mtext("B", cex=1.5, side=3, line=0, at=-0.15)
mtext(bquote(R^2 == .(round(summary(mod)$r.squared,2))), cex=1, side=3, line=-1.3, at=0.7)

# Relationship between geographical distance and resource scarcity
plot(geoDist_rescaled, resources_rescaled, ylim=c(0,1), xlim=c(0,1), xlab="Scaled geographical distance", ylab="Scaled  resource scarcity", main="", xaxt="n", axes=F, pch=20, cex=0.7)
axis(side=2)
axis(side=1, at=c(0,0.2,0.4,0.6,0.8,1), labels = as.character(c(0,0.2,0.4,0.6,0.8,1)))
mod = lm(resources_rescaled ~ geoDist_rescaled + I(geoDist_rescaled^2))
lines(sort(geoDist_rescaled), fitted(mod)[order(geoDist_rescaled)], type="l", col="red")
mtext("C", cex=1.5, side=3, line=0, at=-0.15)
mtext(bquote(R^2 == .(round(summary(mod)$r.squared,2))), cex=1, side=3, line=-2.5, at=0.6)

# Relationship between the three factors
hist(geoDist_rescaled, xlim=c(0,1), xlab="", ylab="", main="", xaxt="n", axes=F, col="light grey", border="grey")
axis(side=4)
axis(side=1, at=c(0,0.2,0.4,0.6,0.8,1), labels = as.character(c(0,0.2,0.4,0.6,0.8,1)))
par(new=T, mar=c(2.5,3,1.5,3), mgp=c(1.5,0.5,0))
plot(geoDist_rescaled, resources_rescaled, ylim=c(0,1), xlim=c(0,1), xlab="Scaled geographical distance", ylab="Scaled  resource scarcity", main="", xaxt="n", axes=F, pch=20, cex=0.7, col="yellow")
axis(side=2)
axis(side=1, at=c(0,0.2,0.4,0.6,0.8,1), labels = as.character(c(0,0.2,0.4,0.6,0.8,1)))
points(geoDist_rescaled[which(nicheDist_rescaled < quantile(nicheDist_rescaled)[4])], resources_rescaled[which(nicheDist_rescaled < quantile(nicheDist_rescaled)[4])], xlim=c(0,1), col="orange", cex=0.7, pch=20)
points(geoDist_rescaled[which(nicheDist_rescaled < quantile(nicheDist_rescaled)[3])], resources_rescaled[which(nicheDist_rescaled < quantile(nicheDist_rescaled)[3])], xlim=c(0,1), col="red", cex=0.7, pch=20)
points(geoDist_rescaled[which(nicheDist_rescaled < quantile(nicheDist_rescaled)[2])], resources_rescaled[which(nicheDist_rescaled < quantile(nicheDist_rescaled)[2])], xlim=c(0,1), col="brown4", cex=0.7, pch=20)
mtext("D", cex=1.5, side=3, line=0, at=-0.15)
mtext("Number of species", side=4, line=1.4, cex=0.85,las=0)
legend("topright", inset=.1, bg="white", box.col="white", title="Scaled niche\ndistance", c(paste("<", round(quantile(nicheDist_rescaled)[2], 2), sep=" "), paste(round(quantile(nicheDist_rescaled)[2], 2), round(quantile(nicheDist_rescaled)[3], 2), sep="–"), paste(round(quantile(nicheDist_rescaled)[3], 2), round(quantile(nicheDist_rescaled)[4], 2), sep="–"), paste(">", round(quantile(nicheDist_rescaled)[4], 2), sep=" ")), fill=c("brown4", "red", "orange", "yellow"), cex=0.8)


ρ




####  Figure 5 - histograms multi-factors analyses

par(mfrow=c(3,3), mar=c(4.5,3.5,1.5,0.5), mgp=c(2,0.5,0))

hist(rank_niche, main="", ylab="Number of species", xlab="Scaled rank for niche distance", xlim=c(0,1), cex.lab=1.3, col="grey", border="white")
mtext("A", cex=1.3, side=3, line=-0.25, at=-0.2)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-2.5, at=0.7)
hist(rank_geo, main="", ylab="", xlab="Scaled rank for geographic distance", xlim=c(0,1), cex.lab=1.3, col="grey", border="white")
mtext("B", cex=1.3, side=3, line=-0.25, at=-0.2)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-2.5, at=0.7)
hist(rank_resources, main="", ylab="", xlab="Scaled rank for resource scarcity", xlim=c(0,1), cex.lab=1.3, col="grey", border="white")
mtext("C", cex=1.3, side=3, line=-0.25, at=-0.2)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-2.5, at=0.7)

par(new=F, mar=c(4.5,3.5,1.5,0.5), mgp=c(3,0.5,0))
hist(rank_niche_geo, main="", ylab="", xlab="Scaled rank for niche distance \n+ geographic distance", xlim=c(0,1), cex.lab=1.3, col="grey", border="white")
mtext("Number of species", cex=0.9, side=2, line=2, at=135)
mtext("D", cex=1.3, side=3, line=-0.25, at=-0.2)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-2.5, at=0.7)
hist(rank_niche_resources, main="", ylab="", xlab="Scaled rank for niche distance \n+ resource scarcity", xlim=c(0,1), cex.lab=1.3, col="grey", border="white")
mtext("E", cex=1.3, side=3, line=-0.25, at=-0.2)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-2.5, at=0.7)
hist(rank_resources_geo, main="", ylab="", xlab="Scaled rank for geographic distance \n+ resource scarcity", xlim=c(0,1), cex.lab=1.3, col="grey", border="white")
mtext("F", cex=1.3, side=3, line=-0.25, at=-0.2)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-2.5, at=0.7)

hist(rank_niche_geo_resources, main="", ylab="", xlab="Scaled rank for niche distance \n+ geographic distance + resource scarcity", xlim=c(0,1), cex.lab=1.3, col="grey", border="white")
mtext("G", cex=1.3, side=3, line=-0.25, at=-0.2)
abline(a=65.2, b=0)
mtext("P < 0.0001", cex=1, side=3, line=-2.5, at=0.7)
mtext("Number of species", cex=0.9, side=2, line=2, at=155)

ks.test(rank_niche, runif(100000,0,1), alternative="greater")$p.value
ks.test(rank_geo, runif(100000,0,1), alternative="greater")$p.value
ks.test(rank_resources, runif(100000,0,1), alternative="greater")$p.value
ks.test(rank_niche_geo, runif(100000,0,1), alternative="greater")$p.value
ks.test(rank_resources_geo, runif(100000,0,1), alternative="greater")$p.value
ks.test(rank_niche_resources, runif(100000,0,1), alternative="greater")$p.value
ks.test(rank_niche_geo_resources, runif(100000,0,1), alternative="greater")$p.value




#### MAPS

PresAbs_BR_NH_WH_compet <- PresAbs_BR_NH_WH[,match(names(selectedSpeciesNH_WH), colnames(PresAbs_BR_NH_WH))]
PresAbs_NB_NH_WH_compet <- PresAbs_NB_NH_WH[,match(names(selectedSpeciesNH_WH), colnames(PresAbs_NB_NH_WH))]
PresAbs_BR_SH_WH_compet <- PresAbs_BR_SH_WH[,match(names(selectedSpeciesSH_WH), colnames(PresAbs_BR_SH_WH))]
PresAbs_NB_SH_WH_compet <- PresAbs_NB_SH_WH[,match(names(selectedSpeciesSH_WH), colnames(PresAbs_NB_SH_WH))]
PresAbs_BR_WH_compet <- cbind(PresAbs_BR_NH_WH_compet, PresAbs_BR_SH_WH_compet)
PresAbs_NB_WH_compet <- cbind(PresAbs_NB_NH_WH_compet, PresAbs_NB_SH_WH_compet)

PresAbs_BR_NH_EH_compet <- PresAbs_BR_NH_EH[,match(names(selectedSpeciesNH_EH), colnames(PresAbs_BR_NH_EH))[-which(is.na(match(names(selectedSpeciesNH_EH), colnames(PresAbs_BR_NH_EH)))==T)]]
PresAbs_NB_NH_EH_compet <- PresAbs_NB_NH_EH[,match(names(selectedSpeciesNH_EH), colnames(PresAbs_NB_NH_EH))[-which(is.na(match(names(selectedSpeciesNH_EH), colnames(PresAbs_NB_NH_EH)))==T)]]
for(j in 380:381){
PresAbs_BR_NH_EH_compet[,j] <- PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesNH_EH[j]))]
PresAbs_NB_NH_EH_compet[,j] <- PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesNH_EH[j]))]
}
PresAbs_BR_SH_EH_compet <- PresAbs_BR_SH_EH[,match(names(selectedSpeciesSH_EH), colnames(PresAbs_BR_SH_EH))[-which(is.na(match(names(selectedSpeciesSH_EH), colnames(PresAbs_BR_SH_EH)))==T)]]
PresAbs_NB_SH_EH_compet <- PresAbs_NB_SH_EH[,match(names(selectedSpeciesSH_EH), colnames(PresAbs_NB_SH_EH))[-which(is.na(match(names(selectedSpeciesSH_EH), colnames(PresAbs_BR_SH_EH)))==T)]]
PresAbs_BR_SH_EH_compet[,13] <- PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesSH_EH[13]))]
PresAbs_NB_SH_EH_compet[,13] <- PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesSH_EH[13]))]
PresAbs_BR_EH_compet <- cbind(PresAbs_BR_NH_EH_compet, PresAbs_BR_SH_EH_compet)
PresAbs_NB_EH_compet <- cbind(PresAbs_NB_NH_EH_compet, PresAbs_NB_SH_EH_compet)


#####   Figure 6 – Weighted Richness Maps

par(mfrow=c(4,2), mar=c(0.1,0.1,0.1,0.1), mgp=c(1.5,0.5,0))

## Weighted richness maps for Niche

weighted.richness.br_WH <- apply(PresAbs_BR_WH_compet, 1, function(x) sum(((c(rank_nicheDistObs_WH, rank_nicheDistObs_SHWH)-1)/200)*x)) / apply(PresAbs_BR_WH_compet, 1, sum)
weighted.richness.nb_WH <- apply(PresAbs_NB_WH_compet, 1, function(x) sum(((c(rank_nicheDistObs_WH, rank_nicheDistObs_SHWH)-1)/200)*x)) / apply(PresAbs_NB_WH_compet, 1, sum)
weighted.richness.br_WH[which(weighted.richness.br_WH == "NaN")] <- 0
weighted.richness.nb_WH[which(weighted.richness.nb_WH == "NaN")] <- 0
weighted.richness.br_EH <- apply(PresAbs_BR_EH_compet, 1, function(x) sum(((c(rank_nicheDistObs_EH, rank_nicheDistObs_SHEH)-1)/200)*x)) / apply(PresAbs_BR_EH_compet, 1, sum)
weighted.richness.nb_EH <- apply(PresAbs_NB_EH_compet, 1, function(x) sum(((c(rank_nicheDistObs_EH, rank_nicheDistObs_SHEH)-1)/200)*x)) / apply(PresAbs_NB_EH_compet, 1, sum)
weighted.richness.br_EH[which(weighted.richness.br_EH == "NaN")] <- 0
weighted.richness.nb_EH[which(weighted.richness.nb_EH == "NaN")] <- 0
weighted.richness.br <- c(weighted.richness.br_WH, weighted.richness.br_EH)
weighted.richness.nb <- c(weighted.richness.nb_WH, weighted.richness.nb_EH)

rbPal <- colorRampPalette(c("yellow", "red3"))

datcol <- rbPal(5)[as.numeric(cut(weighted.richness.br, breaks=c(-0.1,0.1,0.2,0.3,0.4,1.1)))]
plot(hexgrid[match(c(PresAbs_BR_NH_WH[,1], PresAbs_BR_NH_EH[,1]), hexgrid@data[,1]),], col= datcol, border = datcol, bg="grey")
mtext("A", cex=1.3, side=3, line=-2, at=-180)
legend("bottomleft", inset=.04, bg="grey", box.col="grey", title="Average rank\nof migrants", c("> 0.4","0.3–0.4", "0.2–0.3", "0.1–0.2", "0–0.1"), fill=rev(rbPal(5)), cex=0.8)

datcol <- rbPal(5)[as.numeric(cut(weighted.richness.nb, breaks=c(-0.1,0.1,0.2,0.3,0.4,1.1)))]
plot(hexgrid[match(c(PresAbs_NB_NH_WH[,1], PresAbs_NB_NH_EH[,1]), hexgrid@data[,1]),], col= datcol, border = datcol, bg="grey")
mtext("B", cex=1.3, side=3, line=-2, at=-180)




## Weighted richness maps for Resources + GeoDist

weighted.richness.br_WH <- apply(PresAbs_BR_WH_compet, 1, function(x) sum(((c(rank_resources_geoDist_Obs_WH, rank_resources_geoDist_Obs_SHWH)-1)/200)*x)) / apply(PresAbs_BR_WH_compet, 1, sum)
weighted.richness.nb_WH <- apply(PresAbs_NB_WH_compet, 1, function(x) sum(((c(rank_resources_geoDist_Obs_WH, rank_resources_geoDist_Obs_SHWH)-1)/200)*x)) / apply(PresAbs_NB_WH_compet, 1, sum)
weighted.richness.br_WH[which(weighted.richness.br_WH == "NaN")] <- 0
weighted.richness.nb_WH[which(weighted.richness.nb_WH == "NaN")] <- 0
weighted.richness.br_EH <- apply(PresAbs_BR_EH_compet, 1, function(x) sum(((c(rank_resources_geoDist_Obs_EH, rank_resources_geoDist_Obs_SHEH)-1)/200)*x)) / apply(PresAbs_BR_EH_compet, 1, sum)
weighted.richness.nb_EH <- apply(PresAbs_NB_EH_compet, 1, function(x) sum(((c(rank_resources_geoDist_Obs_EH, rank_resources_geoDist_Obs_SHEH)-1)/200)*x)) / apply(PresAbs_NB_EH_compet, 1, sum)
weighted.richness.br_EH[which(weighted.richness.br_EH == "NaN")] <- 0
weighted.richness.nb_EH[which(weighted.richness.nb_EH == "NaN")] <- 0
weighted.richness.br <- c(weighted.richness.br_WH, weighted.richness.br_EH)
weighted.richness.nb <- c(weighted.richness.nb_WH, weighted.richness.nb_EH)

datcol <- rbPal(5)[as.numeric(cut(weighted.richness.br, breaks=c(-0.1,0.1,0.2,0.3,0.4,1.1)))]
plot(hexgrid[match(c(PresAbs_BR_NH_WH[,1], PresAbs_BR_NH_EH[,1]), hexgrid@data[,1]),], col= datcol, border = datcol, bg="grey")
mtext("C", cex=1.3, side=3, line=-2, at=-180)

datcol <- rbPal(5)[as.numeric(cut(weighted.richness.nb, breaks=c(-0.1,0.1,0.2,0.3,0.4,1.1)))]
plot(hexgrid[match(c(PresAbs_NB_NH_WH[,1], PresAbs_NB_NH_EH[,1]), hexgrid@data[,1]),], col= datcol, border = datcol, bg="grey")
mtext("D", cex=1.3, side=3, line=-2, at=-180)




## Weighted richness maps for Niche + Resources + GeoDist

weighted.richness.br_WH <- apply(PresAbs_BR_WH_compet, 1, function(x) sum(((c(rank_nicheDist_resources_geoDist_Obs_WH, rank_nicheDist_resources_geoDist_Obs_SHWH)-1)/200)*x)) / apply(PresAbs_BR_WH_compet, 1, sum)
weighted.richness.nb_WH <- apply(PresAbs_NB_WH_compet, 1, function(x) sum(((c(rank_nicheDist_resources_geoDist_Obs_WH, rank_nicheDist_resources_geoDist_Obs_SHWH)-1)/200)*x)) / apply(PresAbs_NB_WH_compet, 1, sum)
weighted.richness.br_WH[which(weighted.richness.br_WH == "NaN")] <- 0
weighted.richness.nb_WH[which(weighted.richness.nb_WH == "NaN")] <- 0
weighted.richness.br_EH <- apply(PresAbs_BR_EH_compet, 1, function(x) sum(((c(rank_nicheDist_resources_geoDist_Obs_EH, rank_nicheDist_resources_geoDist_Obs_SHEH)-1)/200)*x)) / apply(PresAbs_BR_EH_compet, 1, sum)
weighted.richness.nb_EH <- apply(PresAbs_NB_EH_compet, 1, function(x) sum(((c(rank_nicheDist_resources_geoDist_Obs_EH, rank_nicheDist_resources_geoDist_Obs_SHEH)-1)/200)*x)) / apply(PresAbs_NB_EH_compet, 1, sum)
weighted.richness.br_EH[which(weighted.richness.br_EH == "NaN")] <- 0
weighted.richness.nb_EH[which(weighted.richness.nb_EH == "NaN")] <- 0
weighted.richness.br <- c(weighted.richness.br_WH, weighted.richness.br_EH)
weighted.richness.nb <- c(weighted.richness.nb_WH, weighted.richness.nb_EH)

datcol <- rbPal(5)[as.numeric(cut(weighted.richness.br, breaks=c(-0.1,0.1,0.2,0.3,0.4,1.1)))]
plot(hexgrid[match(c(PresAbs_BR_NH_WH[,1], PresAbs_BR_NH_EH[,1]), hexgrid@data[,1]),], col= datcol, border = datcol, bg="grey")
mtext("E", cex=1.3, side=3, line=-2, at=-180)

datcol <- rbPal(5)[as.numeric(cut(weighted.richness.nb, breaks=c(-0.1,0.1,0.2,0.3,0.4,1.1)))]
plot(hexgrid[match(c(PresAbs_NB_NH_WH[,1], PresAbs_NB_NH_EH[,1]), hexgrid@data[,1]),], col= datcol, border = datcol, bg="grey")
mtext("F", cex=1.3, side=3, line=-2, at=-180)




## richness in migrants

mbr_NHWH <- apply(PresAbs_BR_NH_WH_compet, 1, sum)
mnb_NHWH <- apply(PresAbs_NB_NH_WH_compet, 1, sum)
mbr_SHWH <- apply(PresAbs_BR_SH_WH_compet, 1, sum)
mnb_SHWH <- apply(PresAbs_NB_SH_WH_compet, 1, sum)
mbr_WH <- mbr_NHWH + mbr_SHWH
mnb_WH <- mnb_NHWH + mnb_SHWH
mbr_NHEH <- apply(PresAbs_BR_NH_EH_compet, 1, sum)
mnb_NHEH <- apply(PresAbs_NB_NH_EH_compet, 1, sum)
mbr_SHEH <- apply(PresAbs_BR_SH_EH_compet, 1, sum)
mnb_SHEH <- apply(PresAbs_NB_SH_EH_compet, 1, sum)
mbr_EH <- mbr_NHEH + mbr_SHEH
mnb_EH <- mnb_NHEH + mnb_SHEH

mbr <- c(mbr_WH, mbr_EH)
mnb <- c(mnb_WH, mnb_EH)


rbPal <- colorRampPalette(c("yellow3", "dark green"))

datcol <- rbPal(5)[as.numeric(cut(mbr, breaks=c(-0.1,25,50,75,100,150)))]
plot(hexgrid[match(c(PresAbs_BR_NH_WH[,1], PresAbs_BR_NH_EH[,1]), hexgrid@data[,1]),], col= datcol, border = datcol, bg="grey")
mtext("G", cex=1.3, side=3, line=-2, at=-180)
legend("bottomleft", inset=.04, bg="grey", box.col="grey", title="Richness\nin migrants", c("> 100","75–100", "50–75", "25–50", "0–25"), fill=rev(rbPal(5)), cex=0.8)

datcol <- rbPal(5)[as.numeric(cut(mnb, breaks=c(-0.1,25,50,75,100,150)))]
plot(hexgrid[match(c(PresAbs_NB_NH_WH[,1], PresAbs_NB_NH_EH[,1]), hexgrid@data[,1]),], col= datcol, border = datcol, bg="grey")
mtext("H", cex=1.3, side=3, line=-2, at=-180)






###  Figure 7 – Weighted arrows maps 


centroidsWH <- matrix(nrow=length(selectedSpeciesNH_WH), ncol=4)
centroidsEH <- matrix(nrow=length(selectedSpeciesNH_EH), ncol=4)
centroidsSHWH <- matrix(nrow=length(selectedSpeciesSH_WH), ncol=4)
centroidsSHEH <- matrix(nrow=length(selectedSpeciesSH_EH), ncol=4)
for(j in 1:length(selectedSpeciesNH_WH)){ 
centroidsWH[j,1] <- mean(west_Hem[which(PresAbs_BR_NH_WH[,match(names(selectedSpeciesNH_WH)[j], colnames(PresAbs_BR_NH_WH))]==1),1])
centroidsWH[j,2] <- mean(west_Hem[which(PresAbs_BR_NH_WH[,match(names(selectedSpeciesNH_WH)[j], colnames(PresAbs_BR_NH_WH))]==1),2])
centroidsWH[j,3] <- mean(west_Hem[which(PresAbs_NB_NH_WH[,match(names(selectedSpeciesNH_WH)[j], colnames(PresAbs_NB_NH_WH))]==1),1])
centroidsWH[j,4] <- mean(west_Hem[which(PresAbs_NB_NH_WH[,match(names(selectedSpeciesNH_WH)[j], colnames(PresAbs_NB_NH_WH))]==1),2])
}
for(j in 1:length(selectedSpeciesSH_WH)){ 
centroidsSHWH[j,1] <- mean(west_Hem[which(PresAbs_BR_SH_WH[,match(names(selectedSpeciesSH_WH)[j], colnames(PresAbs_BR_SH_WH))]==1),1])
centroidsSHWH[j,2] <- mean(west_Hem[which(PresAbs_BR_SH_WH[,match(names(selectedSpeciesSH_WH)[j], colnames(PresAbs_BR_SH_WH))]==1),2])
centroidsSHWH[j,3] <- mean(west_Hem[which(PresAbs_NB_SH_WH[,match(names(selectedSpeciesSH_WH)[j], colnames(PresAbs_NB_SH_WH))]==1),1])
centroidsSHWH[j,4] <- mean(west_Hem[which(PresAbs_NB_SH_WH[,match(names(selectedSpeciesSH_WH)[j], colnames(PresAbs_NB_SH_WH))]==1),2])
}
for(j in 1:length(selectedSpeciesNH_EH)){ 
centroidsEH[j,1] <- mean(east_Hem[which(PresAbs_BR_NH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_BR_NH_EH))]==1),1])
centroidsEH[j,2] <- mean(east_Hem[which(PresAbs_BR_NH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_BR_NH_EH))]==1),2])
centroidsEH[j,3] <- mean(east_Hem[which(PresAbs_NB_NH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_NB_NH_EH))]==1),1])
centroidsEH[j,4] <- mean(east_Hem[which(PresAbs_NB_NH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_NB_NH_EH))]==1),2])
}
for(j in 380:381){ 
centroidsEH[j,1] <- mean(east_Hem[which(PresAbs_BR_SH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_BR_SH_EH))]==1),1])
centroidsEH[j,2] <- mean(east_Hem[which(PresAbs_BR_SH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_BR_SH_EH))]==1),2])
centroidsEH[j,3] <- mean(east_Hem[which(PresAbs_NB_SH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_NB_SH_EH))]==1),1])
centroidsEH[j,4] <- mean(east_Hem[which(PresAbs_NB_SH_EH[,match(names(selectedSpeciesNH_EH)[j], colnames(PresAbs_NB_SH_EH))]==1),2])
}
for(j in 1:length(selectedSpeciesSH_EH)){ 
centroidsSHEH[j,1] <- mean(east_Hem[which(PresAbs_BR_SH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_BR_SH_EH))]==1),1])
centroidsSHEH[j,2] <- mean(east_Hem[which(PresAbs_BR_SH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_BR_SH_EH))]==1),2])
centroidsSHEH[j,3] <- mean(east_Hem[which(PresAbs_NB_SH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_NB_SH_EH))]==1),1])
centroidsSHEH[j,4] <- mean(east_Hem[which(PresAbs_NB_SH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_NB_SH_EH))]==1),2])
}
j=13
centroidsSHEH[j,1] <- mean(east_Hem[which(PresAbs_BR_NH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_BR_NH_EH))]==1),1])
centroidsSHEH[j,2] <- mean(east_Hem[which(PresAbs_BR_NH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_BR_NH_EH))]==1),2])
centroidsSHEH[j,3] <- mean(east_Hem[which(PresAbs_NB_NH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_NB_NH_EH))]==1),1])
centroidsSHEH[j,4] <- mean(east_Hem[which(PresAbs_NB_NH_EH[,match(names(selectedSpeciesSH_EH)[j], colnames(PresAbs_NB_NH_EH))]==1),2])

centroids <- rbind(centroidsWH, centroidsSHWH, centroidsEH, centroidsSHEH)

weights.NHWH <- (rank_nicheDistObs_WH-1)/200
weights.SHWH <- (rank_nicheDistObs_SHWH-1)/200
weights.NHEH <- (rank_nicheDistObs_EH-1)/200
weights.SHEH <- (rank_nicheDistObs_SHEH-1)/200
weights.niche <- c(weights.NHWH, weights.SHWH, weights.NHEH, weights.SHEH)
centroids.niche <- centroids[order(weights.niche),]
weights.niche <- weights.niche[order(weights.niche)]
weights.niche1 <- weights.niche[which(weights.niche < 0.1)]
centroids.niche1 <- centroids.niche[which(weights.niche < 0.1),]
weights.niche2 <- weights.niche[which(weights.niche >= 0.1 & weights.niche < 0.5)]
centroids.niche2 <- centroids.niche[which(weights.niche >= 0.1 & weights.niche < 0.5),]
weights.niche3 <- weights.niche[which(weights.niche >= 0.5)]
centroids.niche3 <- centroids.niche[which(weights.niche >= 0.5),]

weights.NHWH <- (rank_resources_geoDist_Obs_WH-1)/200
weights.SHWH <- (rank_resources_geoDist_Obs_SHWH-1)/200
weights.NHEH <- (rank_resources_geoDist_Obs_EH-1)/200
weights.SHEH <- (rank_resources_geoDist_Obs_SHEH-1)/200
weights.ResGeo <- c(weights.NHWH, weights.SHWH, weights.NHEH, weights.SHEH)
centroids.ResGeo <- centroids[order(weights.ResGeo),]
weights.ResGeo <- weights.ResGeo[order(weights.ResGeo)]
weights.ResGeo1 <- weights.ResGeo[which(weights.ResGeo < 0.1)]
centroids.ResGeo1 <- centroids.ResGeo[which(weights.ResGeo < 0.1),]
weights.ResGeo2 <- weights.ResGeo[which(weights.ResGeo >= 0.1 & weights.ResGeo < 0.5)]
centroids.ResGeo2 <- centroids.ResGeo[which(weights.ResGeo >= 0.1 & weights.ResGeo < 0.5),]
weights.ResGeo3 <- weights.ResGeo[which(weights.ResGeo >= 0.5)]
centroids.ResGeo3 <- centroids.ResGeo[which(weights.ResGeo >= 0.5),]

weights.NHWH <- (rank_nicheDist_resources_geoDist_Obs_WH-1)/200
weights.SHWH <- (rank_nicheDist_resources_geoDist_Obs_SHWH-1)/200
weights.NHEH <- (rank_nicheDist_resources_geoDist_Obs_EH-1)/200
weights.SHEH <- (rank_nicheDist_resources_geoDist_Obs_SHEH-1)/200
weights.nicheResGeo <- c(weights.NHWH, weights.SHWH, weights.NHEH, weights.SHEH)
centroids.nicheResGeo <- centroids[order(weights.nicheResGeo),]
weights.nicheResGeo <- weights.nicheResGeo[order(weights.nicheResGeo)]
weights.nicheResGeo1 <- weights.nicheResGeo[which(weights.nicheResGeo < 0.1)]
centroids.nicheResGeo1 <- centroids.nicheResGeo[which(weights.nicheResGeo < 0.1),]
weights.nicheResGeo2 <- weights.nicheResGeo[which(weights.nicheResGeo >= 0.1 & weights.nicheResGeo < 0.5)]
centroids.nicheResGeo2 <- centroids.nicheResGeo[which(weights.nicheResGeo >= 0.1 & weights.nicheResGeo < 0.5),]
weights.nicheResGeo3 <- weights.nicheResGeo[which(weights.nicheResGeo >= 0.5)]
centroids.nicheResGeo3 <- centroids.nicheResGeo[which(weights.nicheResGeo >= 0.5),]

library(maps)
library(geosphere)

#rbPal <- colorRampPalette(c("yellow", "brown4"))
#rbPal1 <- colorRampPalette(c("yellow", "orange"))
#rbPal2 <- colorRampPalette(c("orange", "brown4"))

par(mfrow=c(3,3), mar=c(0.1,0.1,0.1,0.1), mgp=c(1.5,0.5,0))

plot(hexgrid, col= "grey", border = "grey")
for(i in 1:length(weights.niche1)){
	inter <- gcIntermediate(c(centroids.niche1[i,1],centroids.niche1[i,2]), c(centroids.niche1[i,3],centroids.niche1[i,4]), n=50, addStartEnd=T)
	#colcol <- rbPal1(10)[as.numeric(cut(weights.niche1, breaks=10))]
	#lines(inter, lwd=(weights.niche1[i]*1)+1, col=colcol[i])
	lines(inter, lwd=1, col="yellow")
	points(centroids.niche1[i,1], centroids.niche1[i,2], pch=20, col="red", cex=0.3)
	points(centroids.niche1[i,3], centroids.niche1[i,4], pch=20, col="blue", cex=0.3)
}
mtext("A", cex=1.3, side=3, line=-2.2, at=-180)
#legend("bottomleft", inset=c(0, 0.1), bg="white", box.col="white", title="", c("Rank < 20", "Breeding centroid", "Non-breeding centroid"), lty=c(1,NA,NA), pch=c(NA,20,20), col=c("yellow","red", "blue"))

plot(hexgrid, col= "grey", border = "grey")
for(i in 1:length(weights.niche2)){
	inter <- gcIntermediate(c(centroids.niche2[i,1],centroids.niche2[i,2]), c(centroids.niche2[i,3],centroids.niche2[i,4]), n=50, addStartEnd=T)
	#colcol <- rbPal(10)[as.numeric(cut(weights.niche2, breaks=10))]
	#lines(inter, lwd=(weights.niche2[i]*1)+1, col=colcol[i])
	lines(inter, lwd=1, col="orange")
	points(centroids.niche2[i,1], centroids.niche2[i,2], pch=20, col="red", cex=0.3)
	points(centroids.niche2[i,3], centroids.niche2[i,4], pch=20, col="blue", cex=0.3)
}
mtext("B", cex=1.3, side=3, line=-2.2, at=-180)
#legend("bottomleft", inset=c(0, 0.1), bg="white", box.col="white", title="", c("20 <= Rank < 100 ", "Breeding centroid", "Non-breeding centroid"), lty=c(1,NA,NA), pch=c(NA,20,20), col=c("orange","red", "blue"))

plot(hexgrid, col= "grey", border = "grey")
for(i in 1:length(weights.niche3)){
	inter <- gcIntermediate(c(centroids.niche3[i,1],centroids.niche3[i,2]), c(centroids.niche3[i,3],centroids.niche3[i,4]), n=50, addStartEnd=T)
	#colcol <- rbPal(10)[as.numeric(cut(weights.niche3, breaks=10))]
	#lines(inter, lwd=(weights.niche3[i]*1)+1, col=colcol[i])
	lines(inter, lwd=1, col="brown4")
	points(centroids.niche3[i,1], centroids.niche3[i,2], pch=20, col="red", cex=0.3)
	points(centroids.niche3[i,3], centroids.niche3[i,4], pch=20, col="blue", cex=0.3)
}
mtext("C", cex=1.3, side=3, line=-2.2, at=-180)
#legend("bottomleft", inset=c(0, 0.1), bg="white", box.col="white", title="", c("Rank >= 100", "Breeding centroid", "Non-breeding centroid"), lty=c(1,NA,NA), pch=c(NA,20,20), col=c("brown4","red", "blue"))


plot(hexgrid, col= "grey", border = "grey")
for(i in 1:length(weights.ResGeo1)){
	inter <- gcIntermediate(c(centroids.ResGeo1[i,1],centroids.ResGeo1[i,2]), c(centroids.ResGeo1[i,3],centroids.ResGeo1[i,4]), n=50, addStartEnd=T)
	#colcol <- rbPal1(10)[as.numeric(cut(weights.ResGeo1, breaks=10))]
	#lines(inter, lwd=(weights.ResGeo1[i]*1)+1, col=colcol[i])
	lines(inter, lwd=1, col="yellow")
	points(centroids.ResGeo1[i,1], centroids.ResGeo1[i,2], pch=20, col="red", cex=0.3)
	points(centroids.ResGeo1[i,3], centroids.ResGeo1[i,4], pch=20, col="blue", cex=0.3)
}
mtext("D", cex=1.3, side=3, line=-2.2, at=-180)
#legend("bottomleft", inset=c(0, 0.1), bg="white", box.col="white", title="", c("Rank < 20", "Breeding centroid", "Non-breeding centroid"), lty=c(1,NA,NA), pch=c(NA,20,20), col=c("yellow","red", "blue"))

plot(hexgrid, col= "grey", border = "grey")
for(i in 1:length(weights.ResGeo2)){
	inter <- gcIntermediate(c(centroids.ResGeo2[i,1],centroids.ResGeo2[i,2]), c(centroids.ResGeo2[i,3],centroids.ResGeo2[i,4]), n=50, addStartEnd=T)
	#colcol <- rbPal(10)[as.numeric(cut(weights.ResGeo2, breaks=10))]
	#lines(inter, lwd=(weights.ResGeo2[i]*1)+1, col=colcol[i])
	lines(inter, lwd=1, col="orange")
	points(centroids.ResGeo2[i,1], centroids.ResGeo2[i,2], pch=20, col="red", cex=0.3)
	points(centroids.ResGeo2[i,3], centroids.ResGeo2[i,4], pch=20, col="blue", cex=0.3)
}
mtext("E", cex=1.3, side=3, line=-2.2, at=-180)
#legend("bottomleft", inset=c(0, 0.1), bg="white", box.col="white", title="", c("20 <= Rank < 100 ", "Breeding centroid", "Non-breeding centroid"), lty=c(1,NA,NA), pch=c(NA,20,20), col=c("orange","red", "blue"))

plot(hexgrid, col= "grey", border = "grey")
for(i in 1:length(weights.ResGeo3)){
	inter <- gcIntermediate(c(centroids.ResGeo3[i,1],centroids.ResGeo3[i,2]), c(centroids.ResGeo3[i,3],centroids.ResGeo3[i,4]), n=50, addStartEnd=T)
	#colcol <- rbPal(10)[as.numeric(cut(weights.ResGeo3, breaks=10))]
	#lines(inter, lwd=(weights.ResGeo3[i]*1)+1, col=colcol[i])
	lines(inter, lwd=1, col="brown4")
	points(centroids.ResGeo3[i,1], centroids.ResGeo3[i,2], pch=20, col="red", cex=0.3)
	points(centroids.ResGeo3[i,3], centroids.ResGeo3[i,4], pch=20, col="blue", cex=0.3)
}
mtext("F", cex=1.3, side=3, line=-2.2, at=-180)
#legend("bottomleft", inset=c(0, 0.1), bg="white", box.col="white", title="", c("Rank >= 100", "Breeding centroid", "Non-breeding cetroid"), lty=c(1,NA,NA), pch=c(NA,20,20), col=c("brown4","red", "blue"))


plot(hexgrid, col= "grey", border = "grey")
for(i in 1:length(weights.nicheResGeo1)){
	inter <- gcIntermediate(c(centroids.nicheResGeo1[i,1],centroids.nicheResGeo1[i,2]), c(centroids.nicheResGeo1[i,3],centroids.nicheResGeo1[i,4]), n=50, addStartEnd=T)
	#colcol <- rbPal1(10)[as.numeric(cut(weights.nicheResGeo1, breaks=10))]
	#lines(inter, lwd=(weights.nicheResGeo1[i]*1)+1, col=colcol[i])
	lines(inter, lwd=1, col="yellow")
	points(centroids.nicheResGeo1[i,1], centroids.nicheResGeo1[i,2], pch=20, col="red", cex=0.3)
	points(centroids.nicheResGeo1[i,3], centroids.nicheResGeo1[i,4], pch=20, col="blue", cex=0.3)
}
mtext("G", cex=1.3, side=3, line=-2.2, at=-180)
#legend("bottomleft", inset=c(0, 0.1), bg="white", box.col="white", title="", c("Rank < 20", "Breeding centroid", "Non-breeding centroid"), lty=c(1,NA,NA), pch=c(NA,20,20), col=c("yellow","red", "blue"))

plot(hexgrid, col= "grey", border = "grey")
for(i in 1:length(weights.nicheResGeo2)){
	inter <- gcIntermediate(c(centroids.nicheResGeo2[i,1],centroids.nicheResGeo2[i,2]), c(centroids.nicheResGeo2[i,3],centroids.nicheResGeo2[i,4]), n=50, addStartEnd=T)
	#colcol <- rbPal(10)[as.numeric(cut(weights.nicheResGeo2, breaks=10))]
	#lines(inter, lwd=(weights.nicheResGeo2[i]*1)+1, col=colcol[i])
	lines(inter, lwd=1, col="orange")
	points(centroids.nicheResGeo2[i,1], centroids.nicheResGeo2[i,2], pch=20, col="red", cex=0.3)
	points(centroids.nicheResGeo2[i,3], centroids.nicheResGeo2[i,4], pch=20, col="blue", cex=0.3)
}
mtext("H", cex=1.3, side=3, line=-2.2, at=-180)
#legend("bottomleft", inset=c(0, 0.1), bg="white", box.col="white", title="", c("20 <= Rank < 100 ", "Breeding centroid", "Non-breeding centroid"), lty=c(1,NA,NA), pch=c(NA,20,20), col=c("orange","red", "blue"))

plot(hexgrid, col= "grey", border = "grey")
for(i in 1:length(weights.nicheResGeo3)){
	inter <- gcIntermediate(c(centroids.nicheResGeo3[i,1],centroids.nicheResGeo3[i,2]), c(centroids.nicheResGeo3[i,3],centroids.nicheResGeo3[i,4]), n=50, addStartEnd=T)
	#colcol <- rbPal(10)[as.numeric(cut(weights.nicheResGeo3, breaks=10))]
	#lines(inter, lwd=(weights.nicheResGeo3[i]*1)+1, col=colcol[i])
	lines(inter, lwd=1, col="brown4")
	points(centroids.nicheResGeo3[i,1], centroids.nicheResGeo3[i,2], pch=20, col="red", cex=0.3)
	points(centroids.nicheResGeo3[i,3], centroids.nicheResGeo3[i,4], pch=20, col="blue", cex=0.3)
}
mtext("I", cex=1.3, side=3, line=-2.2, at=-180)
#legend("bottomleft", inset=c(0, 0.1), bg="white", box.col="white", title="", c("Rank >= 100", "Breeding centroid", "Non-breeding centroid"), lty=c(1,NA,NA), pch=c(NA,20,20), col=c("brown4","red", "blue"))









####   RANDOMISED NULL MODEL  ####

## NH WH

niche.distance.WH <- matrix(nrow=length(selectedSpeciesNH_WH), ncol=length(selectedSpeciesNH_WH))
resources.scarcity.WH <- matrix(nrow=length(selectedSpeciesNH_WH), ncol=length(selectedSpeciesNH_WH))
geo.distance.WH <- matrix(nrow=length(selectedSpeciesNH_WH), ncol=length(selectedSpeciesNH_WH))
for(k in 1:length(selectedSpeciesNH_WH)){
	for(h in 1:length(selectedSpeciesNH_WH)){
			
		br <- cbind(TempNS_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[h]))]==1)], PrecNS_WH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[h]))]==1)])
		br_kernel <- kde2d(br[,1], br[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
		br_kernel_2_mat <- rep(br_kernel[[1]][1],50)
		for(i in 2:50){
			br_kernel_2_mat <- c(br_kernel_2_mat, rep(br_kernel[[1]][i],50))
		}
		br_kernel_2_mat2 <- rep(br_kernel[[2]],50)
		br_kernel_2_mat3 <- br_kernel[[3]][1,]
		for(i in 2:50){
			br_kernel_2_mat3 <- c(br_kernel_2_mat3, br_kernel[[3]][i,])
		}
		br_kernel_mat <- as.data.frame(cbind(br_kernel_2_mat, br_kernel_2_mat2, br_kernel_2_mat3))
		breeding.niche.raster <- raster(ncol=50, nrow=50)
		extent(breeding.niche.raster) <- extent(c(-3, 3, -3,3))
		breeding.niche.raster <- rasterize(br_kernel_mat[,1:2], breeding.niche.raster, br_kernel_mat[,3])
		breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))
		thres = 0
		i=0
		while(thres <= 0.99){
			i = i+1
			thres = thres + sort(as.vector(breeding.niche.raster), decreasing=T)[i]
		}
		breeding.niche.raster[which(as.vector(breeding.niche.raster) < sort(as.vector(breeding.niche.raster), decreasing=T)[i])] = 0
		breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))

		nb <- cbind(TempNW_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[k]))]==1)], PrecNW_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[k]))]==1)])
		nb_kernel <- kde2d(nb[,1], nb[,2], h=0.5, n=50, lims=c(-3, 3, -3,3))
		nb_kernel_2_mat <- rep(nb_kernel[[1]][1],50)
		for(i in 2:50){
			nb_kernel_2_mat <- c(nb_kernel_2_mat, rep(nb_kernel[[1]][i],50))
		}
		nb_kernel_2_mat2 <- rep(nb_kernel[[2]],50)
		nb_kernel_2_mat3 <- nb_kernel[[3]][1,]
		for(i in 2:50){
			nb_kernel_2_mat3 <- c(nb_kernel_2_mat3, nb_kernel[[3]][i,])
		}
		nb_kernel_mat <- as.data.frame(cbind(nb_kernel_2_mat, nb_kernel_2_mat2, nb_kernel_2_mat3))
		nonbreeding.niche.raster <- raster(ncol=50, nrow=50)
		extent(nonbreeding.niche.raster) <- extent(c(-3, 3, -3,3))
		nonbreeding.niche.raster <- rasterize(nb_kernel_mat[,1:2], nonbreeding.niche.raster, nb_kernel_mat[,3])
		nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))
		thres = 0
		i=0
		while(thres <= 0.99){
			i = i+1
			thres = thres + sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i]
		}
		nonbreeding.niche.raster[which(as.vector(nonbreeding.niche.raster) < sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i])] = 0
		nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))

		niche.distance.WH[k,h] <- emd(breeding.niche.raster, nonbreeding.niche.raster, threshold=2)


		geo.distance.WH[k,h] <- rdist.earth(t(as.matrix(c(mean(west_Hem[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[k]))]==1),1]), mean(west_Hem[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[k]))]==1),2])))),t(as.matrix(c(mean(west_Hem[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[h]))]==1),1]), mean(west_Hem[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[h]))]==1),2])))), miles = F)
		resources.scarcity.WH[k,h] <-  mean(resourcesScarcity_januaryWH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[h]))]==1)]) + mean(resourcesScarcity_julyWH[which(PresAbs_BR_NH_WH[,which(colnames(PresAbs_BR_NH_WH) == names(selectedSpeciesNH_WH[k]))]==1)])
	}
print(k)
}

niche.distance.WH.rescaled <- apply(niche.distance.WH, c(1:2), function(x) (x-min(niche.distance.WH))/(max(niche.distance.WH) - min(niche.distance.WH)))
geo.distance.WH.rescaled <- apply(geo.distance.WH, c(1:2), function(x) (x-min(geo.distance.WH))/(max(geo.distance.WH) - min(geo.distance.WH)))
resources.scarcity.WH.rescaled <- apply(resources.scarcity.WH, c(1:2), function(x) (x-min(resources.scarcity.WH))/(max(resources.scarcity.WH) - min(resources.scarcity.WH)))
niche.geo.mat.WH <- niche.distance.WH.rescaled + geo.distance.WH.rescaled
niche.res.mat.WH <- niche.distance.WH.rescaled + resources.scarcity.WH.rescaled
geo.res.mat.WH <- geo.distance.WH.rescaled + resources.scarcity.WH.rescaled
niche.geo.res.mat.WH <- niche.distance.WH.rescaled + geo.distance.WH.rescaled + resources.scarcity.WH.rescaled

niche.val.obs.WH <- sum(diag(niche.distance.WH.rescaled))
geo.val.obs.WH <- sum(diag(geo.distance.WH.rescaled))
res.val.obs.WH <- sum(diag(resources.scarcity.WH.rescaled))
niche.geo.val.obs.WH <- sum(diag(niche.geo.mat.WH))
niche.res.val.obs.WH <- sum(diag(niche.res.mat.WH))
geo.res.val.obs.WH <- sum(diag(geo.res.mat.WH))
niche.geo.res.val.obs.WH <- sum(diag(niche.geo.res.mat.WH))

niche.val.random.WH <- vector()
geo.val.random.WH <- vector()
res.val.random.WH <- vector()
niche.geo.val.random.WH <- vector()
niche.res.val.random.WH <- vector()
geo.res.val.random.WH <- vector()
niche.geo.res.val.random.WH <- vector()
for(j in 1:1000){
	randomised.breeding.ranges <- sample(1:length(selectedSpeciesNH_WH))
	niche.distance.WH.rescaled.random <- niche.distance.WH.rescaled[,randomised.breeding.ranges]
	geo.distance.WH.rescaled.random <- geo.distance.WH.rescaled[,randomised.breeding.ranges]
	resources.scarcity.WH.rescaled.random <- resources.scarcity.WH.rescaled[,randomised.breeding.ranges]
	niche.geo.mat.random.WH <- niche.geo.mat.WH[,randomised.breeding.ranges]
	niche.res.mat.random.WH <- niche.res.mat.WH[,randomised.breeding.ranges]
	geo.res.mat.random.WH <- geo.res.mat.WH[,randomised.breeding.ranges]
	niche.geo.res.mat.random.WH <- niche.geo.res.mat.WH[,randomised.breeding.ranges]

	niche.val.random.WH[j] <- sum(diag(niche.distance.WH.rescaled.random))
	geo.val.random.WH[j] <- sum(diag(geo.distance.WH.rescaled.random))
	res.val.random.WH[j] <- sum(diag(resources.scarcity.WH.rescaled.random))
	niche.geo.val.random.WH[j] <- sum(diag(niche.geo.mat.random.WH))
	niche.res.val.random.WH[j] <- sum(diag(niche.res.mat.random.WH))
	geo.res.val.random.WH[j] <- sum(diag(geo.res.mat.random.WH))
	niche.geo.res.val.random.WH[j] <- sum(diag(niche.geo.res.mat.random.WH))
}

niche.distance.WH.rescaled <- apply(niche.distance.WH, c(1:2), function(x) (x-min(niche.distance.WH))/(max(niche.distance.WH) - min(niche.distance.WH)))
geo.distance.WH.rescaled <- apply(geo.distance.WH, c(1:2), function(x) (x-min(geo.distance.WH))/(max(geo.distance.WH) - min(geo.distance.WH)))
resources.scarcity.WH.rescaled <- apply(resources.scarcity.WH, c(1:2), function(x) (x-min(resources.scarcity.WH))/(max(resources.scarcity.WH) - min(resources.scarcity.WH)))
niche.geo.mat.WH <- niche.distance.WH.rescaled + geo.distance.WH.rescaled
niche.res.mat.WH <- niche.distance.WH.rescaled + resources.scarcity.WH.rescaled
geo.res.mat.WH <- geo.distance.WH.rescaled + resources.scarcity.WH.rescaled
niche.geo.res.mat.WH <- niche.distance.WH.rescaled + geo.distance.WH.rescaled + resources.scarcity.WH.rescaled

niche.val.obs.WH <- sum(diag(niche.distance.WH.rescaled))
geo.val.obs.WH <- sum(diag(geo.distance.WH.rescaled))
res.val.obs.WH <- sum(diag(resources.scarcity.WH.rescaled))
niche.geo.val.obs.WH <- sum(diag(niche.geo.mat.WH))
niche.res.val.obs.WH <- sum(diag(niche.res.mat.WH))
geo.res.val.obs.WH <- sum(diag(geo.res.mat.WH))
niche.geo.res.val.obs.WH <- sum(diag(niche.geo.res.mat.WH))

niche.val.random.WH <- vector()
geo.val.random.WH <- vector()
res.val.random.WH <- vector()
niche.geo.val.random.WH <- vector()
niche.res.val.random.WH <- vector()
geo.res.val.random.WH <- vector()
niche.geo.res.val.random.WH <- vector()
for(j in 1:1000){
	randomised.breeding.ranges <- sample(1:length(selectedSpeciesNH_WH))
	niche.distance.WH.rescaled.random <- niche.distance.WH.rescaled[,randomised.breeding.ranges]
	geo.distance.WH.rescaled.random <- geo.distance.WH.rescaled[,randomised.breeding.ranges]
	resources.scarcity.WH.rescaled.random <- resources.scarcity.WH.rescaled[,randomised.breeding.ranges]
	niche.geo.mat.random.WH <- niche.geo.mat.WH[,randomised.breeding.ranges]
	niche.res.mat.random.WH <- niche.res.mat.WH[,randomised.breeding.ranges]
	geo.res.mat.random.WH <- geo.res.mat.WH[,randomised.breeding.ranges]
	niche.geo.res.mat.random.WH <- niche.geo.res.mat.WH[,randomised.breeding.ranges]

	niche.val.random.WH[j] <- sum(diag(niche.distance.WH.rescaled.random))
	geo.val.random.WH[j] <- sum(diag(geo.distance.WH.rescaled.random))
	res.val.random.WH[j] <- sum(diag(resources.scarcity.WH.rescaled.random))
	niche.geo.val.random.WH[j] <- sum(diag(niche.geo.mat.random.WH))
	niche.res.val.random.WH[j] <- sum(diag(niche.res.mat.random.WH))
	geo.res.val.random.WH[j] <- sum(diag(geo.res.mat.random.WH))
	niche.geo.res.val.random.WH[j] <- sum(diag(niche.geo.res.mat.random.WH))
}


## SH WH

niche.distance.SHWH <- matrix(nrow=length(selectedSpeciesSH_WH), ncol=length(selectedSpeciesSH_WH))
resources.scarcity.SHWH <- matrix(nrow=length(selectedSpeciesSH_WH), ncol=length(selectedSpeciesSH_WH))
geo.distance.SHWH <- matrix(nrow=length(selectedSpeciesSH_WH), ncol=length(selectedSpeciesSH_WH))
for(k in 1:length(selectedSpeciesSH_WH)){
	for(h in 1:length(selectedSpeciesSH_WH)){
			
		br <- cbind(TempNW_WH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[h]))]==1)], PrecNW_WH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[h]))]==1)])
		br_kernel <- kde2d(br[,1], br[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
		br_kernel_2_mat <- rep(br_kernel[[1]][1],50)
		for(i in 2:50){
			br_kernel_2_mat <- c(br_kernel_2_mat, rep(br_kernel[[1]][i],50))
		}
		br_kernel_2_mat2 <- rep(br_kernel[[2]],50)
		br_kernel_2_mat3 <- br_kernel[[3]][1,]
		for(i in 2:50){
			br_kernel_2_mat3 <- c(br_kernel_2_mat3, br_kernel[[3]][i,])
		}
		br_kernel_mat <- as.data.frame(cbind(br_kernel_2_mat, br_kernel_2_mat2, br_kernel_2_mat3))
		breeding.niche.raster <- raster(ncol=50, nrow=50)
		extent(breeding.niche.raster) <- extent(c(-3, 3, -3,3))
		breeding.niche.raster <- rasterize(br_kernel_mat[,1:2], breeding.niche.raster, br_kernel_mat[,3])
		breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))
		thres = 0
		i=0
		while(thres <= 0.99){
			i = i+1
			thres = thres + sort(as.vector(breeding.niche.raster), decreasing=T)[i]
		}
		breeding.niche.raster[which(as.vector(breeding.niche.raster) < sort(as.vector(breeding.niche.raster), decreasing=T)[i])] = 0
		breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))

		nb <- cbind(TempNW_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[k]))]==1)], PrecNW_WH[which(PresAbs_NB_NH_WH[,which(colnames(PresAbs_NB_NH_WH) == names(selectedSpeciesNH_WH[k]))]==1)])
		nb_kernel <- kde2d(nb[,1], nb[,2], h=0.5, n=50, lims=c(-3, 3, -3,3))
		nb_kernel_2_mat <- rep(nb_kernel[[1]][1],50)
		for(i in 2:50){
			nb_kernel_2_mat <- c(nb_kernel_2_mat, rep(nb_kernel[[1]][i],50))
		}
		nb_kernel_2_mat2 <- rep(nb_kernel[[2]],50)
		nb_kernel_2_mat3 <- nb_kernel[[3]][1,]
		for(i in 2:50){
			nb_kernel_2_mat3 <- c(nb_kernel_2_mat3, nb_kernel[[3]][i,])
		}
		nb_kernel_mat <- as.data.frame(cbind(nb_kernel_2_mat, nb_kernel_2_mat2, nb_kernel_2_mat3))
		nonbreeding.niche.raster <- raster(ncol=50, nrow=50)
		extent(nonbreeding.niche.raster) <- extent(c(-3, 3, -3,3))
		nonbreeding.niche.raster <- rasterize(nb_kernel_mat[,1:2], nonbreeding.niche.raster, nb_kernel_mat[,3])
		nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))
		thres = 0
		i=0
		while(thres <= 0.99){
			i = i+1
			thres = thres + sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i]
		}
		nonbreeding.niche.raster[which(as.vector(nonbreeding.niche.raster) < sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i])] = 0
		nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))

		niche.distance.SHWH[k,h] <- emd(breeding.niche.raster, nonbreeding.niche.raster, threshold=2)


		geo.distance.SHWH[k,h] <- rdist.earth(t(as.matrix(c(mean(west_Hem[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[k]))]==1),1]), mean(west_Hem[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[k]))]==1),2])))),t(as.matrix(c(mean(west_Hem[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[h]))]==1),1]), mean(west_Hem[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[h]))]==1),2])))), miles = F)
		resources.scarcity.SHWH[k,h] <-  mean(resourcesScarcity_julyWH[which(PresAbs_NB_SH_WH[,which(colnames(PresAbs_NB_SH_WH) == names(selectedSpeciesSH_WH[h]))]==1)]) + mean(resourcesScarcity_januaryWH[which(PresAbs_BR_SH_WH[,which(colnames(PresAbs_BR_SH_WH) == names(selectedSpeciesSH_WH[k]))]==1)])
	}
print(k)
}

niche.distance.SHWH.rescaled <- apply(niche.distance.SHWH, c(1:2), function(x) (x-min(niche.distance.SHWH))/(max(niche.distance.SHWH) - min(niche.distance.SHWH)))
geo.distance.SHWH.rescaled <- apply(geo.distance.SHWH, c(1:2), function(x) (x-min(geo.distance.SHWH))/(max(geo.distance.SHWH) - min(geo.distance.SHWH)))
resources.scarcity.SHWH.rescaled <- apply(resources.scarcity.SHWH, c(1:2), function(x) (x-min(resources.scarcity.SHWH))/(max(resources.scarcity.SHWH) - min(resources.scarcity.SHWH)))
niche.geo.mat.SHWH <- niche.distance.SHWH.rescaled + geo.distance.SHWH.rescaled
niche.res.mat.SHWH <- niche.distance.SHWH.rescaled + resources.scarcity.SHWH.rescaled
geo.res.mat.SHWH <- geo.distance.SHWH.rescaled + resources.scarcity.SHWH.rescaled
niche.geo.res.mat.SHWH <- niche.distance.SHWH.rescaled + geo.distance.SHWH.rescaled + resources.scarcity.SHWH.rescaled

niche.val.obs.SHWH <- sum(diag(niche.distance.SHWH.rescaled))
geo.val.obs.SHWH <- sum(diag(geo.distance.SHWH.rescaled))
res.val.obs.SHWH <- sum(diag(resources.scarcity.SHWH.rescaled))
niche.geo.val.obs.SHWH <- sum(diag(niche.geo.mat.SHWH))
niche.res.val.obs.SHWH <- sum(diag(niche.res.mat.SHWH))
geo.res.val.obs.SHWH <- sum(diag(geo.res.mat.SHWH))
niche.geo.res.val.obs.SHWH <- sum(diag(niche.geo.res.mat.SHWH))

niche.val.random.SHWH <- vector()
geo.val.random.SHWH <- vector()
res.val.random.SHWH <- vector()
niche.geo.val.random.SHWH <- vector()
niche.res.val.random.SHWH <- vector()
geo.res.val.random.SHWH <- vector()
niche.geo.res.val.random.SHWH <- vector()
for(j in 1:1000){
	randomised.breeding.ranges <- sample(1:length(selectedSpeciesSH_WH))
	niche.distance.SHWH.rescaled.random <- niche.distance.SHWH.rescaled[,randomised.breeding.ranges]
	geo.distance.SHWH.rescaled.random <- geo.distance.SHWH.rescaled[,randomised.breeding.ranges]
	resources.scarcity.SHWH.rescaled.random <- resources.scarcity.SHWH.rescaled[,randomised.breeding.ranges]
	niche.geo.mat.random.SHWH <- niche.geo.mat.SHWH[,randomised.breeding.ranges]
	niche.res.mat.random.SHWH <- niche.res.mat.SHWH[,randomised.breeding.ranges]
	geo.res.mat.random.SHWH <- geo.res.mat.SHWH[,randomised.breeding.ranges]
	niche.geo.res.mat.random.SHWH <- niche.geo.res.mat.SHWH[,randomised.breeding.ranges]

	niche.val.random.SHWH[j] <- sum(diag(niche.distance.SHWH.rescaled.random))
	geo.val.random.SHWH[j] <- sum(diag(geo.distance.SHWH.rescaled.random))
	res.val.random.SHWH[j] <- sum(diag(resources.scarcity.SHWH.rescaled.random))
	niche.geo.val.random.SHWH[j] <- sum(diag(niche.geo.mat.random.SHWH))
	niche.res.val.random.SHWH[j] <- sum(diag(niche.res.mat.random.SHWH))
	geo.res.val.random.SHWH[j] <- sum(diag(geo.res.mat.random.SHWH))
	niche.geo.res.val.random.SHWH[j] <- sum(diag(niche.geo.res.mat.random.SHWH))
}







## NH EH


niche.distance.EH <- matrix(nrow=length(selectedSpeciesNH_EH)-2, ncol=length(selectedSpeciesNH_EH)-2)
resources.scarcity.EH <- matrix(nrow=length(selectedSpeciesNH_EH)-2, ncol=length(selectedSpeciesNH_EH)-2)
geo.distance.EH <- matrix(nrow=length(selectedSpeciesNH_EH)-2, ncol=length(selectedSpeciesNH_EH)-2)
for(k in 1:(length(selectedSpeciesNH_EH)-2)){
	for(h in 1:(length(selectedSpeciesNH_EH)-2)){
		
		br <- cbind(TempNS_EH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[h]))]==1)], PrecNS_WH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[h]))]==1)])
		br_kernel <- kde2d(br[,1], br[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
		br_kernel_2_mat <- rep(br_kernel[[1]][1],50)
		for(i in 2:50){
			br_kernel_2_mat <- c(br_kernel_2_mat, rep(br_kernel[[1]][i],50))
		}
		br_kernel_2_mat2 <- rep(br_kernel[[2]],50)
		br_kernel_2_mat3 <- br_kernel[[3]][1,]
		for(i in 2:50){
			br_kernel_2_mat3 <- c(br_kernel_2_mat3, br_kernel[[3]][i,])
		}
		br_kernel_mat <- as.data.frame(cbind(br_kernel_2_mat, br_kernel_2_mat2, br_kernel_2_mat3))
		breeding.niche.raster <- raster(ncol=50, nrow=50)
		extent(breeding.niche.raster) <- extent(c(-3, 3, -3,3))
		breeding.niche.raster <- rasterize(br_kernel_mat[,1:2], breeding.niche.raster, br_kernel_mat[,3])
		breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))
		thres = 0
		i=0
		while(thres <= 0.99){
			i = i+1
			thres = thres + sort(as.vector(breeding.niche.raster), decreasing=T)[i]
		}
		breeding.niche.raster[which(as.vector(breeding.niche.raster) < sort(as.vector(breeding.niche.raster), decreasing=T)[i])] = 0
		breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))

		nb <- cbind(TempNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[k]))]==1)], PrecNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[k]))]==1)])
		nb_kernel <- kde2d(nb[,1], nb[,2], h=0.5, n=50, lims=c(-3, 3, -3,3))
		nb_kernel_2_mat <- rep(nb_kernel[[1]][1],50)
		for(i in 2:50){
			nb_kernel_2_mat <- c(nb_kernel_2_mat, rep(nb_kernel[[1]][i],50))
		}
		nb_kernel_2_mat2 <- rep(nb_kernel[[2]],50)
		nb_kernel_2_mat3 <- nb_kernel[[3]][1,]
		for(i in 2:50){
			nb_kernel_2_mat3 <- c(nb_kernel_2_mat3, nb_kernel[[3]][i,])
		}
		nb_kernel_mat <- as.data.frame(cbind(nb_kernel_2_mat, nb_kernel_2_mat2, nb_kernel_2_mat3))
		nonbreeding.niche.raster <- raster(ncol=50, nrow=50)
		extent(nonbreeding.niche.raster) <- extent(c(-3, 3, -3,3))
		nonbreeding.niche.raster <- rasterize(nb_kernel_mat[,1:2], nonbreeding.niche.raster, nb_kernel_mat[,3])
		nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))
		thres = 0
		i=0
		while(thres <= 0.99){
			i = i+1
			thres = thres + sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i]
		}
		nonbreeding.niche.raster[which(as.vector(nonbreeding.niche.raster) < sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i])] = 0
		nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))

		niche.distance.EH[k,h] <- emd(breeding.niche.raster, nonbreeding.niche.raster, threshold=2)


		geo.distance.EH[k,h] <- rdist.earth(t(as.matrix(c(mean(east_Hem[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[k]))]==1),1]), mean(east_Hem[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[k]))]==1),2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[h]))]==1),1]), mean(east_Hem[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[h]))]==1),2])))), miles = F)
		resources.scarcity.EH[k,h] <-  mean(resourcesScarcity_januaryWH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[h]))]==1)]) + mean(resourcesScarcity_julyWH[which(PresAbs_BR_NH_EH[,which(colnames(PresAbs_BR_NH_EH) == names(selectedSpeciesNH_EH[k]))]==1)])
	}
}

niche.distance.EH.rescaled <- apply(niche.distance.EH, c(1:2), function(x) (x-min(niche.distance.EH))/(max(niche.distance.EH) - min(niche.distance.EH)))
geo.distance.EH.rescaled <- apply(geo.distance.EH, c(1:2), function(x) (x-min(geo.distance.EH))/(max(geo.distance.EH) - min(geo.distance.EH)))
resources.scarcity.EH.rescaled <- apply(resources.scarcity.EH, c(1:2), function(x) (x-min(resources.scarcity.EH))/(max(resources.scarcity.EH) - min(resources.scarcity)))
niche.geo.mat.EH <- niche.distance.EH.rescaled + geo.distance.EH.rescaled
niche.res.mat.EH <- niche.distance.EH.rescaled + resources.scarcity.EH.rescaled
geo.res.mat.EH <- geo.distance.EH.rescaled + resources.scarcity.EH.rescaled
niche.geo.res.mat.EH <- niche.distance.EH.rescaled + geo.distance.EH.rescaled + resources.scarcity.EH.rescaled

niche.val.obs.EH <- sum(diag(niche.distance.EH.rescaled))
geo.val.obs.EH <- sum(diag(geo.distance.EH.rescaled))
res.val.obs.EH <- sum(diag(resources.scarcity.EH.rescaled))
niche.geo.val.obs.EH <- sum(diag(niche.geo.mat.EH))
niche.res.val.obs.EH <- sum(diag(niche.res.mat.EH))
geo.res.val.obs.EH <- sum(diag(geo.res.mat.EH))
niche.geo.res.val.obs.EH <- sum(diag(niche.geo.res.mat.EH))

niche.val.random.EH <- vector()
geo.val.random.EH <- vector()
res.val.random.EH <- vector()
niche.geo.val.random.EH <- vector()
niche.res.val.random.EH <- vector()
geo.res.val.random.EH <- vector()
niche.geo.res.val.random.EH <- vector()
for(j in 1:1000){
	randomised.breeding.ranges <- sample(1:(length(selectedSpeciesNH_EH)-2))
	niche.distance.EH.rescaled.random <- niche.distance.EH.rescaled[,randomised.breeding.ranges]
	geo.distance.EH.rescaled.random <- geo.distance.EH.rescaled[,randomised.breeding.ranges]
	resources.scarcity.EH.rescaled.random <- resources.scarcity.EH.rescaled[,randomised.breeding.ranges]
	niche.geo.mat.random.EH <- niche.geo.mat.EH[,randomised.breeding.ranges]
	niche.res.mat.random.EH <- niche.res.mat.EH[,randomised.breeding.ranges]
	geo.res.mat.random.EH <- geo.res.mat.EH[,randomised.breeding.ranges]
	niche.geo.res.mat.random.EH <- niche.geo.res.mat.EH[,randomised.breeding.ranges]

	niche.val.random.EH[j] <- sum(diag(niche.distance.EH.rescaled.random))
	geo.val.random.EH[j] <- sum(diag(geo.distance.EH.rescaled.random))
	res.val.random.EH[j] <- sum(diag(resources.scarcity.EH.rescaled.random))
	niche.geo.val.random.EH[j] <- sum(diag(niche.geo.mat.random.EH))
	niche.res.val.random.EH[j] <- sum(diag(niche.res.mat.random.EH))
	geo.res.val.random.EH[j] <- sum(diag(geo.res.mat.random.EH))
	niche.geo.res.val.random.EH[j] <- sum(diag(niche.geo.res.mat.random.EH))
}


## SH EH

niche.distance.SHEH <- matrix(nrow=length(selectedSpeciesSH_EH), ncol=length(selectedSpeciesSH_EH))
resources.scarcity.SHEH <- matrix(nrow=length(selectedSpeciesSH_EH), ncol=length(selectedSpeciesSH_EH))
geo.distance.SHEH <- matrix(nrow=length(selectedSpeciesSH_EH), ncol=length(selectedSpeciesSH_EH))
for(k in 1:length(selectedSpeciesSH_EH)){
	for(h in 1:length(selectedSpeciesSH_EH)){
			
		br <- cbind(TempNW_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[h]))]==1)], PrecNW_EH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[h]))]==1)])
		br_kernel <- kde2d(br[,1], br[,2], h=0.5, n=50, lims=c(-3, 3, -3,3)) # compute the kernel
		br_kernel_2_mat <- rep(br_kernel[[1]][1],50)
		for(i in 2:50){
			br_kernel_2_mat <- c(br_kernel_2_mat, rep(br_kernel[[1]][i],50))
		}
		br_kernel_2_mat2 <- rep(br_kernel[[2]],50)
		br_kernel_2_mat3 <- br_kernel[[3]][1,]
		for(i in 2:50){
			br_kernel_2_mat3 <- c(br_kernel_2_mat3, br_kernel[[3]][i,])
		}
		br_kernel_mat <- as.data.frame(cbind(br_kernel_2_mat, br_kernel_2_mat2, br_kernel_2_mat3))
		breeding.niche.raster <- raster(ncol=50, nrow=50)
		extent(breeding.niche.raster) <- extent(c(-3, 3, -3,3))
		breeding.niche.raster <- rasterize(br_kernel_mat[,1:2], breeding.niche.raster, br_kernel_mat[,3])
		breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))
		thres = 0
		i=0
		while(thres <= 0.99){
			i = i+1
			thres = thres + sort(as.vector(breeding.niche.raster), decreasing=T)[i]
		}
		breeding.niche.raster[which(as.vector(breeding.niche.raster) < sort(as.vector(breeding.niche.raster), decreasing=T)[i])] = 0
		breeding.niche.raster = breeding.niche.raster / sum(as.vector(breeding.niche.raster))

		nb <- cbind(TempNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[k]))]==1)], PrecNW_EH[which(PresAbs_NB_NH_EH[,which(colnames(PresAbs_NB_NH_EH) == names(selectedSpeciesNH_EH[k]))]==1)])
		nb_kernel <- kde2d(nb[,1], nb[,2], h=0.5, n=50, lims=c(-3, 3, -3,3))
		nb_kernel_2_mat <- rep(nb_kernel[[1]][1],50)
		for(i in 2:50){
			nb_kernel_2_mat <- c(nb_kernel_2_mat, rep(nb_kernel[[1]][i],50))
		}
		nb_kernel_2_mat2 <- rep(nb_kernel[[2]],50)
		nb_kernel_2_mat3 <- nb_kernel[[3]][1,]
		for(i in 2:50){
			nb_kernel_2_mat3 <- c(nb_kernel_2_mat3, nb_kernel[[3]][i,])
		}
		nb_kernel_mat <- as.data.frame(cbind(nb_kernel_2_mat, nb_kernel_2_mat2, nb_kernel_2_mat3))
		nonbreeding.niche.raster <- raster(ncol=50, nrow=50)
		extent(nonbreeding.niche.raster) <- extent(c(-3, 3, -3,3))
		nonbreeding.niche.raster <- rasterize(nb_kernel_mat[,1:2], nonbreeding.niche.raster, nb_kernel_mat[,3])
		nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))
		thres = 0
		i=0
		while(thres <= 0.99){
			i = i+1
			thres = thres + sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i]
		}
		nonbreeding.niche.raster[which(as.vector(nonbreeding.niche.raster) < sort(as.vector(nonbreeding.niche.raster), decreasing=T)[i])] = 0
		nonbreeding.niche.raster = nonbreeding.niche.raster / sum(as.vector(nonbreeding.niche.raster))

		niche.distance.SHEH[k,h] <- emd(breeding.niche.raster, nonbreeding.niche.raster, threshold=2)


		geo.distance.SHEH[k,h] <- rdist.earth(t(as.matrix(c(mean(east_Hem[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[k]))]==1),1]), mean(east_Hem[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[k]))]==1),2])))),t(as.matrix(c(mean(east_Hem[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[h]))]==1),1]), mean(east_Hem[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[h]))]==1),2])))), miles = F)
		resources.scarcity.SHEH[k,h] <-  mean(resourcesScarcity_julyEH[which(PresAbs_NB_SH_EH[,which(colnames(PresAbs_NB_SH_EH) == names(selectedSpeciesSH_EH[h]))]==1)]) + mean(resourcesScarcity_januaryEH[which(PresAbs_BR_SH_EH[,which(colnames(PresAbs_BR_SH_EH) == names(selectedSpeciesSH_EH[k]))]==1)])
	}
print(k)
}

niche.distance.SHEH.rescaled <- apply(niche.distance.SHEH, c(1:2), function(x) (x-min(niche.distance.SHEH))/(max(niche.distance.SHEH) - min(niche.distance.SHEH)))
geo.distance.SHEH.rescaled <- apply(geo.distance.SHEH, c(1:2), function(x) (x-min(geo.distance.SHEH))/(max(geo.distance.SHEH) - min(geo.distance.SHEH)))
resources.scarcity.SHEH.rescaled <- apply(resources.scarcity.SHEH, c(1:2), function(x) (x-min(resources.scarcity.SHEH))/(max(resources.scarcity.SHEH) - min(resources.scarcity.SHEH)))
niche.geo.mat.SHEH <- niche.distance.SHEH.rescaled + geo.distance.SHEH.rescaled
niche.res.mat.SHEH <- niche.distance.SHEH.rescaled + resources.scarcity.SHEH.rescaled
geo.res.mat.SHEH <- geo.distance.SHEH.rescaled + resources.scarcity.SHEH.rescaled
niche.geo.res.mat.SHEH <- niche.distance.SHEH.rescaled + geo.distance.SHEH.rescaled + resources.scarcity.SHEH.rescaled

niche.val.obs.SHEH <- sum(diag(niche.distance.SHEH.rescaled))
geo.val.obs.SHEH <- sum(diag(geo.distance.SHEH.rescaled))
res.val.obs.SHEH <- sum(diag(resources.scarcity.SHEH.rescaled))
niche.geo.val.obs.SHEH <- sum(diag(niche.geo.mat.SHEH))
niche.res.val.obs.SHEH <- sum(diag(niche.res.mat.SHEH))
geo.res.val.obs.SHEH <- sum(diag(geo.res.mat.SHEH))
niche.geo.res.val.obs.SHEH <- sum(diag(niche.geo.res.mat.SHEH))

niche.val.random.SHEH <- vector()
geo.val.random.SHEH <- vector()
res.val.random.SHEH <- vector()
niche.geo.val.random.SHEH <- vector()
niche.res.val.random.SHEH <- vector()
geo.res.val.random.SHEH <- vector()
niche.geo.res.val.random.SHEH <- vector()
for(j in 1:1000){
	randomised.breeding.ranges <- sample(1:length(selectedSpeciesSH_EH))
	niche.distance.SHEH.rescaled.random <- niche.distance.SHEH.rescaled[,randomised.breeding.ranges]
	geo.distance.SHEH.rescaled.random <- geo.distance.SHEH.rescaled[,randomised.breeding.ranges]
	resources.scarcity.SHEH.rescaled.random <- resources.scarcity.SHEH.rescaled[,randomised.breeding.ranges]
	niche.geo.mat.random.SHEH <- niche.geo.mat.SHEH[,randomised.breeding.ranges]
	niche.res.mat.random.SHEH <- niche.res.mat.SHEH[,randomised.breeding.ranges]
	geo.res.mat.random.SHEH <- geo.res.mat.SHEH[,randomised.breeding.ranges]
	niche.geo.res.mat.random.SHEH <- niche.geo.res.mat.SHEH[,randomised.breeding.ranges]

	niche.val.random.SHEH[j] <- sum(diag(niche.distance.SHEH.rescaled.random))
	geo.val.random.SHEH[j] <- sum(diag(geo.distance.SHEH.rescaled.random))
	res.val.random.SHEH[j] <- sum(diag(resources.scarcity.SHEH.rescaled.random))
	niche.geo.val.random.SHEH[j] <- sum(diag(niche.geo.mat.random.SHEH))
	niche.res.val.random.SHEH[j] <- sum(diag(niche.res.mat.random.SHEH))
	geo.res.val.random.SHEH[j] <- sum(diag(geo.res.mat.random.SHEH))
	niche.geo.res.val.random.SHEH[j] <- sum(diag(niche.geo.res.mat.random.SHEH))
}
















###  FIGURES

# Figure 1

par(mfrow=c(1,2), mgp=c(2.2,0.8,0), mar=c(3.5,3,1.5,2))

plot(hexgridWH, col="grey", border = "grey")
plot(hexgridWH[match(PresAbs_BR_NH[which(PresAbs_BR_NH[,which(colnames(PresAbs_NB_NH) == "Dolichonyx.oryzivorus")] == 1),1], hexgridWH@data[,1]),], col="red3", border = "red3", add=T)
plot(hexgridWH[match(PresAbs_NB_NH[which(PresAbs_NB_NH[,which(colnames(PresAbs_NB_NH) == "Dolichonyx.oryzivorus")] == 1),1], hexgridWH@data[,1]),], col="blue2", border = "blue2", add=T)

centroidBR <- c(mean(hexgridWH@data[match(PresAbs_BR_NH[which(PresAbs_BR_NH[,which(colnames(PresAbs_BR_NH) == "Dolichonyx.oryzivorus")] == 1),1], hexgridWH@data[,1]),4]), mean(hexgridWH@data[match(PresAbs_BR_NH[which(PresAbs_BR_NH[,which(colnames(PresAbs_BR_NH) == "Dolichonyx.oryzivorus")] == 1),1], hexgridWH@data[,1]),3]))

centroidNB <- c(mean(hexgridWH@data[match(PresAbs_NB_NH[which(PresAbs_NB_NH[,which(colnames(PresAbs_NB_NH) == "Dolichonyx.oryzivorus")] == 1),1], hexgridWH@data[,1]),4]), mean(hexgridWH@data[match(PresAbs_NB_NH[which(PresAbs_NB_NH[,which(colnames(PresAbs_NB_NH) == "Dolichonyx.oryzivorus")] == 1),1], hexgridWH@data[,1]),3]))

inter <- gcIntermediate(c(centroidBR[1], centroidBR[2]), c(centroidNB[1], centroidNB[2]), n=10, addStartEnd=TRUE)
lines(inter, lwd=2)

mtext("A", cex=1.8, side=3, line=-1.5, at=-190)

temp_bobolink_BR <- TempNS[match(PresAbs_BR_NH[which(PresAbs_BR_NH[,which(colnames(PresAbs_BR_NH) == "Dolichonyx.oryzivorus")] == 1),1], hexid)]
prec_bobolink_BR <- PrecNS[match(PresAbs_BR_NH[which(PresAbs_BR_NH[,which(colnames(PresAbs_BR_NH) == "Dolichonyx.oryzivorus")] == 1),1], hexid)]

temp_bobolink_NB <- TempNW[match(PresAbs_NB_NH[which(PresAbs_NB_NH[,which(colnames(PresAbs_NB_NH) == "Dolichonyx.oryzivorus")] == 1),1], hexid)]
prec_bobolink_NB <- PrecNW[match(PresAbs_NB_NH[which(PresAbs_NB_NH[,which(colnames(PresAbs_NB_NH) == "Dolichonyx.oryzivorus")] == 1),1], hexid)]

temp_bobolink_resBR <- TempNW[match(PresAbs_BR_NH[which(PresAbs_BR_NH[,which(colnames(PresAbs_BR_NH) == "Dolichonyx.oryzivorus")] == 1),1], hexid)]
prec_bobolink_resBR <- PrecNW[match(PresAbs_BR_NH[which(PresAbs_BR_NH[,which(colnames(PresAbs_BR_NH) == "Dolichonyx.oryzivorus")] == 1),1], hexid)]

temp_bobolink_resNB <- TempNS[match(PresAbs_NB_NH[which(PresAbs_NB_NH[,which(colnames(PresAbs_NB_NH) == "Dolichonyx.oryzivorus")] == 1),1], hexid)]
prec_bobolink_resNB <- PrecNS[match(PresAbs_NB_NH[which(PresAbs_NB_NH[,which(colnames(PresAbs_NB_NH) == "Dolichonyx.oryzivorus")] == 1),1], hexid)]


plot(temp_bobolink_resBR, prec_bobolink_resBR, pch=20, xlim=c(-1.8,1.8), ylim=c(-1.5,1.5), col="light blue", xlab = "Temperature", ylab="log(Precipiation+1)")
points(temp_bobolink_resNB, prec_bobolink_resNB, pch=20, xlim=c(-2,2), ylim=c(-2,2), col="orange")
points(temp_bobolink_BR, prec_bobolink_BR, pch=20, xlim=c(-2,2), ylim=c(-2,2), col="red3")
points(temp_bobolink_NB, prec_bobolink_NB, pch=20, xlim=c(-2,2), ylim=c(-2,2), col="blue2")

centroidNicheBR <- c(mean(temp_bobolink_BR), mean(prec_bobolink_BR))
centroidNicheNB <- c(mean(temp_bobolink_NB), mean(prec_bobolink_NB))
centroidNicheResBR <- c(mean(temp_bobolink_resBR), mean(prec_bobolink_resBR))
centroidNicheResNB <- c(mean(temp_bobolink_resNB), mean(prec_bobolink_resNB))

points(centroidNicheBR[1], centroidNicheBR[2], pch=1)
points(centroidNicheBR[1], centroidNicheBR[2], pch=20, col="green")
points(centroidNicheNB[1], centroidNicheNB[2], pch=1)
points(centroidNicheNB[1], centroidNicheNB[2], pch=20, col="green")
points(centroidNicheResBR[1], centroidNicheResBR[2], pch=1)
points(centroidNicheResBR[1], centroidNicheResBR[2], pch=20, col="green")
points(centroidNicheResNB[1], centroidNicheResNB[2], pch=1)
points(centroidNicheResNB[1], centroidNicheResNB[2], pch=20, col="green")

lines(rbind(centroidNicheBR, centroidNicheNB), lwd=2)

mtext("B", cex=1.8, side=3, line=-1.5, at=-2.7)




# Figure 2

par(mfrow=c(1,2), mgp=c(3,1.5,0))

boxplot(nicheDistObs, nicheDistObs_resNB, nicheDistObs_resBR, names=c("Observed \nmigration", "If resident \nin NB", "If resident \nin BR"), ylab = "Niche distance", cex.lab=1.3)

mtext("A", cex=1.8, side=3, line=1, at=0.1)

boxplot((nicheDistObs - nicheDistObs_resNB), (nicheDistObs - nicheDistObs_resBR), names=c("Observed migration \n- if resident in NB", "Observed migration \n- if resident in BR"), ylab="Difference in niche distance", cex.lab=1.3)

mtext("B", cex=1.8, side=3, line=1, at=0.1)

t.test((nicheDistObs - nicheDistObs_resNB), alternative="less", mu=0)
t.test((nicheDistObs - nicheDistObs_resBR), alternative="less", mu=0)


# Figure 3

par(mar=c(2.5,2.5,0,2.5), mgp=c(1.5,0.5,0))
hist(migrationDistanceObs, xlim=c(0,120), xlab="", ylab="", main="", xaxt="n", axes=F, col="light grey", border="grey")
axis(side=4)
axis(side=1, at=c(0,20,40,60,80,100,120), labels = as.character(c(0,20,40,60,80,100,120)^2))
par(new=T, mar=c(2.5,2.5,0,2.5), mgp=c(1.5,0.5,0))
plot(migrationDistanceObs, resourcesObs, pch=20, ylab="Resource gain", xlab="Geographic distance (Km)", xlim=c(0,120), col="yellow", cex=1.2, axes=F)
points(migrationDistanceObs[which(nicheDistObs < 0.9)], resourcesObs[which(nicheDistObs < 0.9)], pch=20, xlim=c(0,120), col="orange", cex=1.2)
points(migrationDistanceObs[which(nicheDistObs < 0.6)], resourcesObs[which(nicheDistObs < 0.6)], pch=20, xlim=c(0,120), col="red", cex=1.2)
points(migrationDistanceObs[which(nicheDistObs < 0.3)], resourcesObs[which(nicheDistObs < 0.3)], pch=20, xlim=c(0,120), col="brown4", cex=1.2)
axis(side=2)
mtext("Number of species", side=4, line=1.5, cex.lab=1,las=0)
legend("topleft", inset=.05, title="Niche distance", c("< 0.3", "0.3–0.6", "0.6–0.9", "> 0.9"), fill=c("brown4", "red", "orange", "yellow"))

arrows(migrationDistanceObs[405]+2, resourcesObs[405]-4, migrationDistanceObs[405], resourcesObs[405], length=0.1, lwd=2, col="blue")


# Figure 4

## re-run model for j=182 (Bobolink)

hexidEH <- PresAbs_BR_NH_EH[,1]

par(mfrow=c(2,2), mar=c(2.5,2.5,0.2,1.5), mgp=c(1.5,0.5,0))

plot(hexgridEH, col="grey", border = "grey")
plot(hexgridEH[match(PresAbs_BR_NH[which(PresAbs_BR_NH[,which(colnames(PresAbs_NB_NH) == "Ixobrychus.eurhythmus")] == 1),1], hexgridEH@data[,1]),], col="red3", border = "red3", add=T)
plot(hexgridEH[match(hexidEH[range.sim_EH[[182]][[3]]], hexgridEH@data[,1]),], col="black", border = "black", add=T)
plot(hexgridEH[match(hexidEH[range.sim_EH[[182]][[17]]], hexgridEH@data[,1]),], col="black", border = "black", add=T)

mtext("A", cex=1.8, side=3, line=-1.5, at=-50)

plot(ggg, nnn, pch=20, axes=F, xlab="Geographic distance", ylab = "Niche distance", cex.lab=1.3, add=F)
 axis(1)
 axis(2)
 abline(a = nobs, b = 0, col="orange")
 abline(a = (nobs + gobs), b = -1, col="blue")
 points(gobs, nobs, pch=20, col="red3", cex=2.5)

mtext("C", cex=1.8, side=3, line=-1.5, at=-0.2)

plot(hexgridEH, col="grey", border = "grey")
plot(hexgridEH[match(PresAbs_NB_NH[which(PresAbs_NB_NH[,which(colnames(PresAbs_NB_NH) == "Ixobrychus.eurhythmus")] == 1),1], hexgridEH@data[,1]),], col="red3", border = "red3", add=T)
plot(hexgridEH[match(hexidEH[range.simNB_EH[[182]][[1]]], hexgridEH@data[,1]),], col="black", border = "black", add=T)
plot(hexgridEH[match(hexidEH[range.simNB_EH[[182]][[4]]], hexgridEH@data[,1]),], col="black", border = "black", add=T)

mtext("B", cex=1.8, side=3, line=-1.5, at=-50)

plot(nnnrrrggg, rnorm(length(nnnrrrggg), 1, 0.1), pch=20, axes=F, xlab="Niche distance + geographic distance + resource scarcity", ylab = "", ylim=c(0,2), xlim = c(0.5, 2.5), cex.lab=1.3)
axis(1)
abline(v=(nobs + robs + gobs), col="green")
points(nobs + robs + gobs, 1, pch=20, ylim=c(0,2), col="red3", cex=2.5)

mtext("D", cex=1.8, side=3, line=-1.5, at=0.1)



# Figure 5

par(mfrow=c(1,3), mar=c(4.5,3,1.5,0.5), mgp=c(1.8,0.5,0))
hist(c(rank_nicheEH_ranges, rank_nicheWH_ranges, rank_nicheSHEH_ranges, rank_nicheSHWH_ranges), main="", ylab="Number of species", xlab="Ranks for niche distance", xlim=c(0,200), cex.lab=1.3)
mtext("A", cex=1.5, side=3, line=-0.5, at=100)
par(new=FALSE, mar=c(4.5,2.5,1.5,0.5), mgp=c(3,0.5,0))
hist(c(rank_niche_geoEH_ranges, rank_niche_geoWH_ranges, rank_niche_geoSHEH_ranges, rank_niche_geoSHWH_ranges), main="", ylab="", xlab="Ranks for niche distance \n+ geographic distance", xlim=c(0,200), cex.lab=1.3)
mtext("B", cex=1.5, side=3, line=-0.5, at=100)
hist(c(rank_niche_resources_geoEH_ranges, rank_niche_resources_geoWH_ranges, rank_niche_resources_geoSHEH_ranges, rank_niche_resources_geoSHWH_ranges), main="", ylab="", xlab="Ranks for niche distance \n+ geographic distance + resource scarcity", xlim=c(0,200), cex.lab=1.3)
mtext("C", cex=1.5, side=3, line=-0.5, at=100)











####  Phylogenetic Signal ####

library(phytools)

outputs <- read.csv("/Users/mariussomveille/Desktop/PhD/Chapter 3 – niche/Submission process/Proceedings b/Revisions/Appendix1.csv")
list_phylo <- birdTrees2[[1]]$tip.label
list_phylo <- unlist(lapply(list_phylo, function(x) paste(strsplit(x, "_")[[1]][1], strsplit(x, "_")[[1]][2], sep=".")))
list_migrant <- as.character(outputs[,1])
list_phylo2 <- list_phylo[match(list_migrant, list_phylo)]
list_phylo3 <- list_phylo[-match(list_migrant, list_phylo)]
list_phylo3 <- unlist(lapply(list_phylo3, function(x) paste(strsplit(x, "[.]")[[1]][1], strsplit(x, "[.]")[[1]][2], sep="_")))

birdConsensusTree2 <- drop.tip(birdConsensusTree, list_phylo3)
birdConsensusTree2 <- compute.brlen(birdConsensusTree2, method="Grafen")

list_phylo4 <- unlist(lapply(birdConsensusTree2$tip.label, function(x) paste(strsplit(x, "_")[[1]][1], strsplit(x, "_")[[1]][2], sep="."))) 


sppdupl <- list_migrant[which(duplicated(list_migrant) == T)]

outputs.both <- as.vector(c(as.character(outputs[which(list_migrant == sppdupl[1]),][1,1]), "both", apply(outputs[which(list_migrant == sppdupl[1]), 3:15], 2, mean)))
outputs <- outputs[-which(list_migrant == sppdupl[1]),]
list_migrant <- list_migrant[-which(list_migrant == sppdupl[1])]
for(i in 2:length(sppdupl)){
outputs.both <- rbind(outputs.both, as.vector(c(as.character(outputs[which(list_migrant == sppdupl[i]),][1,1]), "both", apply(outputs[which(list_migrant == sppdupl[i]), 3:15], 2, mean))))
outputs <- outputs[-which(list_migrant == sppdupl[i]),]
list_migrant <- list_migrant[-which(list_migrant == sppdupl[i])]
}
row.names(outputs.both) <- NULL
outputs.both <- as.data.frame(outputs.both)
colnames(outputs.both) <- colnames(outputs)
outputs <- rbind(outputs, outputs.both)

list_migrant <- as.character(outputs[,1])
outputs <- outputs[match(list_phylo4, list_migrant),]

phylosig(birdConsensusTree2, as.numeric(outputs[,15]), test=T)

















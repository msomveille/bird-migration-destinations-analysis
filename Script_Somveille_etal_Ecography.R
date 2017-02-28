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
#selectedSpeciesNH_both <- names(PresAbs_BR_NH_EH2)[which(is.element(names(PresAbs_BR_NH_EH2), names(PresAbs_BR_NH_WH2)))]
#aa = apply(PresAbs_BR_NH_WH2[,match(selectedSpeciesNH_both, names(PresAbs_BR_NH_WH2))], 2, sum)
#bb = apply(PresAbs_NB_NH_WH2[,match(selectedSpeciesNH_both, names(PresAbs_NB_NH_WH2))], 2, sum)
#cc = apply(PresAbs_BR_NH_EH2[,match(selectedSpeciesNH_both, names(PresAbs_BR_NH_EH2))], 2, sum)
#dd = apply(PresAbs_NB_NH_EH2[,match(selectedSpeciesNH_both, names(PresAbs_NB_NH_EH2))], 2, sum)
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
#Hexagons ID
hexidWH <- PresAbs_BR_NH[match(envData2[which(envData2$LONGITUDE<=-30),1], PresAbs_BR_NH[,1]),1]
hexidEH <- PresAbs_BR_NH[match(envData2[which(envData2$LONGITUDE>-30),1], PresAbs_BR_NH[,1]),1]
hexgrid2 <- rbind(hexgridWH[match(hexidWH, hexgridWH@data[,1]),], hexgridEH[match(hexidEH, hexgridEH@data[,1]),])

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

#Kolmogorov-Smirnov tests of whether the distribution of ranks is left-skewed wompared to the uniform distribution
ks.test(ranks_niche, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_geo, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_resource, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_geo_niche, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_geo_resource, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_niche_resource, runif(100000,0,1), alternative="greater")$p.value
ks.test(ranks_geo_niche_resource, runif(100000,0,1), alternative="greater")$p.value




### Figure 6 - Weighted Richness Maps ###

ranks_geo_niche_resourceNH_WH <- vector()
for(i in 1:length(geo.distances.simulatedNH_WH)){
	geo_niche_resource <- scale(c(geo.distancesNH_WH[i], geo.distances.simulatedNH_WH[[i]])) + scale(c(niche.distancesNH_WH[i], niche.distances.simulatedNH_WH[[i]])) + scale(c(resource.scarcityNH_WH[i], resource.scarcity.simulatedNH_WH[[i]]))
	ranks_geo_niche_resourceNH_WH[i] <- length(which(geo_niche_resource[-1] < geo_niche_resource[1])) / length(geo_niche_resource[-1])	
}
ranks_geo_niche_resourceSH_WH <- vector()
for(i in 1:length(geo.distances.simulatedSH_WH)){
	geo_niche_resource <- scale(c(geo.distancesSH_WH[i], geo.distances.simulatedSH_WH[[i]])) + scale(c(niche.distancesSH_WH[i], niche.distances.simulatedSH_WH[[i]])) + scale(c(resource.scarcitySH_WH[i], resource.scarcity.simulatedSH_WH[[i]]))
	ranks_geo_niche_resourceSH_WH[i] <- length(which(geo_niche_resource[-1] < geo_niche_resource[1])) / length(geo_niche_resource[-1])	
}
ranks_geo_niche_resourceWH <- c(ranks_geo_niche_resourceNH_WH, ranks_geo_niche_resourceSH_WH)
ranks_geo_niche_resourceNH_EH <- vector()
for(i in 1:length(geo.distances.simulatedNH_EH)){
	geo_niche_resource <- scale(c(geo.distancesNH_EH[i], geo.distances.simulatedNH_EH[[i]])) + scale(c(niche.distancesNH_EH[i], niche.distances.simulatedNH_EH[[i]])) + scale(c(resource.scarcityNH_EH[i], resource.scarcity.simulatedNH_EH[[i]]))
	ranks_geo_niche_resourceNH_EH[i] <- length(which(geo_niche_resource[-1] < geo_niche_resource[1])) / length(geo_niche_resource[-1])	
}
ranks_geo_niche_resourceSH_EH <- vector()
for(i in 1:length(geo.distances.simulatedSH_EH)){
	geo_niche_resource <- scale(c(geo.distancesSH_EH[i], geo.distances.simulatedSH_EH[[i]])) + scale(c(niche.distancesSH_EH[i], niche.distances.simulatedSH_EH[[i]])) + scale(c(resource.scarcitySH_EH[i], resource.scarcity.simulatedSH_EH[[i]]))
	ranks_geo_niche_resourceSH_EH[i] <- length(which(geo_niche_resource[-1] < geo_niche_resource[1])) / length(geo_niche_resource[-1])	
}
ranks_geo_niche_resourceEH <- c(ranks_geo_niche_resourceNH_EH, ranks_geo_niche_resourceSH_EH)
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
	hex.weight_BR_WH[i] <- sum(ranks_geo_niche_resourceWH[c(unlist(lapply(pres.list_BR_NH_WH, function(x) is.element(hexidWH[i],x))), unlist(lapply(pres.list_BR_SH_WH, function(x) is.element(hexidWH[i],x))))]) / length(ranks_geo_niche_resourceWH[c(unlist(lapply(pres.list_BR_NH_WH, function(x) is.element(hexidWH[i],x))), unlist(lapply(pres.list_BR_SH_WH, function(x) is.element(hexidWH[i],x))))])
}
hex.weight_NB_WH <- vector()
for(i in 1:length(hexidWH)){
	hex.weight_NB_WH[i] <- sum(ranks_geo_niche_resourceWH[c(unlist(lapply(pres.list_NB_NH_WH, function(x) is.element(hexidWH[i],x))), unlist(lapply(pres.list_NB_SH_WH, function(x) is.element(hexidWH[i],x))))]) / length(ranks_geo_niche_resourceWH[c(unlist(lapply(pres.list_NB_NH_WH, function(x) is.element(hexidWH[i],x))), unlist(lapply(pres.list_NB_SH_WH, function(x) is.element(hexidWH[i],x))))])
}
hex.weight_BR_EH <- vector()
for(i in 1:length(hexidEH)){
	hex.weight_BR_EH[i] <- sum(ranks_geo_niche_resourceEH[c(unlist(lapply(pres.list_BR_NH_EH, function(x) is.element(hexidEH[i],x))), unlist(lapply(pres.list_BR_SH_EH, function(x) is.element(hexidEH[i],x))))]) / length(ranks_geo_niche_resourceEH[c(unlist(lapply(pres.list_BR_NH_EH, function(x) is.element(hexidEH[i],x))), unlist(lapply(pres.list_BR_SH_EH, function(x) is.element(hexidEH[i],x))))])
}
hex.weight_NB_EH <- vector()
for(i in 1:length(hexidEH)){
	hex.weight_NB_EH[i] <- sum(ranks_geo_niche_resourceEH[c(unlist(lapply(pres.list_NB_NH_EH, function(x) is.element(hexidEH[i],x))), unlist(lapply(pres.list_NB_SH_EH, function(x) is.element(hexidEH[i],x))))]) / length(ranks_geo_niche_resourceEH[c(unlist(lapply(pres.list_NB_NH_EH, function(x) is.element(hexidEH[i],x))), unlist(lapply(pres.list_NB_SH_EH, function(x) is.element(hexidEH[i],x))))])
}
hex.weightBR <- c(hex.weight_BR_WH, hex.weight_BR_EH)
hex.weightNB <- c(hex.weight_NB_WH, hex.weight_NB_EH)

#Plot weighted richness
par(mfrow=c(2,2), mar=c(0.1,0.1,1.2,0.1), mgp=c(1.5,0.5,0))
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(5)[as.numeric(cut(hex.weightBR, breaks=c(-0.1,0.1,0.2,0.3,0.4,1.1)))]
datcol[which(is.na(datcol)==T)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey")
mtext("A", cex=1.3, side=3, line=-2, at=-180)
mtext("Breeding season", cex=1.1, side=3, line=0.15, at=0)
legend("bottomleft", inset=.04, bg="grey", box.col="grey", title="Average rank\nof migrants", c("> 0.4","0.3–0.4", "0.2–0.3", "0.1–0.2", "0–0.1", "No species"), fill=c(rev(rbPal(5)),"white"), cex=0.8)
datcol <- rbPal(5)[as.numeric(cut(hex.weightNB, breaks=c(-0.1,0.1,0.2,0.3,0.4,1.1)))]
datcol[which(is.na(datcol)==T)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey")
mtext("B", cex=1.3, side=3, line=-2, at=-180)
mtext("Non-breeding season", cex=1.1, side=3, line=0.15, at=0)

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
mtext("C", cex=1.3, side=3, line=-2, at=-180)
legend("bottomleft", inset=.04, bg="grey", box.col="grey", title="Richness\nin migrants", c("> 100","75–100", "50–75", "25–50", "1–25", "0"), fill=c(rev(rbPal(5)),"white"), cex=0.8)
datcol <- rbPal(5)[as.numeric(cut(mnb, breaks=c(0.9,25,50,75,100,150)))]
datcol[which(is.na(datcol)==T)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey")
mtext("D", cex=1.3, side=3, line=-2, at=-180)



### Figure 7 - Weighted arrows maps  ###

centroids.breeding.grounds <- rbind(centroids.breeding.groundsNH_WH, centroids.breeding.groundsNH_EH, centroids.breeding.groundsSH_WH, centroids.breeding.groundsSH_EH)
centroids.nonbreeding.grounds <- rbind(centroids.nonbreeding.groundsNH_WH, centroids.nonbreeding.groundsNH_EH, centroids.nonbreeding.groundsSH_WH, centroids.nonbreeding.groundsSH_EH)
#Split the species into 3 groups: good rank, intermerdiate rank and bad rank
centroids.breeding.grounds_1 <- centroids.breeding.grounds[which(ranks_geo_niche_resource < 0.1),]
centroids.nonbreeding.grounds_1 <- centroids.nonbreeding.grounds[which(ranks_geo_niche_resource < 0.1),]
centroids.breeding.grounds_2 <- centroids.breeding.grounds[which(ranks_geo_niche_resource >= 0.1 & ranks_geo_niche_resource <= 0.5),]
centroids.nonbreeding.grounds_2 <- centroids.nonbreeding.grounds[which(ranks_geo_niche_resource >= 0.1 & ranks_geo_niche_resource <= 0.5),]
centroids.breeding.grounds_3 <- centroids.breeding.grounds[which(ranks_geo_niche_resource > 0.5),]
centroids.nonbreeding.grounds_3 <- centroids.nonbreeding.grounds[which(ranks_geo_niche_resource > 0.5),]
#Plot
par(mfrow=c(3,1), mar=c(0.1,0.1,0.1,0.1), mgp=c(1.5,0.5,0))
plot(hexgrid2, col= "grey", border = "grey")
points(centroids.breeding.grounds_1, pch=20, col="red", cex=0.3)
points(centroids.nonbreeding.grounds_1, pch=20, col="blue", cex=0.3)
for(i in 1:length(centroids.breeding.grounds_1[,1])){
	inter <- gcIntermediate(centroids.breeding.grounds_1[i,], centroids.nonbreeding.grounds_1[i,], n=50, addStartEnd=T)
	lines(inter, lwd=1, col="yellow")
}
mtext("A", cex=1.3, side=3, line=-2.2, at=-180)
mtext("Scaled rank < 0.1", cex=1, side=1, line=-2.2, at=25)
plot(hexgrid2, col= "grey", border = "grey")
points(centroids.breeding.grounds_2, pch=20, col="red", cex=0.3)
points(centroids.nonbreeding.grounds_2, pch=20, col="blue", cex=0.3)
for(i in 1:length(centroids.breeding.grounds_2[,1])){
	inter <- gcIntermediate(centroids.breeding.grounds_2[i,], centroids.nonbreeding.grounds_2[i,], n=50, addStartEnd=T)
	lines(inter, lwd=1, col="orange")
}
mtext("B", cex=1.3, side=3, line=-2.2, at=-180)
mtext("0.1 <= Scaled rank <= 0.5", cex=1, side=1, line=-2.2, at=25)
plot(hexgrid2, col= "grey", border = "grey")
points(centroids.breeding.grounds_3, pch=20, col="red", cex=0.3)
points(centroids.nonbreeding.grounds_3, pch=20, col="blue", cex=0.3)
for(i in 1:length(centroids.breeding.grounds_3[,1])){
	inter <- gcIntermediate(centroids.breeding.grounds_3[i,], centroids.nonbreeding.grounds_3[i,], n=50, addStartEnd=T)
	lines(inter, lwd=1, col="brown4")
}
mtext("C", cex=1.3, side=3, line=-2.2, at=-180)
mtext("Scaled rank > 0.5", cex=1, side=1, line=-2.2, at=25)




#Export excel table for appendix of publication

species.names <- names(c(selectedSpeciesNH_EH, selectedSpeciesNH_WH, selectedSpeciesSH_EH, selectedSpeciesSH_WH))
hem <- c(rep("EH", length(selectedSpeciesNH_EH)), rep("WH", length(selectedSpeciesNH_WH)), rep("EH", length(selectedSpeciesSH_EH)), rep("WH", length(selectedSpeciesSH_WH)))

xlsTable <- cbind(species.names, hem, geoDist, nicheDist, resources, geoDist_rescaled, nicheDist_rescaled, resources_rescaled, rank_geo, rank_niche, rank_resources, rank_niche_geo, rank_niche_resources, rank_resources_geo, rank_niche_geo_resources)
colnames(xlsTable) <- c("species name", "Longitudinal hemisphere", "Geographical distance", "Niche distance", "Resource scarcity", "Geographical distance rescaled", "Niche distance rescaled", "Resource scarcity rescaled", "Rank geo distance", "Rank niche distance", "Rank resource scarcity", "Rank niche distance + geo distance", "Rank niche distance + resource scarcity", "Rank resource scarcity + geo distance", "Rank niche distance + resource scarcity + geo distance")

xlsTable <- xlsTable[order(species.names),]

write.csv(xlsTable, "/Users/mariussomveille/Desktop/PhD/Chapter 3 – niche/Submission process/Proceedings b/Revisions/Appendix1.csv")














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

















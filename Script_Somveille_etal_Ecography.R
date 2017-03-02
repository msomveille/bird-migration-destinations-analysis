library(geosphere)
library(fields)
library(MASS)
library(raster)
library(move)
library(vioplot)
library(igraph)

setwd("~/Data")



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

#Testing if niche tracking is smaller than if resident
t.test(niche.distances.obs, niche.distances.residentinNB, "less")
t.test(niche.distances.obs, niche.distances.residentinBR, "less")


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
niche.distances.simulatedNH_WH <- list()
for(i in 1:length(nonbreeding.nichesNH_WH)){
	breeding.niches.simulated_NH_WH <- lapply(simulated.ranges_BR_NH_WH[[i]], function(x) nicheDensityRaster(cbind(TempNS_WH[x], PrecNS_WH[x])))
	niche.distances.simulatedBR_NH_WH <- mapply(function(X){ emd(X,nonbreeding.nichesNH_WH[[i]]) }, X=breeding.niches.simulated_NH_WH)
	nonbreeding.niches.simulated_NH_WH <- lapply(simulated.ranges_NB_NH_WH[[i]], function(x) nicheDensityRaster(cbind(TempNW_WH[x], PrecNW_WH[x])))
	niche.distances.simulatedNB_NH_WH <- mapply(function(Y){ emd(breeding.nichesNH_WH[[i]],Y) }, Y=nonbreeding.niches.simulated_NH_WH)
	niche.distances.simulatedNH_WH[[i]] <- c(unlist(niche.distances.simulatedBR_NH_WH), unlist(niche.distances.simulatedNB_NH_WH))
}
niche.distances.simulatedSH_WH <- list()
for(i in 1:length(nonbreeding.nichesSH_WH)){
	breeding.niches.simulated_SH_WH <- lapply(simulated.ranges_BR_SH_WH[[i]], function(x) nicheDensityRaster(cbind(TempNW_WH[x], PrecNW_WH[x])))
	niche.distances.simulatedBR_SH_WH <- mapply(function(X){ emd(X,nonbreeding.nichesSH_WH[[i]]) }, X=breeding.niches.simulated_SH_WH)
	nonbreeding.niches.simulated_SH_WH <- lapply(simulated.ranges_NB_SH_WH[[i]], function(x) nicheDensityRaster(cbind(TempNS_WH[x], PrecNS_WH[x])))
	niche.distances.simulatedNB_SH_WH <- mapply(function(Y){ emd(breeding.nichesSH_WH[[i]],Y) }, Y=nonbreeding.niches.simulated_SH_WH)
	niche.distances.simulatedSH_WH[[i]] <- c(unlist(niche.distances.simulatedBR_SH_WH), unlist(niche.distances.simulatedNB_SH_WH))
}
niche.distances.simulatedNH_EH <- list()
for(i in 1:length(nonbreeding.nichesNH_EH)){
	breeding.niches.simulated_NH_EH <- lapply(simulated.ranges_BR_NH_EH[[i]], function(x) nicheDensityRaster(cbind(TempNS_EH[x], PrecNS_EH[x])))
	niche.distances.simulatedBR_NH_EH <- mapply(function(X){ emd(X,nonbreeding.nichesNH_EH[[i]]) }, X=breeding.niches.simulated_NH_EH)
	nonbreeding.niches.simulated_NH_EH <- lapply(simulated.ranges_NB_NH_EH[[i]], function(x) nicheDensityRaster(cbind(TempNW_EH[x], PrecNW_EH[x])))
	niche.distances.simulatedNB_NH_EH <- mapply(function(Y){ emd(breeding.nichesNH_EH[[i]],Y) }, Y=nonbreeding.niches.simulated_NH_EH)
	niche.distances.simulatedNH_EH[[i]] <- c(unlist(niche.distances.simulatedBR_NH_EH), unlist(niche.distances.simulatedNB_NH_EH))
}
niche.distances.simulatedSH_EH <- list()
for(i in 1:length(nonbreeding.nichesSH_EH)){
	breeding.niches.simulated_SH_EH <- lapply(simulated.ranges_BR_SH_EH[[i]], function(x) nicheDensityRaster(cbind(TempNW_EH[x], PrecNW_EH[x])))
	niche.distances.simulatedBR_SH_EH <- mapply(function(X){ emd(X,nonbreeding.nichesSH_EH[[i]]) }, X=breeding.niches.simulated_SH_EH)
	nonbreeding.niches.simulated_SH_EH <- lapply(simulated.ranges_NB_SH_EH[[i]], function(x) nicheDensityRaster(cbind(TempNS_EH[x], PrecNS_EH[x])))
	niche.distances.simulatedNB_SH_EH <- mapply(function(Y){ emd(breeding.nichesSH_EH[[i]],Y) }, Y=nonbreeding.niches.simulated_SH_EH)
	niche.distances.simulatedSH_EH[[i]] <- c(unlist(niche.distances.simulatedBR_SH_EH), unlist(niche.distances.simulatedNB_SH_EH))
}
niche.distances.simulated <- c(niche.distances.simulatedNH_WH, niche.distances.simulatedNH_EH, niche.distances.simulatedSH_WH, niche.distances.simulatedSH_EH)
#load("simulated_niche_distances.RData")


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
plot(hexgridEH[match(hexidEH[simulated.ranges_BR_NH_EH[[181]][[3]]], hexgridEH@data[,1]),], col="black", border = "black", add=T)
plot(hexgridEH[match(hexidEH[simulated.ranges_BR_NH_EH[[181]][[8]]], hexgridEH@data[,1]),], col="black", border = "black", add=T)
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
plot(hexgridEH[match(hexidEH[simulated.ranges_NB_NH_EH[[181]][[7]]], hexgridEH@data[,1]),], col="black", border = "black", add=T)
plot(hexgridEH[match(hexidEH[simulated.ranges_NB_NH_EH[[181]][[70]]], hexgridEH@data[,1]),], col="black", border = "black", add=T)
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
plot(hexgrid2, col= "dark grey", border = "dark grey", bg="light grey")
points(centroids.breeding.grounds_1, pch=20, col="red", cex=0.3)
points(centroids.nonbreeding.grounds_1, pch=20, col="blue", cex=0.3)
for(i in 1:length(centroids.breeding.grounds_1[,1])){
	inter <- gcIntermediate(centroids.breeding.grounds_1[i,], centroids.nonbreeding.grounds_1[i,], n=50, addStartEnd=T)
	lines(inter, lwd=1, col="yellow")
}
mtext("A", cex=1.3, side=3, line=-2.2, at=-180)
mtext("Scaled rank < 0.1", cex=1, side=1, line=-2.2, at=25)
plot(hexgrid2, col= "dark grey", border = "dark grey", bg="light grey")
points(centroids.breeding.grounds_2, pch=20, col="red", cex=0.3)
points(centroids.nonbreeding.grounds_2, pch=20, col="blue", cex=0.3)
for(i in 1:length(centroids.breeding.grounds_2[,1])){
	inter <- gcIntermediate(centroids.breeding.grounds_2[i,], centroids.nonbreeding.grounds_2[i,], n=50, addStartEnd=T)
	lines(inter, lwd=1, col="orange")
}
mtext("B", cex=1.3, side=3, line=-2.2, at=-180)
mtext("0.1 <= Scaled rank <= 0.5", cex=1, side=1, line=-2.2, at=25)
plot(hexgrid2, col= "dark grey", border = "dark grey", bg="light grey")
points(centroids.breeding.grounds_3, pch=20, col="red", cex=0.3)
points(centroids.nonbreeding.grounds_3, pch=20, col="blue", cex=0.3)
for(i in 1:length(centroids.breeding.grounds_3[,1])){
	inter <- gcIntermediate(centroids.breeding.grounds_3[i,], centroids.nonbreeding.grounds_3[i,], n=50, addStartEnd=T)
	lines(inter, lwd=1, col="brown4")
}
mtext("C", cex=1.3, side=3, line=-2.2, at=-180)
mtext("Scaled rank > 0.5", cex=1, side=1, line=-2.2, at=25)



#Export excel table for appendix of publication
species.names <- c(names(PresAbs_BR_NH_WH2), names(PresAbs_BR_NH_EH2), names(PresAbs_BR_SH_WH2), names(PresAbs_BR_SH_EH2))
longitudinal.hemisphere <- c(rep("WH", length(PresAbs_BR_NH_WH2[1,])), rep("EH", length(PresAbs_BR_NH_EH2[1,])), rep("WH", length(PresAbs_BR_SH_WH2[1,])), rep("EH", length(PresAbs_BR_SH_EH2[1,])))
xlsTable <- cbind(species.names, longitudinal.hemisphere, geo.distances.obs, niche.distances.obs, resource.scarcity.obs, ranks_geo, ranks_niche, ranks_resource, ranks_geo_niche, ranks_niche_resource, ranks_geo_resource, ranks_geo_niche_resource)
colnames(xlsTable) <- c("species name", "Longitudinal hemisphere", "Geographical distance", "Niche distance", "Resource scarcity", "Scaled rank geo distance", "Scaled rank niche distance", "Scaled rank resource scarcity", "Scaled rank geo distance + niche distance", "Scaled rank niche distance + resource scarcity", "Scaled rank geo distance + resource scarcity", "Scaled rank geo distance + niche distance + resource scarcity")
xlsTable <- xlsTable[order(species.names),]
write.csv(xlsTable, "Appendix1.csv")



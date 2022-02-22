### run the megadata file ### 
library(tidyverse)
library(raster)
library(sf)
library(patchwork)

# Load the fish catch distribution data 

megadata <- readRDS('MegaData.rds') 
Cleanmegacell<-readRDS(file = "Cleanmegacell_mollweide.rds")
CleanCoordmegacell<-readRDS(file = "CleanCoordmegacell_mollweide.rds")
KprotectedPerCell_Library<-readRDS(file = "KprotectedPerCell_Library_mollweide.rds")
MPA_coord<-readRDS(file="MPA_coord_mollweide.rds")
land_shp_moll<-readRDS(file = "land_shp_moll.rds")

# This is a dataframe of fish stock assessment output from the North Sea (grabbed somewhere from ICES)
df <- read.csv('TEST.csv')


# Test e.g, cod 
fish <- megadata[grep('Gadus morhua', megadata$SciName),]

# A range of different cod stocks 
stocknames <- unique(fish$stockid)
stocknames
# Print all the species names 


# Plot the distribution 

ShortCoord<-CleanCoordmegacell # This data frame has the location of the fish
ShortCoord$ID<-row.names(ShortCoord)


ShortCoord$rank <- Cleanmegacell[,which(names(Cleanmegacell)== stocknames[20])] # 
   
dist <- ShortCoord %>% ggplot(aes(x=lon,y=lat, fill = rank)) + scale_fill_gradient(low="white", high="blue",name="Abundance")+
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
    geom_raster()+ geom_tile(data=MPA_coord, aes(x=lon, y=lat),fill="#EEA47FFF")+
    geom_sf(data = land_shp_moll, inherit.aes = F)#+coord_sf(ylim = c(2726460.23, 2876760.23), xlim = c(8678495,  8554895), crs = prj)
  
dist

codmsy <- 97777 # from the ices stock assessment 

# Plot the expected MSY in relation to the biomass in the stock assessment 
fish$MSYfin[fish$stockid == stocknames[20]]/codmsy  # The projected MSY is 3x higher than the actual MSY for north sea cod 

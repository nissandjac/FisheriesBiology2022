### run the megadata file ### 
library(tidyverse)
library(raster)
library(sf)

megadata <- readRDS('MegaData.rds') # Run the data from the MPA paper 
Cleanmegacell<-readRDS(file = "Cleanmegacell_mollweide.rds")
CleanCoordmegacell<-readRDS(file = "CleanCoordmegacell_mollweide.rds")
KprotectedPerCell_Library<-readRDS(file = "KprotectedPerCell_Library_mollweide.rds")
MPA_coord<-readRDS(file="MPA_coord_mollweide.rds")
land_shp_moll<-readRDS(file = "land_shp_moll.rds")


df <- read.csv('TEST.csv')
idx <- names(megadata) %in% c('SpeciesID', 'stockid', 'SciName',
  'K', 'B', 'MSYfin', 'Kfin', 'Hopt','Efin', 'E', 'r', 'm', 'Emsy')

megadata <- megadata[,idx] # Clean up the rows


# Test e.g, cod 
cod <- megadata[megadata$SciName == 'Gadus morhua',]

# Find the lat long of the cod species (try spp 1)

#cod.df <- Cleanmegacell[,which(names(Cleanmegacell) == 'CODNS')]

#$cod_MPA<-left_join(cod.df,MPA_coord,by=c("lon","lat"))


# See if it is even close to being realistic 

ShortCoord<-CleanCoordmegacell
ShortCoord$ID<-row.names(ShortCoord)
ShortCoord$rank <- Cleanmegacell$CODNS

dim(ShortCoord)
head(ShortCoord)
# Remember that lat long has to be transformed 
# xmin <- c(50,60)
prj <- "+proj=moll"
# w <- spTransform(xmin, CRS = prj)

GlobalMap<- ShortCoord %>% ggplot(aes(x=lon,y=lat, fill = rank)) + scale_fill_gradient(low="white", high="#00539CFF",name="Rank")+
 # theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  geom_raster()+ geom_tile(data=MPA_coord, aes(x=lon, y=lat),fill="#EEA47FFF")+
  geom_sf(data = land_shp_moll, inherit.aes = F)+coord_sf(ylim = c(2726460.23, 2876760.23), xlim = c(8678495,  8554895), crs = prj)
GlobalMap



# Plot sprat 

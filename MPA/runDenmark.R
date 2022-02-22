### WHere does denmark get their stuff from 
library(countrycode)
library(tidyverse)
library(sf)

#Load MOLLWEIDE projected data
MegaData <- readRDS('MegaData.rds') 
Cleanmegacell<-readRDS(file = "Cleanmegacell_mollweide.rds")
CleanCoordmegacell<-readRDS(file = "CleanCoordmegacell_mollweide.rds")
KprotectedPerCell_Library<-readRDS(file = "KprotectedPerCell_Library_mollweide.rds")
MPA_coord<-readRDS(file="MPA_coord_mollweide.rds")
land_shp_moll<-readRDS(file = "land_shp_moll.rds")
head(MPA_coord)
dim(MPA_coord)

#get MPA positions
CleanCoordmegacell_MPA<-left_join(CleanCoordmegacell,MPA_coord,by=c("lon","lat"))
head(CleanCoordmegacell_MPA)
dim(CleanCoordmegacell_MPA)
sum(CleanCoordmegacell_MPA$MPA,na.rm=T)

#positions of 1s (MPAs)
MPAposition<-which(CleanCoordmegacell_MPA$MPA==1)
head(MPAposition)
length(MPAposition)#2931 --- 2.44% are MPAs
length(MPAposition)*100/dim(Cleanmegacell)[1]

##TRY new approach
numcell<-dim(Cleanmegacell)[1]
celltoiterateFULL<-1:numcell
MPAselect0<-matrix(0, nrow=numcell, ncol=1)
PriorityAreas<-c()
NetworkResult<-vector()

#Make MPAselect0==1 for MPAs
MPAselect0[MPAposition]<-1
head(MPAselect0)
sum(MPAselect0)

#remove MPAs from celltoiterateFULL
celltoiterateFULL<-celltoiterateFULL[-MPAposition]
celltoiterate<-celltoiterateFULL
ncell<-length(celltoiterate)

###Compute spillover---PIXEL-LEVEL spillover 
K<-MegaData$Kfin # K per species (stock)
m<-MegaData$m # mobility per species
r<-MegaData$r

E<-MegaData$Efin_BAU1
ER<-1-E
ER<-1*(ER>1) + ER*(ER<=1)
max(ER)
min(ER)

MPAselect<-MPAselect0
R<-rowSums(KprotectedPerCell_Library[,which(MPAselect==1),drop=FALSE])
ER_redistribute<-1-(1-ER)^(1/(1-R))


b_outside_bau <-
  ((m * K * (1 - R)) / (ER_redistribute * R + m)) * (1 - (ER_redistribute * (1 - R) * m) /
                                                       ((ER_redistribute + m) * r)) 

hbau<-na.omit(ER_redistribute*((m*K*(1-R))/((ER_redistribute*R)+m))*(1-((ER_redistribute*(1-R)*m)/(((ER_redistribute*R)+m)*r))))
hbau<-hbau*(hbau>0)
HBAU<-sum(hbau)
HBAU

PICKSIZE<-100 #number of MPA sites selected

nmax<-floor(length(celltoiterate)/PICKSIZE)
nmax #this is the number of iterations needed for PICKSIZE at a time!

Eevolve<-matrix(nrow=nmax,ncol=dim(MegaData)[1])
# k_per_cell <- t(as.datKprotectedPerCell_Library)

hmpa_stockish <- matrix(NA, nrow = nmax, ncol = nrow(KprotectedPerCell_Library))
# 

dsn <- 'Shapefiles/EEZ_Land_v3_202030.shp'


# pull out case study cells
eezs <- sf::read_sf(dsn) %>% mutate(iso3_code = countrycode(SOVEREIGN1, "country.name",  "iso3c")) %>% 
  mutate(iso3c_name = countrycode(iso3_code, "iso3c",  "country.name")) %>%
  filter(!is.na(iso3_code)) %>%
  sf::st_transform(sf::st_crs(land_shp_moll))


teleporting_eezs <- eezs %>% 
  filter(iso3c_name %in% c("Australia"))


cabral_cells <- CleanCoordmegacell %>%
  sf::st_as_sf(coords = c("lon", "lat"),
               crs = sf::st_crs(land_shp_moll)) %>% 
  mutate(cell = 1:nrow(.))

spcoords <- sf::st_coordinates(cabral_cells) %>% 
  as_tibble()

cabral_cells <- cabral_cells %>% 
  bind_cols(spcoords) 

teleport_cells <- cabral_cells %>% 
  sf::st_intersection(teleporting_eezs)


teleport_cells %>%
  ggplot() +
  geom_sf(
    data = land_shp_moll,
    fill = "black",
    lwd = 0,
    inherit.aes = F
  ) + theme(panel.grid.major = element_line(colour = 'transparent')) +
  geom_sf(color = "red")


celltoiterateFULL <- teleport_cells$cell

celltoiterate <- celltoiterateFULL

PICKSIZE<-100 #number of MPA sites selected

nmax<-floor(length(celltoiterate) <- PICKSIZE)

nmax #this is the number of iterations needed for PICKSIZE at a time!

Eevolve<-matrix(nrow=nmax,ncol=dim(MegaData)[1])
# k_per_cell <- t(as.datKprotectedPerCell_Library)

hmpa_stockish <- matrix(NA, nrow = nmax, ncol = nrow(KprotectedPerCell_Library))

  
NetworkResult <-readRDS(file = "teleport_NetworkResult100_BAU1_mollweide.rds")
PriorityAreas <- readRDS(file = "teleport_PriorityAreas100_BAU1_mollweide.rds")
  

BenefitCurve<-as.data.frame(NetworkResult)/1000000
MPAsize<-(length(MPAposition)+1)*100/dim(Cleanmegacell)[1]

BenefitCurve$MPA <- seq(MPAsize,100, length.out = nrow(BenefitCurve))
zerozero<-data.frame(0,0)
names(zerozero)<-c("NetworkResult","MPA")
zerozero[1,]<-c(-HBAU/(1000000),100)
BenefitCurve<-rbind(BenefitCurve,zerozero)
benefitplot<-ggplot(BenefitCurve, aes(MPA, NetworkResult)) +geom_line(col="blue")+# theme_classic()+
  labs(x="% EEZ protected",y="Change in catch (MMT)",title=paste("Global (max =", round(max(BenefitCurve$NetworkResult),2),"MMT)"))#+geom_hline(yintercept = 0)
benefitplot

Priority<-as.data.frame(PriorityAreas)
Priority$rank <- ((seq.int(nrow(Priority))/ncell))*100

#add dh
NetworkResult_prime<-as.data.frame(NetworkResult)
dH_prime<-NetworkResult_prime %>% mutate(dH = NetworkResult - lag(NetworkResult, default = 0))
plot(dH_prime$dH)
Priority$NetworkResult<-rep(NetworkResult, each=PICKSIZE) ## this is adding delta H
Priority$dH<-rep(dH_prime$dH/PICKSIZE, each=PICKSIZE) ## this is adding delta H
head(Priority)
dim(Priority)

PriorityFrame <- as.data.frame(seq.int(nrow(CleanCoordmegacell)))
names(PriorityFrame) <- "PriorityAreas"
#add the priority
PriorityFrame2<-left_join(PriorityFrame,Priority, by="PriorityAreas")
PriorityFrame2[is.na(PriorityFrame2)] <- 0
head(PriorityFrame2)
dim(PriorityFrame2)

ShortCoord<-CleanCoordmegacell
ShortCoord$rank<-PriorityFrame2$rank
ShortCoord$NetworkResult<- PriorityFrame2$NetworkResult#this is for deriving dh
ShortCoord$dH<- PriorityFrame2$dH#this is for deriving dh
ShortCoord$ID<-row.names(ShortCoord)

dim(ShortCoord)
head(ShortCoord)

ShortCoord_Sort <- ShortCoord[order(-ShortCoord$rank),] %>% filter(rank!=-1000)

ShortCoord_Sort$sign<-ShortCoord_Sort$dH>0
InflectPoint<-min(which(ShortCoord_Sort$sign == TRUE))
InflectMPA<-ShortCoord_Sort$rank[InflectPoint]

MPAcoverage<-ShortCoord_Sort %>% filter(sign == TRUE) %>%  
  ggplot(aes(x=lon,y=lat,fill=dH)) + #scale_fill_gradient2(low="darkred", mid="white",high="#00539CFF",midpoint = 0, space = "Lab",na.value = "grey50", guide = "colourbar", aesthetics = "fill",name=expression(paste(Delta, "H (MT)")))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(), panel.background = element_blank(),legend.position=c(0.63, 0.05), legend.direction = "horizontal")+ #"bottom
  geom_raster()+
  geom_sf(data = land_shp_moll, fill="black", lwd = 0, inherit.aes = F)+ theme(panel.grid.major = element_line(colour = 'transparent'))
MPAcoverage


tunas <- MegaData$SciName %>%  str_detect("Thunnus")


total_dh <- rowSums(hmpa_stockish)

stockish_dh <- hmpa_stockish[total_dh == max(total_dh),]

percent_affected <- mean(abs(stockish_dh) > 0)

percent_benefited <- mean(stockish_dh > 0)

percent_hurt <- mean(stockish_dh < 0)

affected_stockish <- which(stockish_dh > 0 & tunas == FALSE) # to keep it simple just focus on stocks that get benefits from MPAs

check <- MegaData$SciName[affected_stockish]

tmp <- ShortCoord %>% 
  mutate(k_per_cell = as.numeric(KprotectedPerCell_Library[which(MegaData$SciName == check[1])[1],])) %>% 
  sf::st_as_sf(coords = c("lon", "lat"),
               crs = sf::st_crs(land_shp_moll)) %>% 
  filter(k_per_cell > 0) 

tmp %>% 
  ggplot(aes(color = k_per_cell)) + 
  geom_sf(data = land_shp_moll, fill="black", lwd = 0, inherit.aes = F)+ theme(panel.grid.major = element_line(colour = 'transparent')) +
  geom_sf(size = .1) + 
  scale_color_viridis_c()




teleport_k <- as.data.table(KprotectedPerCell_Library[affected_stockish,], keep.rownames = FALSE)

teleport_h <- stockish_dh[stockish_dh > 0 & tunas == FALSE]

test <- rowSums(teleport_h * teleport_k)

teleport_benefits <- colSums(teleport_h * teleport_k)

plot(test, teleport_h)
abline(a = 0, b = 1)

cabral_cells$teleport_benefits <- teleport_benefits

cabral_cells$teleport_benefits[cabral_cells$teleport_benefits == 0] <- NA

cabral_cells <- cabral_cells %>% 
  filter(!is.na(teleport_benefits)) %>% 
  ungroup() %>% 
  mutate(mpa = cell %in% as.numeric(ShortCoord$ID[ShortCoord$dH > 0 & !is.na(ShortCoord$dH)])) 

cabral_mpa <- cabral_cells


cabral_cells <- cabral_cells %>% 
  filter(!mpa) %>% 
  mutate(scaled_benefits = teleport_benefits / max(teleport_benefits))



teleport_map <-cabral_cells %>% filter(scaled_benefits > 0.01) %>% ggplot() + 
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(), panel.background = element_blank(),legend.position= "bottom", legend.direction = "horizontal")+ #"bottom
  geom_sf(data = land_shp_moll, fill="black", lwd = 0, inherit.aes = F)+ theme(panel.grid.major = element_line(colour = 'transparent')) + 
  geom_raster(aes(x = X, y = Y,fill=scaled_benefits))+
  geom_raster(data = cabral_mpa %>% filter(mpa),aes(x = X, y = Y), fill = "red") +
  scale_fill_viridis_c(name = "Benefits / max(Benefits)",guide = guide_colorbar(barwidth = unit(15, "lines")))

compact_teleport_map <-cabral_cells %>% filter(scaled_benefits > 0.01) %>% ggplot() + 
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(), panel.background = element_blank(),legend.position= "bottom", legend.direction = "horizontal")+ #"bottom
  geom_sf(data = land_shp_moll, fill="black", lwd = 0, inherit.aes = F)+ theme(panel.grid.major = element_line(colour = 'transparent')) + 
  geom_raster(aes(x = X, y = Y,fill=scaled_benefits))+
  geom_raster(data = cabral_mpa %>% filter(mpa),aes(x = X, y = Y), fill = "red") +
  scale_fill_viridis_c(name = "Benefits / max(Benefits)")


ggsave(filename = file.path(results_path, "teleport-map.pdf"), teleport_map, height = 10, width = 20)

ggsave(filename = file.path(results_path, "teleport-map.png"), teleport_map, height = 10, width = 20)
```
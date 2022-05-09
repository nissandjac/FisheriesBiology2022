##################
###  R Script: Traits diversity of fish communities
###            Functional Richness and null model
################
# Original code prepared by: Laurene Pecuchet (Thanks a lot!)
#######################################
# Please discuss and answer the 10 questions highlighted throughout the code!
#######################################

## Install packages
# install.packages("FD")
# install.packages("vegan")
# install.packages("reshape2")
# install.packages("FactoMineR")
# install.packages("mapplots")
# install.packages("maps")
# install.packages("mapdata")
# install.packages("mgcv")
# install.packages("qqplot2")

# Load packages
library(FD) ## Calculate functional diversity indices
library(vegan) ## Community ecology package
library(reshape2) ## Rearrange data into matrix
library(FactoMineR)  ## Factor analysis for mixed data

## R library for plotting results
library(mapplots)
library(maps)
library(mapdata)
library(mgcv)
library(ggplot2)
rbPal <- colorRampPalette(c('blue','red'))

##### Community data ####
##### Species presence/abundance per area
##### Example: Mean Species abundance (catch per unit effort) in each ICES squares in the Baltic Sea, 
# for an aggragated period 2003-2014
communities<-read.table("communities.csv", header=T)
head(communities)

#################################
# 1. What are the columns of the community data table? 
# What sources of bias could affect the results (if not standardised as presently done)? 

# 
# Species - mean length (I assumed meaned over hausl. Is selectivity taken into account?) Haul sum: Maybe number of hauls
# Cpue mean average cpue over hauls, biomass = cpue * weight. 


#################################

##### Traits data #### Diet/Mean length/Age at 50% maturity/Fecundity/Shape of caudal fin (proxy for activity)/body shape
traits<-read.table("Species_Traits.csv", header=T)
head(traits)

#################################
# 2. What are the columns of the trait matrix and what do they represent? 
# What sources of bias could one expect in this type of data?

# Large interspecies variety. Lmean and a50 is not a trait it's a life history parameter. Fecundity and shapes are good traits, but 
# may rarely describe stuff. Diet can be a trait.  



#############################

#### Calculate Functional diversity indices ####
### From dbFD function in package FD
### Distance-Based Functional Diversity Indices
# Inputs data: traits data with species rows and traits columns   or   speciesxspecies distance matrix
# Inputs data: abundances of the species in each communities
# !! All species in the communities matrix need to be in the traits data
# Optional inputs data: a vector giving different weights to the traits
# FDiversity common calculation

## First, the communities dataframe have to be in a matrix form: row sites, columns species
communities.mat<-dcast(communities, ICES~Species, fun=mean)  ### ICES, Species and biomass data rearrange in a matrix
head(communities.mat)
rownames(communities.mat) <- communities.mat[,1]  ## sites names (ICES squares) as rowname
communities.mat[,1] <- NULL
communities.mat[is.na(communities.mat)] <- 0  ### When Na == equals 0 abundance
communities.mat<-communities.mat[ , order(names(communities.mat))] ### Need to be in the same order

### Good, now we have traitsxspecies and speciesxsites matrix, we can calculate dbFD
?dbFD  ## look at the function and all its option
#convert to function
traits$fecundity <- as.numeric(traits$fecundity)
traits$Diet <- as.factor(traits$Diet)
traits$caudal <- as.factor(traits$caudal)
traits$body <- as.factor(traits$body)

FD<-dbFD(traits, communities.mat, stand.FRic=T)  #stand.FRic=T, means FRic will be between 0 and 1
    ## traits can also be given weight, adding a w vector
    ## use CWM.type = c("all") if you want for categorical traits the different proportions
    ## CWM.type = c("dom") if you only want the dominant (with highest biomass) category

#### traits data were composed of different types: quantitatives, categorical -> gower distances with sqrt correction
FD  #### All results
FD$nbsp  ### Species Richness              -  presence/absence
FD$FRic  ### Functional Richness           -  presence/absence
FD$FEve  ### Functional Eveness            -  biomass
FD$FDiv  ### Functional Divergence         -  biomass
FD$RaoQ  ### Functional Entropy            -  biomass
FD$CWM   ### Community weighted mean traits -  biomass
FD$CWM$Diet ## When catagorical: return the category with the highest biomass
FD$CWM$Lmean ## When continuous: return the weighted mean of the traits

FD<-do.call(cbind, FD)

#################################
# 3. What are the different indices and what do they mean?
#################################

############# Plotting data ####
## Get longitude and latitude of the ICES squares
FD[,c("lon","lat")]<-ices.rect(rownames(FD))

###### Or get back to your original lon, lat observations data


#### Plot 4 indices: e.g. Species Richness, Functional Richness, Functional Eveness, CWM.Lmean
par(mfrow = c(2, 2))
par(mar = c(2, 2, 2, 1), oma = c(0, 0, 0, 0.5))

### Species Richness
plot(c(9.5, 24), c(53, 61), type="n", xlab="", ylab="", main="Species Richness",  cex.axis=0.6)
map('worldHires', add=T, fill=T, col="grey", xlim=c(9.5, 24), ylim=c(53, 61))
box()
cuts<-levels(cut(FD$nbsp,breaks = 6, dig.lab=0))
cuts<-gsub(","," - ",cuts)
cuts<-gsub("\\(","",cuts)
cuts<-substr(cuts, 1, nchar(cuts)-1)
legend(x=13.25, y=59.25,cuts,col=rbPal(6),pch=16, cex=0.6, text.font = 2,pt.cex=1.1, bty="n")
with(FD, points(lon, lat, pch=21, cex=2.1, bg=rbPal(6)[as.numeric(cut(nbsp,breaks = 6))] ))

### Functional Richness
plot(c(9.5, 24), c(53, 61), type="n", xlab="", ylab="", main="Functional Richness",  cex.axis=0.6)
map('worldHires', add=T, fill=T, col="grey", xlim=c(9.5, 24), ylim=c(53, 61))
box()
cuts<-levels(cut(FD$FRic,breaks = 6, dig.lab=0))
cuts<-gsub(","," - ",cuts)
cuts<-gsub("\\(","",cuts)
cuts<-substr(cuts, 1, nchar(cuts)-1)
legend(x=13.25, y=59.25,cuts,col=rbPal(6),pch=16, cex=0.6, text.font = 2,pt.cex=1.1, bty="n")
with(FD, points(lon, lat, pch=21, cex=2.1, bg=rbPal(6)[as.numeric(cut(FRic,breaks = 6))] ))

### Functional Eveness
plot(c(9.5, 24), c(53, 61), type="n", xlab="", ylab="", main="Functional Evenness",  cex.axis=0.6)
map('worldHires', add=T, fill=T, col="grey", xlim=c(9.5, 24), ylim=c(53, 61))
box()
cuts<-levels(cut(FD$FEve,breaks = 6, dig.lab=0))
cuts<-gsub(","," - ",cuts)
cuts<-gsub("\\(","",cuts)
cuts<-substr(cuts, 1, nchar(cuts)-1)
legend(x=13.25, y=59.25,cuts,col=rbPal(6),pch=16, cex=0.6, text.font = 2,pt.cex=1.1, bty="n")
with(FD, points(lon, lat, pch=21, cex=2.1, bg=rbPal(6)[as.numeric(cut(FEve,breaks = 6))] ))

### CWM Lmean
plot(c(9.5, 24), c(53, 61), type="n", xlab="", ylab="", main="Mean Length",  cex.axis=0.6)
map('worldHires', add=T, fill=T, col="grey", xlim=c(9.5, 24), ylim=c(53, 61))
box()
cuts<-levels(cut(FD$CWM.Lmean,breaks = 6))
cuts<-gsub(","," - ",cuts)
cuts<-gsub("\\(","",cuts)
cuts<-substr(cuts, 1, nchar(cuts)-1)
legend(x=13.25, y=59.25,cuts,col=rbPal(6),pch=16, cex=0.6, text.font = 2,pt.cex=1.1, bty="n")
with(FD, points(lon, lat, pch=21, cex=2.1, bg=rbPal(6)[as.numeric(cut(CWM.Lmean,breaks = 6))] ))

#################################
# 4. What are the main spatial patterns of the different indices?


# Most richness in western baltic in terms of species / functional. 
# Mean length highest in the eastern. Why? Usually they're smnaller in the more saline sea. Is this weighted by abundance? 
#################################


##### Environmental variables #####
Env<-read.table("Environment.csv", header=T, sep=",")
head(Env)

#################################
# 5. What are the environmental factors and how could these potentially affect community composition?
#################################

FD$ICES<-rownames(FD)

FD<-merge(FD, Env, by="ICES")
head(FD)

## Exemple: Explainning SR and FR by the environment

# Backward selection: start with all variables in the model and eliminate non significant variables, 
# or do not increase largely the model fit (AIC)

## Species Richness
gam1<-gam(nbsp~s(Depth, k=3)+s(Oxygen,k=3)+
            s(Depth.sd,k=3)+s(Salinity,k=3)+
            s(Temperature,k=3)+ 
            s(Habitat.dominance, k=3), na.action=na.exclude, data=FD,select=TRUE, gamma=1.4)
summary(gam1)
AIC(gam1)

gam1<-gam(nbsp~s(Depth, k=3)+
            s(Salinity,k=3)+
            s(Temperature,k=3)+ 
            s(Habitat.dominance, k=3), na.action=na.exclude, data=FD,select=TRUE, gamma=1.4)
summary(gam1)
AIC(gam1)

# Salinity alone
gam1<-gam(nbsp~s(Salinity,k=3), na.action=na.exclude, data=FD,select=TRUE, gamma=1.4)
summary(gam1)
AIC(gam1)
plot(gam1)
#################################
# 6. What does the shape of the relationship between species richness and salinity indicate?
#################################
AIC(gam1)

## Functional Richness
gam2<-gam(FRic~s(Depth, k=3)+s(Oxygen,k=3)+s(Depth.sd,k=3)+s(Salinity,k=3)+s(Temperature,k=3)+ s(Habitat.dominance, k=3), na.action=na.exclude, data=FD,select=TRUE, gamma=1.4)
summary(gam2)

gam2<-gam(FRic~s(Depth, k=3)+s(Salinity,k=3)+ s(Habitat.dominance, k=3), na.action=na.exclude, data=FD,select=TRUE, gamma=1.4)
summary(gam2)

# Salinity alone
gam2<-gam(FRic~s(Salinity,k=3), na.action=na.exclude, data=FD,select=TRUE, gamma=1.4)
summary(gam2)
plot(gam2)


### Null model for SR and FRic ####
## Back with the sites x species matrix, this time only presence/absence (1/0)
communities.mat<-dcast(communities, ICES~Species, mean)  ### ICES, Species and biomass data rearrange in a matrix
head(communities.mat)
rownames(communities.mat) <- communities.mat[,1]  ## sites names (ICES squares) as rowname
communities.mat[,1] <- NULL
communities.mat[!is.na(communities.mat)] <- 1
communities.mat[is.na(communities.mat)] <- 0  ### When Na == equals 0 abundance

FR <- dbFD(traits, communities.mat, stand.FRic=T) ## Calculate the "observed" communties Functional Richness
FR<-data.frame(SR=FR$nbsp, FRic=FR$FRic)  ## Store only SR and FR

#### Calculate the Functional richness of each resulting randomization (1000 randomized matrix)
# Let?s use only 10 iterations to save some time....

null<-permatswap(communities.mat, times=10)  ## To create the null model, we randomize the species x site matrix, 
                                     ## by accounting for row (sites) sums and columns (species) sums
for(i in 1:10){
  ex <- dbFD(traits, null$perm[[i]], stand.FRic=T)
  FR[,2+i]<-ex$FRic
  print(i)
}

#################################
# 7. What is the purpose of this null model?
#################################

### Caluculate: the mean value of the null model for each level of Species Richness
### The 50% and 95% quantiles distributons of the null model
FR$null<-rowMeans(FR[,3:12]) ## mean "null" value of FRic for each site, based on 100 randomization at each site
FR$null75<-apply(FR[,3:12], 1, function(x){quantile(x, probs = c(0.75))}) ### 75% quantiles distribution
FR$null25<-apply(FR[,3:12], 1, function(x){quantile(x, probs = c(0.25))}) ### 25% quantiles distribution
FR$null97<-apply(FR[,3:12], 1, function(x){quantile(x, probs = c(0.975))}) ## 97.5% quantiles distribution
FR$null2<-apply(FR[,3:12], 1, function(x){quantile(x, probs = c(0.025))})  ## 2.5% quantiles distribution

### Calculate the difference between the observed FRic and the one expected under random composition of the communities
FR$predas<-FR$FRic-FR$null
FR$predas.bio<-FR$FRic-FR$null75 ## Limiting similarity: observed-75%quantile 
FR[FR$predas.bio<0,]$predas.bio<-0  ## If negative, then observed in the 0-75th distribution, then 0. The rest points towards L.s.
FR$predas.biosign<-FR$FRic-FR$null97 
FR[FR$predas.biosign<0,]$predas.biosign<-0 ## If negative, then observed in the 0-95th distribution, then 0. The rest significant L.s.
FR$predas.env<-FR$FRic-FR$null25 ## Environmental filtering: observed-25%quantile 
FR[FR$predas.env>0,]$predas.env<-0  ## If positive, then observed in the 25-100th distribution, then 0. The rest points towards E.f.
FR$predas.envsign<-FR$FRic-FR$null2
FR[FR$predas.envsign>0,]$predas.envsign<-0 ## If positive, then observed in the 5-100th distribution, then 0. The rest points towards E.f.

#### Return longitude and Latitude of the ICES squares, to be able to map
FR[,c("lon","lat")]<-ices.rect(rownames(FR))

#### Model the Functional richness in function of species richness
gam.as1.sd.neg2<-gam(null2~s(SR, k=3),data=FR) ## 2.5% lowest FRic
summary(gam.as1.sd.neg2)

gam.as1.sd.neg25<-gam(null25~s(SR, k=3),data=FR) ## 25% lowest FRic
summary(gam.as1.sd.neg25)

gam.as.null.sd.pos75<-gam(null75~s(SR, k=3),data=FR) ## 75% highest FRic
summary(gam.as.null.sd.pos75)

gam.as.null.sd.pos97<-gam(null97~s(SR, k=3),data=FR) ## 97.5% highest FRic
summary(gam.as.null.sd.pos97)

gam.as1<-gam(null~s(SR, k=3),data=FR) ## Mean null FRic 
summary(gam.as1)

gam.as2<-gam(FRic~s(SR, k=3),data=FR) ## Observed FRic 
summary(gam.as2)

plot(gam.as2)
#################################
# 8. What does the saturating relationship between SR and FRic indicate?
#################################


######### Figure #############
blue<- colorRampPalette(c('#0072B2',"white"))
red <- colorRampPalette(c("white",'#D55E00'))

ggplot(FR, aes(SR)) + 
  xlab("Species Richness") +
  ylab("Functional Richness") + 
  scale_y_continuous(limits=c(-0.04, 1))+
  theme_bw() +
  theme(axis.text.x = element_text(colour="black",size=6,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=6,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=6,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=6,angle=90,hjust=.5,vjust=.5,face="plain")) +
  geom_line(aes(sort(FR$SR), sort((predict(gam.as1)))), colour="black", cex=0.8) + 
  geom_ribbon(aes(ymin=predict(gam.as1.sd.neg2), ymax=predict(gam.as.null.sd.pos97)), alpha=0.1) +
  geom_ribbon(aes(ymin=predict(gam.as1.sd.neg25), ymax=predict(gam.as.null.sd.pos75)), alpha=0.2) +
  geom_point(data=FR, aes(x=FR$SR, y=FR$FRic), colour="red", size=0.8) + 
  geom_line(data=FR, aes(sort(FR$SR), sort((predict(gam.as2)))), colour="red", cex=0.6)+
  annotate("text", label = "Limiting similarity", size=1.2, x = 18, y = 0.9) +
  annotate("text", label = "Environmental filtering", size=1.2, x = 26, y = 0.45) +
  annotate("text", label = "Neutral", size=1.2, x = 25.5, y = 0.778) +
  annotate("text", label = "a", size=3.3, x = 4, y = 0.95) +
  annotate("text", label = "2.5%", x = 34, y = 0.70, colour="grey30", size = 1) +
  annotate("text", label = "25%", x = 34, y = 0.77, colour="grey30", size = 1) +
  annotate("text", label = "75%", x = 34, y = 0.9, colour="grey30", size = 1) +
  annotate("text", label = "97.5%" ,x = 34, y = 0.97, colour="grey30", size = 1)

par(mfrow=c(1,1))
par(mar=c(2,2,2,1))
plot(c(9.5, 22.5), c(53.9, 59), type="n",  xlab="Longitude", ylab="Latitude", main="", cex.axis=0.5, mgp=c(2,0.5,0))
map('worldHires', add=T, fill=T, col="grey", xlim=c(9.5, 22.5), ylim=c(53.9, 59))
box()
#draw.rect()
with(FR, points(lon, lat, pch=21, cex=3, bg="white" ))
with(subset(FR, predas.bio>0), points(lon, lat, pch=21, cex=3, bg=red(4)[as.numeric(cut(predas,breaks = 4))] ))
with(subset(FR, predas.env<0), points(lon, lat, pch=21, cex=3, bg=blue(4)[as.numeric(cut(predas,breaks = 4))] ))
with(subset(FR, predas.envsign<0), points(lon, lat, pch=3, cex=1,col="black"))
with(subset(FR, predas.biosign>0), points(lon, lat, pch=3, cex=1,col="black"))
legend(x=11.26, y=59.2,c("Environmental filtering","Neutral","Limiting similarity"),col=c("#0072B2","white","#D55E00"),pch=16, cex=0.35, text.font = 2,pt.cex=1, bty="n")
text(10,58.75, "b", col="black",font = 1, cex=0.8)

#################################
# 9. In what areas do we find evidence of environmental filtering and limiting similarity?
# 10. What leads us to these conclusions and what are the underlying mechanisms?
#################################

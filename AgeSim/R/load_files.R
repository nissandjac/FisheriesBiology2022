# Run MSE from the RAM stocks # 
library(reshape2)
library(gtools)
library(patchwork)
library(tidyverse)
library(TMB)
library(FishLife)

# Load the life history parameters 
source('R/getParameters_RAM.R')
df <- load('data/RAM/processed/SurplusProdStocks.RData')
# TMB::compile('src/runAssessment.cpp')
# dyn.load(dynlib("src/runAssessment"))
source('R/load_functions.R')
source('R/calcMSY.R')
source('R/getParameters_RAM.R')
load('Results/RAMdata.Rdata')

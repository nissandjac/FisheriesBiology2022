# Run MSE from the RAM stocks # 
library(reshape2)
library(gtools)
library(patchwork)
library(tidyverse)
library(FishLife)

# Load the life history parameters 
source('R/getParameters_RAM.R')
# TMB::compile('src/runAssessment.cpp')
# dyn.load(dynlib("src/runAssessment"))
source('R/load_functions.R')
source('R/calcMSY.R')
source('R/getParameters_RAM.R')
load('Results/RAMdata.Rdata')

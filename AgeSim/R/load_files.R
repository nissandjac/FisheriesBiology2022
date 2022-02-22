# Run MSE from the RAM stocks # 
# Load the life history parameters 
source('R/getParameters_RAM.R')
# TMB::compile('src/runAssessment.cpp')
# dyn.load(dynlib("src/runAssessment"))
source('R/load_functions.R')
source('R/calcMSY.R')
source('R/getParameters_RAM.R')
load('Results/RAMdata.Rdata')

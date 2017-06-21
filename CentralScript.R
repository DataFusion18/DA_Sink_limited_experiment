#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#- Central analysis script for DA with Sink limited experiment manuscript, entitled
#  "Inferring the effects of sink strength on carbon balance from experimental measurements".

#  The idea is to keep this script nice and tidy, but reproducibly do all the
#  analysis and make all of the figures for the manuscript. Raw and processed data will be 
#  placed in the "raw_data" and "processed_data" folders respectively, while figures 
#  and tables will be placed in the "output" folder.
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# Clear the workspace (if needed)
rm(list=ls())
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- Load the packages and custom functions that do all the work.
#  This will install a bunch of required libraries, including some non-standard stuff.
source("R/load_packages_CBM.R")

#- load the custom analysis and plotting functions that do all of the actual work
source("R/functions_CBM.R")
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
#- Download all raw data files. This downloads the zipfile from figshare
download.file("https://ndownloader.figshare.com/files/8724376", "raw_data.zip", mode="wb")
# Extract data to a folder named "raw_data".
unzip("raw_data.zip",overwrite=F)
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
#- This script imports and processes the raw Sink limited container volume experiment data 
#  to model the carbon pools and fluxes using MCMC
source("R/initial_data_processing.R")
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
#- Read the processed data and clean the workspace from data pre-processing
#  Processed data are placed in "raw_data" folder
source("R/read_data_CBM.R")
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
#- Make figure 1. Model representation of storage, allocation and autotrophic respiration processes and 
# pathways in the CBM with storage pool, separate growth and maintenance respiration components.
# source("R/plot.model.R")
plot.model()
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
#- Carbon Dynamic Model (CDM) comparison with and without storage
# Model run without storage
result.without.storage = CBM.grouping(chainLength=3000, no.param.par.var=c(1,2,3), vol.group=list(c(1,2,3),c(4,5,6),7), with.storage=F, model.comparison=T, model.optimization=F)
bic.without.storage = result.without.storage[[6]]

# Model run with storage
result.with.storage = CBM.grouping(chainLength=3000, no.param.par.var=c(1,2,3), vol.group=list(c(1,2,3),c(4,5,6),7), with.storage=T, model.comparison=T, model.optimization=F)
bic.with.storage = result.with.storage[[6]]

# Make figure 2. BIC for grouped treatments with and without storage pool
# source("R/plot.with.without.storage.R")
plot.with.without.storage(bic.with.storage, bic.without.storage)
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
#- Optimize both Treatment groupings and Parameter settings
#- This can take a long time (~5 minutes) to download and interpolate all of the climate data
bic.group1 = CBM.grouping(chainLength=3000, no.param.par.var=c(1,2,3), vol.group=list(c(1,2,3,4,5,6,7)), with.storage=T, model.comparison=F, model.optimization=T)
bic.group2 = CBM.grouping(chainLength=3000, no.param.par.var=c(1,2,3), vol.group=list(c(1,2,3,4,5,6),7), with.storage=T, model.comparison=F, model.optimization=T)
bic.group3 = CBM.grouping(chainLength=3000, no.param.par.var=c(1,2,3), vol.group=list(c(1,2,3),c(4,5,6),7), with.storage=T, model.comparison=F, model.optimization=T)
bic.group4 = CBM.grouping(chainLength=3000, no.param.par.var=c(1,2,3), vol.group=list(1,2,3,4,5,6,7), with.storage=T, model.comparison=F, model.optimization=T)

# Make figure 3. BIC values for different treatment groupings and parameter settings
# source("R/plot.parameter.settings.R")
plot.parameter.settings(bic.group1, bic.group2, bic.group3, bic.group4)
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
#- Parameter estimates and Analysis of carbon stock dynamics
# Model run without LA feedback
result = CBM.grouping(chainLength = 3000, no.param.par.var=c(3), vol.group=list(c(1,2,3),c(4,5,6),7), with.storage=T, model.comparison=F, model.optimization=F)

#- Make figure 4 and 5
# This script creates the figures and saves those in "output" folder
# source("R/generate_figures_CBM_grouping.R")
plot.Modelled.parameters(result)
plot.Modelled.biomass(result)
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
#- Perform sensitivity analysis and make figure 6
source("R/Parameter_sensitivity_shifting_LA.R")
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------

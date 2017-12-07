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
#  Processed data are placed in "processed_data" folder
source("R/read_data_CBM.R")
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
#- Make figure 1. Model representation of storage, allocation and autotrophic respiration processes and 
# pathways in the CBM with storage pool, separate growth and maintenance respiration components.
# source("R/plot.model.R")
plot.model()

# # Model run without LA feedback
# result = CBM.grouping(chainLength = 300, no.param.par.var=c(3), vol.group=list(c(1,2,3),c(4,5,6),7), with.storage=T, model.comparison=F, model.optimization=F)
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
#- Carbon Dynamic Model (CDM) simulation
# 1. Compare the model performance with and without storage - 2 simulations
# 2. Optimize the model settings: Try 3 different parameter numbers (variation over time) and 4 different treatment groupings (parameter variation for treatments) - 4 simulations
# 3. Final parameter estimates and Analysis of carbon stock dynamics with optimized parameter setting - 1 simulation
#- This can take long time (~30 minutes) to run all 7 simulations in parallel using 7 CPU cores

cluster <- makeCluster(detectCores()-1)
# clusterEvalQ(cluster, library(xts))
clusterExport(cl=cluster, list("Cday.data.processed","GPP.data.processed","Rd.data.processed","Mleaf.data.processed",
                               "Mstem.data.processed","Mroot.data.processed","tnc.data.processed","tnc"))
ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
clusterExport(cluster, ex)
result.cluster = list()
bic.cluster = list()

start <- proc.time() # Start clock
result.cluster <- clusterMap(cluster, CBM.grouping, with.storage=c(F,T,T,T,T,T,T), model.comparison=c(T,T,F,F,F,F,F), model.optimization=c(F,F,T,T,T,T,F), 
                             no.param.par.var=list(c(1,2,3),c(1,2,3),c(1,2,3),c(1,2,3),c(1,2,3),c(1,2,3),3),
         vol.group=c(list(list(c(1,2,3),c(4,5,6),7)),list(list(c(1,2,3),c(4,5,6),7)),list(list(c(1,2,3,4,5,6,7))),list(list(c(1,2,3,4,5,6),7)),
                     list(list(c(1,2,3),c(4,5,6),7)),list(c(1,2,3,4,5,6,7)),list(list(c(1,2,3),c(4,5,6),7))),
         MoreArgs=list(chainLength=3000))

time_elapsed_series <- proc.time() - start # End clock
bic.without.storage = result.cluster[[1]][[6]]
bic.with.storage = result.cluster[[2]][[6]]
bic.group1 = result.cluster[[3]]
bic.group2 = result.cluster[[4]]
bic.group3 = result.cluster[[5]]
bic.group4 = result.cluster[[6]]
result = result.cluster[[7]]
stopCluster(cluster)

#-------------------------------------------------------------------------------------
# #- Carbon Dynamic Model (CDM) comparison with and without storage
# #- This can take significant time (~25 minutes) to run both the analyses with and without storage options, 
# # 3 different parameter numbers and 3 different treatment groupings to generate figure 2
# 
# # Model run without storage
# result.without.storage = CBM.grouping(chainLength=3000, no.param.par.var=c(1,2,3), vol.group=list(c(1,2,3),c(4,5,6),7), with.storage=F, model.comparison=T, model.optimization=F)
# bic.without.storage = result.without.storage[[6]]
# 
# # Model run with storage
# result.with.storage = CBM.grouping(chainLength=3000, no.param.par.var=c(1,2,3), vol.group=list(c(1,2,3),c(4,5,6),7), with.storage=T, model.comparison=T, model.optimization=F)
# bic.with.storage = result.with.storage[[6]]
# 
# # Make figure 2. BIC for grouped treatments with and without storage pool
# # source("R/plot.with.without.storage.R")
# plot.with.without.storage(bic.with.storage, bic.without.storage)
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
# #- Optimize both Treatment groupings and Parameter settings
# #- This can take long time (~50 minutes) to run all 4 analyses with different treatment groupings and
# # 3 different parameter numbers to generate figure 3
# bic.group1 = CBM.grouping(chainLength=3000, no.param.par.var=c(1,2,3), vol.group=list(c(1,2,3,4,5,6,7)), with.storage=T, model.comparison=F, model.optimization=T)
# bic.group2 = CBM.grouping(chainLength=3000, no.param.par.var=c(1,2,3), vol.group=list(c(1,2,3,4,5,6),7), with.storage=T, model.comparison=F, model.optimization=T)
# bic.group3 = CBM.grouping(chainLength=3000, no.param.par.var=c(1,2,3), vol.group=list(c(1,2,3),c(4,5,6),7), with.storage=T, model.comparison=F, model.optimization=T)
# bic.group4 = CBM.grouping(chainLength=3000, no.param.par.var=c(1,2,3), vol.group=list(1,2,3,4,5,6,7), with.storage=T, model.comparison=F, model.optimization=T)
# 
# # Make figure 3. BIC values for different treatment groupings and parameter settings
# # source("R/plot.parameter.settings.R")
# plot.parameter.settings(bic.group1, bic.group2, bic.group3, bic.group4)
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
# #- Parameter estimates and Analysis of carbon stock dynamics
# # Model run without LA feedback
# result = CBM.grouping(chainLength = 3000, no.param.par.var=c(3), vol.group=list(c(1,2,3),c(4,5,6),7), with.storage=T, model.comparison=F, model.optimization=F)
# # result = CBM.grouping(chainLength = 3000, no.param.par.var=c(3), vol.group=list(7), with.storage=T, model.comparison=F, model.optimization=F)
# 
# #- Make figure 4 and 5
# # This script creates the figures and saves those in "output" folder
# # source("R/generate_figures_CBM_grouping.R")
# # source("R/functions_CBM.R")
plot.Modelled.parameters(result)
plot.Modelled.biomass(result)
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
#- Perform sensitivity analysis and make figure 6
# source("R/functions_CBM.R")
source("R/Parameter_sensitivity_shifting_LA.R")
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
# Calculate total C partitioning for individual treatments 
# and make figure 7 and Table S1
# source("R/functions_CBM.R")
source("R/C_partitioning.R")
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------

# mean(Cday.data.processed[Cday.data.processed$volume==10,3]) / mean(Cday.data.processed[Cday.data.processed$volume==5,3])
# mean(Cday.data.processed[Cday.data.processed$volume==20,3]) / mean(Cday.data.processed[Cday.data.processed$volume==10,3])
# mean(Cday.data.processed[Cday.data.processed$volume==35,3]) / mean(Cday.data.processed[Cday.data.processed$volume==15,3])

#-------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------




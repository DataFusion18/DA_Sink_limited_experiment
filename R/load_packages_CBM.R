#- This script cleans the workspace, loads necessary Rfunctions and packages
#---------------------------------------------------------------------
#- function to load a package, and install it if necessary
Library <- function(pkg, ...){
  
  PACK <- .packages(all.available=TRUE)
  pkgc <- deparse(substitute(pkg))
  
  if(pkgc %in% PACK){
    library(pkgc, character.only=TRUE)
  } else {
    install.packages(pkgc, ...)
    library(pkgc, character.only=TRUE)
  }
  
}
#---------------------------------------------------------------------


#---------------------------------------------------------------------
#- load all the libraries (and install them if needed)
Library(mvtnorm) # Creates candidate parameter vector as a multivariate normal jump away from the current candidate
Library(reshape2)
Library(ggplot2)
Library(lubridate)
Library(rio)
Library(dplyr)
Library(zoo)
Library(doBy)
Library(corrplot)
Library(png)
Library(grid)
Library(gridExtra)
Library(parallel)
Library(snow)
# Library(plyr)


# Packages from Court's script
Library(RColorBrewer)
Library(sciplot)
Library(reshape)
Library(nlme)
Library(visreg)
Library(multcomp)
Library(lme4)
Library(lmerTest)
Library(plotrix)
Library(plantecophys)
Library(scales)
Library(effects)

# install.packages("lubridate")
# library(lubridate)

# # Load model packages from Court's script
# Library(devtools)
# install_bitbucket("remkoduursma/plantecophys")
# library(plantecophys)
# # install.packages("YplantQMC")
# # install.packages("plantecophys")
# install_bitbucket("remkoduursma/HIEv")
# library(HIEv)
# # install.packages("plantecophys")
# # install_bitbucket("remkoduursma/plantecophys")
# # library(plantecophys)
# install_bitbucket("remkoduursma/optipipe")
# library(OptiPipe)
# # install_bitbucket("yplantqmc","remkoduursma", quick=FALSE)
# install_bitbucket("remkoduursma/yplantqmc")
# library(YplantQMC)
# # install.packages("YplantQMC")

#---------------------------------------------------------------------
#- check if the data and output directories exist. If they don't, create them.
dir.create(file.path("raw_data"),showWarnings=F)
dir.create(file.path("processed_data"),showWarnings=F)
dir.create(file.path("output"),showWarnings=F)
#---------------------------------------------------------------------


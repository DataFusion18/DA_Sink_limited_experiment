#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-- Script to analysis the quantum yield with pot size effect and recalculate Cday and shading factor 
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

###### R script to process Cday and shading factor considering the lower values of quantun yield
#-----------------------------------------------------------------------------------------

#apply photosynthesis model to generate Gross Photosynthesis (Cday) over the experiment
#apply values to leaf area that has been generated
#compare to actual measurements
source("R/load_packages_CBM.R")
source("R/functions_Court.R")
# source("functions and packages/plot objects.R")

#READ DATA, calculated data for A parameters, met data, and leaf area

# read jmax and vcmax from Court's Github
jmax_vcmax <- read.csv("raw_data/jmax_vcmax_clean.csv")
jmax_vcmax$SF = jmax_vcmax$Jmax.mean / max(jmax_vcmax$Jmax.mean)

jmax_vcmax$alpha = 0.24 * jmax_vcmax$SF

#   
#   phys_agg <- summaryBy(Jmax.mean+Vcmax.mean ~ volume, data=jmax_vcmax2, FUN=mean)
#   names(phys_agg)[2:5]<- c("Jmax", "Vcmax", "Jmax_se", "Vcmax_se")

#site weather data, rename variables, format date stuff
eucpve_met <- read.csv("raw_data/eucpve_met.csv")
names(eucpve_met)[2:5] <- c("record", "PPFD", "temp", "RH")
eucpve_met$Date <- as.Date(eucpve_met$DateTime15)
#need to turn the datetime 15 into hms
eucpve_met$DateTime15 <- ymd_hms(eucpve_met$DateTime15)
eucpve_met$time15 <- format(eucpve_met$DateTime15, format='%H:%M:%S')

#subset by Date range of experiment
eucpve_met1 <- subset(eucpve_met[,3:7], Date  >= "2013-01-21" & Date  <= "2013-05-21")

#Rdark Q10 equations by volume
rdarkq10 <- read.csv("raw_data/rdarkq10.csv")
rd25_clean <- read.csv("raw_data/rdark_clean.csv")

#---------------------------------------------------------------------------------------------
#A_model, use met data from site to generate parameters to enter into the model
#enter a vector of parameters over the period of the study to model PS with volume specific vcmax and jmax

#merge ps parameters to met data
A_model <- merge(eucpve_met1, jmax_vcmax)

# #need to calculate Rdark through time using rdarkq10 equation by volume
# A_model <- merge(A_model, rdarkq10[,1:2], by="volume")
# A_model <- merge(A_model, rd25_clean[,1:2], by="volume")
# 
# ####this q10 is too high, use crous for salinga from wtc1
# q10_crous <- 1.95
# q25_drake <- 1.86
# 
# #rdark_eq$rd25_euct <- with(rdark_eq, rd12.3*(q25_drake^((abs(CTleaf-25))/10)))
# 
# #refit RD
# A_model$Rd_pred2 <- with(A_model, rd12.3 * q25_drake^((temp-12.3)/10))
# #with(A_model, plot(temp,Rd_pred2, col=volume))

#------------------------------------------------------------------------------------------------------
#input parameters from optimal conductance model (using nls)

g1 <- read.csv("raw_data/g1_pred.csv")
g1_agg <- summaryBy(g1_date ~ volume, data=g1, keep.names=TRUE)


A_model <- merge(A_model, g1_agg, by="volume")

#-------------------------------------------------------------------------------------------------
#now run the model (includes pred Rdark from crouseq10, and gs parameters modelled from spot measurements)

#convert RH to VPD
A_model$VPD <- RHtoVPD(A_model$RH, A_model$temp, Pa=101)

#model, should return the aleaf for every 15 minutes. (will retrun Aleaf at 15min interval)
A_pred <- Photosyn(VPD=A_model$VPD,Ca=400, PPFD=A_model$PPFD, Tleaf=A_model$temp, 
                   Jmax=A_model$Jmax.mean, Vcmax=A_model$Vcmax.mean, alpha=A_model$alpha, Rd=0, g1=A_model$g1_date)

#need a new dfr with Aleaf and Anet across the day
Aleaf <- A_pred[,c(1:4, 7:11)]
Aleaf_15min <- cbind(Aleaf, A_model[,c(1, 5:6)])
# write.csv(Aleaf_15min, "processed_data/Aleaf_pred_15min_gross.csv", row.names=FALSE)

Aleaf_15min$Date <- as.Date(Aleaf_15min$Date)
Aleaf_15min$volume <- as.factor(Aleaf_15min$volume)
Aleaf_15min$photo15gc <- with(Aleaf_15min, ALEAF*15*60*10^-6*12)

Aleaf <- summaryBy(photo15gc ~ Date+volume, data=Aleaf_15min, FUN=sum, keep.names=TRUE )
names(Aleaf)[3] <- "carbon_day"
write.csv(Aleaf, "processed_data/cday_120_clean_gross.csv", row.names=FALSE)

# Aleaf_agg <- summaryBy(carbon_day ~ volume, data=Aleaf, FUN=mean, keep.names=TRUE )
# write.csv(Aleaf_agg, "processed_data/gCday_means_clean_gross.csv", row.names=FALSE)


#------------------------------------------------------------------------------------------------------
# plot the changes in Cday
Cday.prev <- read.csv("raw_data/cday_120_clean_gross.csv")
Cday.prev$alpha = as.factor("Fixed")
Cday.prev$Date = as.Date(Cday.prev$Date, format = "%Y-%m-%d")
Cday.prev$volume = as.factor(Cday.prev$volume)

Aleaf$alpha = as.factor("Variable")
Cday = rbind(Cday.prev, Aleaf)

font.size=12
pd <- position_dodge(0.2) 
p = ggplot(Cday, aes(x=Date, y=carbon_day, group = interaction(alpha,volume), shape=as.factor(alpha), colour=as.factor(volume))) + 
  geom_point(position=pd) +
  geom_line(position=pd,data = Cday, aes(x = Date, y = carbon_day, group = interaction(alpha,volume), linetype=as.factor(alpha), colour=as.factor(volume))) +
  ylab(expression(C[day]~"(g C "*m^"-2"*" "*d^"-1"*")")) + xlab("") +
labs(shape="alpha", linetype="alpha", colour="Soil Volume") +
theme_bw() +
theme(legend.title = element_text(colour="black", size=font.size)) +
theme(legend.text = element_text(colour="black", size = font.size)) +
theme(legend.key.height=unit(0.7,"line")) +
theme(legend.position = c(0.8,0.8), legend.direction = "vertical", legend.box = "horizontal") +
theme(legend.key = element_blank()) +
theme(text = element_text(size=font.size)) +
theme(axis.title.x = element_blank()) +
theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png("output/Figure_extra_Cday_prev_vs_modified.png", units="px", width=2000, height=1000, res=110)
p
dev.off()
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
# Run "pve_yplant.R" to get the shading factor M
#First need to run the sample trees through construct plant and get summary for each to use with raw data
#use readplantlist and summary.plant3d

#use the csv file to create the correct list of files and add the path so that readplant can use them
euckey <- read.csv("yplant/euc_plfiles/euc_key.csv")
euckey <- as.data.frame(euckey)
euckey$pfile <- paste("yplant/euc_plfiles/", euckey$pfile, sep = "")
euckey$lfile <- paste("yplant/euc_plfiles/", euckey$lfile, sep = "")

# #test
# test <- constructplant("yplant/euc_plfiles/Eletr5.p", "yplant/euc_plfiles/Elelf5.l")
# summary(test)
# # plot(test)

#use readplant list in order to construct multiple plants, complete by species
euc3d <- readplantlist(pfiles=euckey$pfile, lfiles=euckey$lfile)
#look at one random tree
# plot(euc3d[[15]])

#summarise
eucs3d <- summary(euc3d, writefile=TRUE)
# write.csv(eucs3d, "yplant/eucs_constructplant.csv", row.names=FALSE)


#--------------------------------------------------------------------------------------------------------
#run simulations for a sunny data from met data

#1. set location with long and latitude (can plot() to see if it is correct)
richmond <- setLocation(lat=-33.6, long=150.75, tzlong=150)

#2. get one met day 
#met data on 15 minute average from jan to june from HIEV (ros shelters)
weather <- read.csv("raw_data/eucpve_met.csv")
weather$DateTime15 <- ymd_hms(weather$DateTime15)
weather$Date <- as.Date(weather$DateTime15)
weather$time <- format(weather$DateTime15, format='%H:%M')

#search for some sunny days
sunmax <- subset(weather, PPFD_Avg.mean >= 1800)
unique(sunmax$Date)
#search for some cloudy days
sunmin <- subset(weather,  weather$time== "12:00:00" & weather$PPFD_Avg.mean <= 1000)

#subset one sunny day
metday <- subset(weather, Date=="2013-02-07")
plot(PPFD_Avg.mean~DateTime15, data=metday, ylim=c(0,2250))
metday$timeofday <- as.numeric(gsub(":", ".", metday$time))

#calculate total par
metday$par15_mol <- metday$PPFD_Avg.mean/1000000
plot(par15_mol~DateTime15, data=metday)
metday$par15_mol_s <- metday$par15_mol*15*60
daypar <- sum(metday$par15_mol_s)#molsm2d
daypar_mj <- daypar/(4.57/2)


#object with weather data from chosen date
sunnyday <- setMet(richmond, month=2, day=07, year=2013, nsteps=12, Tmin=14.83, Tmax=31.25, PARday=daypar_mj)
plot(sunnyday)

parpred <- as.data.frame(sunnyday[1][[1]])
#compare PAR from weather station and from setMET
plot(PAR~timeofday, data=parpred, ylim=c(0,2500), xlim=c(0,24))
points(PPFD_Avg.mean~timeofday, data=metday, pch=16)


# Test direct vs. diffuse--------------------------------------------------------------------
# sunny_fbeam0 <- setMet(richmond, month=2, day=16, nsteps=12, Tmin=14.44, Tmax=25.83, PARday=daypar,
#                      fbeam=0, fbeammethod="constant")
# plot(sunny_fbeam0)
# sunny_fbeam1 <- setMet(richmond, month=2, day=16, nsteps=12, Tmin=14.44, Tmax=25.83, PARday=daypar,
#                      fbeam=1, fbeammethod="constant")
# plot(sunny_fbeam1)

#------------------------------------------------------------------------------------------------
#3. measured and predicted photosynthetic parameters (from my data)
phys <- read.csv("raw_data/jmax_vcmax.csv")
phys_agg <- summaryBy(Jmax.mean+Vcmax.mean~ volume, data=phys, keep.names=TRUE)
names(phys_agg)[2:3]<- c("Jmax", "Vcmax")
phys_agg$SF = phys_agg$Jmax / max(phys_agg$Jmax)
phys_agg$alpha = 0.24 * phys_agg$SF
phys_agg = phys_agg[,-4]

g1 <- read.csv("raw_data/g1_pred.csv")
g1.court <- read.csv("https://raw.githubusercontent.com/kashifmahmud/EucPVE/master/calculated%20data/g1_pred.csv")

g1_mean <- mean(g1$g1_date)
g1_agg <- summaryBy(g1_vol ~ volume, data=g1, keep.names=TRUE)
g1_agg.c <- summaryBy(g1_vol ~ volume, data=g1.court, keep.names=TRUE)

# rd <- read.csv("calculated data/Rd_leaf.csv")
rd2 <- read.csv("raw_data/rdark_clean.csv")
# rd_agg <- rd[,c(7,9,15)]
# rd_agg <- summaryBy(.~ volume, data=rd_agg, keep.names=TRUE)
names(rd2)[2]<- "respdark"

Aparam <- merge(rd2[,1:2], g1_agg)
Aparam <- merge(Aparam, phys_agg)
write.csv(Aparam, "processed_data/A_parameters.csv", row.names=FALSE)

# Aparam.1 <- read.csv("https://raw.githubusercontent.com/kashifmahmud/EucPVE/master/calculated%20data/A_parameters.csv")
# Aparam$g1_vol = Aparam.1$g1_vol
A_free.1 <- subset(Aparam.1, volume=="1000")

#A parameters for each volume
A_free <- subset(Aparam, volume=="1000")
A_5 <- subset(Aparam, volume=="5")
A_10 <- subset(Aparam, volume=="10")
A_15 <- subset(Aparam, volume=="15")
A_20 <- subset(Aparam, volume=="20")
A_25 <- subset(Aparam, volume=="25")
A_35 <- subset(Aparam, volume=="35")

#5.setPhy = Constructs an object of class 'ypphy', which calculates A and E from weather data and PAR.
# eucphy_free <- setPhy("Photosyn",leafpars=list(Vcmax=A_free[1,5], Jmax=A_free[1,4], g1=A_free[1,3], Rd0=-(A_free[1,2]), alpha=A_free[1,6]))
# eucphy5 <- setPhy("Photosyn",leafpars=list(Vcmax=A_5[1,5], Jmax=A_5[1,4], g1=A_5[1,3], Rd0=-(A_5[1,2]), alpha=A_5[1,6]))
# eucphy10 <- setPhy("Photosyn",leafpars=list(Vcmax=A_10[1,5], Jmax=A_10[1,4], g1=A_10[1,3], Rd0=-(A_10[1,2]), alpha=A_10[1,6]))
# eucphy15 <- setPhy("Photosyn",leafpars=list(Vcmax=A_15[1,5], Jmax=A_15[1,4], g1=A_15[1,3], Rd0=-(A_15[1,2]), alpha=A_15[1,6]))
# eucphy20 <- setPhy("Photosyn",leafpars=list(Vcmax=A_20[1,5], Jmax=A_20[1,4], g1=A_20[1,3], Rd0=-(A_20[1,2]), alpha=A_20[1,6]))
# eucphy25 <- setPhy("Photosyn",leafpars=list(Vcmax=A_25[1,5], Jmax=A_25[1,4], g1=A_25[1,3], Rd0=-(A_25[1,2]), alpha=A_25[1,6]))
# eucphy35 <- setPhy("Photosyn",leafpars=list(Vcmax=A_35[1,5], Jmax=A_35[1,4], g1=A_35[1,3], Rd0=-(A_35[1,2]), alpha=A_35[1,6]))

eucphy_free <- setPhy("Photosyn",leafpars=list(Vcmax=A_free[1,5], Jmax=A_free[1,4], g1=A_free[1,3], Rd0=-(A_free[1,2])))
eucphy5 <- setPhy("Photosyn",leafpars=list(Vcmax=A_5[1,5], Jmax=A_5[1,4], g1=A_5[1,3], Rd0=-(A_5[1,2])))
eucphy10 <- setPhy("Photosyn",leafpars=list(Vcmax=A_10[1,5], Jmax=A_10[1,4], g1=A_10[1,3], Rd0=-(A_10[1,2])))
eucphy15 <- setPhy("Photosyn",leafpars=list(Vcmax=A_15[1,5], Jmax=A_15[1,4], g1=A_15[1,3], Rd0=-(A_15[1,2])))
eucphy20 <- setPhy("Photosyn",leafpars=list(Vcmax=A_20[1,5], Jmax=A_20[1,4], g1=A_20[1,3], Rd0=-(A_20[1,2])))
eucphy25 <- setPhy("Photosyn",leafpars=list(Vcmax=A_25[1,5], Jmax=A_25[1,4], g1=A_25[1,3], Rd0=-(A_25[1,2])))
eucphy35 <- setPhy("Photosyn",leafpars=list(Vcmax=A_35[1,5], Jmax=A_35[1,4], g1=A_35[1,3], Rd0=-(A_35[1,2])))


# eucphy_free.1 <- setPhy("Farquhar",leafpars=list(Vcmax=A_free[1,5], Jmax=A_free[1,4], g1=A_free[1,3], Rd0=-(A_free[1,2])))
# euc3d.test = euc3d[[1]]
# eucphyList.test <- list(eucphy_free = eucphy_free,
#                    eucs5 = eucphy5)
# euc_list.test <- list()
# for(i in 1:length(eucphyList.test)){
#   euc_list.test[[i]] <- YplantDay(euc3d.test, phy = eucphyList.test[[i]], met = sunnyday, PSRsuffix=names(eucphyList.test)[i])
# }
# psrdata(euc_list[[1]])
# psrdata(euc_list[[2]])
# 
# eucsumm <- lapply(euc_list, summary)
# listnames <- c("eucs_free", "eucs_free.1")
# names(eucsumm) <- listnames
# eucsumm[[1]]
# l_ply(names(eucsumm), function(x) write.csv(eucsumm[[x]], file = paste("yplant/sim2_summary/", x, ".csv", sep = "")))

# #test one plant
# testday <- YplantDay(test, phy=eucphy20, met=sunnyday)
# summary(testday)
# plot(testday)
# testdata<-psrdata(testday)


# #direct v. diffuse
# run2 <- YplantDay(plant, phy=eucphy, met=sunny_fbeam0)
# run3 <- YplantDay(plant, phy=eucphy, met=sunny_fbeam1)
# 
# with(psrdata(run2), plot(timeofday, A/A0, type='l'))
# with(psrdata(run3), points(timeofday, A/A0, type='l', col="red"))



####diurnal simulation of all euc plants---------------------------------------------------------
#from summary of these plants  get variables in which to correlate to self shading multiplier

#create list of setPhy objects
eucphyList <- list(eucphy_free = eucphy_free, 
                   eucs5 = eucphy5, 
                   eucs10 = eucphy10, 
                   eucs15 = eucphy15, 
                   eucs20 = eucphy20,
                   eucs25 = eucphy25, 
                   eucs35 = eucphy35)

#run yplantday on the eucphylist, with richmond sunnday, and euc3d plants (61) 
# euc_test <- lapply(list(eucphyList[[1]], eucphyList[[2]]), function(x) YplantDay(test, phy = x, met = sunnyday))
# summary(euc_test[1])
# plot(euc_test[[2]])


# euc_list_try <- list()
# euc_list_try[[1]] <- YplantDay(euc3d[[1]], phy = eucphyList[[1]], met = sunnyday, PSRsuffix=names(eucphyList)[1])
# summary(euc_list_try[[1]])

# sum_fun <- function(x){
#   p <- psrdata(x)
#   vec <- c(ALEAF=mean(p$ALEAF), ALEAF0=mean(p$ALEAF0))
#   
#   return(vec)
# }

#run simulation, make empty list, run each phy=volume, and output psr files for each tree, then lapply for summary
euc_list <- list()
for(i in 1:length(eucphyList)){
  euc_list[[i]] <- YplantDay(euc3d, phy = eucphyList[[i]], met = sunnyday, PSRsuffix=names(eucphyList)[i])
}
# summary(euc_list[[1]])
# psrdata(euc_list[[i]])

# sum_fun <- function(x){
#   p <- psrdata(x)
#   vec <- c(ALEAF=mean(p$ALEAF), ALEAF0=mean(p$ALEAF0))
# 
#   return(vec)
# }
# 
# yrunsa <- as.data.frame(do.call(rbind, lapply(yruns1, sum_fun)))

#-----------------------------------------------------------------------------------------
# # perform a test analysis
# eucphy5.test <- setPhy("Photosyn",leafpars=list(Vcmax=A_5[1,5], Jmax=A_5[1,4], g1=A_5[1,3], Rd0=-(A_5[1,2]), alpha=A_5[1,6]))
# eucphy10.test <- setPhy("Photosyn",leafpars=list(Vcmax=A_10[1,5], Jmax=A_10[1,4], g1=A_10[1,3], Rd0=-(A_10[1,2]), alpha=A_10[1,6]))
# 
# eucphyList.test <- list(eucs5.test = eucphy5.test, 
#                    eucs10.test = eucphy10.test)
# euckey.test <- read.csv("/Users/kashifmahmud/WSU/Final_projects/yplantupscale/plantfiles/euc_key.csv", stringsAsFactors = FALSE)
# euckey.test <- as.data.frame(euckey.test[61,])
# # euckey.test$pfile <- paste("yplant/euc_plfiles/", euckey.test$pfile, sep = "")
# # euckey.test$lfile <- paste("yplant/euc_plfiles/", euckey.test$lfile, sep = "")
# 
# # if(!require(pacman))install.packages("pacman")
# # pacman::p_load(YplantQMC, withr, lubridate, plantecophys)
# # euckey.test <- read.csv("/Users/kashifmahmud/WSU/Final_projects/yplantupscale/plantfiles/euc_key.csv", stringsAsFactors = FALSE)
# euc3d.test <- with_dir("/Users/kashifmahmud/WSU/Final_projects/yplantupscale/plantfiles", 
#                  readplantlist(pfiles=euckey.test$pfile, lfiles=euckey.test$lfile))
# # euc3d.test <- readplantlist(pfiles=euckey.test$pfile, lfiles=euckey.test$lfile)
# summary(euc3d.test)
# # euc3d.test <- constructplant(pfiles=euckey.test$pfile, lfiles=euckey.test$lfile)
# eucs3d.test <- summary(euc3d.test, writefile=TRUE)
# 
# euc_list.test <- list()
# for(i in 1:length(eucphyList.test)){
#   euc_list.test[[i]] <- YplantDay(euc3d.test, phy = eucphyList.test[[i]], met = sunnyday, PSRsuffix=names(eucphyList.test)[i])
# }
# summary(euc_list.test[[1]])

#-----------------------------------------------------------------------------------------


#saveRDS(euc_list, "yplant/euc_sim2.rds")
#euc_list <- readRDS("somefilename.rds")
eucsumm <- lapply(euc_list, summary)
#need to add names the list by their volume
listnames <- c("eucs_5", "eucs_10", "eucs_15", "eucs_20", "eucs_25", "eucs_35", "eucs_free")
names(eucsumm) <- listnames

#save each list as a dfr with the name (##use names of list for apply function, see use of [[x]] for func arg)

#lapply(names(eucsumm), function(x) write.csv(eucsumm[[x]], file = paste(x, ".csv", sep = "")))
l_ply(names(eucsumm), function(x) write.csv(eucsumm[[x]], file = paste("yplant/sim_prev_summary/", x, ".csv", sep = "")))

#run on a cloudy day
# eucs_cloud <- YplantDay(euc3d, phy=eucphy, met=cloudyday)
# plot(eucs_cloud)
# # add summary variables:
# cloudy_summary <- summary(eucs_cloud) 

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#calculate M from yplant sim and then regression against leaf number for each plant ()
#use simulations from sunny day and A parameters for each pot volume to calc M
#M = Aplant/Aobs

# source("functions and packages/plot objects.R")

#read eucs3d summary from construct plant
eucs3d <- read.csv("yplant/eucs_constructplant.csv")

#read in summary files for all volumes into a list
eucs_summs <- list.files(path = "yplant/sim_prev_summary/", pattern="*.csv", full.names = TRUE)
eucs_list <- lapply(eucs_summs, function(x) read.csv(x))

#name list elements
listnames <- c("10", "15", "20", "25", "35", "5", "free")
names(eucs_list) <- listnames

#function to calculate M (twp arguments:1st is dfr, 2nd is the name of each list element)
Mcalc <- function(dfr, vol){
  dfr$plant_id <- gsub("yplant/euc_plfiles/", "", dfr$pfile)
  dfr$plant_id <- gsub(".p", "", dfr$plant_id)
  dfr$M <- dfr$totALEAF/dfr$totALEAF0
  M_dfr <- subset(dfr, select = c("plant_id", "M", "totPARleaf"))
  M_dfr$volume <- as.factor(vol)
  return(M_dfr)
}

#call the names of the list into plyr, run the function on each element (eucs_list[[x]]) where vol argument = names
M_eucs <- llply(names(eucs_list), function(x) Mcalc(dfr = eucs_list[[x]], vol = x))
#write.csv(M_eucs, "processed_data/M_eucs.csv", row.names=FALSE)


####now regress M for each plant against leaf------------------------------------------------------------------------ 

#1: first merge M for each plant (7 volumes) with leaf# for each plant id
eucs3d$plant_id <- gsub("yplant/euc_plfiles/", "", eucs3d$pfile)
eucs3d$plant_id <- gsub(".p", "", eucs3d$plant_id)
eucs_allom <- eucs3d[,c("plant_id", "LA")]

M_eucs2<- lapply(M_eucs, function(x) merge(x, eucs_allom, by="plant_id", all=TRUE))
#write.csv(M_eucs3d, "processed_data/M_eucs3d.csv", row.names=FALSE)


#2. for each volume trt (list) run model of M vs leafN, 

M_regress <- lapply(M_eucs2, function(x) lm(M ~ LA, data=x))
# 
#  test <- M_eucs2[[7]]
#  test_lm <- lm(M ~ LA, data=test)
#  windows(8,10)
#  visreg(test_lm)
#  dev.copy2pdf(file= "gC_day_model/model_output/Mleafarea.pdf")
#  dev.off()


#3. extract coefs for each model by trt

M_coefs <- lapply(M_regress, function(x) extract_func(x))
##make a dfr but order of list these come from was 5-free, need to reorder carefully
M_coefs2 <-rbind.fill(M_coefs)
volorder <- names(eucs_list)
M_coefs2$volume <- as.factor(volorder)
M_coefs2$volume <- gsub("free", 1000, M_coefs2$volume)
#reorder
M_coefs3 <- M_coefs2[c(6, 1:5, 7),]

write.csv(M_coefs3, "processed_data/M_leafarea_model_prev.csv", row.names=FALSE)

#------------------------------------------------------------------------------------------------------
# # plot the changes in self shading
# sigma_prev <- read.csv("raw_data/M_leafarea_model.csv")
# # sigma_new <- M_coefs3
# sigma_new <- read.csv("processed_data/M_leafarea_model.csv")
# sigma.test = 
# plot()

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------





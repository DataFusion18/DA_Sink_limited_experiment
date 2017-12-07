#----------------------------------------------------------------------------------------------------------------
#- This is a collection of many functions that do the actual work of data analysis, run the model and figure plotting. 
# These functions are called by just a few lines of code in "CentralScript.R" to recreate the analyses and figures.

# Analysing Plant Storage (TNC) partitioning for Cstorage pool prediction
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
standard.error <- function(dat,na.rm=F,...){
  if(na.rm==T){
    dat <- subset(dat,is.na(dat)==F)
  }
  std <- sd(dat)
  n <- length(dat)
  se <- std/sqrt(n)
  return(se)
}

#----------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
##### This function calculates TNC partitioning to tree organs (without considering organ biomass)
tnc.analysis <- function(carbohydrates,harvest) {
  # carbohydrates = read.csv("raw_data/Duan_carbohydrates.csv")
  carbohydrates = subset(carbohydrates,CO2 == 400 & Temp == "Amb" & Water == "Well watered")
  carbohydrates$tnc = carbohydrates$StarchW + carbohydrates$SolSugW # Unit = mg of tnc per g of dry weight biomass
  carbohydrates$tnc = carbohydrates$tnc / 10 # Unit = % of dry weight biomass
  
  # unit conversion from g of tnc per g of dry weight biomass to gC in tnc per gC in plant biomass
  # 1 g of tnc has 0.4 gC and 1 g of dry weight biomass has 0.48 gC
  carbohydrates$tnc = carbohydrates$tnc * (0.4/c1) # Unit = gC in tnc / gC in plant biomass
  #-----------------------------------------------------------------------------------------
  ##### Total TNC calculation considering tree organ biomass partitioning
  # harvest = read.csv("raw_data/Duan_harvest.csv")
  harvest = subset(harvest,CO2 == 400 & Temp == "Amb" & Water == "Well watered")
  
  leaf.tnc = subset(carbohydrates,Organ == "Leaf") # Unit = % of dry weight leafmass
  stem.tnc = subset(carbohydrates,Organ == "Stem") # Unit = % of dry weight stemmass
  root.tnc = subset(carbohydrates,Organ == "Root") # Unit = % of dry weight rootmass
  
  tnc = data.frame(harvest$LeafDW,leaf.tnc$tnc,harvest$StemDW,stem.tnc$tnc,harvest$RootDW,root.tnc$tnc)
  names(tnc) <- c("leaf.C","leaf.tnc.C","stem.C","stem.tnc.C","root.C","root.tnc.C") 
  
  tnc$total.leaf.tnc.C = tnc$leaf.tnc.C * tnc$leaf.C / 100 # Unit = gC
  tnc$total.stem.tnc.C = tnc$stem.tnc.C * tnc$stem.C / 100 # Unit = gC
  tnc$total.root.tnc.C = tnc$root.tnc.C * tnc$root.C / 100 # Unit = gC
  
  tnc$leaf_to_all = tnc$total.leaf.tnc.C / (tnc$total.leaf.tnc.C + tnc$total.stem.tnc.C + tnc$total.root.tnc.C) * 100 # Unit = %
  tnc$stem_to_all = tnc$total.stem.tnc.C / (tnc$total.leaf.tnc.C + tnc$total.stem.tnc.C + tnc$total.root.tnc.C) * 100 # Unit = %
  tnc$root_to_all = tnc$total.root.tnc.C / (tnc$total.leaf.tnc.C + tnc$total.stem.tnc.C + tnc$total.root.tnc.C) * 100 # Unit = %
  
  tnc[nrow(tnc)+1, ] = colMeans(tnc, na.rm = TRUE) # R7 = Average of data
  tnc[nrow(tnc)+1, ] = apply(tnc, 2, sd) # R8 = Standard deviation of data
  dimnames(tnc)[[1]] <- c(1:6, "mean", "SD")
  
  return(tnc)
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# This script calcualtes LogLikelihood to find the most accurate model
logLikelihood <- function (data,output,model.comparison) {
  logLi <- matrix(0, nrow=nrow(data), ncol = 1) # Initialising the logLi
  for (i in 1:nrow(data)) {
    if (!is.na(data$Mleaf[i])) {
      logLi[i] = - 0.5*((output$Mleaf[i] - data$Mleaf[i])/data$Mleaf_SD[i])^2 - log(data$Mleaf_SD[i]) - log(2*pi)^0.5
    }
    if (!is.na(data$Mstem[i])) {
      logLi[i] = logLi[i] - 0.5*((output$Mstem[i] - data$Mstem[i])/data$Mstem_SD[i])^2 - log(data$Mstem_SD[i]) - log(2*pi)^0.5
    }
    if (!is.na(data$Mroot[i])) {
      logLi[i] = logLi[i] - 0.5*((output$Mroot[i] - data$Mroot[i])/data$Mroot_SD[i])^2 - log(data$Mroot_SD[i]) - log(2*pi)^0.5
    }
    if (model.comparison==F) {
      if (!is.null(data$Sleaf)) {
        if (!is.na(data$Sleaf[i])) {
          logLi[i] = logLi[i] - 0.5*((output$Sleaf[i] - data$Sleaf[i])/data$Sleaf_SD[i])^2 - log(data$Sleaf_SD[i]) - log(2*pi)^0.5
        }
      }
    }
  }
  return(sum(logLi))
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Carbon balance model (CBM)
# Developed by Kashif Mahmud and Belinda Medlyn (March 2017)
# k.mahmud@westernsydney.edu.au

# This version tries to group various treatments according to their similarities to have a trend in paramter settings

# This code carries out Bayesian calibration for 5 variables (allocation fractions: "k","Y",af","as","sf") on 
# various temporal scales (e.g. 1,2,...,121 days) to estimate Carbon pools (Cstorage,Cleaf,Cstem,Croot) and fluxes
#-------------------------------------------------------------------------------------
CBM.grouping <- function(chainLength, no.param.par.var, vol.group, with.storage, model.comparison, model.optimization) {
  
  source("R/load_packages_CBM.R")
  
  # Assign inputs for MCMC
  bunr_in = chainLength * 0.1 # Discard the first 10% iterations for Burn-IN in MCMC (According to Oijen, 2008)
  if (with.storage==T) {
    no.var = 5 # variables to be modelled are: k,Y,af,as,sf
  } else {
    no.var = 4 # variables to be modelled are: k,Y,af,as
  }
  
  # Assign pot volumes and number of parameters per varible in temporal scale
  vol = unique(Cday.data.processed$volume)[order(unique(Cday.data.processed$volume))] # Assign all treatment pot volumes
  
  param.mean = data.frame(matrix(ncol = no.var+1, nrow = length(no.param.par.var)*length(vol.group)))
  if (with.storage==T) {
    names(param.mean) = c("k","Y","af","as","ar","sf")
  } else {
    names(param.mean) = c("Y","af","as","ar","sf")
  }
  aic.bic = data.frame(matrix(ncol = 7, nrow = length(no.param.par.var)*length(vol.group)))
  names(aic.bic) <- c("logLi","aic","bic","time","volume.group","no.param","volume")
  time = data.frame(no.param=rep(no.param.par.var,length(vol.group)),
                    start.time=numeric(length(no.param.par.var)*length(vol.group)),
                    end.time=numeric(length(no.param.par.var)*length(vol.group)),
                    time.taken=numeric(length(no.param.par.var)*length(vol.group)))
  
  q = 0 # Indicates the iteration number
  # set.seed(3)
  set.seed(15) # final seed for reproducible results
  # set.seed(18) 

  # Start the iteration for different treatment group and number of parameters
  for (v1 in 1:length(vol.group)) {
    v = unlist(vol.group[v1])
    # This script take the subset of processed data for particular treatment group
    source("R/data_processing_CBM.R", local=TRUE)
    
    for (z in 1:length(no.param.par.var)) {
      # Initialize few output data files
      q = q + 1
      time$start.time[q] <- Sys.time()
      param.vary = ceiling(nrow(data)/no.param.par.var[z]) # How many days the parameter set remain unchanged (weekly = 7; monthly = 30; just one parameter = nrow(data))
      no.param = ceiling(nrow(data)/param.vary) # number of parameter set for the whole duration of experiment (121 days)
      
      # This script initializes the parameter setting
      source("R/parameter_setting.R", local=TRUE) # initialize 'sf' prior differently for grouped treatments
      
      # Defining the variance-covariance matrix for proposal generation
      vcov = (0.005*(pMaxima-pMinima))^2
      vcovProposal =  vcov # The higher the coefficient, the higher the deviations in parameter time series
      
      
      browser()
      # Find the Prior probability density
      prior.dist = vector("list", no.var)
      for (i in 1:no.var) {
        prior.dist[i] = list(log(dnorm(pValues[ , i], (pMinima[ , i] + pMaxima[ , i])/2, (pMaxima[ , i] - pMinima[ , i])/3))) # Prior normal gaussian distribution
      }
      logPrior0 <- sum(unlist(prior.dist))
      
      # Calculating model outputs for the starting point of the chain
      for (j in 1:length(v)) {
        data.set = subset(data,(volume %in% vol[v[j]]))
        Mleaf = Mstem = Mroot = LA = c()
        Mleaf[1] <- data.set$Mleaf[1]
        Mstem[1] <- data.set$Mstem[1]
        Mroot[1] <- data.set$Mroot[1]
        
        if (with.storage==T) {
          output.set = model(data.set$GPP,data.set$Rd,no.param,Mleaf,Mstem,Mroot,tnc,pValues$Y,pValues$k,pValues$af,pValues$as,pValues$sf)
        } else {
          output.set = model.without.storage(data.set$GPP,data.set$Rd,no.param,Mleaf,Mstem,Mroot,pValues$Y,pValues$af,pValues$as,pValues$sf)
        }
        
        output.set$volume = as.factor(vol[v[j]])
        if (j == 1) {
          output = output.set
        }
        if (j > 1) {
          output = rbind(output,output.set)
        }
      }
      
      data = data[order(data$volume),]
      logL0 <- logLikelihood(data,output,model.comparison) # Calculate log likelihood of starting point of the chain
      
      if (with.storage==T) {
        pChain[1,] <- c(pValues$k,pValues$Y,pValues$af,pValues$as,pValues$sf,logL0) # Assign the first parameter set with log likelihood
      } else {
        pChain[1,] <- c(pValues$Y,pValues$af,pValues$as,pValues$sf,logL0) # Assign the first parameter set with log likelihood
      }
      
      
      # Calculating the next candidate parameter vector, as a multivariate normal jump away from the current point
      for (c in (2 : chainLength)) {
        candidatepValues = matrix(ncol = no.var, nrow = no.param)
        for (i in 1:no.var) {
          candidatepValues[,i] = rmvnorm(n=1, mean=pValues[,i],
                                         sigma=diag(vcovProposal[,i],no.param)) 
        }
        candidatepValues = data.frame(candidatepValues)
        if (with.storage==T) {
          names(candidatepValues) <- c("k","Y","af","as","sf")
        } else {
          names(candidatepValues) <- c("Y","af","as","sf")
        }
        
        # Reflected back to generate another candidate value
        reflectionFromMin = pmin( unlist(matrix(0,nrow=no.param,ncol=no.var)), unlist(candidatepValues-pMinima) )
        reflectionFromMax = pmax( unlist(list(rep(0, no.param))), unlist(candidatepValues-pMaxima) )
        candidatepValues = candidatepValues - 2 * reflectionFromMin - 2 * reflectionFromMax 
        
        
        # Calculating the prior probability density for the candidate parameter vector
        if (all(candidatepValues>pMinima) && all(candidatepValues<pMaxima)){
          uni.dist = vector("list", no.var)
          for (i in 1:no.var) {
            uni.dist[i] = list(log(dnorm(candidatepValues[ , i], (pMinima[ , i] + pMaxima[ , i])/2, (pMaxima[ , i] - pMinima[ , i])/3))) # Prior normal gaussian distribution
          }
          logPrior1 <- sum(unlist(uni.dist))
          Prior1 = 1
        } else {
          Prior1 <- 0
        }
        
        
        # Calculating the outputs for the candidate parameter vector and then log likelihood
        if (Prior1 > 0) {
          for (j in 1:length(v)) {
            data.set = subset(data,(volume %in% vol[v[j]]))
            Mleaf = Mstem = Mroot = c()
            Mleaf[1] <- data.set$Mleaf[1]
            Mstem[1] <- data.set$Mstem[1]
            Mroot[1] <- data.set$Mroot[1]
            
            
            if (with.storage==T) {
              out.cand.set = model(data.set$GPP,data.set$Rd,no.param,Mleaf,Mstem,Mroot,tnc,candidatepValues$Y,
                                   candidatepValues$k,candidatepValues$af,candidatepValues$as,candidatepValues$sf)
            } else {
              out.cand.set = model.without.storage(data.set$GPP,data.set$Rd,no.param,Mleaf,Mstem,Mroot,candidatepValues$Y,
                                                   candidatepValues$af,candidatepValues$as,candidatepValues$sf)
            }
            
            out.cand.set$volume = as.factor(vol[v[j]])
            if (j == 1) {
              out.cand = out.cand.set
            }
            if (j > 1) {
              out.cand = rbind(out.cand,out.cand.set)
            }
          }
          
          data = data[order(data$volume),]
          logL1 <- logLikelihood(data,out.cand,model.comparison) # Calculate log likelihood
          
          # Calculating the logarithm of the Metropolis ratio
          logalpha <- (logPrior1+logL1) - (logPrior0+logL0) 
          
          # Accepting or rejecting the candidate vector
          if ( log(runif(1, min = 0, max =1)) < logalpha && candidatepValues$af[1] + candidatepValues$as[1] <= 1
               && candidatepValues$as[1] >= 0 && candidatepValues$af[1] >= 0) {
            pValues <- candidatepValues
            logPrior0 <- logPrior1
            logL0 <- logL1
          }
        }
        if (with.storage==T) {
          pChain[c,] <- c(pValues$k,pValues$Y,pValues$af,pValues$as,pValues$sf,logL0)
        } else {
          pChain[c,] <- c(pValues$Y,pValues$af,pValues$as,pValues$sf,logL0)
        }
      }
      
      # Discard the first 500 iterations for Burn-IN in MCMC
      pChain <- pChain[(bunr_in+1):nrow(pChain),]
      pChain = as.data.frame(pChain)
      if (with.storage==T) {
        if (no.param.par.var[z]==1) {
          names(pChain) <- c("k1","Y1","af1","as1","sf1","logli")
        } else if (no.param.par.var[z]==2) {
          names(pChain) <- c("k1","k2","Y1","Y2","af1","af2","as1","as2","sf1","sf2","logli")
        } else if (no.param.par.var[z]==3) {
          names(pChain) <- c("k1","k2","k3","Y1","Y2","Y3","af1","af2","af3","as1","as2","as3","sf1","sf2","sf3","logli")
        }
      } else {
        if (no.param.par.var[z]==1) {
          names(pChain) <- c("Y1","af1","as1","sf1","logli")
        } else if (no.param.par.var[z]==2) {
          names(pChain) <- c("Y1","Y2","af1","af2","as1","as2","sf1","sf2","logli")
        } else if (no.param.par.var[z]==3) {
          names(pChain) <- c("Y1","Y2","Y3","af1","af2","af3","as1","as2","as3","sf1","sf2","sf3","logli")
        }
      }
      
      
      # Store the final parameter set values
      param.set = colMeans(pChain[ , 1:(no.param*no.var)])
      param.SD = apply(pChain[ , 1:(no.param*no.var)], 2, sd)
      param.final = data.frame(matrix(ncol = (no.var)*2, nrow = no.param))
      if (with.storage==T) {
        names(param.final) <- c("k","Y","af","as","sf","k_SD","Y_SD","af_SD","as_SD","sf_SD")
        param.final$k = param.set[1:no.param]
        param.final$Y = param.set[(1+no.param):(2*no.param)]
        param.final$af = param.set[(1+2*no.param):(3*no.param)]
        param.final$as = param.set[(1+3*no.param):(4*no.param)]
        param.final$sf = param.set[(1+4*no.param):(5*no.param)]
        
        param.final$k_SD = param.SD[1:no.param]
        param.final$Y_SD = param.SD[(1+no.param):(2*no.param)]
        param.final$af_SD = param.SD[(1+2*no.param):(3*no.param)]
        param.final$as_SD = param.SD[(1+3*no.param):(4*no.param)]
        param.final$sf_SD = param.SD[(1+4*no.param):(5*no.param)]
      } else {
        names(param.final) <- c("Y","af","as","sf","Y_SD","af_SD","as_SD","sf_SD")
        param.final$Y = param.set[1:no.param]
        param.final$af = param.set[(1+no.param):(2*no.param)]
        param.final$as = param.set[(1+2*no.param):(3*no.param)]
        param.final$sf = param.set[(1+3*no.param):(4*no.param)]
        
        param.final$Y_SD = param.SD[1:no.param]
        param.final$af_SD = param.SD[(1+no.param):(2*no.param)]
        param.final$as_SD = param.SD[(1+2*no.param):(3*no.param)]
        param.final$sf_SD = param.SD[(1+3*no.param):(4*no.param)]
      }
      
      # Calculate final output set from the predicted parameter set
      for (j in 1:length(v)) {
        data.set = subset(data,(volume %in% vol[v[j]]))
        Mleaf = Mstem = Mroot = c()
        Mleaf[1] <- data.set$Mleaf[1]
        Mstem[1] <- data.set$Mstem[1]
        Mroot[1] <- data.set$Mroot[1]
        
        if (with.storage==T) {
          output.final.set = model(data.set$GPP,data.set$Rd,no.param,Mleaf,Mstem,Mroot,tnc,param.final$Y,
                                   param.final$k,param.final$af,param.final$as,param.final$sf)
        } else {
          output.final.set = model.without.storage(data.set$GPP,data.set$Rd,no.param,Mleaf,Mstem,Mroot,param.final$Y,
                                                   param.final$af,param.final$as,param.final$sf)
        } 
        
        output.final.set$volume = as.factor(vol[v[j]])
        if (j == 1) {
          output.final = output.final.set
        }
        if (j > 1) {
          output.final = rbind(output.final,output.final.set)
        }
      }
      
      # #----------------------------------------------------------------------------------------------------------------
      # if (with.storage==T) {
      #   output.final$Sleaf = output.final$Sleaf / output.final$Mleaf * 100
      # }
      # #----------------------------------------------------------------------------------------------------------------
      
      # Calculate daily parameter values with SD
      Days <- seq(1,nrow(data.set), length.out=nrow(data.set))
      param.daily = param.final[1,]
      
      if (no.param == 1) {
        for (i in 2:length(Days)) {
          param.daily[i,] = param.final[1,]
        }
      }
      if (no.param == 2) {
        for (i in 2:length(Days)) {
          param.daily[i,1:no.var] = param.final[1,1:no.var] + param.final[2,1:no.var] * i
        }
        for (i in (no.var+1):(2*no.var)) {
          param.daily[,i] = ((param.final[1,i]^2 + param.final[2,i]^2)/2)^0.5
        }
      }
      if (no.param == 3) {
        for (i in 2:length(Days)) {
          param.daily[i,1:no.var] = param.final[1,1:no.var] + param.final[2,1:no.var] * i + param.final[3,1:no.var] * i^2
        }
        for (i in (no.var+1):(2*no.var)) {
          param.daily[,i] = ((param.final[1,i]^2 + param.final[2,i]^2 + param.final[3,i]^2)/3)^0.5
        }
      }
      param.daily$ar = 1 - param.daily$af - param.daily$as
      param.daily$ar_SD = with(param.daily, ((af_SD*af_SD + as_SD*as_SD)/2)^0.5)
      param.daily$Date = as.Date(data.set$Date)
      
      
      # Plotting the parameter sets over time
      if (with.storage==T) { 
        melted.param1 = melt(param.daily[,c("k","Y","af","as","ar","sf","Date")], id.vars="Date")
        melted.param2 = melt(param.daily[,c("k_SD","Y_SD","af_SD","as_SD","ar_SD","sf_SD","Date")], id.vars="Date")
      } else {
        melted.param1 = melt(param.daily[,c("Y","af","as","ar","sf","Date")], id.vars="Date")
        melted.param2 = melt(param.daily[,c("Y_SD","af_SD","as_SD","ar_SD","sf_SD","Date")], id.vars="Date")
      }
      melted.param = data.frame(melted.param1$Date, melted.param1$variable, melted.param1$value, melted.param2$value)
      names(melted.param) = c("Date","variable","Parameter","Parameter_SD")
      melted.param$Date = as.Date(melted.param$Date)
      melted.param$volume = list(vol[unlist(vol.group[v1])])
      melted.param$volume.group = as.factor(v1)
      melted.param$no.param = as.factor(no.param.par.var[z])
      
      
      # Plotting the Measured (data) vs Modelled Plant Carbon pools for plotting and comparison
      #----------------------------------------------------------------------------------------------------------------
      # lm.daily = read.csv("processed_data/Cleaf_daily_data.csv") # Unit gC per gC plant
      #----------------------------------------------------------------------------------------------------------------
      
      for (j in 1:length(v)) {
        data.set = subset(data,(volume %in% vol[v[j]]))
        
        #----------------------------------------------------------------------------------------------------------------
        # leafmass.daily = subset(lm.daily,(volume %in% vol[v[j]]))
        # leafmass.daily$Date = as.Date(leafmass.daily$Date)
        # # leafmass.daily.gas = leafmass.daily[leafmass.daily$Date %in% c(unique(as.Date(tnc.data.processed$Date))), ]
        # # browser()
        # data.set$Sleaf = data.set$Sleaf / leafmass.daily$leafmass * 100
        # data.set$Sleaf_SD = ((data.set$Sleaf_SD*data.set$Sleaf_SD + leafmass.daily$leafmass_SE*leafmass.daily$leafmass_SE)/2)^0.5 / leafmass.daily$leafmass * 100
        #----------------------------------------------------------------------------------------------------------------
        
        output.final.set = subset(output.final,(volume %in% vol[v[j]]))
        output.final.set$Date = data.set$Date
        
        if (with.storage==T) { 
          names(output.final.set) = c("Cstorage.modelled","Mleaf.modelled","Mstem.modelled","Mroot.modelled","Sleaf.modelled","volume","Date")
          melted.output = melt(output.final.set[,c("Mleaf.modelled","Mstem.modelled","Mroot.modelled","Sleaf.modelled","Date")], id.vars="Date")
          melted.data = melt(data.set[ , c("Mleaf","Mstem","Mroot","Sleaf","Date")], id.vars="Date")
          melted.error = melt(data.set[ , c("Mleaf_SD","Mstem_SD","Mroot_SD","Sleaf_SD","Date")], id.vars="Date")
        } else {
          names(output.final.set) = c("Mleaf.modelled","Mstem.modelled","Mroot.modelled","volume","Date")
          melted.output = melt(output.final.set[,c("Mleaf.modelled","Mstem.modelled","Mroot.modelled","Date")], id.vars="Date")
          melted.data = melt(data.set[ , c("Mleaf","Mstem","Mroot","Date")], id.vars="Date")
          melted.error = melt(data.set[ , c("Mleaf_SD","Mstem_SD","Mroot_SD","Date")], id.vars="Date")
        }
        melted.output$Date = as.Date(melted.output$Date)
        melted.output$volume = as.factor(vol[v[j]])
        melted.output$no.param = as.factor(no.param.par.var[z])
        
        if (with.storage==T) { 
          melted.Cstorage = output.final.set[,c("Cstorage.modelled","Date")]
          melted.Cstorage$Date = as.Date(melted.Cstorage$Date)
          melted.Cstorage$volume = as.factor(vol[v[j]])
          melted.Cstorage$no.param = as.factor(no.param.par.var[z])
        }
        
        melted.data$Date = as.Date(melted.data$Date)
        melted.data$volume = as.factor(vol[v[j]])
        
        melted.error$Date = as.Date(melted.error$Date)
        melted.error$volume = as.factor(vol[v[j]])
        melted.error$parameter = melted.data$value
        melted.error$no.param = as.factor(no.param.par.var[z])
        
        if (v1 < 8){
          melted.output$volume.group = as.factor(1)
          if (with.storage==T) { 
            melted.Cstorage$volume.group = as.factor(1)
          }
          melted.error$volume.group = as.factor(1)
        }
        if (v1 == 8){
          melted.output$volume.group = as.factor(2)
          if (with.storage==T) { 
            melted.Cstorage$volume.group = as.factor(2)
          }
          melted.error$volume.group = as.factor(2)
        }
        
        
        # Storing the summary of this volume group of data, outputs, Cstorage (Parameter is same for the group, will be stored later)
        if (j == 1) {
          summary.data.set = melted.data
          summary.error.set = melted.error
          summary.output.set = melted.output
          if (with.storage==T) { 
            summary.Cstorage.set = melted.Cstorage
          }
        }
        if (j > 1) {
          summary.output.set = rbind(summary.output.set,melted.output)
          if (with.storage==T) { 
            summary.Cstorage.set = rbind(summary.Cstorage.set,melted.Cstorage)
          }
          summary.error.set = rbind(summary.error.set,melted.error)
          if (z == 1) {
            summary.data.set = rbind(summary.data.set,melted.data)
          }
        }
      }
      
      
      # Storing the summary of all volume group's data, outputs, Cstorage, parameters
      if (q == 1) {
        summary.data = summary.data.set
        summary.error = summary.error.set
        summary.output = summary.output.set
        if (with.storage==T) { 
          summary.Cstorage = summary.Cstorage.set
        }
        summary.param = melted.param
      }
      if (q > 1) {
        summary.output = rbind(summary.output,summary.output.set)
        if (with.storage==T) { 
          summary.Cstorage = rbind(summary.Cstorage,summary.Cstorage.set)
        }
        summary.param = rbind(summary.param,melted.param)
        summary.error = rbind(summary.error,summary.error.set)
        if (z == 1) {
          summary.data = rbind(summary.data,summary.data.set)
        }
      }
      
      # Display the Acceptance rate of the chain
      nAccepted = length(unique(pChain[,1]))
      acceptance = (paste("Volume =",vol[v],", Total Parameter number =",no.param.par.var[z],": ", nAccepted, "out of ", chainLength-bunr_in, "candidates accepted ( = ",
                          round(100*nAccepted/chainLength), "%)"))
      print(acceptance)
      
      
      # # Plotting all parameter whole iterations for Day 1 only to check the convergance
      # png(file = paste("output/Parameter_iterations_day1_",v1,"_vol_",vol[v],"_par_",no.param.par.var[z], ".png", sep = ""))
      # par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
      # plot(pChain[,1],col="red",main="Utilization coefficient at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="k",ylim=c(param.k[1,1],param.k[1,3]))
      # plot(pChain[,1+no.param],col="green",main="Alloc frac to Biomass at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="Y",ylim=c(param.Y[1,1],param.Y[1,3]))
      # plot(pChain[,1+2*no.param],col="magenta",main="Alloc frac to foliage at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="af",ylim=c(param.af[1,1],param.af[1,3]))
      # plot(pChain[,1+3*no.param],col="blue",main="Alloc frac to stem at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="as",ylim=c(param.as[1,1],param.as[1,3]))
      # plot(pChain[,1+4*no.param],col="green",main="Foliage turnover at Day 1",cex.lab = 1.5,xlab="Iterations",ylab="sf",ylim=c(param.sf[1,1],param.sf[1,3]))
      # plot(pChain[,1+5*no.param],col="magenta",main="Log-likelihood",cex.lab = 1.5,xlab="Iterations",ylab="Log-likelihood")
      # title(main = paste("First day Parameter iterations for volume group",v1,"with par",no.param.par.var[z]), outer=TRUE, cex = 1.5)
      # dev.off()
      
      # Calculate LogLi, AIC, BIC, Time to find the most accurate model for best balance between model fit and complexity
      output.final1 = output.final
      if (with.storage==T) { 
        names(output.final1) = c("Cstorage","Mleaf","Mstem","Mroot","Sleaf","volume") # Rename for the logLikelihood function
      } else {
        names(output.final1) = c("Mleaf","Mstem","Mroot","volume")
      }
      
      data = data[with(data, order(volume)), ]
      row.names(data) = c(1:nrow(data))
      aic.bic[q,1] <- logLikelihood(data,output.final1,model.comparison) # Calculate logLikelihood
      
      k1 = 2 # k = 2 for the usual AIC
      npar = no.param*no.var # npar = total number of parameters in the fitted model
      aic.bic[q,2] = -2*aic.bic[q,1] + k1*npar
      
      if (model.comparison==F) {
        n = sum(!is.na(data$Sleaf)) + sum(!is.na(data$Mleaf)) + sum(!is.na(data$Mstem)) + sum(!is.na(data$Mroot))
      } else {
        n = sum(!is.na(data$Mleaf)) + sum(!is.na(data$Mstem)) + sum(!is.na(data$Mroot))
      }
      k2 = log(n) # n being the number of observations for the so-called BIC
      aic.bic[q,3] = -2*aic.bic[q,1] + k2*npar
      
      time$end.time[q] <- Sys.time()
      time$time.taken[q] <- time$end.time[q] - time$start.time[q]
      aic.bic[q,4] = time$time.taken[q]
      aic.bic$volume.group[q] = v1
      aic.bic[q,6] = no.param.par.var[z]
      aic.bic$volume[q] = list(vol[unlist(vol.group[v1])])
    }
  }
  bic = data.frame(aic.bic[,c("bic","volume.group","no.param","volume")])
  melted.aic.bic = melt(aic.bic[,c(1:6)], id.vars=c("no.param","volume.group"))
  
  # if (model.comparison==T | model.optimization==T) {
  if (model.optimization==T) {
    return(bic)
  } else if (model.optimization==F & with.storage==T) {
    result = list(no.param.par.var,summary.param,summary.data,summary.output,summary.error,bic,summary.Cstorage)
    return(result)
  } else {
    result = list(no.param.par.var,summary.param,summary.data,summary.output,summary.error,bic)
    return(result)
  }
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#-- Function to run Carbon Balance Model (CBM). 
#-----------------------------------------------------------------------------------------
# This script define the model equations to carry out Bayesian calibration for 
# 5 variables (allocation fractions: "k","Y",af","as","sf") on 
# various temporal scales to estimate Carbon pools (Cstorage,Cleaf,Cstem,Croot)
#-----------------------------------------------------------------------------------------
# Defining the model to iteratively calculate Cstorage, Cleaf, Cstem, Croot, Sleaf, Sstem, Sroot
model <- function (GPP,Rd,no.param,Mleaf,Mstem,Mroot,tnc,Y,k,af,as,sf) {
  Cstorage = Sleaf = Sstem = Sroot = c()
  
  # From Duan's experiment for TNC partitioning to tree organs
  # Leaf TNC C / Leaf C =  0.1167851; Stem TNC C / Stem C =  0.03782242; Root TNC C / Root C =  0.01795031
  Sleaf[1] = Mleaf[1] * tnc$leaf.tnc.C[7]/100
  Sstem[1] = Mstem[1] * tnc$stem.tnc.C[7]/100
  Sroot[1] = Mroot[1] * tnc$root.tnc.C[7]/100
  Cstorage[1] <- Sleaf[1] + Sstem[1] + Sroot[1] 
  
  Cleaf <- Croot <- Cstem <- c()
  Cleaf[1] <- Mleaf[1] - Sleaf[1]
  Cstem[1] <- Mstem[1] - Sstem[1]
  Croot[1] <- Mroot[1] - Sroot[1]
  
  if (no.param == 1) {
    k.i = k[1]; Y.i = Y[1]; af.i = af[1]; as.i = as[1]; sf.i = sf[1]
  }
  for (i in 2:length(GPP)) {
    if (no.param == 2) {
      k.i = k[1] + k[2]*i; Y.i = Y[1]+ Y[2]*i; af.i = af[1]+ af[2]*i; as.i = as[1]+ as[2]*i; sf.i = sf[1]+ sf[2]*i
    }
    if (no.param == 3) {
      k.i = k[1] + k[2]*i + k[3]*i*i; Y.i = Y[1]+ Y[2]*i + Y[3]*i*i; af.i = af[1]+ af[2]*i + af[3]*i*i; 
      as.i = as[1]+ as[2]*i + as[3]*i*i; 
      sf.i = sf[1]+ sf[2]*i + sf[3]*i*i
    }
    Cstorage[i] <- Cstorage[i-1] + GPP[i-1] - Rd[i-1]*(Mleaf[i-1] + Mroot[i-1] + Mstem[i-1]) - k.i*Cstorage[i-1]
    Sleaf[i] <- Cstorage[i] * tnc$leaf_to_all[7]/100 # 75% of storage goes to leaf (Duan's experiment)
    Sstem[i] <- Cstorage[i] * tnc$stem_to_all[7]/100 # 16% of storage goes to stem (Duan's experiment)
    Sroot[i] <- Cstorage[i] * tnc$root_to_all[7]/100 # 9% of storage goes to root (Duan's experiment)
    
    Cleaf[i] <- Cleaf[i-1] + k.i*Cstorage[i-1]*af.i*(1-Y.i) - sf.i*Mleaf[i-1]
    Cstem[i] <- Cstem[i-1] + k.i*Cstorage[i-1]*as.i*(1-Y.i)
    Croot[i] <- Croot[i-1] + k.i*Cstorage[i-1]*(1-af.i-as.i)*(1-Y.i)
    
    Mleaf[i] <- Cleaf[i] + Sleaf[i]
    Mstem[i] <- Cstem[i] + Sstem[i]
    Mroot[i] <- Croot[i] + Sroot[i]
  }
  output = data.frame(Cstorage,Mleaf,Mstem,Mroot,Sleaf)
  
  return(output)
}
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

#-- Function to run Carbon Balance Model (CBM) without considering storage pool. 
#-----------------------------------------------------------------------------------------
# This script define the model equations to carry out Bayesian calibration for 
# 4 variables (allocation fractions: "Y",af","as","sf") on 
# various temporal scales to estimate Carbon pools (Cleaf,Cstem,Croot)
#-----------------------------------------------------------------------------------------
# This version does not consider storage pool

# Defining the model to iteratively calculate Cleaf, Cstem, Croot
model.without.storage <- function (GPP,Rd,no.param,Mleaf,Mstem,Mroot,Y,af,as,sf) {
  
  if (no.param == 1) {
    Y.i = Y[1]; af.i = af[1]; as.i = as[1]; sf.i = sf[1]
  }
  
  for (i in 2:length(GPP)) {
    if (no.param == 2) {
      Y.i = Y[1]+ Y[2]*i; af.i = af[1]+ af[2]*i; as.i = as[1]+ as[2]*i; sf.i = sf[1]+ sf[2]*i
    }
    if (no.param == 3) {
      Y.i = Y[1]+ Y[2]*i + Y[3]*i*i; af.i = af[1]+ af[2]*i + af[3]*i*i; 
      as.i = as[1]+ as[2]*i + as[3]*i*i; 
      sf.i = sf[1]+ sf[2]*i + sf[3]*i*i
    }
    
    Mleaf[i] <- Mleaf[i-1] + (GPP[i-1] - Rd[i-1]*(Mleaf[i-1] + Mroot[i-1] + Mstem[i-1])) * af.i*(1-Y.i) - sf.i*Mleaf[i-1]
    Mstem[i] <- Mstem[i-1] + (GPP[i-1] - Rd[i-1]*(Mleaf[i-1] + Mroot[i-1] + Mstem[i-1])) * as.i*(1-Y.i)
    Mroot[i] <- Mroot[i-1] + (GPP[i-1] - Rd[i-1]*(Mleaf[i-1] + Mroot[i-1] + Mstem[i-1])) * (1-af.i-as.i)*(1-Y.i)
  }
  output = data.frame(Mleaf,Mstem,Mroot)
  return(output)
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
################ Figure 1 #####################
# Script to read and gerenate the model diagram
plot.model <- function() { 
  img <- readPNG("raw_data/Figure_1.png")
  
  #get size
  h<-dim(img)[1]
  w<-dim(img)[2]
  
  #open new file for saving the image in "output" folder
  png("output/Figure_1_CBM.png", width=w, height=h)
  par(mar=c(0,0,0,0), xpd=NA, mgp=c(0,0,0), oma=c(0,0,0,0), ann=F)
  plot.new()
  plot.window(0:1, 0:1)
  
  #fill plot with image
  usr<-par("usr")    
  rasterImage(img, usr[1], usr[3], usr[2], usr[4])
  
  #close image
  dev.off()
}
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

################ Figure 2 #####################
# Plot Model Measures ("bic") against "models with and without storage pool" 
# - 3 Grouped treatments with various number of parameters (Constant, Linear and Quadratic)
#-------------------------------------------------------------------------------------
plot.with.without.storage <- function(bic.with.storage, bic.without.storage) { 
  keeps = c("bic", "volume.group", "no.param")
  bic.with.storage = bic.with.storage[ , keeps, drop = FALSE]
  bic.without.storage = bic.without.storage[ , keeps, drop = FALSE]
  bic.with.storage.sub = subset(bic.with.storage, no.param %in% 3)
  bic.without.storage.sub = subset(bic.without.storage, no.param %in% 3)
  bic = data.frame(bic.without.storage.sub$bic, bic.with.storage.sub$bic, bic.with.storage.sub$volume.group)
  names(bic)[1:3] <- c("without storage", "with storage", "Treatment")
  bic$Group = ifelse(bic$Treatment==1,"Small",ifelse(bic$Treatment==2,"Large","Free"))
  bic$Group = factor(bic$Group, levels = bic$Group[c(1,2,3)])
  
  keeps = c("without storage", "with storage", "Group")
  bic = bic[ , keeps, drop = FALSE]
  bic.melt <- melt(bic, id.vars = "Group")
  names(bic.melt)[3] <- "BIC"
  
  p1 = ggplot(data = bic.melt, aes(x = Group, y = BIC, group = variable, fill = variable)) +
    # ggplot(data = bic.melt, aes(x=factor(bic.melt$Treatment, levels=unique(as.character(bic.melt$Treatment))), y=BIC, fill = Model_setting)) +
    geom_bar(stat="identity", width = 0.5, position = "dodge") +
    scale_fill_brewer(palette = "Set1", name = "") +
    ylab("BIC") +
    theme_bw() +
    theme(legend.title = element_text(colour="black", size=10)) +
    theme(legend.text = element_text(colour="black", size=10)) +
    theme(legend.position = c(0.75,0.82)) +
    theme(legend.key.height=unit(1,"line")) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=12)) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = 12, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  png("output/Figure_2_bic.with.without.storage.png", units="px", width=800, height=600, res=220)
  print (p1)
  dev.off()
  
  # bic = data.frame(bic.without.storage$bic, bic.with.storage$bic, bic.with.storage$volume.group, bic.with.storage$no.param)
  # names(bic)[1:4] <- c("without storage", "with storage", "Treatment", "no.param")
  # bic$Group = ifelse(bic$Treatment==1,"Group 1",ifelse(bic$Treatment==2,"Group 2","Group 3"))
  # bic$Group = as.factor(bic$Group)
  # bic$Parameter.setting = ifelse(bic$no.param==1,"Constant",ifelse(bic$no.param==2,"Linear","Quadratic"))
  # bic$Parameter.setting = as.factor(bic$Parameter.setting)
  
  # keeps = c("without storage", "with storage", "Group", "Parameter.setting")
  # bic = bic[ , keeps, drop = FALSE]
  # bic.melt <- melt(bic, id.vars = c("Group", "Parameter.setting"))
  # names(bic.melt)[3:4] <- c("Model_setting", "BIC")
  #-------------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------------
  # p1 = ggplot(data = bic.melt, aes(x = Parameter.setting, y = BIC, group = Model_setting, fill = Model_setting)) +
  #   # ggplot(data = bic.melt, aes(x=factor(bic.melt$Treatment, levels=unique(as.character(bic.melt$Treatment))), y=BIC, fill = Model_setting)) +
  #   geom_bar(stat="identity", width = 0.5, position = "dodge") +
  #   facet_grid(. ~ Group) +
  #   scale_fill_brewer(palette = "Set1", name = "Model setting") +
  #   xlab("Parameter setting") +
  #   ylab("BIC") +
  #   # ggtitle("BIC for various model settings") +
  #   theme_bw() +
  #   # theme(plot.title = element_text(size = 12)) +
  #   theme(legend.title = element_text(colour="black", size=10)) +
  #   theme(legend.text = element_text(colour="black", size=10)) +
  #   theme(legend.position = c(0.85,0.8)) +
  #   theme(legend.key.height=unit(1,"line")) +
  #   theme(legend.key = element_blank()) +
  #   theme(text = element_text(size=12)) +
  #   theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
  #   theme(axis.title.y = element_text(size = 12, vjust=0.3)) +
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  # 
  # png("output/Figure_2_bic.with.without.storage.png", units="px", width=2000, height=1600, res=220)
  # print (p1)
  # dev.off()
}
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

################ Figure 3 #####################
# Plot Model Measures ("bic") against "Treatment groupings" and "number of parameters"
#-------------------------------------------------------------------------------------
# cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
plot.parameter.settings <- function(bic.group1, bic.group2, bic.group3, bic.group4) { 
  cbPalette = c("gray", "orange", "skyblue", "green3", "yellow3", "#0072B2", "#D55E00")
  bic.group = bic.group1[rep(rownames(bic.group1), length(bic.group1$volume[[1]])), ]
  bic.group = bic.group[with(bic.group, order(no.param)), ]
  bic.group$volume = unlist(bic.group1$volume)
  
  bic.group2.order = bic.group2[rep(rownames(bic.group2[1:3,]), length(bic.group2$volume[[1]])), ]
  bic.group2.order = bic.group2.order[with(bic.group2.order, order(no.param)), ]
  bic.group2.order = rbind(bic.group2.order, bic.group2[bic.group2$volume.group==2, ])
  bic.group2.order$volume = unlist(bic.group2$volume)
  bic.group2.order$volume.group = 2
  bic.group = rbind(bic.group, bic.group2.order)
  
  bic.group3.order = bic.group3[rep(rownames(bic.group3[1:6,]), length(bic.group3$volume[[1]])), ]
  bic.group3.order = bic.group3.order[with(bic.group3.order, order(volume.group,no.param)), ]
  bic.group3.order = rbind(bic.group3.order, bic.group3[bic.group3$volume.group==3, ])
  bic.group3.order$volume = unlist(bic.group3$volume)
  bic.group3.order$volume.group = 3
  bic.group = rbind(bic.group, bic.group3.order)
  
  bic.group4$volume.group = 4
  bic.group = rbind(bic.group, bic.group4)
  bic.group$volume = unlist(bic.group$volume)
  bic.group$volume = as.factor(bic.group$volume)
  #-------------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------------
  pd <- position_dodge(0.3)
  p2 = ggplot(data = bic.group, aes(x = factor(bic.group$volume, levels=unique(as.character(bic.group$volume)) ), y = bic, group = interaction(volume.group,no.param), colour=factor(volume.group), shape=factor(no.param))) +
    geom_line(position=pd, data = bic.group, aes(x = factor(bic.group$volume, levels=unique(as.character(bic.group$volume)) ), y = bic, group = interaction(volume.group,no.param), colour=factor(volume.group), linetype=factor(no.param))) +
    geom_point(position=pd, size=2) +
    xlab("Treatments") +
    labs(colour="Grouping options", shape="Parameters", linetype="Parameters") +
    scale_y_continuous(name="BIC",limits = c(0, round_any(max(bic.group$bic), 50, f = ceiling)), breaks=seq(0,round_any(max(bic.group$bic), 50, f = ceiling),250)) +
    scale_x_discrete(name="Treatment", breaks=c("5", "10", "15", "20", "25", "35", "1000"),
                     labels=c("5L", "10L", "15L", "20L", "25L", "35L", "1000L")) +
    scale_color_manual(breaks = c("1", "2", "3", "4"), values=cbPalette[1:4]) +
    theme_bw() +
    theme(legend.title = element_text(colour="black", size=10)) +
    theme(legend.text = element_text(colour="black", size=10)) +
    theme(legend.key.height=unit(1,"line")) +
    theme(legend.key.width=unit(2,"line")) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=12)) +
    theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
    theme(axis.title.y = element_text(size = 12, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  png("output/Figure_3_bic_treat_group.png", units="px", width=1500, height=1200, res=220)
  print (p2)
  dev.off()
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

################ Figure 4 #####################
# Plot modelled parameters with 3 Grouped treatments and quadratic parameter setting
#-------------------------------------------------------------------------------------
plot.Modelled.parameters <- function(result) { 
  cbPalette = c("gray", "orange", "skyblue", "green3", "yellow3", "#0072B2", "#D55E00")
  i = 0
  font.size = 12
  plot = list() 
  var = as.factor(c("k","Y","af","as","ar","sf"))
  # var = as.factor(c("Y","af","as","ar","sf"))
  title = as.character(c("A","B","C","D","E","F"))
  pd <- position_dodge(0.5)
  no.param.par.var = result[[1]]
  summary.param = result[[2]]
  summary.data = result[[3]]
  summary.output = result[[4]]
  summary.error = result[[5]]
  
  for (p in 1:length(var)) {
    summary.param.set.limit = subset(summary.param, variable %in% var[p])
    for (z in 1:length(no.param.par.var)) {
      summary.param.set = subset(summary.param, variable %in% var[p] & no.param %in% no.param.par.var[z])
      i = i + 1
      plot[[i]] = ggplot(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume.group, colour=factor(volume.group))) +
        geom_ribbon(data = summary.param.set, aes(ymin=Parameter-Parameter_SD, ymax=Parameter+Parameter_SD), linetype=2, alpha=0.1,size=0.1) +
        geom_point(position=pd,size=0.01) +
        geom_line(position=pd,data = summary.param.set, aes(x = Date, y = Parameter,  group = volume.group, colour=factor(volume.group)),size=1) +
        # ylab(paste(as.character(var[p]),"(fraction)")) +
        ylab(paste(as.character(var[p]))) +
        labs(colour="Treatment Group") +
        scale_color_manual(breaks = c("1", "2", "3"), values=cbPalette[2:4]) +
        scale_y_continuous(limits = c(min(summary.param.set.limit$Parameter)-2*max(summary.param.set.limit$Parameter_SD),
                                      max(summary.param.set.limit$Parameter)+2*max(summary.param.set.limit$Parameter_SD))) +
        annotate("text", x = min(summary.param.set$Date), y = max(summary.param.set$Parameter) + 2*max(summary.param.set$Parameter_SD), size = font.size-7, label = paste(title[p])) +
        theme_bw() +
        theme(legend.title = element_text(colour="black", size=font.size)) +
        theme(legend.text = element_text(colour="black", size=font.size-1)) +
        theme(legend.position = c(0.22,0.18)) +
        theme(legend.key = element_blank()) +
        theme(text = element_text(size=font.size)) +
        theme(axis.title.x = element_blank()) +
        theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      
      if (p==1) {
        plot[[i]] = plot[[i]] + scale_colour_manual(name="", breaks=c("1", "2", "3"),
                                                      labels=c("Small", "Large", "Free"), values=cbPalette[2:4]) +
          ylab(expression(k~"(g C "*g^"-1"*" C "*d^"-1"*")"))
        plot[[i]] = plot[[i]] + theme(legend.key.height=unit(0.7,"line"))
      } else if (p>1) {
        plot[[i]] = plot[[i]] + guides(colour=FALSE)
      } 
      if (p==2) {
        plot[[i]] = plot[[i]] + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 0.8), units="line"))
      }
      if (p==3) {
        # plot[[i]] = plot[[i]] + ylab(expression(a[f]~"(fraction)")) + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 1), units="line"))
        plot[[i]] = plot[[i]] + ylab(expression(a[f])) + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 1), units="line"))
      }
      if (p==4) {
        plot[[i]] = plot[[i]] + ylab(expression(a[w])) + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 1), units="line"))
      }
      if (p==5) {
        plot[[i]] = plot[[i]] + ylab(expression(a[r])) + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 1), units="line"))
      }
      if (p==6) {
        plot[[i]] = plot[[i]] + ylab(expression(s[f]~"(g C "*g^"-1"*" C "*d^"-1"*")"))
      }
    }
  }
  
  png("output/Figure_4_modelled_parameters.png", units="px", width=2000, height=2000, res=250)
  print (do.call(grid.arrange,  plot))
  dev.off()
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

################ Figure 5 #####################
# Plot Daily analysis (lines) with optimum parameter setting and intermittent observations (symbols) of selected carbon stocks
#-------------------------------------------------------------------------------------
plot.Modelled.biomass <- function(result) { 
  cbPalette = c("gray", "orange", "skyblue", "green3", "yellow3", "#0072B2", "#D55E00")
  i = 0
  font.size = 12
  plot = list() 
  no.param.par.var = result[[1]]
  summary.param = result[[2]]
  summary.data = result[[3]]
  summary.output = result[[4]]
  summary.error = result[[5]]
  meas = as.factor(c("Mleaf","Mstem","Mroot","Sleaf"))
  res = as.factor(c("Mleaf.modelled","Mstem.modelled","Mroot.modelled","Sleaf.modelled"))
  error = as.factor(c("Mleaf_SD","Mstem_SD","Mroot_SD","Sleaf_SD"))
  title = as.character(c("A","B","C","D"))
  pd <- position_dodge(2) # move the overlapped errorbars horizontally
  for (p in 1:length(meas)) {
    summary.data.Cpool = subset(summary.data,variable %in% meas[p])
    summary.output.Cpool = subset(summary.output,variable %in% res[p])
    summary.error.Cpool = subset(summary.error,variable %in% error[p])
    
    i = i + 1
    plot[[i]] = ggplot(summary.error.Cpool, aes(x=Date, y=parameter, group = volume, colour=volume)) + 
      geom_point(position=pd) +
      geom_errorbar(position=pd,aes(ymin=parameter-value, ymax=parameter+value), colour="grey", width=2) +
      # geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = interaction(volume,volume.group,no.param), linetype=volume.group, colour=volume, size=no.param)) +
      # geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = interaction(volume,no.param), linetype=no.param, colour=volume)) +
      geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = volume, colour=volume)) +
      ylab(paste(as.character(meas[p]),"(g C)")) + xlab("Month") +
      # ggtitle("C pools - Measured (points) vs Modelled (lines)") +
      # labs(colour="Soil Volume", linetype="Grouping treatment", size="Total No of Parameter") +
      # labs(colour="Pot Volume (L)", linetype="No. of Parameters") +
      labs(colour="Pot Volume (L)") +
      scale_color_manual(breaks=c("5","10","15","20","25","35","1000"), values=cbPalette[1:7]) +
      # scale_color_manual(labels = c("Individuals", "One Group"), values = c("blue", "red")) +
      # coord_trans(y = "log10") + ylab(paste(as.character(meas[p]),"(g C plant-1)")) +
      theme_bw() +
      annotate("text", x = max(summary.output.Cpool$Date), y = min(summary.output.Cpool$value), size = font.size-7, label = paste(title[p])) +
      # theme(plot.title = element_text(size = 20, face = "bold")) +
      theme(legend.title = element_text(colour="black", size=font.size)) +
      theme(legend.text = element_text(colour="black", size = font.size-1)) +
      # theme(legend.key.height=unit(0.9,"line")) +
      theme(legend.position = c(0.18,0.73)) +
      theme(legend.key = element_blank()) +
      theme(text = element_text(size=font.size)) +
      theme(axis.title.x = element_blank()) +
      theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
      # theme(plot.title = element_text(hjust = 0)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    
    if (p==1) {
      plot[[i]] = plot[[i]] + scale_y_log10(breaks=c(.5,1,2,5,10,20),labels=c(.5,1,2,5,10,20)) + ylab(expression(C["t,f"]~"(g C "*plant^"-1"*")")) +
        scale_colour_manual(name="Treatments", breaks=c("5","10","15","20","25","35","1000"),
                              labels=c("5 L", "10 L", "15 L", "20 L", "25 L", "35 L", "FS"), values=cbPalette[1:7])
      plot[[i]] = plot[[i]]  + theme(legend.key.height=unit(0.7,"line"))
      
    } else if (p==2) {
      plot[[i]] = plot[[i]] + scale_y_log10(breaks=c(.5,1,2,5,10,20),labels=c(.5,1,2,5,10,20)) + ylab(expression(C["t,w"]~"(g C "*plant^"-1"*")"))
      plot[[i]] = plot[[i]] + theme(plot.margin=unit(c(0.4, 0.4, 0.4, 0.75), units="line"))
    } else if (p==3) {
      plot[[i]] = plot[[i]] + scale_y_log10(breaks=c(.5,1,2,5,10,20,40),labels=c(.5,1,2,5,10,20,40)) + ylab(expression(C["t,r"]~"(g C "*plant^"-1"*")"))
    } else {
      plot[[i]] = plot[[i]] + scale_y_log10(breaks=c(.05,.1,.2,.5,1,2),labels=c(.05,.1,.2,.5,1,2)) + ylab(expression(C["n,f"]~"(g C "*plant^"-1"*")"))
    }
    if (p>1) {
      plot[[i]] = plot[[i]] + guides(colour=FALSE)
    }
    
    # #----------------------------------------------------------------------------------------------------------------
    # # keeps <- c("Date", "volume", "tnc.conc", "tnc.conc_SE")
    # # tnc.data = tnc.data.processed[ , keeps, drop = FALSE]
    # 
    # if (p == 4) {
    #   plot[[i]] = ggplot(summary.error.Cpool, aes(x=Date, y=parameter, group = volume, colour=volume)) +
    #     geom_point(position=pd) +
    #     geom_errorbar(position=pd,aes(ymin=parameter-value, ymax=parameter+value), colour="grey", width=2) +
    #     geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = volume, colour=volume)) +
    #     ylab(paste(as.character(meas[p]),"(g C)")) + xlab("Month") +
    #     labs(colour="Pot Volume (L)") +
    #     theme_bw() +
    #     annotate("text", x = min(summary.output.Cpool$Date), y = max(summary.output.Cpool$value), size = font.size-7, label = paste(title[p])) +
    #     theme(legend.title = element_text(colour="black", size=font.size)) +
    #     theme(legend.text = element_text(colour="black", size = font.size)) +
    #     theme(legend.position = c(0.17,0.7)) +
    #     theme(legend.key = element_blank()) +
    #     theme(text = element_text(size=font.size)) +
    #     theme(axis.title.x = element_blank()) +
    #     theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
    #     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    #     ylab(expression(S[leaf]~"(% of"~M[leaf]~")")) + guides(colour=FALSE)
    # }
    # #----------------------------------------------------------------------------------------------------------------
    
  }
  
  png("output/Figure_5_modelled_biomass.png", units="px", width=1600, height=1300, res=220)
  print (do.call(grid.arrange,  plot))
  dev.off()
  
  # # #----------------------------------------------------------------------------------------------------------------
  # # # Represent Sleaf as a concentration of Mleaf instead of total mass
  # if (p == 4) {
  #   summary.output.Mleaf = subset(summary.output,variable %in% "Mleaf.modelled")
  #   summary.output.Sleaf = subset(summary.output,variable %in% "Sleaf.modelled")
  #   summary.error.Sleaf = subset(summary.error,variable %in% "Sleaf_SD")
  #   summary.output.Sleaf$value = summary.output.Sleaf$value / summary.output.Mleaf$value * 100
  #   summary.output.Sleaf = summary.output.Sleaf[,-c(5,6)]
  #   
  #   # summary.error.Sleaf$value = summary.error.Sleaf$value / lm.daily.m$leafmass * 100
  #   leafmass.daily = read.csv("processed_data/Cleaf_daily_data.csv") # Unit gC per gC plant
  #   leafmass.daily = leafmass.daily[with(leafmass.daily, order(volume,Date)), ]
  #   summary.error.Sleaf$value = ((summary.error.Sleaf$value*summary.error.Sleaf$value + leafmass.daily$leafmass_SE*leafmass.daily$leafmass_SE)/2)^0.5 / lm.daily.m$leafmass * 100
  #   summary.error.Sleaf$parameter = summary.error.Sleaf$parameter / lm.daily.m$leafmass * 100
  #   summary.error.Sleaf = summary.error.Sleaf[,-c(6,7)]
  #   
  #   pd <- position_dodge(4) # move the overlapped errorbars horizontally
  #   plot[[i]] = ggplot(summary.error.Sleaf, aes(x=Date, y=parameter, group = volume, colour=volume)) +
  #     geom_errorbar(position=pd,aes(ymin=parameter-value, ymax=parameter+value), colour="grey", width=0.2) +
  #     geom_line(position=pd,data = summary.output.Sleaf, aes(x = Date, y = value, group = volume, colour=volume)) +
  #     geom_point(position=pd) +
  #     # ylab("Sleaf (g C)") + xlab("Month") +
  #     ylab(paste(as.character(meas[p]),"(g C)")) + xlab("Month") +
  #     labs(colour="Pot Volume (L)") +
  #     theme_bw() +
  #     annotate("text", x = min(summary.output.Sleaf$Date), y = max(summary.output.Sleaf$value), size = font.size-7, label = paste(title[p])) +
  #     theme(legend.title = element_text(colour="black", size=font.size)) +
  #     theme(legend.text = element_text(colour="black", size = font.size)) +
  #     theme(legend.position = c(0.17,0.7)) +
  #     theme(legend.key = element_blank()) +
  #     theme(text = element_text(size=font.size)) +
  #     theme(axis.title.x = element_blank()) +
  #     theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  #     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #     ylab(expression(S[leaf]~"(% of"~M[leaf]~")")) + guides(colour=FALSE)
  # }
  # 
  # png("output/Figure_5_modelled_biomass_Sleaf_conc.png", units="px", width=2200, height=1600, res=220)
  # print (do.call(grid.arrange,  plot))
  # dev.off()
  # # #----------------------------------------------------------------------------------------------------------------
  
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
################ Figure 6A #####################
# Function to gerenate plot for Daily net C assimilation per unit leaf area (Cday) 
# cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
# cbPalette = c("gray", "orange", "skyblue", "black", "yellow", "vermilion", "reddishpurple")
plot.Cday <- function(Cday.set, iteration) { 
  ggplot(data = Cday.set, aes(x = Date, y = carbon_day,  group = volume, colour=factor(volume))) +
    geom_point(size=0.01) +
    geom_line(data = Cday.set, aes(x = Date, y = carbon_day,  group = volume, colour=factor(volume))) +
    ylab(expression(C[day]~"(g C "*d^"-1"*" "*leafarea^"-1"*")")) +
    scale_colour_manual(name="", breaks=c("5", "1000"), labels=c("5L", "FS"), values=cbPalette[1:2]) +
    annotate("text", x = min(Cday.set$Date), y = max(Cday.set$carbon_day)*0.98, size = font.size-7, label = paste(title[1])) +
    theme_bw() +
    theme(legend.position = c(0.85,0.85)) +
    theme(legend.title = element_blank()) +
    theme(legend.key = element_blank(), plot.margin=unit(c(0.25, 0.25, 0, 0.55), units="line")) +
    theme(text = element_text(size=font.size)) +
    theme(legend.key.height=unit(0.65,"line")) +
    # theme(legend.key.width=unit(2,"line")) +
    theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(plot.title = element_text(vjust=-2))
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
################ Figure 6B #####################
# Function to gerenate plot for Day respiration rate (Rd)
plot.Rd <- function(Rd.set, iteration) { 
  plot.shift[[2]] = ggplot(data = Rd.set, aes(x = Date, y = Rd_daily,  group = volume, colour=factor(volume))) +
    geom_point(size=0.01) +
    geom_line(data = Rd.set, aes(x = Date, y = Rd_daily,  group = volume, colour=factor(volume))) +
    # xlab("Month") +
    ylab(expression(R[d]~"(g C "*g^"-1"*" plant "*d^"-1"*")")) + 
    # ggtitle("B - Case 2: Rd (5L pot -> Free)") +
    scale_colour_manual(name="", breaks=c("5", "1000"), labels=c("5L", "FS"), values=cbPalette[2:3]) +
    # scale_y_continuous(limits = c(min(summary.param.set$Parameter), max(summary.param.set$Parameter)), breaks=c(.005,.01,.015),labels=c(.005,.01,.015)) +
    annotate("text", x = min(Rd.set$Date), y = max(Rd.set$Rd_daily)*0.98, size = font.size-7, label = paste(title[2])) +
    theme_bw() +
    theme(legend.position = c(0.85,0.85)) +
    # theme(plot.title = element_text(size = font.size)) +
    theme(legend.title = element_blank()) +
    theme(legend.key = element_blank(), plot.margin=unit(c(0.25, 0.25, 0, 0), units="line")) +
    theme(text = element_text(size=font.size)) +
    theme(legend.key.height=unit(0.65,"line")) +
    # theme(legend.key.width=unit(2,"line")) +
    theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
################ Figure 6C #####################
# Function to gerenate plot biomass allocation fractions (af, as, ar)
plot.allocation.fractions <- function(summary.param.set, iteration) { 
  ggplot(data = summary.param.set, aes(x = Date, y = Parameter,  group = interaction(volume.group,variable), colour=factor(volume.group))) +
    # p0 = ggplot(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume.group, colour=factor(volume.group))) +
    # geom_point(size=0.01) +
    # geom_line(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume.group, colour=factor(volume.group))) +
    geom_line(data = summary.param.set, aes(x = Date, y = Parameter,  group = interaction(volume.group,variable), colour=factor(volume.group), linetype=factor(variable))) +
    # xlab("Month") +
    # ggtitle(paste(title[p],"- Case",iteration,":",as.character(var[p]),"(5L pot -> Free)")) +
    scale_colour_manual(breaks=c("1", "3"), labels=c("5L", "FS"), values=cbPalette[iteration:(iteration+1)]) +
    # scale_linetype_manual(values=c("solid","dashes", "dotted")) +
    scale_linetype_manual(breaks=c("af","as","ar"), labels=c(expression(a[f]),expression(a[w]),expression(a[r])),values=c(19,17,15)) +
    scale_y_continuous(name=expression(Allocations~"(g C "*g^"-1"*" C "*d^"-1"*")"),limits = c(0,0.9), breaks=seq(0,1,0.2)) +
    annotate("text", x = min(summary.param.set$Date), y = 0.86, size = font.size-7, label = paste(title[iteration])) +
    theme_bw() +
    theme(legend.position = c(0.47,0.83),legend.direction = "vertical",legend.box = "horizontal") +
    # theme(plot.title = element_text(size = font.size)) +
    theme(legend.title = element_blank()) +
    theme(legend.key = element_blank(), plot.margin=unit(c(0.25, 0.25, 0, 1), units="line")) +
    theme(text = element_text(size=font.size)) +
    theme(legend.key.height=unit(0.75,"line")) +
    theme(legend.key.width=unit(2,"line")) +
    theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
################ Figure 6D #####################
# Function to gerenate plot for growth respiration rate (Y) 
plot.Y <- function(summary.param.set, iteration) { 
  ggplot(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume.group, colour=factor(volume.group))) +
    geom_point(size=0.01) +
    geom_line(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume.group, colour=factor(volume.group))) +
    # geom_line(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume, colour=factor(volume))) +
    # xlab("Month") +
    # ggtitle(paste(title[p],"- Case",iteration,":",as.character(var[p]),"(5L pot -> Free)")) +
    scale_colour_manual(breaks=c("1", "3"), labels=c("5L", "FS"), values=cbPalette[iteration:(iteration+1)]) +
    scale_y_continuous(limits = c(min(summary.param.set$Parameter)*0.95, max(summary.param.set$Parameter)*1.05)) +
    annotate("text", x = min(summary.param.set$Date), y = max(summary.param.set$Parameter)*1.04, size = font.size-7, label = paste(title[iteration])) +
    ylab(expression(Y~"(g C "*g^"-1"*" C "*d^"-1"*")")) +
    theme_bw() +
    theme(legend.position = c(0.85,0.55)) +
    # theme(plot.title = element_text(size = font.size)) +
    theme(legend.title = element_blank()) +
    theme(legend.key = element_blank(), plot.margin=unit(c(0.25, 0.25, 0, 0.65), units="line")) +
    theme(text = element_text(size=font.size)) +
    theme(legend.key.height=unit(0.75,"line")) +
    # theme(legend.key.width=unit(3,"line")) +
    theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
################ Figure 6E #####################
# Function to gerenate plot for leaf turnover rate (sf)
plot.sf <- function(summary.param.set, iteration) { 
  ggplot(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume.group, colour=factor(volume.group))) +
    geom_point(size=0.01) +
    geom_line(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume.group, colour=factor(volume.group))) +
    # geom_line(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume, colour=factor(volume))) +
    # xlab("Month") +
    # ggtitle(paste(title[p],"- Case",iteration,":",as.character(var[p]),"(5L pot -> Free)")) +
    scale_colour_manual(breaks=c("1", "3"), labels=c("5L", "FS"), values=cbPalette[iteration:(iteration+1)]) +
    # scale_y_continuous(limits = c(min(summary.param.set$Parameter)*0.95, max(summary.param.set$Parameter)*1.05), breaks=c(.005,.007,.009),labels=c(.005,.007,.009)) +
    scale_y_continuous(limits = c(min(summary.param.set$Parameter)*0.95, max(summary.param.set$Parameter)*1.05)) +
    annotate("text", x = min(summary.param.set$Date), y = max(summary.param.set$Parameter)*1.04, size = font.size-7, label = paste(title[iteration])) +
    ylab(expression(s[f]~"(g C "*g^"-1"*" C "*d^"-1"*")")) +
    theme_bw() +
    theme(legend.position = c(0.85,0.2)) +
    # theme(plot.title = element_text(size = font.size)) +
    theme(legend.title = element_blank()) +
    theme(legend.key = element_blank(), plot.margin=unit(c(0.25, 0.25, 0, 0.3), units="line")) +
    theme(text = element_text(size=font.size)) +
    theme(legend.key.height=unit(0.75,"line")) +
    # theme(legend.key.width=unit(3,"line")) +
    theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
################ Figure 6F #####################
# Function to gerenate plot for utilization coefficient (K)
plot.k <- function(summary.param.set, iteration) { 
  ggplot(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume.group, colour=factor(volume.group))) +
    geom_point(size=0.01) +
    geom_line(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume.group, colour=factor(volume.group))) +
    # geom_line(data = summary.param.set, aes(x = Date, y = Parameter,  group = volume, colour=factor(volume))) +
    # xlab("Month") +
    # ggtitle(paste(title[p],"- Case",iteration,":",as.character(var[p]),"(5L pot -> Free)")) +
    scale_colour_manual(breaks=c("1", "3"), labels=c("5L", "FS"), values=cbPalette[iteration:(iteration+1)]) +
    scale_y_continuous(limits = c(min(summary.param.set$Parameter)*0.95, max(summary.param.set$Parameter)*1.05)) +
    annotate("text", x = min(summary.param.set$Date), y = max(summary.param.set$Parameter)*1.04, size = font.size-7, label = paste(title[iteration])) +
    ylab(expression(k~"(g C "*g^"-1"*" C "*d^"-1"*")")) +
    theme_bw() +
    theme(legend.position = c(0.85,0.85)) +
    # theme(plot.title = element_text(size = font.size)) +
    theme(legend.title = element_blank()) +
    theme(legend.key = element_blank(), plot.margin=unit(c(0.25, 0.25, 0.25, 0.65), units="line")) +
    theme(text = element_text(size=font.size)) +
    theme(legend.key.height=unit(0.75,"line")) +
    # theme(legend.key.width=unit(3,"line")) +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# This sript plots leaf biomass pools (Figure 6G) for various test cases 
# with parameter shifted from potted seedling to free seedling
plot.Mleaf <- function(shift.output.Mleaf) { 
  ggplot() +
    geom_line(data = shift.output.Mleaf, aes(x = Date, y = value, group = Case, colour=Case), size=1) +
    geom_point(size=2) +
    ylab(expression("Foliage mass"~"(g C "*plant^"-1"*")")) + 
    annotate("text", x = min(shift.output.Mleaf$Date), y = max(shift.output.Mleaf$value), size = font.size-7, label = paste(title[7])) +
    scale_colour_manual(breaks=c("0","1","2","3","4","5","6"), labels=c("Baseline (5L)","+ Cday of FS",expression(+ R[d]~"of FS"),
                                                                        expression(+ (a[f] + a[w] + a[r])~"of FS"),"+ Y of FS",expression(+ s[f]~"of FS"),"+ k of FS (Complete FS)"), 
                        values=cbPalette) +
    theme_bw() +
    theme(legend.position = c(0.2,0.7),legend.text.align = 0) +
    theme(legend.title = element_blank()) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=font.size+2)) +
    theme(legend.key.height=unit(1,"line")) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# This sript plots stem biomass pools (Figure 6H) for various test cases 
# with parameter shifted from potted seedling to free seedling
plot.Mstem <- function(shift.output.Mstem) { 
  ggplot() +
    geom_line(data = shift.output.Mstem, aes(x = Date, y = value, group = Case, colour=Case), size=1) +
    geom_point(size=2) +
    ylab(expression("Wood mass"~"(g C "*plant^"-1"*")")) + 
    scale_colour_manual(breaks=c("0","1","2","3","4","5","6"), labels=c("Baseline (5L)","+ Cday of FS",expression(+ R[d]~"of FS"),
                                                                        expression(+ (a[f] + a[w] + a[r])~"of FS"),"+ Y of FS",expression(+ s[f]~"of FS"),"+ k of FS (Complete FS)"), 
                        values=cbPalette) +
    annotate("text", x = min(shift.output.Mstem$Date), y = max(shift.output.Mstem$value), size = font.size-7, label = paste(title[8])) +
    theme_bw() +
    theme(legend.position = c(0.2,0.7),legend.text.align = 0) +
    theme(legend.title = element_blank()) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=font.size+2)) +
    theme(legend.key.height=unit(1,"line")) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# This sript plots root biomass pools (Figure 6I) for various test cases 
# with parameter shifted from potted seedling to free seedling
plot.Mroot <- function(shift.output.Mroot) { 
  ggplot() +
    geom_line(data = shift.output.Mroot, aes(x = Date, y = value, group = Case, colour=Case), size=1) +
    geom_point(size=2) +
    ylab(expression("Root mass"~"(g C "*plant^"-1"*")")) + 
    scale_colour_manual(breaks=c("0","1","2","3","4","5","6"), labels=c("Baseline (5L)","+ Cday of FS",expression(+ R[d]~"of FS"),
                                                                        expression(+ (a[f] + a[w] + a[r])~"of FS"),"+ Y of FS",expression(+ s[f]~"of FS"),"+ k of FS (Complete FS)"), 
                        values=cbPalette) +
    annotate("text", x = min(shift.output.Mroot$Date), y = max(shift.output.Mroot$value), size = font.size-7, label = paste(title[9])) +
    theme_bw() +
    theme(legend.position = c(0.2,0.7),legend.text.align = 0) +
    theme(legend.title = element_blank()) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=font.size+2)) +
    theme(legend.key.height=unit(1,"line")) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# This sript plots root biomass pools (Figure 6I) for various test cases 
# with parameter shifted from potted seedling to free seedling
plot.biomass <- function(shift.output.biomass) { 
  ggplot() +
    geom_line(data = shift.output.biomass, aes(x = Date, y = value, group = Case, colour=Case), size=1) +
    geom_point(size=2) +
    ylab(expression("Total Biomass"~"(g C "*plant^"-1"*")")) + 
    scale_colour_manual(breaks=c("0","1","2","3","4","5","6"), labels=c("Baseline (5L)","+ Cday of FS",expression(+ R[d]~"of FS"),
                                                                        expression(+ (a[f] + a[w] + a[r])~"of FS"),"+ Y of FS",expression(+ s[f]~"of FS"),"+ k of FS (Complete FS)"), 
                        values=cbPalette) +
    # annotate("text", x = min(shift.output.biomass$Date), y = max(shift.output.biomass$value), size = font.size-7, label = paste(title[9])) +
    theme_bw() +
    theme(legend.position = c(0.2,0.7),legend.text.align = 0) +
    theme(legend.title = element_blank()) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=font.size+2)) +
    theme(legend.key.height=unit(1,"line")) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

# ################ Figure 7 #####################
# # Calculate total C partitioning for individual treatments
# # # Set working directory for saving figures
# # setwd("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/archive/processeddata")
# vol = c(5,10,15,20,25,35,1000)
# vol_group <- list(c(5,10,15),c(20,25,35),1000)
# 
# Ct.group = data.frame(matrix(ncol = 8, nrow = length(vol)))
# names(Ct.group) = c("GPP","Rgrowth","Rd","Cstorage","Croot","Cstem","Cleaf","Clit")
# 
# GPP.data = read.csv("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/rawdata/GPP.csv")
# GPP.data = GPP.data[with(GPP.data, order(Date,volume)), ]
# names(GPP.data)[3] = "GPP"
# Rd.data = read.csv("/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/rawdata/Rd.csv")
# Rd.data = Rd.data[with(Rd.data, order(Date,volume)), ]
# names(Rd.data)[3] = "Rd"
# Cstorage.data = read.csv("Cstorage_summary_group.csv")
# Cstorage.data = Cstorage.data[with(Cstorage.data, order(Date,volume)), 1:3]
# names(Cstorage.data)[1] = "Cstorage"
# data = merge(GPP.data, Rd.data, by=c("Date","volume"))
# data = merge(data, Cstorage.data, by=c("Date","volume"))
# param.summary = read.csv("param_summary_group.csv")
# output.summary = read.csv("output_summary_group.csv")
# 
# cpool = as.factor(c("Mleaf.modelled","Mstem.modelled","Mroot.modelled","Sleaf.modelled"))
# for (i in 1:length(cpool)) {
#   cpool.data = subset(output.summary, variable==cpool[i])
#   cpool.data = cpool.data[, c("Date","value","volume")]
#   cpool.data = cpool.data[with(cpool.data, order(Date,volume)), ]
#   names(cpool.data)[2] = as.character(cpool[i])
#   data = merge(data,cpool.data, all = TRUE)
# }
# 
# # Mleaf.data = subset(output.summary, variable==res[i])
# # Mleaf.data = Mleaf.data[, c("Date","value","volume")]
# # Mleaf.data = Mleaf.data[with(Mleaf.data, order(Date,volume)), ]
# # names(Mleaf.data)[2] = "Mleaf"
# # data = merge(data,Mleaf.data, all = TRUE)
# # 
# # Mstem.data = subset(output.summary, variable==res[2])
# # Mstem.data = Mstem.data[, c("Date","value","volume")]
# # Mstem.data = Mstem.data[with(Mstem.data, order(Date,volume)), ]
# # names(Mstem.data)[2] = "Mstem"
# # 
# # Mroot.data = subset(output.summary, variable==res[3])
# # Mroot.data = Mroot.data[, c("Date","value","volume")]
# # Mroot.data = Mroot.data[with(Mroot.data, order(Date,volume)), ]
# # names(Mroot.data)[2] = "Mroot"
# # 
# # Sleaf.data = subset(output.summary, variable==res[4])
# # Sleaf.data = Sleaf.data[, c("Date","value","volume")]
# # Sleaf.data = Sleaf.data[with(Sleaf.data, order(Date,volume)), ]
# # names(Sleaf.data)[2] = "Sleaf"
# 
# # data = merge(data,Mstem.data, all = TRUE)
# # data = merge(data,Mroot.data, all = TRUE)
# # data = merge(data,Sleaf.data, all = TRUE)
# # data = data[with(data, order(volume,Date)), ]
# 
# param = data.frame(matrix(ncol = 3, nrow = nrow(data)))
# names(param) = c("Date","volume","volume.group")
# param[,c(1,2)] = data[,c(1,2)]
# param$volume = as.factor(param$volume)
# for (i in 1:length(vol_group)) {
#   for (j in 1:length(vol_group[[i]])) {
#     param <- within(param, volume.group[volume == vol_group[[i]][j]] <- i)
#   }
# }
# 
# var = as.factor(c("k","Y","af","as","ar","sf"))
# for (i in 1:length(var)) {
#   ind.param = subset(param.summary, variable==var[i])
#   ind.param = ind.param[, c("Date","Parameter","volume.group")]
#   names(ind.param)[2] = as.character(var[i])
#   param <- merge(param,ind.param,by=c("Date","volume.group"))
# }
# 
# # combine data and parameters together
# data <- merge(data,param,by=c("Date","volume"))
# data = data[with(data, order(volume,Date)), ]
# 
# # v = 1
# # res = as.factor(c("Mleaf.modelled","Mstem.modelled","Mroot.modelled","Sleaf.modelled"))
# # # Data processing for different pot volumes (e.g. 5L, 10L, ....., 1000L)
# # for (v in 1:length(vol)) {
# #   GPP = subset(GPP.data,(volume %in% vol[v])) # Consider one pot volume at a time to run MCMC on CBM
# #   names(GPP)[3] = "GPP"
# #   Rd = subset(Rd.data,(volume %in% vol[v]))
# #   data = subset(output.summary,(volume %in% vol[v]))
# #   Mleaf = subset(data, variable==res[1])
# #   Mstem = subset(data, variable==res[2])
# #   Mroot = subset(data, variable==res[3])
# #   Sleaf = subset(data, variable==res[4])
# # }
# 
# dates = as.Date(c("2013-02-21","2013-03-21","2013-04-21","2013-05-21"))
# data$Date = as.Date(data$Date)
# for (i in 1:length(dates)) {
#   data.set = data[data$Date <= dates[i], ]
#   for (v in 1:length(vol)) {
#     Ct.group$GPP[v] = sum ( data.set$GPP[which(data.set$volume == vol[v])] )
#     Ct.group$Rd[v] = sum ( data.set$Rd[which(data.set$volume == vol[v])] * (data.set$Mleaf.modelled[which(data.set$volume == vol[v])] + 
#                                                                               data.set$Mstem.modelled[which(data.set$volume == vol[v])] + data.set$Mroot.modelled[which(data.set$volume == vol[v])]) )
#     Ct.group$Cstorage[v] = data.set$Cstorage[which(data.set$volume == vol[v] & data.set$Date == dates[i])]
#     Ct.group[v, c(5:7)] = data.set[which(data.set$volume == vol[v] & data.set$Date == dates[i]), 8:6] - data.set[which(data.set$volume == vol[v] & data.set$Date == as.Date("2013-01-21")), 6:8]
#     Ct.group$Clit[v] = sum ( data.set$sf [which(data.set$volume == vol[v])] * data.set$Mleaf.modelled [which(data.set$volume == vol[v])])
#     Ct.group$Rgrowth[v] = Ct.group$GPP[v] - sum(Ct.group[v,c(3:8)])
#     # Ct.group$Rgrowth[v] = sum ( (data.set$GPP[which(data.set$volume == vol[v])]) - (data.set$Rd[which(data.set$volume == vol[v])] * (data.set$Mleaf.modelled[which(data.set$volume == vol[v])] + 
#     #                                                                                                                                 data.set$Mstem.modelled[which(data.set$volume == vol[v])] + data.set$Mroot.modelled[which(data.set$volume == vol[v])])) 
#     # * (data.set$k [which(data.set$volume == vol[v])] * data.set$Y [which(data.set$volume == vol[v])]) )
#   }
#   
#   Ct.fraction.group = Ct.group[, c(2:8)]
#   Ct.fraction.group[,] = Ct.fraction.group[,] / Ct.group[, 1] * 100
#   row.names(Ct.fraction.group) <- c("5L","10L","15L","20L","25L","35L","Free")
#   
#   # png(file = "C Partitioning.png")
#   png(file = paste("C Partitioning_month_",i,".png",sep=""))
#   par(mfrow = c(1, 1), mar=c(5, 5, 5, 6))
#   barplot(as.matrix(t(Ct.fraction.group)), main = paste("Accumulated C from",as.Date("2013-01-21"),"to",dates[i]), ylab = "C Partitioning", xlab = "Treatments (Pot size)",  
#           col = rainbow(20),legend = colnames(Ct.fraction.group),
#           args.legend = list(x = "topright", bty = "n", inset=c(-0.25, 0)))
#   dev.off()
# }
# plots0 <- lapply(ll <- list.files(patt='.*[.]png'),function(x){
#   img <- as.raster(readPNG(x))
#   rasterGrob(img, interpolate = FALSE)
# })
# ggsave("Barplot_C_accumulation.pdf", marrangeGrob(grobs=plots0, nrow=2, ncol=2))
# write.csv(Ct.group, file = "/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/archive/processeddata/carbon_allocations.csv", row.names = FALSE)
# write.csv(Ct.fraction.group, file = "/Users/kashifmahmud/WSU/ARC_project/CBM_Kashif/archive/processeddata/carbon_allocation_fractions.csv", row.names = FALSE)


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------






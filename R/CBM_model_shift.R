# This sript runs the model equations for parameter shifting from potted seedling to free seedling

# Calculating model outputs
Mleaf = Mstem = Mroot = LA = c()
Mleaf[1] <- Mleaf.data.set$leafmass[1]
Mstem[1] <- Mstem.data.set$stemmass[1]
Mroot[1] <- Mroot.data.set$rootmass[1]
LA[1] <- data.set$LA[1]

k=param.casted$k; Y=param.casted$Y; af=param.casted$af; as=param.casted$as; ar=param.casted$ar; sf=param.casted$sf
GPP = Cstorage = Sleaf = Sstem = Sroot = M = c()

##########
# M <- sigma.data$b * LA.data$LA + sigma.data$intercept
# GPP <- LA.data$LA * Cday.data$Cday * M # calculate total daily C gain with self shading
##########

# From Duan's experiment for TNC partitioning to tree organs
# Leaf TNC/Leaf DW =  0.1401421; Stem TNC/Stem DW =  0.0453869; Root TNC/Root DW =  0.02154037
# Sleaf[1] = Mleaf[1] / 0.65 * 0.1401421
# Sstem[1] = Mstem[1] / 0.65 * 0.0453869
# Sroot[1] = Mroot[1] / 0.65 * 0.02154037
Sleaf[1] = Mleaf[1] * 0.1167851
Sstem[1] = Mstem[1] * 0.03782242
Sroot[1] = Mroot[1] * 0.01795031
Cstorage[1] <- Sleaf[1] + Sstem[1] + Sroot[1] 

Cleaf <- Croot <- Cstem <- c()
Cleaf[1] <- Mleaf[1] - Sleaf[1]
Cstem[1] <- Mstem[1] - Sstem[1]
Croot[1] <- Mroot[1] - Sroot[1]

# LA[1] <- LA.data$LA[1]

for (i in 2:nrow(Cday.data.set)) {
  M[i-1] <- sigma.data.set$b * LA[i-1] + sigma.data.set$intercept
  GPP[i-1] <- LA[i-1] * Cday.data.set$carbon_day[i-1] * M[i-1] # calculate total daily C gain with self shading
  Cstorage[i] <- Cstorage[i-1] + GPP[i-1] - Rd.data.set$Rd_daily[i-1]*(Mleaf[i-1] + Mroot[i-1] + Mstem[i-1]) - k[i-1]*Cstorage[i-1]
  
  # Cstorage[i] <- Cstorage[i-1] + GPP.data$GPP[i-1] - Rd[i-1]*(Mleaf[i-1] + Mroot[i-1] + Mstem[i-1]) - k[i-1]*Cstorage[i-1]
  Sleaf[i] <- Cstorage[i] * 0.75 # 75% of storage goes to leaf (Duan's experiment)
  Sstem[i] <- Cstorage[i] * 0.16 # 16% of storage goes to stem (Duan's experiment)
  Sroot[i] <- Cstorage[i] * 0.09 # 9% of storage goes to root (Duan's experiment)
  
  Cleaf[i] <- Cleaf[i-1] + k[i-1]*Cstorage[i-1]*af[i-1]*(1-Y[i-1]) - sf[i-1]*Mleaf[i-1]
  Cstem[i] <- Cstem[i-1] + k[i-1]*Cstorage[i-1]*as[i-1]*(1-Y[i-1])
  Croot[i] <- Croot[i-1] + k[i-1]*Cstorage[i-1]*(1-af[i-1]-as[i-1])*(1-Y[i-1])
  
  Mleaf[i] <- Cleaf[i] + Sleaf[i]
  Mstem[i] <- Cstem[i] + Sstem[i]
  Mroot[i] <- Croot[i] + Sroot[i]
  
  # Leaf area (t) = Leaf area (T) * Leaf count (t) / Leaf count (T); t = time, T = time of harvest
  # LA[i] <- ((leaf.data$final_LA / leaf.data$final_LM) * Mleaf[i] + (leaf.data$initial_LA / leaf.data$initial_LM) * Mleaf[i])/2 
  # LA[i] <- hd.final$SLA[which(hd.final$volume == 5)] * Mleaf[i]
  LA[i] <- sigma.data.set$SLA * Cleaf[i]
}
output.final = data.frame(Cstorage,Mleaf,Mstem,Mroot,Sleaf)

# Plant Carbon pools for various parameter sensitivity
output.final$Date = Cday.data.set$Date
names(output.final) = c("Cstorage","Mleaf","Mstem","Mroot","Sleaf","Date")
melted.output = melt(output.final[,c("Mleaf","Mstem","Mroot","Cstorage","Sleaf","Date")], id.vars="Date")
melted.Cstorage = output.final[,c("Cstorage","Date")]
melted.output$Date = as.Date(melted.output$Date)
melted.Cstorage$Date = as.Date(melted.Cstorage$Date)
melted.output$Case = as.factor(q)
melted.Cstorage$Case = as.factor(q)

# Storing the summary of data, outputs, Cstorage, parameters
if (q == 0) {
  shift.output = melted.output
  shift.Cstorage = melted.Cstorage
}
if (q > 0) {
  shift.output = rbind(shift.output,melted.output)
  shift.Cstorage = rbind(shift.Cstorage,melted.Cstorage)
}

################ Figure S2 #####################
# Calculate temporal evolution of CUE 
#-------------------------------------------------------------------------------------
cbPalette = c("gray", "orange", "skyblue", "#009E73", "yellow3", "#0072B2", "#D55E00")
vol_group <- list(c(5,10,15),c(20,25,35),1000)

GPP.data = GPP.data.processed[with(GPP.data.processed, order(Date,volume)), ]
names(GPP.data)[3] = "GPP"
GPP.data$Date = as.Date(GPP.data$Date)
Rd.data = Rd.data.processed[with(Rd.data.processed, order(Date,volume)), ]
names(Rd.data)[3] = "Rm"

param.summary = result[[2]]
output.summary = result[[4]]
Cstorage.data = result[[7]]
Cstorage.data = Cstorage.data[with(Cstorage.data, order(Date,volume)), 1:3]
names(Cstorage.data)[1] = "Cs"
data.part = merge(GPP.data, Rd.data, by=c("Date","volume"))
data.part = merge(data.part, Cstorage.data, by=c("Date","volume"))

cpool = as.factor(c("Mleaf.modelled","Mstem.modelled","Mroot.modelled","Sleaf.modelled"))
for (i in 1:length(cpool)) {
  cpool.data = subset(output.summary, variable==cpool[i])
  cpool.data = cpool.data[, c("Date","value","volume")]
  cpool.data = cpool.data[with(cpool.data, order(Date,volume)), ]
  names(cpool.data)[2] = as.character(cpool[i])
  data.part = merge(data.part,cpool.data, all = TRUE)
}

param = data.frame(matrix(ncol = 3, nrow = nrow(data.part)))
names(param) = c("Date","volume","volume.group")
param[,c(1,2)] = data.part[,c(1,2)]
# param$volume = as.integer(param$volume)
for (i in 1:length(vol_group)) {
  for (j in 1:length(vol_group[[i]])) {
    param <- within(param, volume.group[volume == vol_group[[i]][j]] <- i)
  }
}
# param$Date = as.Date(param$Date)
# param$volume.group = as.numeric(param$volume.group)

var = as.factor(c("k","Y","af","as","ar","sf"))
for (i in 1:length(var)) {
  ind.param = subset(param.summary, variable==var[i])
  ind.param = ind.param[, c("Date","Parameter","volume.group")]
  names(ind.param)[2] = as.character(var[i])
  param <- merge(param,ind.param,by=c("Date","volume.group"))
}

# combine data and parameters together
data.part <- merge(data.part,param,by=c("Date","volume"))
data.part = data.part[with(data.part, order(volume,Date)), ]

data.set = data.part[data.part$Date <= as.Date("2013-05-21"), ]
data.set$Rm.tot = data.set$Rm * (data.set$Mleaf.modelled + data.set$Mstem.modelled + data.set$Mroot.modelled)
data.set$total.biomass = data.set$Mleaf.modelled + data.set$Mstem.modelled + data.set$Mroot.modelled

keeps = c("Date", "volume", "GPP", "Rm.tot", "Y", "total.biomass")
data.set = data.set[ , keeps, drop = FALSE]

for (i in 1:length(vols)) {
  data.set.vol = subset(data.set, volume %in% vols[i])
  data.set.vol$biomass.growth.initial = 0
  data.set.vol$biomass.growth.initial[2:nrow(data.set.vol)] = diff(as.matrix(data.set.vol$total.biomass))
  data.set.vol$biomass.growth = data.set.vol$biomass.growth.initial
  for (j in 1:(nrow(data.set.vol)-1)) {
    if (data.set.vol$biomass.growth[j] < 0 ) {
      data.set.vol$biomass.growth[j+1] = data.set.vol$biomass.growth.initial[j+1] + data.set.vol$biomass.growth.initial[j]
      data.set.vol$biomass.growth[j] = 0
    }
  }
  data.set.vol$Rg = data.set.vol$biomass.growth * data.set.vol$Y 
  data.set.vol$cue = 1 - ((data.set.vol$Rm.tot + data.set.vol$Rg) / data.set.vol$GPP)
  data.set.vol$cue[data.set.vol$cue < 0] <- 0
  if (i == 1) {
    cue = data.set.vol
  } else {
    cue = rbind(cue,data.set.vol)
  }
}

cbPalette = c("gray", "orange", "skyblue", "green3", "yellow3", "#0072B2", "#D55E00")
font.size = 11
pd <- position_dodge(0)
p.cue = ggplot(data = cue, aes(x = Date, y = cue,  group = volume, colour=factor(volume))) +
  geom_line(position=pd,data = cue, aes(x = Date, y = cue,  group = volume, colour=factor(volume)),size=0.4) +
  ylab("CUE") + ylim(0,0.9) +
  scale_colour_manual(name="Treatments", breaks=c("5","10","15","20","25","35","1000"),
                      labels=c("5 L", "10 L", "15 L", "20 L", "25 L", "35 L", "FS"), values=cbPalette[1:7]) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=font.size)) +
  theme(legend.text = element_text(colour="black", size = font.size-1)) +
  # theme(legend.key.height=unit(0.9,"line")) +
  theme(legend.position = c(0.65,0.93),legend.direction = "horizontal") +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) + theme(legend.key.height=unit(0.75,"line")) + theme(legend.key.width=unit(0.75,"line")) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  # theme(plot.title = element_text(hjust = 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

png("output/Figure_S2_cue.png", units="px", width=600, height=500, res=120)
p.cue
dev.off()

################ Figure 7 #####################
# Calculate total C partitioning for individual treatments
#-------------------------------------------------------------------------------------
# cbPalette = c("gray", "orange", "skyblue", "green3", "yellow3", "#0072B2", "#D55E00")
cbPalette = c("gray", "orange", "skyblue", "#009E73", "yellow3", "#0072B2", "#D55E00")
vol_group <- list(c(5,10,15),c(20,25,35),1000)
Ct.group = data.frame(matrix(ncol = 8, nrow = length(vols)))
names(Ct.group) = c("GPP","Rg","Rm","Cs","Cr","Cw","Cf","Clit")

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
for (v in 1:length(vols)) {
  Ct.group$GPP[v] = sum ( data.set$GPP[which(data.set$volume == vols[v])] )
  Ct.group$Rm[v] = sum ( data.set$Rm[which(data.set$volume == vols[v])] * (data.set$Mleaf.modelled[which(data.set$volume == vols[v])] + 
                                                                             data.set$Mstem.modelled[which(data.set$volume == vols[v])] + data.set$Mroot.modelled[which(data.set$volume == vols[v])]) )
  Ct.group$Cs[v] = data.set$Cs[which(data.set$volume == vols[v] & data.set$Date == as.Date("2013-05-21"))]
  Ct.group[v, c(5:7)] = data.set[which(data.set$volume == vols[v] & data.set$Date == as.Date("2013-05-21")), 8:6] - data.set[which(data.set$volume == vols[v] & data.set$Date == as.Date("2013-01-21")), 6:8]
  Ct.group$Clit[v] = sum ( data.set$sf [which(data.set$volume == vols[v])] * data.set$Mleaf.modelled [which(data.set$volume == vols[v])])
  Ct.group$Rg[v] = Ct.group$GPP[v] - sum(Ct.group[v,c(3:8)])
}

Ct.fraction.group = Ct.group[, c(2:8)]
Ct.fraction.group[,] = Ct.fraction.group[,] / Ct.group[, 1] * 100
row.names(Ct.fraction.group) <- c("5L","10L","15L","20L","25L","35L","Free")

Ct.group$Volume = vols
colnames(Ct.group) <- c("GPP (g C)", "Rg (g C)", "Rm (g C)", "Cs (g C)", "Cr (g C)", "Cw (g C)", "Cf (g C)", "Clit (g C)", "Volume (L)")
Ct.group = Ct.group[,c(9,1,2,3,4,7,5,6,8)]
write.csv(Ct.group, file = "output/C_partitioning.csv", row.names = FALSE)

png("output/Figure_7_C_partitioning.png", units="px", width=1200, height=1000, res=200)
par(mfrow = c(1, 1), mar=c(5, 4, 2, 6))
# bb = barplot(as.matrix(t(Ct.fraction.group)), ylim=c(0, 107), ylab = "C Partitioning (%)", xlab = "Treatments (Container size)",  
#         col = rainbow(20),legend = colnames(Ct.fraction.group), 
#         args.legend = list(x = "topright", bty = "n", inset=c(-0.15, 0)))
bb = barplot(as.matrix(t(Ct.fraction.group)), ylim=c(0, 107), ylab = "C Partitioning (%)", xlab = "Container size (L))",  
             col = cbPalette[1:7],legend = c(expression(R[g]),expression(R["m,tot"]),expression(C[n]),expression(C["s,r"]),expression(C["s,w"]),expression(C["s,f"]),expression(C["t,lit"])), 
             args.legend = list(x = "topright", bty = "n", inset=c(-0.18, 0)))
# text( bb, Ct.fraction.group[,1]-3, labels = round(Ct.group[,3],1), cex=.9)
# text( bb, Ct.fraction.group[,1]+Ct.fraction.group[,2]-4, labels = round(Ct.group[,4],1), cex=.9)
# text( bb, Ct.fraction.group[,1]+Ct.fraction.group[,2]+Ct.fraction.group[,3]-1, labels = round(Ct.group[,5],1), cex=.9)
# text( bb, Ct.fraction.group[,1]+Ct.fraction.group[,2]+Ct.fraction.group[,3]+Ct.fraction.group[,4]-3, labels = round(Ct.group[,6],1), cex=.9)
# text( bb, Ct.fraction.group[,1]+Ct.fraction.group[,2]+Ct.fraction.group[,3]+Ct.fraction.group[,4]+Ct.fraction.group[,5]-3, labels = round(Ct.group[,7],1), cex=.9)
# text( bb, Ct.fraction.group[,1]+Ct.fraction.group[,2]+Ct.fraction.group[,3]+Ct.fraction.group[,4]+Ct.fraction.group[,5]+Ct.fraction.group[,6]-2, labels = round(Ct.group[,8],1), cex=.9)
# text( bb, Ct.fraction.group[,1]+Ct.fraction.group[,2]+Ct.fraction.group[,3]+Ct.fraction.group[,4]+Ct.fraction.group[,5]+Ct.fraction.group[,6]+Ct.fraction.group[,7]-1, labels = round(Ct.group[,9],1), cex=.9)
text( bb, rowSums(Ct.fraction.group)+0.5, labels = round(Ct.group[,2],1), pos = 3, cex=1, col="red")

dev.off()


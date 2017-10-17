# Test g1 values over time for various pot treatments
# Read g1 data directly from Court's Github folder
g1.data <- read.csv("https://raw.githubusercontent.com/kashifmahmud/EucPVE/master/calculated%20data/g1_pred.csv")
g1.data$Date = as.Date(g1.data$Date)
vols = c(5,10,15,20,25,35,1000)

# Import weekly height diameter data for 3 months (Diameter is in mm; Height is in cm)
height.dia <- read.csv("raw_data/height_diameter.csv")
height.dia$Date <- parse_date_time(height.dia$Date,"d m y")
height.dia$Date = as.Date(height.dia$Date, format = "%d/%m/%Y")
plot.summary = read.csv("raw_data/plot_summary.csv")
height.dia <- merge(plot.summary,height.dia,by=c("pot","plot"))
height.dia = height.dia[ ,c(-1,-2)]
height.dia = height.dia[with(height.dia, order(volume,Date)), ]
height.dia <- height.dia[c("Date", "volume", "height", "diameter")]

for(i in 1:length(vols)) {
  height.dia.idn = subset(height.dia,volume==vols[i]) 
  for(j in 1:length(unique(height.dia.idn$Date))) {
    height.dia.idn.date = subset(height.dia.idn, Date == unique(height.dia.idn$Date)[j])
    height.dia.idn.date[nrow(height.dia.idn.date)+1, 2:ncol(height.dia.idn.date)] = colMeans(height.dia.idn.date[2:ncol(height.dia.idn.date)], na.rm = TRUE) # R7 = Average of leaf data
    height.dia.idn.date[nrow(height.dia.idn.date)+1, 2:ncol(height.dia.idn.date)] = (apply(height.dia.idn.date[2:ncol(height.dia.idn.date)], 2, sd, na.rm = TRUE))/(nrow(height.dia.idn.date)-1)^0.5 # R8 = Standard error of leaf counts
    # height.dia.idn.date[nrow(height.dia.idn.date)+1, 2:ncol(height.dia.idn.date)] = apply(height.dia.idn.date[2:ncol(height.dia.idn.date)], 2, sd, na.rm = TRUE) # R8 = Standard deviation of leaf counts
    height.dia.idn.date$Date = height.dia.idn.date$Date[1]
    dimnames(height.dia.idn.date)[[1]] <- c(1:(nrow(height.dia.idn.date)-2), "Mean", "SE")
    if (i == 1 && j == 1) {
      height.dia.final <- height.dia.idn.date[0,]
    }
    height.dia.final[j+(i-1)*length(unique(height.dia.idn$Date)), ] <- height.dia.idn.date["Mean", ]
    height.dia.final$height_SE[j+(i-1)*length(unique(height.dia.idn$Date))] <- height.dia.idn.date["SE", 3]
    height.dia.final$dia_SE[j+(i-1)*length(unique(height.dia.idn$Date))] <- height.dia.idn.date["SE", 4]
  }
}

# interpolate heights on the dates of g1 measurements from fortnightly data
height.dia.final.g1 = data.frame(Date = as.Date(min(height.dia.final$Date):max(height.dia.final$Date)))
for(i in 1:length(vols)) {
  height.dia.final.sub = subset(height.dia.final, volume %in% vols[i])
  height.dia.final.g1$volume = vols[i]
  ptsLin   <- approx(height.dia.final.sub$Date, height.dia.final.sub$height, method="linear", n=(max(height.dia.final$Date)-min(height.dia.final$Date)))
  height.dia.final.g1$height = c(ptsLin$y[1],ptsLin$y)
  # height.dia.final.g1$height[(1+(i-1)*nrow(height.dia.final.sub)) : (i*nrow(height.dia.final.sub))] = c(ptsLin$y[1],ptsLin$y)
  ptsLin   <- approx(height.dia.final.sub$Date, height.dia.final.sub$height_SE, method="linear", n=(max(height.dia.final$Date)-min(height.dia.final$Date)))
  height.dia.final.g1$height_SE = c(ptsLin$y[1],ptsLin$y)
  if (i == 1) {
    g1.height = height.dia.final.g1
  }
  if (i > 1) {
    g1.height = rbind(g1.height,height.dia.final.g1)
  }
}
g1.data = merge(g1.height, g1.data, by = c("Date", "volume"))
keeps = c("Date","volume","height","g1_date")
g1.data = g1.data[ , keeps, drop = FALSE]
g1.data = unique(g1.data)

# Import daily LA data (unit = m^2)
LA <- read.csv("processed_data/LA_daily_data.csv")
g1.data = merge(g1.data, LA, by = c("Date", "volume"))


# png("output/Figure_g1_with_height.png", units="px", width=2000, height=1600, res=220)
# par(mfrow=c(4,2), mar=c(4, 4, 0.5, 0.2))
# # produce plot for each of the first six species
# for(i in 1:length(vols)) {
#   g1.data.sub = unique(subset(g1.data, volume %in% vols[i]))
#   plot(g1.data.sub$height, g1.data.sub$g1_date, xlab="height", ylab="g1")
#   abline(lm(g1_date ~ height, g1.data.sub))
#   legend('bottomleft', paste("volume =",as.character(vols[i]), "L"), bty='n')
# }
# dev.off()

plots = list()
plots[[1]] = ggplot(data = g1.data, aes(x = Leaf_area, y = g1_date, group = volume, colour=factor(volume))) +
  geom_point(size=1) +
  stat_smooth(method = "lm") +
  ylab("g1") + xlab(expression(Leafarea~"("*m^"2"*")")) +
  theme_bw() + labs(colour="Pot Volume") +
  theme(legend.title = element_text(colour="black", size=10)) +
  theme(legend.text = element_text(colour="black", size=10)) +
  theme(legend.key.height=unit(1,"line")) +
  theme(legend.key.width=unit(2,"line")) +
  theme(legend.key = element_blank()) + theme(legend.position = c(0.7,0.8),legend.direction = "horizontal",legend.box = "vertical") +
  theme(text = element_text(size=12)) +
  theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
  theme(axis.title.y = element_text(size = 12, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plots[[2]] = ggplot(data = g1.data, aes(x = height, y = g1_date, group = volume, colour=factor(volume))) +
  geom_point(size=1) +
  stat_smooth(method = "lm") +
  ylab("g1") + xlab("Height (cm)") + 
  theme_bw() + labs(colour="Pot Volume") +
  theme(legend.title = element_text(colour="black", size=10)) +
  theme(legend.text = element_text(colour="black", size=10)) +
  theme(legend.key.height=unit(1,"line")) +
  theme(legend.key.width=unit(2,"line")) +
  theme(legend.key = element_blank()) + theme(legend.position = c(0.7,0.8),legend.direction = "horizontal",legend.box = "vertical") +
  theme(text = element_text(size=12)) +
  theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
  theme(axis.title.y = element_text(size = 12, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plots[[3]] = ggplot(data = g1.data, aes(x = Date, y = g1_date, group = volume, colour=factor(volume))) +
  geom_point(size=1) +
  stat_smooth(se=FALSE,method = "lm") +
  ylab("g1") + xlab("") + 
  theme_bw() + labs(colour="Pot Volume") +
  theme(legend.title = element_text(colour="black", size=10)) +
  theme(legend.text = element_text(colour="black", size=10)) +
  theme(legend.key.height=unit(1,"line")) +
  theme(legend.key.width=unit(2,"line")) +
  theme(legend.key = element_blank()) + 
  theme(text = element_text(size=12)) +
  theme(axis.title.x = element_text(size = 12, vjust=-.2)) +
  theme(axis.title.y = element_text(size = 12, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png("output/Figure_g1_with_LA_height_time.png", units="px", width=2000, height=1600, res=220)
print (do.call(grid.arrange,  plots))
dev.off()
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
# # Read g1 data directly from Court's Github folder
# Vcmax.Jmax <- read.csv("https://raw.githubusercontent.com/kashifmahmud/EucPVE/master/calculated%20data/jmax_vcmax.csv")
# Vcmax.Jmax$Date = as.Date(Vcmax.Jmax$Date)
# keeps = c("Date","volume","Vcmax.mean","Jmax.mean")
# Vcmax.Jmax.mean = Vcmax.Jmax[ , keeps, drop = FALSE]
# keeps = c("Date","volume","Vcmax.se","Jmax.se")
# Vcmax.Jmax.se = Vcmax.Jmax[ , keeps, drop = FALSE]
# 
# 
# Vcmax.Jmax.mean.melt <- melt(Vcmax.Jmax.mean, id.vars = c("Date","volume"))
# Vcmax.Jmax.se.melt <- melt(Vcmax.Jmax.se, id.vars = c("Date","volume"))
# Vcmax.Jmax.mean.melt$se <- Vcmax.Jmax.se.melt$value
# 
# plots = list()
# plots[[1]] = ggplot(Vcmax.Jmax.mean.melt, aes(x=Date, y=value, group = interaction(volume,variable), colour=as.factor(volume), shape=as.factor(variable))) + 
#   geom_point(position=pd) +
#   geom_errorbar(position=pd,aes(ymin=value-se, ymax=value+se), colour="grey", width=2) +
#   geom_line(position=pd,data = Vcmax.Jmax.mean.melt, aes(x = Date, y = value, group = interaction(volume,variable), linetype=as.factor(variable), colour=as.factor(volume))) 
#   # ylab(paste(as.character(meas[p]),"(g C)")) + xlab("Month") +
#   # # ggtitle("C pools - Measured (points) vs Modelled (lines)") +
#   # # labs(colour="Soil Volume", linetype="Grouping treatment", size="Total No of Parameter") +
#   # # labs(colour="Pot Volume (L)", linetype="No. of Parameters") +
#   # labs(colour="Pot Volume (L)") +
#   # # scale_color_manual(labels = c("Individuals", "One Group"), values = c("blue", "red")) +
#   # # coord_trans(y = "log10") + ylab(paste(as.character(meas[p]),"(g C plant-1)")) +
#   # theme_bw() +
#   # annotate("text", x = min(summary.output.Cpool$Date), y = max(summary.output.Cpool$value), size = font.size-7, label = paste(title[p])) +
#   # # theme(plot.title = element_text(size = 20, face = "bold")) +
#   # theme(legend.title = element_text(colour="black", size=font.size)) +
#   # theme(legend.text = element_text(colour="black", size = font.size)) +
#   # # theme(legend.key.height=unit(0.9,"line")) +
#   # theme(legend.position = c(0.17,0.7)) +
#   # theme(legend.key = element_blank()) +
#   # theme(text = element_text(size=font.size)) +
#   # theme(axis.title.x = element_blank()) +
#   # theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
#   # # theme(plot.title = element_text(hjust = 0)) +
#   # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# #-------------------------------------------------------------------------------------
# 
# 
# #-------------------------------------------------------------------------------------
# # Read g1 data directly from Court's Github folder
# Vcmax.Jmax <- read.csv("https://raw.githubusercontent.com/kashifmahmud/EucPVE/master/calculated%20data/jmax_vcmax_clean.csv")
# Vcmax.Jmax$SF = Vcmax.Jmax$Jmax.mean / max(Vcmax.Jmax$Jmax.mean)
# 
# Vcmax.Jmax$alpha = 0.24 * Vcmax.Jmax$SF



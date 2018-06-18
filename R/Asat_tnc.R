#- Plot photosynthesis and storage relationship from the gas exchange measurements

Asat = read.csv("raw_data/Asat_obs.csv")
keeps = c("Date","volume","Photo")
Asat.set = Asat[ , keeps, drop = FALSE]
Asat.set <- summaryBy(Photo ~ Date+volume, data=Asat.set, FUN=c(mean,standard.error))
names(Asat.set)[3:4] = c("Asat", "Asat_SE")
Asat.set$Date = as.Date(Asat.set$Date, format = "%Y-%m-%d")

#-----------------------------------------------------------------------------------------
leaf.data = read.csv("raw_data/leaf_data.csv")
leaf.data$Date = as.Date(leaf.data$Date, format = "%d/%m/%Y")

keeps = c("Date","volume","starch","sugars","tnc")
tnc.set = leaf.data[ , keeps, drop = FALSE]
tnc.set <- summaryBy(starch+sugars+tnc ~ Date+volume, data=tnc.set, FUN=c(mean,standard.error))
names(tnc.set)[3:8] = c("starch","sugars","tnc","starch_SE","sugars_SE","tnc_SE")
Asat.set = merge(Asat.set,tnc.set,by=c("Date","volume"))


#-----------------------------------------------------------------------------------------
plot = list() 
plot[[1]] = ggplot(data = Asat.set, aes(x = Asat, y = tnc, group = as.factor(volume), colour=as.factor(volume))) +
  geom_point(shape=17, size=1) +
  stat_smooth(method=lm, se = TRUE) +
  labs(colour="Pot Volume (L)") + xlab (expression(paste("Asat (",mu, mol, " ", m^-2, " ", s^-1,")", sep=""))) + ylab ("TNC (%)") + 
  # ggtitle("Asat vs TNC") +
  theme_bw() + theme(legend.position = c(0.85,0.86),legend.key.height=unit(0.5,"line"),legend.key.width=unit(0.5,"line"),legend.direction = "horizontal", legend.text.align = 0) +
  theme(legend.text = element_text(colour="black", size = 6)) +
  theme(legend.title = element_text(colour="black", size = 7)) +
  theme(axis.title.x = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot[[2]] = ggplot(data = Asat.set, aes(x = Asat, y = starch, group = as.factor(volume), colour=as.factor(volume))) +
  geom_point(shape=17, size=1) +
  stat_smooth(method=lm, se = TRUE) +
  labs(colour="Pot Volume (L)") + xlab (expression(paste("Asat (",mu, mol, " ", m^-2, " ", s^-1,")", sep=""))) + ylab ("Starch (%)") + 
  # ggtitle("Asat vs Starch") +
  theme_bw() + theme(legend.position = c(0.85,0.86),legend.key.height=unit(0.5,"line"),legend.key.width=unit(0.5,"line"),legend.direction = "horizontal", legend.text.align = 0) +
  theme(legend.text = element_text(colour="black", size = 6)) +
  theme(legend.title = element_text(colour="black", size = 7)) +
  theme(axis.title.x = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot[[3]] = ggplot(data = Asat.set, aes(x = Asat, y = sugars, group = as.factor(volume), colour=as.factor(volume))) +
  geom_point(shape=17, size=1) +
  stat_smooth(method=lm, se = TRUE) +
  labs(colour="Pot Volume (L)") + xlab (expression(paste("Asat (",mu, mol, " ", m^-2, " ", s^-1,")", sep=""))) + ylab ("Sugars (%)") + 
  # ggtitle("Asat vs Sugars") +
  theme_bw() + theme(legend.position = c(0.85,0.86),legend.key.height=unit(0.5,"line"),legend.key.width=unit(0.5,"line"),legend.direction = "horizontal", legend.text.align = 0) +
  theme(legend.text = element_text(colour="black", size = 6)) +
  theme(legend.title = element_text(colour="black", size = 7)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot[[2]] = plot[[2]] + guides(colour=FALSE)
plot[[3]] = plot[[3]] + guides(colour=FALSE)

png("output/Figure_Asat_tnc_starch_sugars_SE.png", units="px", width=1600, height=1300, res=220)
print (do.call(grid.arrange,  plot))
dev.off()

#-----------------------------------------------------------------------------------------
plot = list()
plot[[1]] = ggplot(Asat.set, aes(x = Date)) +
  geom_line(aes(y = Asat, colour = "Asat")) + geom_errorbar(aes(ymin=Asat-Asat_SE, ymax=Asat+Asat_SE), colour="grey", width=1) +
  geom_line(aes(y = tnc*100, colour = "TNC")) + geom_errorbar(aes(ymin=100*(tnc-tnc_SE), ymax=100*(tnc+tnc_SE)), colour="grey", width=1) +
  scale_y_continuous(sec.axis = sec_axis(~./100, name = "TNC, Starch, Sugars (%)")) +
  geom_line(aes(y = starch*100, colour = "starch")) + geom_errorbar(aes(ymin=100*(starch-starch_SE), ymax=100*(starch+starch_SE)), colour="grey", width=1) +
  geom_line(aes(y = sugars*100, colour = "sugars")) + geom_errorbar(aes(ymin=100*(sugars-sugars_SE), ymax=100*(sugars+sugars_SE)), colour="grey", width=1) +
  facet_wrap( ~ volume) +
  
  scale_colour_manual(values = c("red","blue","green","orange")) +
  labs(y = expression(paste("Asat (",mu, mol, " ", m^-2, " ", s^-1,")", sep="")), x = "Date", colour = "Parameter") +
  
  theme_bw() + theme(legend.position = c(0.7,0.15),legend.key.height=unit(0.8,"line"),legend.key.width=unit(0.8,"line"),legend.direction = "vertical", legend.text.align = 0) +
  theme(legend.text = element_text(colour="black", size = 9)) +
  theme(legend.title = element_text(colour="black", size = 10)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(size=7, angle=0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png("output/Figure_Asat_tnc_starch_sugars_individual.png", units="px", width=1600, height=1300, res=220)
print (do.call(grid.arrange,  plot))
dev.off()


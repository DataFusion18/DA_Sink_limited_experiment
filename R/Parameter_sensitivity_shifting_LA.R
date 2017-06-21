# This code performs sensitivity analysis with parameter shifting from 5L pot to free seedling
# and make Figure 6
#-------------------------------------------------------------------------------------
############## Combining all data in one dataframe
# Merge all Cday, Rd, Cleaf, Cstem, Croot data
data = merge(subset(Cday.data.raw,volume %in% c(5,1000)),subset(Rd.data.processed,volume %in% c(5,1000)), all = TRUE)
data = merge(data,subset(LA.data.processed,volume %in% c(5,1000)), all = TRUE)
data = merge(data,subset(tnc.data.processed,volume %in% c(5,1000)), all = TRUE)
data = merge(data,subset(Mleaf.data.processed,volume %in% c(5,1000)), all = TRUE)
data = merge(data,subset(Mstem.data.processed,volume %in% c(5,1000)), all = TRUE)
data = merge(data,subset(Mroot.data.processed,volume %in% c(5,1000)), all = TRUE)
# names(data)[4:ncol(data)] = c("Rd","Sleaf","Sleaf_SD","Mleaf","Mleaf_SD","Mstem","Mstem_SD","Mroot","Mroot_SD")
names(data)[3:ncol(data)] = c("Cday","Rd","LA","LA_SD","Sleaf","Sleaf_SD","Mleaf","Mleaf_SD","Mstem","Mstem_SD","Mroot","Mroot_SD")
# data[ , c(9:ncol(data))] = data[ , c(9:ncol(data))] * 0.65 # Unit conversion: gDM to gC


# summary.param.set = subset(summary.param, variable %in% var[p] & volume %in% c(5,1000))
# keeps = c("Date", "variable", "Parameter", "volume")
# param.set = summary.param.set[ , keeps, drop = FALSE]
# param.set.casted = dcast( param.set , Date ~ variable )

##################------------------------------
# Consider everything (Cday, LA, Rd, sigma, parameters) for potted seedling 5L (group 1)
q=0 # Case 0
Cday.data = subset(Cday.data.processed,volume==5) # Consider the free seedling to test the parameter sensitivity
Rd.data = subset(Rd.data.processed,volume==5)
Mleaf.data = subset(Mleaf.data.processed,volume==5)
Mstem.data = subset(Mstem.data.processed,volume==5)
Mroot.data = subset(Mroot.data.processed,volume==5)
Sleaf.data = tnc.data = subset(tnc.data.processed,volume==5)
LA.data = subset(LA.data.processed,volume==5)
sigma.data = subset(sigma.data.processed,volume==5)

sla.harvest.data = subset(sla.harvest.processed,volume %in% 5)
sigma.data$SLA = sla.harvest.data$sla_no_tnc

summary.param = result[[2]]
param = subset(summary.param,(volume.group %in% 1)) # volume.group = 1 is for potted seedling 5L
# param = subset(summary.param,(volume==5)) # volume.group = 1 is for potted seedling 5L
keeps = c("Date", "variable", "Parameter")
param = param[ , keeps, drop = FALSE]
param.casted = dcast( param , Date ~ variable )

# Set the colours for the graph (colourblind friendly palette)
# cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift.R")

##################------------------------------
# Take the Cday for free seedling
q=1 # Case 1
# Raw data processing for free seedling only (1000L)
Cday.data = subset(Cday.data.processed,volume==1000) # Consider the free seedling to test the parameter sensitivity

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift.R")

##################------------------------------
# Take the Cday, Rd for free seedling
q=2 # Case 2
Rd.data = subset(Rd.data.processed,volume==1000) # Consider the free seedling to test the parameter sensitivity

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift.R")

############----------------------------------------
# Take the parameters af, as, ar for free seedling
q=3 # Case 4
param.pot = subset(summary.param,(volume.group %in% 1 & variable %in% c("k","Y","sf")))
param.free = subset(summary.param,(volume.group %in% 3 & variable %in% c("af","as","ar")))
param = rbind(param.pot, param.free)
keeps = c("Date", "variable", "Parameter")
param = param[ , keeps, drop = FALSE]
param.casted = dcast( param , Date ~ variable )

sla.harvest.data = subset(sla.harvest.processed,volume %in% 1000)
sigma.data$SLA = sla.harvest.data$sla_no_tnc

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift.R")

############----------------------------------------
# Take the parameters Y, af, as, ar for free seedling
q=4 # Case 5
param.pot = subset(summary.param,(volume.group %in% 1 & variable %in% c("k","sf")))
param.free = subset(summary.param,(volume.group %in% 3 & variable %in% c("Y","af","as","ar")))
param = rbind(param.pot, param.free)
keeps = c("Date", "variable", "Parameter")
param = param[ , keeps, drop = FALSE]
param.casted = dcast( param , Date ~ variable )

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift.R")

############----------------------------------------
# Take the parameters Y, af, as, ar, sf for free seedling
q=5 # Case 6
param.pot = subset(summary.param,(volume.group %in% 1 & variable %in% c("k")))
param.free = subset(summary.param,(volume.group %in% 3 & variable %in% c("Y","af","as","ar","sf")))
param = rbind(param.pot, param.free)
keeps = c("Date", "variable", "Parameter")
param = param[ , keeps, drop = FALSE]
param.casted = dcast( param , Date ~ variable )

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift.R")

############----------------------------------------
# Take the parameters Y, k, af, as, ar, sf for free seedling
q=6 # Case 7
param = subset(summary.param,(volume.group %in% 3)) # volume.group = 3 is for free seedling from the 3 groups
keeps = c("Date", "variable", "Parameter")
param = param[ , keeps, drop = FALSE]
param.casted = dcast( param , Date ~ variable )

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift.R")

############ Summarize the C pools plots
######## Plot both set of Cday
plot.shift = list() 
font.size = 12
title = as.character(c("A","B","C","D","E","F","G","H","I"))

Cday.data.processed$Date = as.Date(Cday.data.processed$Date)
Cday.set = subset(Cday.data.processed, volume %in% c(5,1000))
plot.shift[[1]] = plot.Cday(Cday.set, 1)
# plot.shift[[1]] = ggplot(data = Cday.set, aes(x = Date, y = carbon_day,  group = volume, colour=factor(volume))) +
#   geom_point(size=0.01) +
#   geom_line(data = Cday.set, aes(x = Date, y = carbon_day,  group = volume, colour=factor(volume))) +
#   # xlab("Month") +
#   ylab(expression(C[day]~"(g C "*d^"-1"*" "*leafarea^"-1"*")")) + 
#   # ylab("Cday (g C d-1 leafarea-1)") +
#   # ggtitle("A - Case 1: An (5L pot -> Free)") +
#   # scale_colour_discrete(name="Treatments", breaks=c("5", "1000"), labels=c("5L Pot", "Free seedling")) +
#   scale_colour_manual(name="", breaks=c("5", "1000"), labels=c("5L", "FS"), values=cbPalette[1:2]) +
#   annotate("text", x = min(Cday.set$Date), y = max(Cday.set$carbon_day)*0.98, size = font.size-7, label = paste(title[1])) +
#   theme_bw() +
#   theme(legend.position = c(0.85,0.85)) +
#   # theme(plot.title = element_text(size = font.size)) +
#   theme(legend.title = element_blank()) +
#   theme(legend.key = element_blank(), plot.margin=unit(c(0.25, 0.25, 0, 0.45), units="line")) +
#   theme(text = element_text(size=font.size)) +
#   theme(legend.key.height=unit(0.65,"line")) +
#   # theme(legend.key.width=unit(2,"line")) +
#   theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y = element_text(size = font.size, vjust=0.3)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(plot.title = element_text(vjust=-2))

######## Plot both set of Rd
Rd.data.processed$Date = as.Date(Rd.data.processed$Date)
Rd.data = subset(Rd.data.processed, volume %in% c(5,1000))
plot.shift[[2]] = plot.Rd(Rd.data, 2)
# plot.shift[[2]] = ggplot(data = Rd.data, aes(x = Date, y = Rd_daily,  group = volume, colour=factor(volume))) +
#   geom_point(size=0.01) +
#   geom_line(data = Rd.data, aes(x = Date, y = Rd_daily,  group = volume, colour=factor(volume))) +
#   # xlab("Month") +
#   ylab(expression(R[d]~"(g C "*g^"-1"*" plant "*d^"-1"*")")) + 
#   # ggtitle("B - Case 2: Rd (5L pot -> Free)") +
#   scale_colour_manual(name="", breaks=c("5", "1000"), labels=c("5L", "FS"), values=cbPalette[2:3]) +
#   annotate("text", x = min(Rd.data$Date), y = max(Rd.data$Rd_daily)*0.98, size = font.size-7, label = paste(title[2])) +
#   theme_bw() +
#   theme(legend.position = c(0.85,0.85)) +
#   # theme(plot.title = element_text(size = font.size)) +
#   theme(legend.title = element_blank()) +
#   theme(legend.key = element_blank(), plot.margin=unit(c(0.25, 0.25, 0, 0.2), units="line")) +
#   theme(text = element_text(size=font.size)) +
#   theme(legend.key.height=unit(0.65,"line")) +
#   # theme(legend.key.width=unit(2,"line")) +
#   theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y = element_text(size = font.size, vjust=0.3)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Plot individual modelled parameters ("k","Y","af","sf") against "volume"
summary.param.set = subset(summary.param, variable %in% as.factor(c("af","as","ar")) & volume.group %in% c(1,3))
plot.shift[[3]] = plot.allocation.fractions(summary.param.set, 3)

summary.param.set = subset(summary.param, variable %in% as.factor("Y") & volume.group %in% c(1,3))
plot.shift[[4]] = plot.Y(summary.param.set, 4)
# plot.shift[[4]] = plot.shift[[4]] + ylab(expression(Y~"(g C "*g^"-1"*" C "*d^"-1"*")"))
# + theme(axis.text.x=element_blank())

summary.param.set = subset(summary.param, variable %in% as.factor("sf") & volume.group %in% c(1,3))
plot.shift[[5]] = plot.sf(summary.param.set, 5)
# plot.shift[[5]] = plot.shift[[5]] + ylab(expression(s[f]~"(g C "*g^"-1"*" C "*d^"-1"*")")) 
#   + theme(axis.text.x=element_blank())

summary.param.set = subset(summary.param, variable %in% as.factor("k") & volume.group %in% c(1,3))
plot.shift[[6]] = plot.k(summary.param.set, 6)
# plot.shift[[6]] = plot.shift[[6]] + ylab(expression(k~"(g C "*g^"-1"*" C "*d^"-1"*")"))
#   # + theme(axis.text.x = element_text(size = font.size, vjust=0))
#   + theme(plot.margin=unit(c(0.25, 0.5, 0.25, 0.5), units="line"))

# # This sript plots the C pools for various test cases with parameter shifted from potted seedling to free seedling
# source("R/generate_figures_param_shift.R")
shift.output.Mleaf = subset(shift.output,(variable %in% "Mleaf"))
plot.shift[[7]] = plot.Mleaf(shift.output.Mleaf)

shift.output.Mstem = subset(shift.output,(variable %in% "Mstem"))
plot.shift[[8]] = plot.Mstem(shift.output.Mstem)

shift.output.Mroot = subset(shift.output,(variable %in% "Mroot"))
plot.shift[[9]] = plot.Mroot(shift.output.Mroot)

png("output/parameter_shifting.png", units="px", width=2000, height=3000, res=220)
lay <- rbind(c(1,7,7),c(2,7,7),c(3,8,8),c(4,8,8),c(5,9,9),c(6,9,9))
grid.arrange(grobs = plot.shift, layout_matrix = lay)
dev.off()


# Make figure for AGU presentation

shift.output.biomass <- dplyr::summarize(group_by(subset(shift.output, variable %in% c("Mleaf","Mstem","Mroot")),Date,Case),
                                         value = sum(value,na.rm=T))

shift.output.biomass.set = subset(shift.output.biomass, Case %in% c("0"))
plot.shift[[7]] = ggplot() +
  geom_line(data = shift.output.biomass.set, aes(x = Date, y = value, group = Case, colour=Case), size=1) +
  geom_point(size=2) +
  ylab(expression("Total Biomass"~"(g C "*plant^"-1"*")")) + 
  scale_colour_manual(breaks=c("0"), labels=c("Baseline (5L)"),values=cbPalette) +
  scale_y_continuous(limits = c(min(shift.output.biomass$value), max(shift.output.biomass$value)), breaks=c(10,20,30,40,50),labels=c(10,20,30,40,50)) +
  theme_bw() +
  theme(legend.position = c(0.2,0.7),legend.text.align = 0) +
  theme(legend.title = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size+2)) +
  theme(legend.key.height=unit(1,"line")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

shift.output.biomass.set = subset(shift.output.biomass, Case %in% c("0","1"))
plot.shift[[8]] = ggplot() +
  geom_line(data = shift.output.biomass.set, aes(x = Date, y = value, group = Case, colour=Case), size=1) +
  geom_point(size=2) +
  ylab(expression("Total Biomass"~"(g C "*plant^"-1"*")")) + 
  scale_colour_manual(breaks=c("0","1"), labels=c("Baseline (5L)","+ Cday of FS"),values=cbPalette) +
  scale_y_continuous(limits = c(min(shift.output.biomass$value), max(shift.output.biomass$value)), breaks=c(10,20,30,40,50),labels=c(10,20,30,40,50)) +
  theme_bw() +
  theme(legend.position = c(0.2,0.7),legend.text.align = 0) +
  theme(legend.title = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size+2)) +
  theme(legend.key.height=unit(1,"line")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

shift.output.biomass.set = subset(shift.output.biomass, Case %in% c("0","1","2"))
plot.shift[[9]] = ggplot() +
  geom_line(data = shift.output.biomass.set, aes(x = Date, y = value, group = Case, colour=Case), size=1) +
  geom_point(size=2) +
  ylab(expression("Total Biomass"~"(g C "*plant^"-1"*")")) + 
  scale_colour_manual(breaks=c("0","1","2"), labels=c("Baseline (5L)","+ Cday of FS",expression(+ R[d]~"of FS")), 
                      values=cbPalette) +
  scale_y_continuous(limits = c(min(shift.output.biomass$value), max(shift.output.biomass$value)), breaks=c(10,20,30,40,50),labels=c(10,20,30,40,50)) +
  theme_bw() +
  theme(legend.position = c(0.2,0.7),legend.text.align = 0) +
  theme(legend.title = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size+2)) +
  theme(legend.key.height=unit(1,"line")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

shift.output.biomass.set = subset(shift.output.biomass, Case %in% c("0","1","2","3"))
plot.shift[[10]] = ggplot() +
  geom_line(data = shift.output.biomass.set, aes(x = Date, y = value, group = Case, colour=Case), size=1) +
  geom_point(size=2) +
  ylab(expression("Total Biomass"~"(g C "*plant^"-1"*")")) + 
  scale_colour_manual(breaks=c("0","1","2","3"), labels=c("Baseline (5L)","+ Cday of FS",expression(+ R[d]~"of FS"),
                       expression(+ (a[f] + a[w] + a[r])~"of FS")),values=cbPalette) +
  scale_y_continuous(limits = c(min(shift.output.biomass$value), max(shift.output.biomass$value)), breaks=c(10,20,30,40,50),labels=c(10,20,30,40,50)) +
  theme_bw() +
  theme(legend.position = c(0.2,0.7),legend.text.align = 0) +
  theme(legend.title = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size+2)) +
  theme(legend.key.height=unit(1,"line")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

shift.output.biomass.set = subset(shift.output.biomass, Case %in% c("0","1","2","3","4"))
plot.shift[[11]] = ggplot() +
  geom_line(data = shift.output.biomass.set, aes(x = Date, y = value, group = Case, colour=Case), size=1) +
  geom_point(size=2) +
  ylab(expression("Total Biomass"~"(g C "*plant^"-1"*")")) + 
  scale_colour_manual(breaks=c("0","1","2","3","4"), labels=c("Baseline (5L)","+ Cday of FS",expression(+ R[d]~"of FS"),
                     expression(+ (a[f] + a[w] + a[r])~"of FS"),"+ Y of FS"),values=cbPalette) +
  scale_y_continuous(limits = c(min(shift.output.biomass$value), max(shift.output.biomass$value)), breaks=c(10,20,30,40,50),labels=c(10,20,30,40,50)) +
  theme_bw() +
  theme(legend.position = c(0.2,0.7),legend.text.align = 0) +
  theme(legend.title = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size+2)) +
  theme(legend.key.height=unit(1,"line")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

shift.output.biomass.set = subset(shift.output.biomass, Case %in% c("0","1","2","3","4","5"))
plot.shift[[12]] = ggplot() +
  geom_line(data = shift.output.biomass.set, aes(x = Date, y = value, group = Case, colour=Case), size=1) +
  geom_point(size=2) +
  ylab(expression("Total Biomass"~"(g C "*plant^"-1"*")")) + 
  scale_colour_manual(breaks=c("0","1","2","3","4","5"), labels=c("Baseline (5L)","+ Cday of FS",expression(+ R[d]~"of FS"),
                                                                      expression(+ (a[f] + a[w] + a[r])~"of FS"),"+ Y of FS",expression(+ s[f]~"of FS")), 
                      values=cbPalette) +
  scale_y_continuous(limits = c(min(shift.output.biomass$value), max(shift.output.biomass$value)), breaks=c(10,20,30,40,50),labels=c(10,20,30,40,50)) +
  theme_bw() +
  theme(legend.position = c(0.2,0.7),legend.text.align = 0) +
  theme(legend.title = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size+2)) +
  theme(legend.key.height=unit(1,"line")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot.shift[[13]] = ggplot() +
  geom_line(data = shift.output.biomass, aes(x = Date, y = value, group = Case, colour=Case), size=1) +
  geom_point(size=2) +
  ylab(expression("Total Biomass"~"(g C "*plant^"-1"*")")) + 
  scale_colour_manual(breaks=c("0","1","2","3","4","5","6"), labels=c("Baseline (5L)","+ Cday of FS",expression(+ R[d]~"of FS"),
                                                                      expression(+ (a[f] + a[w] + a[r])~"of FS"),"+ Y of FS",expression(+ s[f]~"of FS"),"+ k of FS (Complete FS)"), 
                      values=cbPalette) +
  scale_y_continuous(limits = c(min(shift.output.biomass$value), max(shift.output.biomass$value)), breaks=c(10,20,30,40,50),labels=c(10,20,30,40,50)) +
  theme_bw() +
  theme(legend.position = c(0.2,0.7),legend.text.align = 0) +
  theme(legend.title = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size+2)) +
  theme(legend.key.height=unit(1,"line")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot.shift[[1]] = plot.shift[[1]] + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = font.size, vjust=0.3))
plot.shift[[2]] = plot.shift[[2]] + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = font.size, vjust=0.3))
plot.shift[[3]] = plot.shift[[3]] + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = font.size, vjust=0.3))
plot.shift[[4]] = plot.shift[[4]] + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = font.size, vjust=0.3))
plot.shift[[5]] = plot.shift[[5]] + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = font.size, vjust=0.3))
plot.shift[[6]] = plot.shift[[6]] + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = font.size, vjust=0.3))

select_grobs <- function(lay) {
  id <- unique(c(t(lay))) 
  id[!is.na(id)]
}

png("output/AGU_Figure/AGU_Figure_1_parameter_shifting.png", units="px", width=2000, height=1000, res=220)
hlay <- rbind(c(NA,7,7),c(NA,7,7))
grid.arrange(grobs = plot.shift[select_grobs(hlay)], layout_matrix = hlay)
dev.off()

png("output/AGU_Figure/AGU_Figure_2_parameter_shifting.png", units="px", width=2000, height=1000, res=220)
hlay <- rbind(c(1,8,8),c(NA,8,8))
grid.arrange(grobs = plot.shift[select_grobs(hlay)], layout_matrix = hlay)
dev.off()

png("output/AGU_Figure/AGU_Figure_3_parameter_shifting.png", units="px", width=2000, height=1000, res=220)
hlay <- rbind(c(2,9,9),c(NA,9,9))
grid.arrange(grobs = plot.shift[select_grobs(hlay)], layout_matrix = hlay)
dev.off()

png("output/AGU_Figure/AGU_Figure_4_parameter_shifting.png", units="px", width=2000, height=1000, res=220)
hlay <- rbind(c(3,10,10),c(NA,10,10))
grid.arrange(grobs = plot.shift[select_grobs(hlay)], layout_matrix = hlay)
dev.off()

png("output/AGU_Figure/AGU_Figure_5_parameter_shifting.png", units="px", width=2000, height=1000, res=220)
hlay <- rbind(c(4,11,11),c(NA,11,11))
grid.arrange(grobs = plot.shift[select_grobs(hlay)], layout_matrix = hlay)
dev.off()

png("output/AGU_Figure/AGU_Figure_6_parameter_shifting.png", units="px", width=2000, height=1000, res=220)
hlay <- rbind(c(5,12,12),c(NA,12,12))
grid.arrange(grobs = plot.shift[select_grobs(hlay)], layout_matrix = hlay)
dev.off()

png("output/AGU_Figure/AGU_Figure_7_parameter_shifting.png", units="px", width=2000, height=1000, res=220)
hlay <- rbind(c(6,13,13),c(NA,13,13))
grid.arrange(grobs = plot.shift[select_grobs(hlay)], layout_matrix = hlay)
dev.off()

shift.output.biomass.set = subset(shift.output.biomass, Date %in% as.Date("2013-05-21"))
shift.output.biomass.set$attribution = (c("Baseline (5L)", "+ Cday of FS", "+ Rd of FS",
                                                      "+ (af + aw + ar) of FS", "+ Y of FS", "+ sf of FS", "+ k of FS (Complete FS)"))
write.csv(shift.output.biomass.set, file = "output/AGU_Figure/final_mass_changes.csv", row.names = FALSE)


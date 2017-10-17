

library(devtools)
install_bitbucket("yplantqmc","remkoduursma")
# do whatever the message says:
library(YplantQMC)


# make a random plant:
plant <- constructplant(randomplant(LA=1,   # leaf area (m2)
                                    leaflen=70, # leaf length mm
                                    height=1000,  # plant height mm
                                    radius=500,   # crown radius mm
                                    crownbase=100, # height to crown base mm
                                    crownshape="CONE"))  # crown shape, see ?randomplant


eucphy <- setPhy("Farquhar",leafpars=list(Vcmax=50, Jmax=100, G1=8, G0=0.01, Rd0=1))
sydney <- setLocation(lat=-33.5, long=152, tzlong=150)

sunnyday <- setMet(sydney, month=3, day=15, nsteps=12, Tmin=15, Tmax=29, PARday=22)

run1 <- YplantDay(plant, phy=eucphy, met=sunnyday)



# Test direct vs. diffuse
day_fbeam0 <- setMet(sydney, month=3, day=15, nsteps=12, Tmin=15, Tmax=29, PARday=22,
                     fbeam=0, fbeammethod="constant")
day_fbeam1 <- setMet(sydney, month=3, day=15, nsteps=12, Tmin=15, Tmax=29, PARday=22,
                     fbeam=1, fbeammethod="constant")
run2 <- YplantDay(plant, phy=eucphy, met=day_fbeam0)
run3 <- YplantDay(plant, phy=eucphy, met=day_fbeam1)

with(psrdata(run2), plot(timeofday, A/A0, type='l'))
with(psrdata(run3), points(timeofday, A/A0, type='l', col="red"))


# 
#euckey <- read.csv(....)
pfiles <- as.character(euckey$pfile)
lfiles <- as.character(euckey$lfile)


p <- constructplant(pfiles[2], lfiles[2])

plants <- readplantlist(pfiles,lfiles)

# send this to YplantDay :
# untested
runall <- YplantDay(plants, phy=eucphy, met=sunnyday)

# add summary variables:
plants_summary <- summary(plants)  # slow









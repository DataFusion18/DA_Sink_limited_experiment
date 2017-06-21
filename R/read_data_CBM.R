#- This sript reads the processed data from Sink limited container volume experiment 
#-----------------------------------------------------------------------------------------
# Import daily GPP, daily Rd data
GPP.data.processed = read.csv("processed_data/GPP.csv") # Units gC d-1
Cday.data.processed <- read.csv("raw_data/cday_120_clean_gross.csv") # Unit gC d-1
Rd.data.processed = read.csv("processed_data/Rd.csv") # Unit gC per gC plant per day
tnc.data.processed = read.csv("processed_data/tnc_fortnightly_data.csv") # Unit gC per gC plant

# Import weekly Cleaf, weekly Cstem, initial/harvest Croot data with Mean and SD
Mleaf.data.processed = read.csv("processed_data/Cleaf_weekly_data.csv") # Unit gC
Mstem.data.processed = read.csv("processed_data/Cstem_weekly_data.csv") # Unit gC
Mroot.data.processed = read.csv("processed_data/Croot_twice_data.csv") # Unit gC
LA.data.processed = read.csv("processed_data/LA_daily_data.csv") # Unit m^2

# Import harvest data
sla.harvest.processed = read.csv("processed_data/sla.harvest.csv") # Unit of SLA = m2 leaf area per gC of leaf biomass

# Import the self shading factors
sigma.data.processed <- read.csv("raw_data/M_leafarea_model.csv") # Unitless


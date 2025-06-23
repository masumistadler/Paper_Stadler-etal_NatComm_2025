## -------------------------------------------------------------------------
##
## Script name: 1_velocity.wrt.R
##
## Purpose of script: We will establish models to predict velocity and discharge
##                    across the entire watershed. Extract WRT for each pixel
##                    in the watershed, which reflects different ecosystems.
##                    E.g. reseroirs versus fluvial.
##
## Author: Masumi Stadler
##
## Date Finalized: 11/03/2025
##
## Copyright (c) Masumi Stadler, 2025
## Email: m.stadler.jp.at@gmail.com
##
## -------------------------------------------------------------------------
##
## Notes: None of the data are made available due to confidentiality
##        request from Hydro-Quebec. Hence, the script is provided for the sake
##        of transparency. Some data may be made available upon request.
##
##
## -------------------------------------------------------------------------

## Use R project with regular scripts, all paths are relative 

# Server set-up -----------------------------------------------------------
## Working directory is set from where the job is submitted
## Load library path, if on a server
# .libPaths( c( .libPaths(), "/home/mstadler/projects/def-pauldel/R/x86_64-pc-linux-gnu-library/4.2") )

# R-setup -----------------------------------------------------------------
## Load Packages -----------------------------------------------------------
pckgs <- list("readxl", # read excel
              "data.table", "stringr", "plyr","tidyverse", # wrangling
              "plotly", "ggpubr", 
              "mgcv", # multivariate
              "raster"
             ) # change as needed

## Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)

## Load
invisible(lapply(pckgs, library, character.only = T))
rm(pckgs, miss.pckgs)

## Load custom functions --------------------------------------------------
funs <- list.files("./Functions", full.names = T)
invisible(lapply(funs, source))

## Other set-up -----------------------------------------------------------
options(scipen = 6, digits = 4) # view outputs in non-scientific notation
 
## Parallel environment ---------------------------------------------------
## Server version
# cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
# registerDoMC(cores = cores)

## Personal version
# detectCores(); registerDoMC(); getDoParWorkers()
# numCores <- detectCores()
# cl <- makeCluster(numCores, type = "FORK")

# Steps ---------------------------------------------------------------------------
# In this script we will gather the data needed to estimate travel time / water age
# We gather:
#
# Reservoir ------------------------------------------------------------------------
# - Discharge (Q) into reservoir
# - Q out from reservoir
# - Reservoir volume
# - Calculate overall WRT of reservoirs
# - Get reservoir pixel number
# - Calculate WRT in each pixel of main channel in reservoir
#
# Streams/Rivers --------------------------------------------------------------------
# - Find models to accurately predict discharge
# - Estimate discharge for each pixel
# - Estimate velocity for each pixel
#
 
# Read in data -----------------------------------------------------------

# Split samples by category for snapping to different stream networks
# Reservoirs and River samples should snap to the main channel of La Romaine
# Rest, lakes, tributaries etc should snap to finer stream network = flow acc. thres.: 20000
meta <- read.csv("./Data/LaRomaine_MasterFile_withDist_2021-02-22.csv", sep = ",", stringsAsFactors = F)
setDT(meta)

main <- rbind(meta[str_detect(meta$sample.name, "LR"),],
      meta[str_detect(meta$sample.name, "RO"),],
      meta[str_detect(meta$sample.name, "RH"),])

rest <- meta[!(sample.name %in% main$sample.name),]

# write.table(main, "./Data/Summary/mother_mainsamples.csv", sep = ",", row.names = F)
# write.table(rest, "./Data/Summary/mother_restsamples.csv", sep = ",", row.names = F)

# after samples were separated by sample type
# main samples were snapped to main channel of La Romaine
# rest samples were snapped to finer stream network of water in ArcGIS
# new coordinates of snapped locations were extracted

# Reservoir residence time ----------------------------------------------------------------------------------------
# Now we need to estimate the Water Residence Time within the reservoirs
# Warning: Data not made available

# Volume from water level (HQ)
# We have a water level to volume and surface area conversion table from HQ
wat.key <- rbind(read_excel("./Data/Hydrology/Water_level_HQ_21July2020.xlsx", sheet = 3, skip = 5) %>% mutate(Reservoir = "RO1"),
                 read_excel("./Data/Hydrology/Water_level_HQ_21July2020.xlsx", sheet = 4, skip = 5) %>% mutate(Reservoir = "RO2"),
                 read_excel("./Data/Hydrology/Water_level_HQ_21July2020.xlsx", sheet = 5, skip = 5) %>% mutate(Reservoir = "RO3"))

# rename column names
colnames(wat.key) <- c("niveau_m","volume_hm3", "area_km2","comments","Reservoir")

# In order to convert the high resolution (daily) water level data we have from each reservoir to volume numbers,
# we have to generate models

# plot water level ~ volume relationship (needed for residence time)
ggplot(wat.key, aes(x = niveau_m, y = volume_hm3)) +
  geom_point() +
  facet_wrap(.~Reservoir, scale = "free") +
  geom_smooth(method = "gam")

# plot water level ~ surface area relationship
ggplot(wat.key, aes(x = niveau_m, y = area_km2)) +
  geom_point() +
  facet_wrap(.~Reservoir, scale = "free") +
  geom_smooth(method = "gam")

# Relationships are non linear for all reservoirs for both volume and surface area

# make regression for water level ~ volume (GAM) for each reservoir
wat.gams <- dlply(wat.key, .(Reservoir), function(x){
  gam.mod <- gam(formula = volume_hm3 ~ s(niveau_m, bs = "cs"), data = x)
})

# make regression for water level ~ surface area (GAM) for each reservoir
area.gams <- dlply(wat.key, .(Reservoir), function(x){
  gam.mod <- gam(formula = area_km2 ~ s(niveau_m, bs = "cs"), data = x)
})

plot(wat.gams[[1]]);plot(wat.gams[[2]]);plot(wat.gams[[3]])
plot(area.gams[[1]]);plot(area.gams[[2]]);plot(area.gams[[3]])

# read in daily measurements of water level across all reservoirs
wat.lev <- read_excel("./Data/Hydrology/Water_level_HQ_21July2020.xlsx", sheet = 7, skip = 5,
                      col_types = c("date","numeric","numeric","numeric","numeric"))
colnames(wat.lev) <- c("Date","RO1","RO2","RO3","RO4")
setDT(wat.lev)

# melt into long format
wat.lev <- melt.data.table(wat.lev, id.vars = "Date", value.name = "niveau_m", variable.name = "Reservoir")
wat.lev <- wat.lev[!is.na(niveau_m),]

# predict volume
wat.lev[Reservoir == "RO1", volume_hm3 := predict.gam(wat.gams[[1]], newdata = wat.lev[Reservoir == "RO1",])]
wat.lev[Reservoir == "RO2", volume_hm3 := predict.gam(wat.gams[[2]], newdata = wat.lev[Reservoir == "RO2",])]
wat.lev[Reservoir == "RO3", volume_hm3 := predict.gam(wat.gams[[3]], newdata = wat.lev[Reservoir == "RO3",])]

# predict area
wat.lev[Reservoir == "RO1", area_km2 := predict.gam(area.gams[[1]], newdata = wat.lev[Reservoir == "RO1",])]
wat.lev[Reservoir == "RO2", area_km2 := predict.gam(area.gams[[2]], newdata = wat.lev[Reservoir == "RO2",])]
wat.lev[Reservoir == "RO3", area_km2 := predict.gam(area.gams[[3]], newdata = wat.lev[Reservoir == "RO3",])]

# Plot predicted volume and surface area from water level
ggplot(wat.lev[Reservoir != "RO4",], aes(x = niveau_m, y = volume_hm3)) +
  theme_bw() +
  geom_point(colour = 'grey20') +
  facet_wrap(.~Reservoir, scale = "free") +
  labs(x = "Water level (m)", y = expression(paste("Volume (hm"^3,")")))

ggplot(wat.lev, aes(x = niveau_m, y = area_km2)) +
  geom_point() +
  facet_wrap(.~Reservoir, scale = "free")

# some negative values in RO1 (close to minimum 52 m), overwrite with 0
wat.lev[Reservoir == "RO1" & volume_hm3 < 0, volume_hm3 := 0]

# calculate mean volume by month, for years when the reservoir is actually flooded
# 2015: RO2
# 2016: RO2 > RO1
# 2017 & 2018: RO3 > RO2 > RO1
wat.lev[, year := year(Date)]
wat.lev[, month := month(Date)]

wat.lev <- rbind(wat.lev[Reservoir == "RO2" & year >=2015,],
      wat.lev[Reservoir == "RO1" & year >=2016,],
      wat.lev[Reservoir == "RO3" & year >=2017,])

# Remove first month of when RO3 was being flooded (too low values)
# other reservoirs' flooding period where not recorded
wat.lev <- wat.lev[!(Reservoir == "RO3" & year == 2017 & month ==5),]

wat.lev[, month.f := factor(month, levels = 1:12)]
(p <- ggplot(wat.lev, aes(x = Date, y = niveau_m, colour = month.f)) +
  theme_bw() +
  geom_point() +
  facet_wrap(.~Reservoir, scale = "free_y") +
  scale_colour_viridis_d(name = "Month") +
  labs(y = "Water level (m)"))

#ggsave("./Figures/Methods/wat.lev_RO_timeseries.png", p, height = 9, width =25, units = "cm")

# calculate means for each year-month
sum.wat <- wat.lev[, .(mean.vol_hm3 = mean(volume_hm3, na.rm = T),
                       mean.area_km2 = mean(area_km2, na.rm = T)), by = .(Reservoir, year, month)]

sum.wat[, .(mean.vol_m3 = mean(mean.vol_hm3) * 1000000), by = .(Reservoir, month)]

# 1. Read bathymetry data -----------------------------------------------------------------------------------------
# Files generated in ArcMap using the bathymetry raster files provided by HQ
# Volume was derived via "Surface volume" tool in 3D Analyst toolbox
# Other statistics such as number of pixels in reservoir were calculated using "Zonal Statistics as Table" using the
# Reservoir shape files provided by HQ

b.stats <- list.files("./Data/Hydrology/", pattern = "bathy_stats")
vol.files <- list.files("./Data/Hydrology/", pattern = "bathy_vol")

# ignore warnings
b.df<-lapply(paste0("./Data/Hydrology/",b.stats), read.table, header = T, sep = ",", stringsAsFactors = F) # gives warnings but that's ok
b.df <- bind_rows(b.df)
b.df$Reservoir <- sapply(strsplit(b.stats, split = "\\_|\\."),"[",3)

vol.df<-lapply(paste0("./Data/Hydrology/",vol.files), read.table, header = T, sep = ",", stringsAsFactors = F)
vol.df <- bind_rows(vol.df)
vol.df$Reservoir <- sapply(strsplit(vol.files, split = "\\_|\\."),"[",3)

# Rename columns
vol.df <- vol.df %>% dplyr::select(Reservoir,
                         water.lev_m = Plane_Height,
                         volume_m3 = Volume)

b.df <- b.df %>% dplyr::select(Reservoir,
                               pxl_no = COUNT,
                               area_m2 = AREA,
                               max.depth_m = MAX,
                               min.depth_m = MIN,
                               mean.depth_m = MEAN,
                               std.depth_m = STD)

setDT(vol.df); setDT(b.df)
# merge
res.stats <- b.df[vol.df, , on = .(Reservoir)]

# convert to hm3 to compare with predicted data
res.stats[, volume_hm3 := volume_m3 * 0.000001]
# convert to km2 to compare with predicted data
res.stats[, area_km2 := area_m2 * 0.000001]
# remove unnecessary objects
rm(vol.df, b.df, vol.files, b.stats)

sum.wat[res.stats, c("gis.volume_hm3","gis.area_km2","gis.pxl_no") := list(i.volume_hm3,
                                                              i.area_km2, i.pxl_no), on = .(Reservoir)]
rm(area.gams, wat.gams, wat.key, wat.lev, p, res.stats, main, meta, rest)
#sum.wat[, pxl_no := mean.area_km2 / (0.018*0.018)] # Pixel size is 18 m x 18 m

# 2. Read in discharge data -------------------------------------------------------------------------
# Code used to summarise discharge data by HQ:
#
# Hydro-Quebec discharge data 
# load following libraries when dealing with summarising the raw data
# library(plyr)
# library(doMC)
# register number of cores
# detectCores()
# registerDoMC()
# getDoParWorkers() # 12

# The following code is to make the summary files, we will use the summary files for downstream analyses

#-- Turbines --#
#col.n <- names(read_excel("./Data/Débits_Turb_Diverse_Historique_Romaine 2014-fev 2020.xlsx", n_max = 0))
#hydro <- read_excel("./Data/Débits_Turb_Diverse_Historique_Romaine 2014-fev 2020.xlsx",
#                    skip = 2, col_names = F, na = c("", "NA"),
#                    col_types = c("date","numeric", "numeric", "numeric",
#                                  "numeric", "numeric", "numeric"))
#colnames(hydro) <- c("Date", "RO1.into.turb", "RO1.overflow",
#                     "RO2.into.turb", "RO2.overflow",
#                     "RO3.into.turb", "RO3.overflow")
# get years and months from date column
#hydro$year <- year(hydro$Date)
#hydro$month <- month(hydro$Date)

# define function for coefficient of variation
#cv <- function(x, na.rm = F){
#  if(na.rm == T){
#    sd(x, na.rm = T) / mean(x, na.rm = T)
#  } else {
#    sd(x) / mean(x)
#  }
#}

# summarise by year and month for downstream analyses
#annual <- hydro %>% dplyr::group_by(year) %>%
#  dplyr::summarise_at(.vars = c("RO1.into.turb", "RO1.overflow",
#                         "RO2.into.turb", "RO2.overflow",
#                         "RO3.into.turb", "RO3.overflow"),
#               .funs = list(~mean(., na.rm = T),
#                            ~sd(., na.rm = T),
#                            ~cv(., na.rm = T),
#                            ~sum(., na.rm = T)))

#monthly <- hydro %>% dplyr::group_by(year, month) %>%
#  dplyr::summarise_at(.vars = c("RO1.into.turb", "RO1.overflow",
#                         "RO2.into.turb", "RO2.overflow",
#                         "RO3.into.turb", "RO3.overflow"),
#               .funs = list(~mean(., na.rm = T),
#                            ~sd(., na.rm = T),
#                            ~cv(., na.rm = T),
#                            ~sum(., na.rm = T)))
# Export
#write.table(annual, file = "./Data/Summary/RO_discharge_annual.csv",
#            sep = ";", dec = ".", row.names = F)
#write.table(monthly, file = "./Data/Summary/RO_discharge_monthly.csv",
#            sep = ";", dec = ".", row.names = F)

#-- River --
#col.n <- excel_sheets("./Data/Raw/Debits_Stations_Romaine_FH_20200525_PourAlainTremblay.xlsx")
# we want: 
# ROMA0832 (inflow of RO3)
# ROMA0977 (outflow of RO3 / inflow of RO2)
# ROMA0945 (outflow of RO1)
#stations <- list("ROMA0832", "ROMA0977", "ROMA0945")
#hydro <- llply(stations, function(x){
#  df <- read_excel("./Data/Raw/Débits_Stations_Romaine_FH_20200525_PourAlainTremblay.xlsx",
#             sheet = x,
#             skip = 7, col_names = T, na = c("", "NA"),
#             col_types = c("date","numeric"))
#  colnames(df) <- c("Date", "Discharge")
#  return(df)
#}, .parallel = T)
#names(hydro) <- unlist(stations)
# Discharge in m3/s

#annual <- ldply(hydro, function(x){
# get years and months from date column
#  x$year <- year(x$Date)
#  x$month <- month(x$Date)

#  out <- x %>% dplyr::group_by(year) %>%
#      dplyr::summarise_at(.vars = "Discharge",
#                   .funs = list(~mean(., na.rm = T),
#                                ~sd(., na.rm = T),
#                                ~cv(., na.rm = T),
#                                ~sum(., na.rm = T)))
#  return(out)
#}, .parallel = T)

#monthly <- ldply(hydro, function(x){
# get years and months from date column
#  x$year <- year(x$Date)
#  x$month <- month(x$Date)

#  out <- x %>% dplyr::group_by(year, month) %>%
#    dplyr::summarise_at(.vars = "Discharge",
#                        .funs = list(~mean(., na.rm = T),
#                                     ~sd(., na.rm = T),
#                                     ~cv(., na.rm = T),
#                                     ~sum(., na.rm = T)))
#  return(out)
#}, .parallel = T)

# Export
#write.table(annual, file = "./Data/Summary/river_discharge_annual.csv",
#            sep = ";", dec = ".", row.names = F)
#write.table(monthly, file = "./Data/Summary/river_discharge_monthly.csv",
#            sep = ";", dec = ".", row.names = F)

#########################################################################

# Inflows to reservoirs are either:
# 1. Sampled river section before the reservoir OR
# 2. Outflow of the prior reservoir within continuum
# The sequence is:
# 2015: RO2
# 2016: RO2 > RO1
# 2017 & 2018: RO3 > RO2 > RO1
# Thus, the outflow of RO3 is the inflow of RO2

# read in summary files of HQ reservoir discharge data
mon.dis <- read.csv("./Data/Summary/RO_discharge_monthly.csv",
                    sep = ";", dec = ".", stringsAsFactors = F)
# read in summary files of HQ river discharge data
mon.riv <- read.csv("./Data/Summary/river_discharge_monthly.csv",
                    sep = ";", dec = ".", stringsAsFactors = F)

# Hydro-Quebec Discharge is given in m3 per second
setDT(mon.dis)
setDT(mon.riv)

# melt HQ data into long format
molten.dis <- melt.data.table(
  mon.dis,
  id.vars = c("year","month"),
  measure.vars = patterns(".into.turb_mean", ".overflow_mean"),
  variable.name = "Reservoir",
  value.name = c("dis.into.turb", "dis.overflow")
)
# rename reservoir factors
molten.dis[, Reservoir := paste0("RO", Reservoir)]
# add turbine and overflow discharge
molten.dis[, dis.sum := dis.into.turb + dis.overflow]

# calculate mean across years
#molten.dis <- molten.dis[dis.sum > 0, .(mean.dis_m3s = mean(dis.sum, na.rm = T)), by = .(Reservoir, month)]

# Get river data, we need:
# RO3 inflow (river) 2015-2018 (ROMA0832)
# RO2 inflow (river) 2015-2016 (ROMA0977)
# calculate means across years (2014-20)
#mon.riv <- mon.riv[mean > 0 & year >= 2014, .(mean.dis_m3s = mean(mean, na.rm = T)), by = .(.id, month)]
# Assign reservoir correspondence
mon.riv[.id == "ROMA0832", Reservoir := "RO3"]
mon.riv[.id == "ROMA0977", Reservoir := "RO2"]

# create data.frame to fill
# we will do this for all three reservoirs
res.sum <- data.frame(Reservoir = c(rep("RO1", times = 12*4),
                         rep("RO2", times = 12*4),
                         rep("RO3", times = 12*4)),
                      year = rep(2015:2018, each = 12, times = 3),
           month = rep(1:12, times = 3*4))

setDT(res.sum); setDT(molten.dis); setDT(mon.riv)

#---- Populate data.frame ----#
#
#- RO2 inflow -#
## 2015, 2016: river
## 2017 onwards: RO3 outflow
#- RO2 outflow -#
## 2015 onwards: RO2 outflow
#
#- RO1 inflow -#
## 2016 onwards: RO2 outflow
#- RO1 outflow -#
## 2016 onwards: RO1 outflow
#
#- RO3 inflow -#
## 2017-2020 river
#- RO3 outflow -#
# 2017 onwards: RO3 outflow

# Initiate discharge columns out and into reservoir
# Outflows from reservoirs are constant over the years, simple join
# Outflow from RO1 (outflow)
res.sum[molten.dis[Reservoir == "RO1"], 
       out.dis_m3s := i.dis.sum, on = .(year, month, Reservoir)]
# Outflow from RO2 (outflow)
res.sum[molten.dis[Reservoir == "RO2"], 
        out.dis_m3s := i.dis.sum, on = .(year, month, Reservoir)]
# Outflow from RO3 (outflow)
res.sum[molten.dis[Reservoir == "RO3"], 
       out.dis_m3s := i.dis.sum, on = .(year, month, Reservoir)]


# Inflows depend on the year, sometimes river, other times the reservoir before in the continuum
# Inflow into RO1 (constant, RO2 outflow)
res.sum[molten.dis[Reservoir == "RO2", ][, Reservoir := "RO1"], 
        in.dis_m3s := i.dis.sum, on = .(year, month, Reservoir)]

# Inflow into RO2
## River
res.sum[mon.riv[Reservoir == "RO2" & (year == 2015 | year == 2016),], in.dis_m3s := i.mean,
       on = .(year, month, Reservoir)]
## RO3 outflow
res.sum[molten.dis[Reservoir == "RO3" & !(year == 2015 | year == 2016), ][, Reservoir := "RO2"], 
        in.dis_m3s := i.dis.sum, on = .(year, month, Reservoir)]

# Inflow into RO3
## River
res.sum[mon.riv[Reservoir == "RO3" & (year < 2020),], in.dis_m3s := i.mean,
       on = .(year, month, Reservoir)]
## RO4 outflow
#master[molten.dis[Reservoir == "RO4" & !(year < 2020), ][, Reservoir := "RO3"], 
#       dis.res.in := i.dis.out, on = .(year, sample.type.year == Reservoir)]

# get relevant months
res.sum <- res.sum[month == 6 | month == 8 | month == 10,]
# omit any NAs (= dam not built yet)
res.sum <- res.sum[!is.na(out.dis_m3s),][out.dis_m3s != 0,]

# merge with reservoir stats
res.sum <- res.sum[sum.wat, , on = .(Reservoir, year, month), nomatch=NULL]

# save this
#write.table(res.sum, "./Data/Summary/mean_res.vol_dis.inout.csv", sep = ',', dec = ".", row.names = F)

# Calculate residence time:
# 1. Transform units
res.sum[,c("mean.vol_m3",
           "gis.vol_m3") := list(mean.vol_hm3 * 1000000,
                                 gis.volume_hm3 * 1000000)]

# 2. Calculate residence time in s and convert to hours and days
# decided to take estimated volume based on models rather than GIS
res.sum[, wrt_s := mean.vol_m3 / out.dis_m3s]
res.sum[, wrt_h := wrt_s / (60 * 60)]
res.sum[, wrt_d := wrt_h / 24]

# Number of pixels in each reservoir
# Pixel area in km2
pixel.area2 <- (18*18) * 0.000001

ggplot(res.sum, aes(x = out.dis_m3s, y = mean.vol_hm3, colour = Reservoir)) +
  geom_point()

# Add main channel area from GIS
# Estimated by taking the La Romaine river shape, creating a buffer of 300m and 200m.
# Then, subtract the buffer area  from the Reservoir shapes
# 300m was chosen for RO2 as the river size before flooding ranged between 150-250m
# 200m was chosen for RO1, since the reservoir itself is very small and resembles greatly the river

# res.sum[Reservoir == "RO1", buffer.area_km2 := 10.139482]
# res.sum[Reservoir == "RO2", buffer.area_km2 := 38.303566]
# res.sum[Reservoir == "RO3", buffer.area_km2 := 20.922903]

res.sum[Reservoir == "RO1", pxl.no_18 := 1267] # number of pixels along main channel, times pixel area
res.sum[Reservoir == "RO2", pxl.no_18 := 2932]
res.sum[Reservoir == "RO3", pxl.no_18 := 1520]

# WRT in each pixel
res.sum[, wrt.px18_d := wrt_d / pxl.no_18]

setcolorder(res.sum, neworder = c("Reservoir","year","month","out.dis_m3s","in.dis_m3s",
                                  "mean.vol_m3", "gis.vol_m3", "mean.area_km2", "gis.area_km2",
                                  "wrt_d", "pxl.no_18","wrt.px18_d"))
out <- res.sum

out[, wrt_s := wrt_d *86400] # convert days to seconds
out[, wrt.px18_s := wrt_s / pxl.no_18] # calculate the time in seconds within each pixel
out[, veloc.px18_ms := 18 / wrt.px18_s ] # get velocity by dividing the length of pixel

mean.wrt <- out[, .(mean.wrt_d = mean(wrt_d, na.rm = T),
                    mean.veloc_ms = mean(veloc.px18_ms, na.rm = T)), by = .(Reservoir)]

# export another discharge summary

(p <- ggplot(out, aes(x = as.factor(month), y = wrt_d)) +
  ggpubr::theme_pubr() +
  geom_boxplot(aes(fill = Reservoir)) +
  labs(x = "Month", y = "WRT (d)") +
  scale_fill_viridis_d())
#ggsave("./Figures/Methods/WRT_byseason_RO.png", p, height = 10, width = 10, units = "cm")

# write.table(mean.wrt, "./Output/reservoir_wrt_seasonal_means.csv", sep = ",", row.names = F)
# write.table(out, paste0("./Output/reservoir_wrt_seasonal_perpixel_", Sys.Date(), ".csv"), sep = ",", row.names = F)

rm(mean.wrt, molten.dis, mon.dis, mon.riv, p, res.sum)

# At this point we have the water residence time and
# velocity in each pixel along the reservoir main channel

# Fluvial systems --------------------------------------------------------------------------------------------------------------

# We had done many variations to come up with discharge and velocity models
# A discharge model for the entire watershed using precipitation, flow accumulation and difference in snow depth between months.
# A cross-sectional area model from sampled rivers and streams to predict log10(cross-sectional area) from flow accumulation by flow condition.
# We calculated velocity using these two models with equations v = Q / A
# this only captured streams and rivers with a catchment area above 550 km2
#
# For streams with a catchment area less than 550 km2, we had created a separate model to directly predict
# velocity from the slope in the pixel of where the stream was sampled
#
# These solutions were acceptable, albeit not perfect
#
# In 2021-2022 campaigns, an additional sub-watershed Bernard was sampled that covered more stream orders
# This additional dataset allowed us to create a universal discharge and velocity model by flow condition

# Old scripts for models can be found in 'old_hydrology_models.R'

# Read in velocity + morphology data from Hydro-Quebec
vel <- read.table("./GIS/veloc_stations_flacc.txt", sep = ",", stringsAsFactors = F, header = T) %>% setDT()
# Two stations were wrongly snapped
resnap <- read.table("./GIS/resnap_veloc_stations_flacc.txt", sep = ",", 
                     stringsAsFactors = F, header = T) %>% setDT()

# remove unnecessary columns
vel <- vel[,c("FID","FID_") := list(NULL, NULL)]
resnap <- resnap[,c("FID","FID1","FID_") := list(NULL, NULL,NULL)]
# combine
vel <- bind_rows(vel[!(station %in% resnap$station),], resnap)
# Transform to date
vel[, date := as.Date(date, format = "%Y-%m-%d")]
vel[, month := month(date)]

# Rename column names, were cut off by ArcMap
cols <- c("station","date", "lat", "long", "flow.condition", "flow.condition.fr", "location",
          "width_m", "area_m2", "velocity_ms", "discharge_m3s", "max.vel_ms", "max.depth_m",
          "mean.depth_m", "water.temp_C", "snapped_id", "snapped_dist","snapped_x", "snapped_y",
          "elev_m", "flacc_npx", "strahler.order", "slope_deg", "month")
colnames(vel) <- cols 
rm(cols)
# ROMA0944, ROMA0943 wrongly snapped, strahler 0
# check ROMA0832, too small width for flacc - don't bother, it's an NA field that was converted to 0 by ArcMap

# Correct 0 values to NA
vel[is.na(date), c("width_m","area_m2","velocity_m", "discharge_", "max_vel_ms", "max_depth_",
                   "mean_depth", "water_temp") := list(NA, NA, NA, NA, NA, NA, NA, NA)]

# Convert flow accumulation in numbers of pixels to km2 = catchment area
# (easier to be implemented by other people using a DEM with a different resolution)
pixel.area <- 18*18 # meters
vel[, flacc_km2 := (flacc_npx * pixel.area) * 0.000001] # m to km

# Streams ----------------------------------------------------------------------------------------------------
# Small streams sampled by us range from 0.97 to 230

# Read Marie's Petite-Romaine (PR) data ----------------------------------------------------------------------
pr <- read_excel("./Data/Hydrology/Velocity/PR_WRT_results_V20190507.xlsx", sheet = 1)
# petite romaine data
rest <- read.table("./GIS/restsamples_flacc.txt", sep = ",", header = T)
# samples not non main channel
sub.ws <- read.table("./GIS/prsamples_subws_avg_slope.txt", sep = ",", header = T)
# avg slope of samples in Petite Romaine
prsamp <- read.table("./GIS/prsamples.txt", sep = ",", header = T) %>%
  dplyr::select(match.ws = FID_, sample.name = sample_nam)
prsamp$match.ws <- as.character(prsamp$match.ws)
sub.ws$match.ws <- as.character(sub.ws$VALUE) # slope is in VALUE
sub.ws <- merge(prsamp, sub.ws, by = "match.ws")

pr$match %in% rest$sample_nam # ok
pr$match %in% sub.ws$sample.name # ok

rest <- rest %>% dplyr::select(sample.name = sample_nam, flacc_npx = flacc_px, 
                               strahler.order = strahler_o, slope_deg = slope_deg)
sub.ws <- sub.ws %>% dplyr::select(sample.name, slope_subws_deg = MEAN)

pr <- merge(pr, rest, by.x = "match", by.y = "sample.name") %>%
  mutate(date = as.Date(pr$sampling_date, origin = "1904-01-01")) %>% setDT()
pr <- merge(pr, sub.ws, by.x = "match", by.y = "sample.name")
pr[, area_m2 := width_m * depth_m]
pr[, c("Year", "month") := list(year(date), month(date))]
pr[, flow.condition := ifelse(month == 6, "high", "low")]
pr[, dis_m3s := as.numeric(meas_disch_m3s)]
pixel.area <- 18*18 # meters
pr[, flacc_km2 := (flacc_npx * pixel.area) * 0.000001] # m to km
pr[, c("velocity_ms", "velocity_m.s1") := list(velocity_m.s1, NULL)]

# Read in Gabriel's Bernard data (2021-2022) ---------------------------------------
ber <- read.csv("./GIS/bernard_streams_stats.txt", sep = ";", dec = ',', stringsAsFactors = F)
ber <- ber %>% dplyr:: select(sample.name = site, flacc_npx, strahler.order = strah_ord, slope_deg = slope_dec)

ber.vel <- read.csv("./Data/Hydrology/bernard_discharge_sonde.csv", sep = ",", stringsAsFactors = F)
ber.vel <- ber.vel %>% dplyr::select(sample.name = site, date, transect_distance_m, area_m2,
                                     velocity_ms = velocity_m_s, avg_velocity, dis_m3s = discharge_m3_s) %>% setDT()
# summarize
# clean zero velocity
ber.vel <- ber.vel[avg_velocity <= 0, avg_velocity := NA]
ber.vel <- ber.vel[, .(width_m = max(transect_distance_m, na.rm = T),
                       vel_ms = mean(avg_velocity, na.rm = T),
                       area_m2 = sum(area_m2, na.rm = T),
                       dis_m3s = unique(dis_m3s)), by = list(sample.name, date)]
ber.vel[, date := as.Date(date, format = "%m/%e/%Y")] # make date

ber.vel <- merge(ber, ber.vel, by = "sample.name") %>% setDT()
ber.vel[, c("Year","Month") := list(year(date), month(date))]
ber.vel[, flow.condition := ifelse(Month == 6, "high", "low")]
ber.vel[, flacc_km2 := (flacc_npx * pixel.area) * 0.000001] # m to km
ber.vel[, depth_m := area_m2 / width_m]

# merge PR and Bernard
st <- bind_rows(pr %>% dplyr::select(sample.name = sample, date, Year, month, flow.condition,
                     strah_ord = strahler.order, flacc_km2,
                     slope_deg, width_m, depth_m, area_m2, dis_m3s, velocity_ms) %>% mutate(project = "PR"),
          ber.vel %>% dplyr::select(sample.name, date, Year, month = Month, flow.condition,
                          strah_ord = strahler.order, flacc_km2,
                          slope_deg, width_m, depth_m, area_m2, dis_m3s, velocity_ms = vel_ms) %>% mutate(project = "BR"))

# combine with Hydro-Quebec
fluv <- bind_rows(vel %>% dplyr::select(sample.name = station, date, flow.condition, flacc_km2, strah_ord = strahler.order,
                                       width_m, depth_m = mean.depth_m, area_m2, dis_m3s = discharge_m3s, velocity_ms) %>% 
                   mutate(project = "HQ"),
                 st %>% dplyr::select(project, sample.name, date, flow.condition, flacc_km2, strah_ord,
                                      width_m, depth_m, area_m2, dis_m3s, velocity_ms))

fluv <- fluv[!is.na(date),]

# plot some potential relationships
ggplot(fluv, aes(y = log10(dis_m3s), x = log10(flacc_km2), colour = flow.condition)) +
  theme_bw() +
  geom_smooth(method = "lm", formula = "y~poly(x,2)") +
  geom_point(aes( shape = project))

ggplot(fluv, aes(y = log10(velocity_ms), x = log10(flacc_km2), colour = flow.condition)) +
  theme_bw() +
  geom_smooth(method = "lm", formula = "y~poly(x,2)") +
  geom_point(aes( shape = project))

ggplot(fluv, aes(y = velocity_ms, x = as.factor(strah_ord), colour = flow.condition)) +
  theme_bw() +
  geom_boxplot()

ggplot(fluv, aes(y = log10(dis_m3s), x = as.factor(strah_ord), colour = flow.condition)) +
  theme_bw() +
  geom_boxplot()

# Velocity model ---------------------------------------------------------
lm.one <- lm(log10(velocity_ms) ~ log10(flacc_km2), data = fluv)
lm.fc <- lm(log10(velocity_ms) ~ log10(flacc_km2) + flow.condition, data = fluv)
poly.one <- lm(log10(velocity_ms) ~ poly(log10(flacc_km2),2), data = fluv)
poly.fc <- lm(log10(velocity_ms) ~ poly(log10(flacc_km2),2) + flow.condition, data = fluv)

anova(lm.one, lm.fc) # lm.fc is better
anova(poly.one, poly.fc) # poly.fc is better
anova(lm.fc, poly.fc) #poly.fc is best

summary(poly.fc) #R2 = 0.44, p < 0.0001

fluv[, est.vel_ms := 10^predict(poly.fc, newdata = fluv)]

ggplot(fluv, aes(x = velocity_ms, y = est.vel_ms)) +
  geom_point(aes(colour = flow.condition)) +
  scale_colour_viridis_d(option = "cividis", end = 0.7, name = "Flow condition") +
  geom_abline(slope = 1, intercept = 0) 

# Discharge model -------------------------------------------------------------
q.lm.one <- lm(log10(dis_m3s) ~ log10(flacc_km2), data = fluv)
q.lm.fc <- lm(log10(dis_m3s) ~ log10(flacc_km2) + flow.condition, data = fluv)
q.poly.one <- lm(log10(dis_m3s) ~ poly(log10(flacc_km2),2), data = fluv)
q.poly.fc <- lm(log10(dis_m3s) ~ poly(log10(flacc_km2),2) + flow.condition, data = fluv)

anova(q.lm.one, q.lm.fc) # lm.fc is better
anova(q.poly.one, q.poly.fc) # poly.fc is better
anova(q.lm.fc, q.poly.fc) #poly.fc is best

summary(q.poly.fc) #R2 = 0.87, p < 0.0001

# Prepare for plotting --------------------------------------------------------
# Significance function
abbrev.p <- function(x){
  if(x < 0.0001){
    out <- c("< 0.0001","****")
  } else if(x <= 0.001){
    out <- c("< 0.001","***")
  } else if(x <= 0.01){
    out <- c("< 0.01","**")
  } else if(x <= 0.05){
    out <- c("< 0.05","*")
  } else {
    out <- c(as.character(paste("=",round(x, 2))),"n.s.", round(x, 2))
  }
  return(out)
}

fluv[, project := factor(project, levels = c("PR", "BR","HQ"), labels = c("Petite Romaine", "Bernard", "Hydro-Québec"))]

lm_coef <- list(a = as.numeric(round(summary(q.poly.fc)$coefficients[1], digits = 2)),
                b = as.numeric(round(summary(q.poly.fc)$coefficients[2], digits = 5)),
                b1 = as.numeric(round(summary(q.poly.fc)$coefficients[3], digits = 2)),# average difference between high and low flow
                r2 = round(summary(q.poly.fc)$r.squared, digits = 2),
                pval = abbrev.p(anova(q.poly.fc)$`Pr(>F)`[1])[1])

lm_coef$b <- abs(lm_coef$b)
lm_coef$b1 <- abs(lm_coef$b1)
lm_coef$log <- "log"
lm.eq <- c(substitute(log ~ italic(y) == a + b %.% log ~ italic(x), lm_coef), # high flow
           substitute(log ~ italic(y) == a + b - b1 %.% log ~ italic(x), lm_coef), # low flow
           substitute(~~italic(R)^2~"="~r2*","~~italic(p)~pval, lm_coef))# Just b is high flow

(qp <- ggplot(fluv, aes(y = log10(dis_m3s), x = log10(flacc_km2), colour = flow.condition)) +
  theme_bw() +
  scale_colour_viridis_d(option = "cividis", end = 0.7, name = "Flow condition") +
  geom_point(aes(shape = project), size = 2.5) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)") +
  annotate(geom = "text",
           y = -2,
           x = 2.3,
           label = lm.eq[1],
           parse = T, colour = viridis::cividis(n = 2, end = 0.7)[1]) +
  annotate(geom = "text",
           y = -2.3,
           x = 2.3,
           label = lm.eq[2],
           parse = T, colour = viridis::cividis(n = 2, end = 0.7)[2]) +
  annotate(geom = "text",
           y = -2.6,
           x = 2.3,
           label = lm.eq[3],
           parse = T) +
  scale_shape_manual(name = "Project", values = c(25,24,21)) +
  labs(x = expression(paste("Catchment area (km"^2,", log)")),
       y = expression(paste("Discharge (m"^3," s"^-1,", log)"))))

lm_coef <- list(a = as.numeric(round(summary(poly.fc)$coefficients[1], digits = 2)),
                b = as.numeric(round(summary(poly.fc)$coefficients[2], digits = 5)),
                b1 = as.numeric(round(summary(poly.fc)$coefficients[3], digits = 2)),# average difference between high and low flow
                r2 = round(summary(poly.fc)$r.squared, digits = 2),
                pval = abbrev.p(anova(poly.fc)$`Pr(>F)`[1])[1])

lm_coef$b <- abs(lm_coef$b)
lm_coef$b1 <- abs(lm_coef$b1)
lm_coef$log <- "log"
lm.eq <- c(substitute(log ~ italic(y) == a + b %.% log ~ italic(x), lm_coef), # high flow
           substitute(log ~ italic(y) == a + b - b1 %.% log ~ italic(x), lm_coef), # low flow
           substitute(~~italic(R)^2~"="~r2*","~~italic(p)~pval, lm_coef))# Just b is high flow

(vp <- ggplot(fluv, aes(y = log10(velocity_ms), x = log10(flacc_km2), colour = flow.condition)) +
  theme_bw() +
  scale_colour_viridis_d(option = "cividis", end = 0.7, name = "Flow condition") +
  geom_point(aes(shape = project), size = 2.5) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)") +
  scale_shape_manual(name = "Project", values = c(25,24,21)) +
  annotate(geom = "text",
           y = -2,
           x = 2.3,
           label = lm.eq[1],
           parse = T, colour = viridis::cividis(n = 2, end = 0.7)[1]) +
  annotate(geom = "text",
           y = -2.2,
           x = 2.3,
           label = lm.eq[2],
           parse = T, colour = viridis::cividis(n = 2, end = 0.7)[2]) +
  annotate(geom = "text",
           y = -2.4,
           x = 2.3,
           label = lm.eq[3],
           parse = T) +
  labs(x = expression(paste("Catchment area (km"^2,", log)")),
       y = expression(paste("Velocity (m s"^-1,", log)"))))


(p <- ggarrange(qp, vp, common.legend = T, ncol = 2, legend = "right"))
#ggsave("./Figures/Methods/models_flacc_Q&v_05.02.2024.png", width = 25, height = 10, units = "cm")

# saveRDS(q.poly.fc, "./Objects/poly.model_flacc_discharge_flowcond_05.02.2024.rds")
# saveRDS(poly.fc, "./Objects/poly.model_flacc_velocity_flowcond_05.02.2024.rds")


# Bring all models together ------------------------------------------------------------------------------
streamnet <- read.csv("./Data/Traveltime/streamnet_pnts.txt", sep = ";", dec = ",",
                      stringsAsFactors = F) %>% setDT()

# Calculate catchment area = drainage area
streamnet[, flacc_km2 := flacc_px * (18*18) * 0.000001]

unique(streamnet[flacc_km2 < 550, ]$strah_ord)
length(unique(streamnet[flacc_km2 <= 550 & strah_ord == 4, ]$reach_id))
length(unique(streamnet[flacc_km2 > 550 & strah_ord == 4, ]$reach_id)) # only 1 reach

# Create seasonal data
# 2/3 months, 4 years
# 2015 - 2016 - 2017 -2018
# 6 - 8 - 10
newdat <- rbind(streamnet, streamnet, streamnet, streamnet)

# Add Year and month identifier
newdat[, Year := rep(2015:2016, each = nrow(streamnet), times = 2)] # set-up: 2015, 2016, 2015, 2016
newdat[, Month := rep(c(6,8), each = nrow(streamnet)*2)]

# sanity check
newdat[,.(n = .N), by = list(Year, Month)]

newdat.plus <- rbind(streamnet, streamnet, streamnet, streamnet, streamnet, streamnet)

# Add Year and month identifier
newdat.plus[, Year := rep(2017:2018, each = nrow(streamnet), times = 3)] # set-up: 2017, 2018, 2017, 2018,  2017, 2018
newdat.plus[, Month := rep(c(6,8,10), each = nrow(streamnet)*2)]

# sanity check
newdat.plus[,.(n = .N), by = list(Year, Month)]

# Add flow condition
newdat <- rbind(newdat, newdat.plus)
newdat[, flow.condition := ifelse(Month == 6, "high","low")]

# arrange
newdat <- newdat[order(Year, Month, pointid),]

# Apply models -------------------------------------------------------------------------------------------------------------------
dis_model <- readRDS("./Objects/poly.model_flacc_discharge_flowcond_05.02.2024.rds")
vel_model <- readRDS("./Objects/poly.model_flacc_velocity_flowcond_05.02.2024.rds")

newdat[, dis_m3s := 10^predict(dis_model, newdat = newdat)]
newdat[, velocity_ms := 10^predict(vel_model, newdat = newdat)]

# newdat[flacc_km2 < 550, width_m := predict(width.streams, newdat = newdat[flacc_km2 < 550, ])]
ggplot(newdat[,.(vel = mean(velocity_ms)), by = list(Year, Month, flow.condition, strah_ord)],
       aes(x = strah_ord, y = vel)) +
  geom_point(aes(colour = flow.condition))

# Boxplot

(bq <- ggplot() +
    theme_bw() +
    geom_boxplot(data = fluv, aes(y = log10(dis_m3s), x = as.factor(strah_ord), colour = flow.condition)) +
    geom_point(data = newdat[,.(dis = median(dis_m3s)), by = list(flow.condition, strah_ord)],
               aes(x = strah_ord, y = log10(dis), group = flow.condition), 
               position = position_dodge(0.8), shape = 18, size = 3, colour = "grey20") +
    scale_colour_viridis_d(option = "cividis", end = 0.7, name = "Flow\ncondition") +
    labs(x = "Strahler order",
         y = expression(paste("Discharge (m"^3," s"^-1,", log)"))) +
    theme(legend.box.margin = margin(2)))

(bv <- ggplot() +
    theme_bw() +
    geom_boxplot(data = fluv, aes(y = log10(velocity_ms), x = as.factor(strah_ord), colour = flow.condition)) +
    geom_point(data = newdat[,.(vel = median(velocity_ms)), by = list(flow.condition, strah_ord)],
               aes(x = strah_ord, y = log10(vel), group = flow.condition), 
               position = position_dodge(0.8), shape = 18, size = 3, colour = "grey20") +
    scale_colour_viridis_d(option = "cividis", end = 0.7, name = "Flow\ncondition") +
    labs(x = "Strahler order",
         y = expression(paste("Velocity (m s"^-1,", log)"))) +
    theme(legend.box.margin = margin(2)))

(p <- ggarrange(bq, bv, ncol = 2, common.legend = T, legend = "right"))
#ggsave("./Figures/Methods/dis_vel_est.vs.meas_byStreamorder_05.02.2024.png", p, width = 21, height = 10, units = "cm")

# Sanity check
any(is.na(newdat$velocity_ms)) # all pixels have an estimate? = should be FALSE
#newdat[velocity_ms < min.vel,] # all pixels are above minimum velocity?

# calculate water residence time in pixel --------------------------------------------------------------------
newdat[, wrt_s := 18 / velocity_ms] # pixel has 18 m
# Convert to minutes
newdat[, wrt_min := wrt_s / 60]

# Extract only necessary columns and save as table
out <- newdat %>% dplyr::select(pointid, coord_x, coord_y,
                  Year, Month, flow.condition,
                  reach_id, strah_ord,
                  #elev_m, slope_deg,
                  #mmean.mean.temp, delta.snow,
                  flacc_px, flacc_km2,
                  dis_m3s, velocity_ms, wrt_min)
                  #area_m2)

# Extract a minimum data frame (without pointids)
mean.wrt <- out[,.(mean.wrt = mean(wrt_min), mean.vel = mean(velocity_ms),
                   #mean.slope = mean(slope_deg), mean.area = mean(area_m2), 
                   mean.dis = mean(dis_m3s),
                   mean.flacc = mean(flacc_km2)) , by = list(Year, Month, flow.condition, strah_ord)]
# seasonality is expressed in strahler order 5-7 rivers, order 1-4 are seasonality less

ggplot(mean.wrt, aes(x = mean.dis, y = mean.vel, colour = strah_ord)) +
  theme_bw() +
  facet_grid(Month~Year) +
  geom_point() +
  scale_colour_viridis_c(option = "cividis")

# export data
out <- readRDS("./Data/Traveltime/streamnet_pnts_time_05.02.2024.rds")

# write.table(out, "./Data/Traveltime/streamnet_pnts_time_05.02.2024.csv", sep = ",", dec = ".", row.names = F)
# saveRDS(out, "./Data/Traveltime/streamnet_pnts_time_05.02.2024.rds")

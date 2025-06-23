## -------------------------------------------------------------------------
##
## Script name: 2_find_downstream.pixel.R
##
## Purpose of script: We track down the pixels along the flow path and identify
##                    which pixels are inlets, outlets, tributaries, and which
##                    flow paths are within standing water bodies.
##
## Author: Masumi Stadler
##
## Date Finalized: 11/3/2025
##
## Copyright (c) Masumi Stadler, 2025
## Email: m.stadler.jp.at@gmail.com
##
## -------------------------------------------------------------------------
##
## Notes: Be careful, depending on the year being processed, uncomment or
##        re-comment lines of code.
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
pckgs <- list("data.table", "tidyverse", # wrangling
              "plyr", "ggpubr",
              "doParallel", "foreach", "doSNOW"
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

# Read in data -----------------------------------------------------------

# get stream network point dataset
streamnet_pnts <- readRDS("./Data/Traveltime/streamnet_pnts_time_05.02.2024.rds")  %>% setDT()

# flow length of each pixel
streamnet_fllength <- read.csv("./Data/Traveltime/flacc3000_streamnet_pnts_withfllength.txt", sep = ";", dec = ",",
                           stringsAsFactors = F)  %>% setDT()

# merge flow length info to streamnet point
streamnet_pnts[streamnet_fllength, "fllength_dec" := i.cdem_fllen, on = .(pointid)]
rm(streamnet_fllength)

# Get polygon lines
streamnet_lines <- read.csv("./Data/Traveltime/streamnetwork_poly.txt",
                            sep = ";", dec = ",", stringsAsFactors = F)  %>% setDT()

streamnet_stats <- read.csv("./Data/Traveltime/streamnetwork_poly_reachstats.txt",
                            sep = ";", dec = ",", stringsAsFactors = F) %>% setDT()

# streamnet_lines has all the lines in the stream network
# Nodes indicate which lines are connected to each other

# Identify the reach ID of each polyline
streamnet_lines[!(streamnet_lines$FID %in% streamnet_stats$FID_),]
# Some lines were too short, and connected to separate stream reaches
# Dropped during "Zonal statistics as Table"

# Add manually
# Numbers were checked manually on ArcGIS

short_lines <- data.frame(OID_ = rep(NA, times = 4),
                          FID_ = c(1004, 1973, 2408, 4731),
                          COUNT = rep(NA, times = 4),
                          AREA = rep(NA, times = 4),
                          MIN = rep(NA, times = 4),
                          MAX = rep(NA, times = 4),
                          RANGE = rep(NA, times = 4),
                          MEAN = rep(NA, times = 4),
                          STD = rep(NA, times = 4),
                          SUM = rep(NA, times = 4),
                          VARIETY = rep(NA, times = 4),
                          MAJORITY = c(1395, 2810, 3001, 4761),
                          MINORITY = c(1396, 2811, 3002, 4762),
                          MEDIAN = rep(NA, times = 4)) %>% setDT()

# Merge
streamnet_stats <- bind_rows(streamnet_stats, short_lines)

# Check if we have all streams now
streamnet_lines[!(streamnet_lines$FID %in% streamnet_stats$FID_),] # Empty, so yes

# Identify sources and confluences ----------------------------------------------
# Node numbers are unique
# All "From" nodes that do not appear in "To" nodes, should be sources
streamnet_lines[!(streamnet_lines$from_node %in% streamnet_lines$to_node), flag_source := "yes"]
streamnet_lines[streamnet_lines$from_node %in% streamnet_lines$to_node, flag_conf := "yes"]

# Add reach ID information
# We need majority and minority
streamnet_lines[streamnet_stats, c("reach_id.maj", "reach_id.min") :=
                  list(i.MAJORITY, i.MINORITY), on = c("FID" = "FID_")]

# Export source and junction coordinate starts to check in GIS
# write.table(streamnet_lines[flag_source == "yes",] %>% dplyr::select(FID, start_x, start_y),
#             "./Data/Traveltime/streamlines_sources.txt", sep = ";", dec = ",", row.names = F)
# 
# write.table(streamnet_lines[flag_conf == "yes",] %>% dplyr::select(FID, start_x, start_y),
#             "./Data/Traveltime/streamlines_conf.txt", sep = ";", dec = ",", row.names = F)


# Organize stream network -------------------------------------------------------------------------------------------------------------
# turn coordinates into character
streamnet_pnts[, c("coord_x", "coord_y") := list(as.character(coord_x),
                                                 as.character(coord_y))]

streamnet_lines[, c("start_x","start_y", "end_x", "end_y") := 
                  list(as.character(start_x),
                       as.character(start_y),
                       as.character(end_x),
                       as.character(end_y))]

# Merge coordinates
streamnet_pnts[, coord := paste(coord_x, coord_y, sep = " ")]
streamnet_lines[, c("start_coord", "end_coord") := list(paste(start_x, start_y, sep = " "),
                                                        paste(end_x, end_y, sep = " "))]

# Calculate elevation difference -------------------------------------------------------------------------------------------------
# Order by reach ID and flow accumulation
streamnet_pnts <- streamnet_pnts[order(reach_id, flacc_px),]

pixel.size <- 18
#streamnet_pnts[, slope_m := (elev_m - shift(elev_m, 1)) / pixel.size, by = .(reach_id)]

# Add node information to point df ------------------------------------------------------------------------------------------------
# Merge by reach ID
streamnet_pnts[streamnet_lines, c("from_node","to_node") := list(i.from_node,
                                                                 i.to_node), on = c("reach_id==reach_id.maj")]
# Which reach IDs were not identified using zonal statistics on polylines?
not.found <- unique(streamnet_pnts[!(reach_id %in% streamnet_lines$reach_id.maj),]$reach_id)
out <- data.frame(reach_id = not.found[order(not.found)])

# Export to csv and manually check via...
# Streamnet_pnts -> Search by Attributes ['reach_id' = xx] -> Zoom to -> Identify polyline -> Extract nodes from and to
#write.table(out, "./Data/Traveltime/missing_reaches.csv", sep = ",", dec = ".", row.names = F)

# Read in manually added Node info
missing.reaches <- read.csv("./Data/Traveltime/missing_reaches_man.filled.csv",
                            sep =",", dec = ".", stringsAsFactors = F) %>% setDT()

# Merge with points
streamnet_pnts[missing.reaches, c("from_node", "to_node") := list(i.from_node, i.to_node),
               on = .(reach_id)]
# Sanity check
any(is.na(streamnet_pnts$from_node)) # all should have a node

# All points have a node now
# Most of these points were on reaches that only had one or 
# max two points on the polyline

# Sanity check
# Do we have a start node for each end node?
out <- data.frame(from_node = unique(streamnet_pnts[!(to_node %in% from_node),]$to_node))
# save lines without a start 
# write.table(out, "./Data/Traveltime/missing_starts.csv", sep = ",", dec = ".", row.names = F)

# identified manually in ArcGIS
missing.starts <- read.csv("./Data/Traveltime/missing_starts_man.filled.csv", sep = ",", dec = ".", stringsAsFactors = F) %>% setDT()

# Two reasons for missing starts
wrong.extract <- missing.starts[is.na(pointid),] # Somehow, reach ID and node info was wrongly extracted
incomplete.reach <- missing.starts[!is.na(pointid),] # Somehow, a reach was not divided into two although there were two lines

# Update points data frame
streamnet_pnts[wrong.extract, c("from_node", "to_node") := list(i.from_node, i.to_node), on = .(reach_id)]
streamnet_pnts[incomplete.reach, c("from_node", "to_node") := list(i.from_node, i.to_node), on = .(pointid)]

unique(streamnet_pnts[!(to_node %in% from_node),]$to_node)
# 7869 is the mouth

# Add flag
mouth.dt <- streamnet_pnts[to_node == 7869, .SD[which.max(fllength_dec)]][, flag_mouth := 1]
# Merge flag info back to main point dataframe
streamnet_pnts[mouth.dt, flag_mouth := i.flag_mouth, on = .(pointid)]
streamnet_pnts[is.na(flag_mouth), flag_mouth := 0]
rm(mouth.dt)

# Get source information ---------------------------------------------------------------------------------------------------------
# within each reach of strahler order 1, the smallest flow accumulation pixel is the source
source.dt <- streamnet_pnts[strah_ord == 1, .SD[which.min(fllength_dec)], by = .(reach_id)][, flag_source := 1]

# Sanity check
# Do all reaches with strahler order 1 have a source?
any(!(streamnet_pnts[strah_ord == 1,]$reach_id %in% source.dt$reach_id)) # should be FALSE

# Merge flag info back to main point dataframe
streamnet_pnts[source.dt, flag_source := i.flag_source, on = .(pointid)]
streamnet_pnts[is.na(flag_source), flag_source := 0]
rm(source.dt)

# At the same time, the first pixel in a reach of strahler order >= 2 should be a confluence
conf.dt <- streamnet_pnts[strah_ord >= 2, .SD[which.min(fllength_dec)], by = .(reach_id)][, flag_conf := 1]

# Sanity check
# Do all reaches with strahler order >2 have a confluence?
any(!(streamnet_pnts[strah_ord >= 2,]$reach_id %in% conf.dt$reach_id)) # should be FALSE

# Merge flag info back to main point dataframe
streamnet_pnts[conf.dt, flag_conf := i.flag_conf, on = .(pointid)]
streamnet_pnts[is.na(flag_conf), flag_conf := 0]
rm(conf.dt)

# Identify reach end --------------------------------------------------------------------------------------------------------------
# point within reach with largest flow legnth is the end point
end.dt <- streamnet_pnts[, .SD[which.max(fllength_dec)], by = .(reach_id)][, flag_end := 1]

# Merge flag info back to main point dataframe
streamnet_pnts[end.dt, flag_end := i.flag_end, on = .(pointid)]
streamnet_pnts[is.na(flag_end), flag_end := 0]
rm(end.dt)

# Reconstruct downstream pixel ----------------------------------------------------------------------------------------------------
# Let's try and reconstruct where the pixels flow into
# Order by node and flow accumulation
streamnet_pnts <- streamnet_pnts[order(from_node, fllength_dec),]

saveRDS(streamnet_pnts, "./Data/Traveltime/flacc3000_streamnet_pnts_cleaned_05.02.2024.rds")

# Extract a sub data set since downstream flow is the same for all years and months
dt <- streamnet_pnts[Month == 8 & Year == 2015,]
# get rid of duplicates
dt <- dt[!duplicated(pointid),]

# Count number of pixels in reach
dt[, px_count := .N, by = .(reach_id)]

# Warning: This loop takes a while
# Loop: find the next downstream pixel ---------------------------------------------------------------------------------------------
# for(i in 1:nrow(dt)){
#     if(dt[i,]$flag_mouth == 1L){ # if the loop round is at the end of the network
#       dt[i, next.pixel := NA]
#     } else if(dt[i,]$flag_end == 1L){ # if the loop round is at the end of the reach
#       
#       if(dt[i,]$fllength_dec <= max(dt[from_node == dt[i,]$from_node,]$fllength_dec))
#       end <- dt[i,]$to_node # extract the node of the next reach
#       next_reach <- dt[dt$from_node == end, ][which.min(fllength_dec),]$pointid # get first pixel of next reach and it's ID
#       dt[i, next.pixel := next_reach] # fill row of next pixel
#     } else { # if we're within a reach, take the next point in the reach (matrix is ordered by flow length and reach ID)
#       next_px <- dt[i+1,]$pointid
#       dt[i, next.pixel := next_px]
#     }
#   cat(paste0("\rWorking on row ", i, " of ", nrow(dt), "."))
#   #print(paste("Working on row...", i))
# }

# out <- read.csv("./Data/Traveltime/streamnet_pnts_downstream.px.csv", stringsAsFactors = F)
# next.df <- out %>% dplyr::select(pointid, next.pixel)

# dt <- dt[next.df, , on = .(pointid)]
# dt[is.na(next.pixel),] # just the mouth

# Sanity check
dups <- dt[duplicated(next.pixel),]$next.pixel
dups <- dt[pointid %in% dups,]
dups[flag_conf != 1,] # should be empty, duplicates are only confluences

# Save next downstream pixel as data frame
out <- dt %>% 
  dplyr::select(pointid, next.pixel, reach_id, px_count, from_node, to_node, flag_source, flag_conf, flag_end, flag_mouth)

write.table(out, "./Data/Traveltime/streamnet_pnts_downstream.px.csv", sep = ",", dec = ".", row.names = F)

# # Add water body residence times --------------------------------------------------------------------------------------------------
# # Read in pixels that are inside lakes
# out <- read.csv("./Data/Traveltime/streamnet_pnts_downstream.px.csv", sep = ",", stringsAsFactors = F) %>% setDT()
# hl <- read.csv("./Data/Traveltime/flacc3000_streamnet_pnts_hydrolakes.txt",
#                sep = ";", dec = ",", stringsAsFactors = F) %>% setDT()
# buf <- read.csv("./Data/Traveltime/flacc3000_streamnet_pnts_withbuflakes.txt",
#                sep = ";", dec = ",", stringsAsFactors = F) %>% setDT()
# colnames(buf)[7] <- "buf.lake_id"
# 
# # merge hydrolakes lake ID and estimated average discharge to next pixel output
# temp <- temp %>% dplyr::select(pointid, fllength_dec)
# out <- merge(out, temp, by = 'pointid')
# out <- merge(out, buf %>% dplyr::select(pointid, lake_id, buf.lake_id, hl_avg_dis, hl_wrt_d), by = "pointid") %>% distinct() %>% setDT()
# 
# # replace -9999 = Nodata to 0
# out[, lake_id := ifelse(lake_id == -9999, 0, lake_id)]
# out[, buf.lake_id := ifelse(buf.lake_id == -9999, 0, buf.lake_id)]
# out[, hl_avg_dis := ifelse(hl_avg_dis == -9999, 0, hl_avg_dis)]
# 
# # Add identifier whether points are in reservoirs/lakes
# ro1o <- read.csv("./Data/Traveltime/flacc3000_streamnet_pnts_RO1.txt", sep = ";", stringsAsFactors = F)
# ro2o <- read.csv("./Data/Traveltime/flacc3000_streamnet_pnts_RO2.txt", sep = ";", stringsAsFactors = F)
# ro3o <- read.csv("./Data/Traveltime/flacc3000_streamnet_pnts_RO3.txt", sep = ";", stringsAsFactors = F)
# 
# # Identify points within reservoirs and lakes
# out[, lake_id := as.character(lake_id)]
# out[pointid %in% ro1o$pointid, lake_id := "RO1"]
# out[pointid %in% ro2o$pointid, lake_id := "RO2"]
# out[pointid %in% ro3o$pointid, lake_id := "RO3"]
# 
# write.table(out, "./Data/Traveltime/flacc3000_streamnet_wHL.csv", sep = ",", dec = ".", row.names = F)

out <- read.csv("./Data/Traveltime/flacc3000_streamnet_wHL.csv", sep = ",", stringsAsFactors = F) %>% setDT()
out[, buf.lake_id := as.character(buf.lake_id)]

# Update reservoir identifier with buffer (manually change) ----------------------------------------------------------------------------
# Added a buffer to standing water bodies

# Add identifier whether points are in reservoirs/lakes
ro1o <- read.csv("./Data/Traveltime/flacc3000_streamnet_pnts_RO1_buf.0004.txt", sep = ";", stringsAsFactors = F)
ro2o <- read.csv("./Data/Traveltime/flacc3000_streamnet_pnts_RO2_buf.0004.txt", sep = ";", stringsAsFactors = F)
ro3o <- read.csv("./Data/Traveltime/flacc3000_streamnet_pnts_RO3_buf.0004.txt", sep = ";", stringsAsFactors = F)

# Identify points within reservoirs and lakes
out[, lake_id := as.character(lake_id)]

#--- change me depending on the year ---#
# 2015
 out[pointid %in% ro2o$pointid, buf.lake_id := "RO2"]
# # 2016
# out[pointid %in% ro2o$pointid, buf.lake_id := "RO2"]
# out[pointid %in% ro1o$pointid, buf.lake_id := "RO1"]
# # 2017+2018
# out[pointid %in% ro2o$pointid, buf.lake_id := "RO2"]
# out[pointid %in% ro1o$pointid, buf.lake_id := "RO1"]
# out[pointid %in% ro3o$pointid, buf.lake_id := "RO3"]
rm(ro1o, ro2o, ro3o)
#--- change me depending on the year ---#

# Fill in pixels by flooded reservoir
# here, we added a buffer around lakes to identify
# a few outlets that were just outside of the lake
# Buffer size = 0.0004 decimal degrees
unique(out[grep("RO", buf.lake_id),]$buf.lake_id)

# Find lake outlets ---------------------------------------------------------------------------------
# (take buffered lakes)
all.wb <- unique(out[buf.lake_id != "0",]$buf.lake_id)
length(all.wb) # 1464 - 1465 - 1465
for(i in 1:length(all.wb)){
  #identify lake outlet, lake outlet is the point within lake that has the longest flow length
  outlet <- out[buf.lake_id == all.wb[i],][which.max(fllength_dec), ]$pointid
  # flag as outlet in main file
  out[pointid %in% outlet, flag_lake.out := 1]
}
# count number of outlets = sanity check
nrow(out[flag_lake.out == 1,]) == length(all.wb) # should be TRUE

# rest is no flag
out[is.na(flag_lake.out), flag_lake.out := 0]

# manually correct a few lakes where the outlet is out of the lake or inlet was misclassified

# Lake corrections applied to all years ---------------------------------------------------------------
out[pointid == 40310, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"7141","7141")]
out[pointid == 40034, flag_lake.out := 0]
out[pointid == 139200, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"93169","93169")]
out[pointid == 139201, c("flag_lake.out","buf.lake_id","lake_id") := list(0,"0","0")]
out[pointid == 157494, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"93376","93376")]
out[pointid == 156874, c("flag_lake.out","buf.lake_id","lake_id") := list(0,"0","0")]

out[pointid == 172570, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"910069","910069")]
out[pointid == 172320, c("flag_lake.out","buf.lake_id","lake_id") := list(0,"0","0")]
out[pointid == 172223, c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"910069","910069")] # two streams part of upstream lake
out[pointid == 172381, c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"910069","910069")]
out[pointid == 172632, c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"93524","93524")] # new inlet to the downstream lake

out[pointid == 186558, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"911456","911456")]
out[pointid == 186551, c("flag_lake.out","buf.lake_id","lake_id") := list(0,"0","0")]

out[pointid == 220362, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"93934","93934")] # new outlet
out[pointid == 217806, c("flag_lake.out","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 3777 & buf.lake_id == "93934", c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]

out[pointid == 216355, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"913609","913609")] # new outlet upstream lake
out[pointid == 217295, c("flag_lake.out","buf.lake_id","lake_id") := list(0,"0","0")] # remove old outlet upstream lake
out[pointid == 216265, c("flag_lake.trib","buf.lake_id","lake_id") := list(1,"913421","913421")] # new inlet downstream lake
out[reach_id == 3179 & buf.lake_id == "913421", c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")] # remove previous tribs to downstream lake
out[reach_id == 3170 & buf.lake_id == "913421", c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]

out[pointid == 235174, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"7358","7358")]
out[pointid == 236263, c("flag_lake.out","buf.lake_id","lake_id") := list(0,"0","0")]

out[pointid == 242776, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"94157","94157")]
out[pointid == 242927, c("flag_lake.out","buf.lake_id","lake_id") := list(0,"0","0")]

out[pointid == 309406, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"922993","922993")]
out[pointid == 309322, c("flag_lake.out","buf.lake_id","lake_id") := list(0,"0","0")]

out[pointid == 310568, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"95056","95056")]
out[pointid == 310570, c("flag_lake.out","buf.lake_id","lake_id") := list(0,"0","0")]

out[pointid == 354561, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"929077","929077")]
out[pointid == 354416, c("flag_lake.out","buf.lake_id","lake_id") := list(0,"0","0")]

out[pointid == 361270, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"95757","95757")]
out[pointid == 361103, c("flag_lake.out","buf.lake_id","lake_id") := list(0,"0","0")]

out[pointid == 402846, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"933686","933686")]
out[pointid == 402949, c("flag_lake.out","buf.lake_id","lake_id") := list(0,"0","0")]

out[pointid == 409583, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"934468","934468")]
out[pointid == 409177, c("flag_lake.out","buf.lake_id","lake_id") := list(0,"0","0")]

out[pointid == 414904, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"935101","935101")]
out[pointid == 413745, c("flag_lake.out","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 6184 & buf.lake_id == "935101", c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]

out[pointid == 476060, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"946145","946145")]
out[pointid == 476137, c("flag_lake.out","buf.lake_id","lake_id") := list(0,"0","0")]

out[pointid == 505473, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"98382","98382")] # new outlet upstream
out[pointid == 505116, "flag_lake.out" := 0] # old outlet upstream lake
out[pointid == 505642, c("flag_lake.in", "flag_lake.trib") := list(1,1)] # new inlet downstream lake
out[pointid == 505371, c("flag_lake.trib", "buf.lake_id","lake_id") := list(0,"0","0")]  # old tribs of downstream lake
out[pointid == 505511, c("flag_lake.trib", "buf.lake_id","lake_id") := list(0,"0","0")]  # old tribs of downstream lake

out[pointid == 521263, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"98780","98780")]
out[pointid == 521064, c("flag_lake.out","buf.lake_id","lake_id") := list(0,"0","0")]

out[pointid == 524663, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"98837","98837")]
out[pointid == 524590, c("flag_lake.out","buf.lake_id","lake_id") := list(0,"0","0")]

out[pointid == 4260, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"90555","90555")]
out[pointid == 4210, flag_lake.out := 0]
out[pointid == 19785, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"889690","889690")]
out[pointid == 19741, flag_lake.out := 0]
out[pointid == 56071, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"893830","893830")]
out[pointid == 56072, flag_lake.out := 0]
out[pointid == 376078, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"931079","931079")]
out[pointid == 375954, flag_lake.out := 0]
out[pointid == 378323, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"95947","95947")]
out[pointid == 378249, flag_lake.out := 0]
out[pointid == 508973, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"952464","952464")]
out[pointid == 508935, flag_lake.out := 0]
out[pointid == 526769, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"98878","98878")]
out[pointid == 526733, flag_lake.out := 0]

# lakes with two outlets
out[reach_id == 853, c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 1646, c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 1929, c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 2059, c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 2172, c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 3650, c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 3120, c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 3152, c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 2447 & buf.lake_id == "94064", c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 2446 & buf.lake_id == "915281", c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 2566, c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 3392 & buf.lake_id == "917326", c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 7457 & buf.lake_id == "947107", c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 7446 & buf.lake_id == "947107", c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 7505 & buf.lake_id == "947881", c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 3727 & buf.lake_id == "912416", c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 7520 & buf.lake_id == "RO2", c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 7605 & buf.lake_id == "RO2", c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 7607 & buf.lake_id == "RO2", c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 7610 & buf.lake_id == "RO2", c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]

# wrongly identified as tributary
out[pointid == 203494, c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]
out[reach_id == 6588, c("flag_lake.trib","buf.lake_id","lake_id") := list(0,"0","0")]

# Next pixel wrong
out[pointid == 138795, next.pixel := 138796]
out[pointid == 282185, next.pixel := 282088]

# misclassified as not lake
out[pointid == 433764, buf.lake_id := "937654"]

#--- change me depending on the year ---#
# Lake corrections depending on year -------------------------------------------------------------------------
# in all years RO1 exists (2016-2018)
# out[pointid == 530898, c("flag_lake.out","buf.lake_id","lake_id") := list(1,"RO1","RO1")]
# out[pointid == 530680, flag_lake.out := 0]
# out[pointid == 530706, flag_lake.out := 0]
#--- change me depending on the year ---#

# Clean lakes that don't have an inlet = Headwater lakes --------------------------------------------------------
# and any lakes that have a headwater stream submerged
lakes.with.source <- unique(out[flag_source == 1 & buf.lake_id != "0",]$buf.lake_id)
# all lakes with one source stream, do not have a tributary/inlet = identify first point as inlet
for(i in 1:length(lakes.with.source)){
  reaches <- unique(out[buf.lake_id == lakes.with.source[i] & flag_source == 1,]$reach_id) # & flag_remove == 0
  # number of other reaches in lake
  other.reaches <- unique(out[buf.lake_id == lakes.with.source[i],]$reach_id)
  
  temp <- out[reach_id %in% reaches,]
  
  if(length(other.reaches) > length(reaches)){
    # if there are several other tributaries that flow into lake except headwater streams
    # flag headwater streams as no count
    no.count <- temp[reach_id %in% reaches,]
    out[buf.lake_id == lakes.with.source[i] & reach_id %in% no.count & flag_source == 1,
        flag_source := 0]
    out[buf.lake_id == lakes.with.source[i] & reach_id %in% no.count,
        flag_skip := 1]
  } else if(length(reaches) > 1L & all(reaches %in% other.reaches)) {
    # if there are multiple headwater streams submerged in lake
    # pick the longest headwater stream = number of pixels
    temp <- out[buf.lake_id == lakes.with.source[i] & reach_id %in% reaches,]
    temp <- temp[, .(length.px = .N), by = .(reach_id)] # calculate how many pixels in headwater stream
    # flag the main headwater stream as such
    main.hw <- temp[which.max(length.px),]$reach_id
    out[buf.lake_id == lakes.with.source[i] & reach_id == main.hw, flag_main.chan := 1]
    out[buf.lake_id == lakes.with.source[i] & reach_id == main.hw & flag_source == 1,
        flag_lake.in := 1]
    # flag other streams as no count
    no.count <- temp[!(reach_id %in% main.hw),]
    out[buf.lake_id == lakes.with.source[i] & reach_id %in% no.count & flag_source == 1,
        flag_source := 0]
    out[buf.lake_id == lakes.with.source[i] & reach_id %in% no.count,
        flag_skip := 1]
  } else if(length(reaches) == 1L & reaches == other.reaches){
    # if there is a single stream in the lake do...
    out[buf.lake_id == lakes.with.source[i] & reach_id == reaches, flag_main.chan := 1] # define main channel in lake
    out[buf.lake_id == lakes.with.source[i] & flag_source == 1, flag_lake.in := 1] # define stream source as lake inlet
  } else {
    print(paste("Iteration",i,"has no condition match", sep = " "))
  }
}

# get number of headwater lakes
length(unique(out[flag_lake.in == 1,]$buf.lake_id))
# all 407
# number of lakes that had a source stream
length(unique(out[flag_skip == 1,]$buf.lake_id)) #12

# Identify tributary confluences to lakes -----------------------------------------------------------
# we have 1469 lakes in the watershed, of which, 407 are lakes that do not have tributaries
# give me all lakes that don't have a tributary assigned yet = non headwater lakes
lakes.trib <- unique(out[flag_lake.in == 1,]$buf.lake_id)
length(lakes.trib) # lakes with trib identified = 407
no.trib <- unique(out[buf.lake_id != "0" ,][!(buf.lake_id %in% lakes.trib),]$buf.lake_id)
length(no.trib) #1058

#setup parallel backend to use many processors
cores <- detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer # outfile=""
#registerDoParallel(cl)
registerDoSNOW(cl)
# prep progressbar
pb <- txtProgressBar(max = length(no.trib), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

#length(no.trib)
# for all lakes that do not have a tributary yet, do...
tribs <-foreach(i = 1:length(no.trib), .combine=c, .options.snow = opts) %dopar% {
  # get the reaches that flow into lake
  reaches <- unique(out[out$buf.lake_id == no.trib[i] & is.na(out$flag_skip),]$reach_id)
  
  tribs <- c()
  for(j in 1:length(reaches)){ # for all reaches in lake...
    reach.path <- out[out$reach_id %in% reaches[j],] # get the flow path of the reach
    # keep reaches that have switches in water bodies (stream -> lake, lake -> lake)
    if(length(unique(reach.path$buf.lake_id)) > 1){
      k <- reach.path[which.min(reach.path$fllength_dec),]$pointid # get upper most point on reach
      
      # if the upper point is not a lake
      if(reach.path[reach.path$pointid == k, ]$buf.lake_id != no.trib[i]){
        # search until it hits the lake of interest
        while(reach.path[reach.path$pointid == k,]$buf.lake_id != no.trib[i]){
          if(reach.path[reach.path$pointid == k,]$flag_mouth != 1){
            k <- reach.path[reach.path$pointid == k,]$next.pixel
          } else {
            # if it hits the mouth, stop the loop and save an NA
            k<-NA
            break
          }}
        tribs <- c(tribs,k)
      } else if(reach.path[reach.path$pointid == k, ]$buf.lake_id == no.trib[i] &
                reach.path[reach.path$pointid == k, ]$flag_conf != 1){
        # if it already hit the lake, and it's not a stream confluence
        tribs <- c(tribs,k)
      }
    }

  }

  return(tribs)
}
close(pb)
stopCluster(cl)

out[pointid %in% tribs, flag_lake.trib := 1]

# correct a few lakes that have not true tribs because of their shapes
not.tribs <- c(28739,31885,32424,36550,
               1185, 87353, 172280, 172438, 307313, 317352,
               361621, 392455, 430489, 480761, 514599, 508973)
# 29412 should be inlet
out[pointid %in% not.tribs, flag_lake.trib := 0]

# for each tributary, fill the pixels until outlet as lake -------------------------------------------------
#rm.reaches <- readRDS("./Objects/remove.reaches.rds")
#out[reach_id %in% rm.reaches, flag_remove := 1]
#out[is.na(flag_remove), flag_remove := 0]

# get all tributaries
out[buf.lake_id != "0" & is.na(flag_lake.trib), flag_lake.trib := 0]
out[buf.lake_id != "0" & is.na(flag_skip), flag_skip := 0]

tribs <- unique(out[flag_lake.trib == 1 & flag_skip != 1,]$pointid)

# save and submit to cluster
# saveRDS(tribs, "./Data/Traveltime/tribs_2016.rds")
# saveRDS(out, "./Data/Traveltime/out_2016.rds")

# until 500 ok
#setup parallel backend to use many processors
cores <- detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer # outfile=""
#registerDoParallel(cl)
registerDoSNOW(cl)
# prep progressbar
pb <- txtProgressBar(max = length(tribs), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# 2112
# time.out <- c()
# fill.paths <- data.frame()
#length(no.trib)
# for all lakes that do not have a tributary yet, do...
fill.tribs <-foreach(i = 1:length(tribs),
                     .combine='rbind', .multicombine=TRUE,
                     #.init=list(list(), list()), # return a list with two bins
                     .options.snow = opts) %dopar% {
# fill.tribs <- list()
# for(i in 1:length(tribs)){
#   print(paste0("Working on...", i))
  require(R.utils)
  # get the outlet of lake
  lake <- out[out$pointid == tribs[i],]$buf.lake_id
  outlet <- out[out$buf.lake_id == lake & out$flag_lake.out == 1,]$pointid
  # often the actual outlet is the next pixel after the outlet
  next.out <- out[out$pointid == outlet,]$next.pixel 
  
  #print(i)
  pixels <- c() # vector to fill with flow path
  k <- tribs[i] # starting point
  
  # try catching the lakes that do not work
  # tryCatch(
  #   expr = {withTimeout({while(!(k %in% c(outlet,next.out))){
  #     # get next pixel in flow path
  #     k <- out[out$pointid == k,]$next.pixel
  #     if(out[out$pointid == k, ]$flag_mouth == 1){
  #       # if we reach the watershed mouth, stop the iteration
  #       break
  #     }
  #     # combine the flow path into a vector during while loop
  #     pixels <- c(pixels, k)
  #   }}, timeout = 20)}, # stop while loop if it takes longer than 20s
  #   TimeoutException = function(ex) cat("Timeout. Skipping.\n"))
  
  while(!(k %in% c(outlet,next.out))){
    print(paste(k))
        # get next pixel in flow path
        k <- out[out$pointid == k,]$next.pixel
        
        if(out[out$pointid == k,]$fllength_dec > out[out$pointid == outlet,]$fllength_dec |
           out[out$pointid == k,]$fllength_dec > out[out$pointid == next.out,]$fllength_dec){
          # if we surpass the outlet, stop the iteration
          print("Outlet not found")
          break
        } else if(out[out$pointid == k, ]$flag_mouth == 1){
          # if we reach the watershed mouth, stop the iteration
          break
        }
        # combine the flow path into a vector during while loop
        pixels <- c(pixels, k)
      }
  
  if(is.null(pixels) == T){
    pixels <- k
  }
  # if the pixel where the loop finished is neither outlet nor the pixel after the outlet, save iteration
  # if(!(k %in% c(outlet,next.out)) == T){
  #   time.out <- c(time.out, i)
  # }
  #if the pixel where the loop finished is not the outlet but the pixel after outlet
  if(k == next.out){
    # update the new outlet info
    out[out$pointid == next.out,"flag_lake.out"] <- 1
    out[out$pointid == next.out,"buf.lake_id"] <- lake
    out[out$pointid == next.out,"lake_id"] <- lake

    out[out$buf.lake_id == lake & out$pointid == outlet, "flag_lake.out"] <- 0
  }
  
  # save the flow path as a dataframe to merge back to out later 
  if(!(k %in% c(outlet,next.out)) == T){ # if it timed out,
  df <- data.frame(buf.lake_id = NA,
                     pointid = NA, 
                     outlet = NA) 
  time.out <- data.frame(buf.lake_id = lake,
                         trib = tribs[i])
  } else {
    df <- data.frame(buf.lake_id = lake, # processed lake
                     pointid = pixels, # flow path from trib to outlet as pointids
                     outlet = k) # final outlet pointid
    time.out <- data.frame(buf.lake_id = NA,
               trib = NA)
  }
  return(list(df, time.out))
  #fill.tribs[[i]] <- list(df, time.out)
}
close(pb)
stopCluster(cl)

saveRDS(fill.tribs, "./Objects/tributary_fill_2015.rds")

#fill.tribs <- readRDS("./Objects/tributary_fill_2017.rds")
# separate the two different data frames
flow.paths <- fill.tribs[1:length(tribs)]
time.out <- fill.tribs[(length(tribs)+1):length(fill.tribs)]

# bind list into data frame
flow.paths <- bind_rows(flow.paths)
time.out <- bind_rows(time.out)

# clean time out lakes
time.out <- time.out[!is.na(time.out$buf.lake_id),]

flow.paths <- flow.paths %>% distinct()
dups <- flow.paths[duplicated(flow.paths$pointid),]$pointid
#View(flow.paths[flow.paths$pointid %in% dups,])
setDT(flow.paths)
# correct
flow.paths[buf.lake_id == 93524 & pointid %in% dups, c("buf.lake_id", "outlet") := list(910069, 172570)]
flow.paths[buf.lake_id == 97279 & pointid %in% dups,  c("buf.lake_id", "outlet") := list(941620, 453940)]
flow.paths[buf.lake_id == 93399, outlet := 158249]

out[pointid == 453981, flag_lake.trib := 0]
out[pointid == 454036, flag_lake.trib := 1]
out[pointid == 158148, flag_lake.out := 0]
out[pointid == 158249, c("buf.lake_id","flag_lake.out") := list(93399, 1)]

flow.paths <- flow.paths %>% distinct()
saveRDS(flow.paths, "./Objects/lake_flow.paths_2015.rds")
saveRDS(out, "./Objects/streamnet_out_2015.rds")

# out <- readRDS("./Objects/streamnet_out_2017.rds")
# flow.paths <- readRDS( "./Objects/lake_flow.paths_2017.rds")
# RO1 outlet
#out[pointid == 530898, c("flag_lake.out","buf.lake_id","lake_id") := list(0,"0","0")]

# make sure that the flow paths within lakes are labelled as such
out[flow.paths, c("buf.lake_id","outlet_pointid") := list(i.buf.lake_id, i.outlet), on = .(pointid)]
# make sure outlet is correct
t <- out[flag_lake.out == 1,][!is.na(outlet_pointid),]
t[pointid != outlet_pointid,]
# one lake below RO2

# clean overwritten pixels by flow paths
out[lake_id != buf.lake_id, lake_id := buf.lake_id]

# sanity check
lakes.trib <- unique(out[out$flag_lake.trib == 1,]$buf.lake_id) #1020
lakes.in <- unique(out[out$flag_lake.in == 1,]$buf.lake_id) # 407
length(unique(out$buf.lake_id)) - (length(lakes.trib)+length(lakes.in)) # lakes with trib identified = 1101, 39 missing
no.trib <- unique(out[buf.lake_id != "0" ,][!(buf.lake_id %in% c(lakes.trib,lakes.in)),]$buf.lake_id)

# for all lakes, identify the biggest tributary as inlet
# clean reservoirs
#out[pointid == 429846, flag_lake.in := 1] # RO2
#out[pointid == 495826, flag_lake.in := 1] # RO1
#out[pointid == 368881, flag_lake.in := 1] # RO3

# prep progressbar
with.in <- unique(out[flag_lake.in == 1,]$buf.lake_id)
no.in <- unique(out[buf.lake_id != "0",]$buf.lake_id)[!(unique(out[buf.lake_id != "0",]$buf.lake_id) %in% with.in)]
# these are actually lakes with two "source" streams, should be ok for the next step

inlets <- c()
for(i in 1:length(no.in)){
  temp <- out[out$buf.lake_id == no.in[i],]
  
  if(nrow(temp[!is.na(temp$flag_lake.in),]) == 0L){
    sub <- temp[temp$flag_lake.trib == 1,]
    if(nrow(sub) == 0L){ # multiple entirely covered streams
      sub <- temp[temp$flag_source == 1,]
    }
    inlets <- c(inlets,sub[which.max(sub$fllength_dec),]$pointid)
  }
}

out[pointid %in% inlets, flag_lake.in := 1]

# Are there as many outlets as inlets?
outlets <- unique(out[flag_lake.out == 1,]$pointid)
inlets <- unique(out[flag_lake.in == 1, ]$pointid)
length(outlets) == length(inlets) # should be true
all.wb <- unique(out[buf.lake_id != "0",]$buf.lake_id)
length(all.wb)

# sanity check, all outlets identified in loop are not duplicated
outlets <- out[buf.lake_id != "0", .(outlet = unique(outlet_pointid)), by = .(buf.lake_id)]
outlets <- outlets[!is.na(outlet),]
nrow(outlets) # 1020
nrow(outlets[duplicated(buf.lake_id)]) #0 none

# sanity check, 
out[pointid == 476060, flag_lake.out := 0]
# for each lake, get the outlet
outlet.init <- out[flag_lake.out == 1, .(out.init = pointid), by = .(buf.lake_id)]
# are there any outlets in the final dataset that are misclassified?
temp <- merge(outlet.init, outlets, by = "buf.lake_id", all = T)
temp[!is.na(outlet),][out.init != outlet,] # none
# save
out[temp, outlet_pointid := i.out.init, on = .(buf.lake_id)]

# length of outlets and inlets are same?
outlets <- unique(out[flag_lake.out == 1,]$pointid)
inlets <- unique(out[flag_lake.in == 1, ]$pointid)
length(outlets) == length(inlets) # should be true
all.wb <- unique(out[buf.lake_id != "0",]$buf.lake_id)
length(all.wb)
# all good

# Do all water bodies have an inlet and outlet?
all.wb[!(all.wb %in% out[pointid %in% outlets,]$buf.lake_id)]
all.wb[!(all.wb %in% out[pointid %in% inlets,]$buf.lake_id)]

# for all inlets, find points until outlet = main channel
ins <- unique(out[flag_lake.in == 1,]$pointid)

cl <- makeCluster(cores[1]-1) #not to overload your computer # outfile=""
registerDoSNOW(cl)

pb <- txtProgressBar(max = length(ins), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

t <- foreach(i = 1:length(ins), .combine=rbind, .options.snow = opts) %dopar% {
 #for(i in 1:length(ins)){

  # get identified inlet
  lake <- out[out$pointid %in% ins[i],]$buf.lake_id
  # get outlet
  outlet <- out[out$buf.lake_id == lake & out$flag_lake.out == 1,]$pointid
  k <- out[out$pointid == ins[i],]$pointid
  #print(paste("I'm at ...",i))
  if(ins[i] != outlet){
    path <- c()
    while(k != outlet){
      k <- out[out$pointid == k,]$next.pixel
      path <- c(path, k)
    }
   df <- data.frame(buf.lake_id = lake, path = path, n = length(path))
  } else {
    df <- data.frame(buf.lake_id = lake, path = ins[i], n = 1)
  }
  #t <- rbind(t, df)
 return(df)
}
close(pb)
stopCluster(cl)

# merge, make sure that the main channel path has the same ID as the corresponding lake
out[t, c("buf.lake_id", "flag_main.chan", "main.chan.n") := list(i.buf.lake_id, 1, i.n), on = c('pointid==path')]

# Skip pixels of tributary of lakes that are within the lake ------------------------------------------------
# for all tributaries, find confluence pixel to main channel
out[is.na(flag_lake.in), flag_lake.in := 0]
tribs <- unique(out[flag_lake.trib == 1 & flag_lake.in != 1,]$pointid)

#t <- foreach(i = 1:length(tribs), .combine=rbind, .options.snow = opts) %dopar% {
# trib.df <- data.frame()
for(i in 1:length(tribs)){
  #print(paste("Working on...",i))
  # get identified inlet
  lake <- out[out$pointid %in% tribs[i],]$buf.lake_id
  # get main channel
  outlet <- out[out$buf.lake_id == lake & out$flag_lake.out == 1,]$pointid
  main.chan <- out[out$buf.lake_id == lake & !is.na(out$main.chan.n),]$pointid
  k <- out[out$pointid == tribs[i],]$pointid
  
  path <- c() # track the flow path until main channel
  while(!(any(k %in% main.chan) | k == outlet)){
    #print(k)
      k <- out[out$pointid == k,]$next.pixel
      path <- c(path, k)
  }
  if(!is.null(path)){
    if(out[out$pointid == tail(path,1),]$flag_conf != 1){
      stop(paste0("Iteration: ",i," - Last point in path not a confluence"))
    } else {
      out[out$pointid == tribs[i], next.pixel := tail(path,1)]
    }
    }
}

#out[t, c("buf.lake_id", "flag_skip") := list(i.buf.lake_id, 1), on = c('pointid==path')]
# all main channel pixels should not be skipped
#out[flag_main.chan == 1 & flag_skip == 1, flag_skip := 0]
out[is.na(flag_skip), flag_skip := 0]

saveRDS(out, "./Objects/cleaned_streamnet_2015.rds")

#--- Until here, run by year ---#

#--- From here, we process all years at the same time ---#

# sanity check
st15 <- readRDS("./Objects/cleaned_streamnet_2015.rds")
st16 <- readRDS("./Objects/cleaned_streamnet_2016.rds")
st17 <- readRDS("./Objects/cleaned_streamnet_2017.rds")

# merge
temp <- merge(st15 %>% dplyr::select(pointid, next.pixel, buf.lake_id),
      st16 %>% dplyr::select(pointid, next.pixel, buf.lake_id), by = "pointid")

View(temp[next.pixel.x != next.pixel.y,]) # only RO1

# merge
temp <- merge(st15 %>% dplyr::select(pointid, next.pixel, buf.lake_id),
              st17 %>% dplyr::select(pointid, next.pixel, buf.lake_id), by = "pointid")

View(temp[next.pixel.x != next.pixel.y,]) # only RO1, RO3

rm(st15,st16,st17, temp)

for(y in 2015:2018){
  for(mon in c("6","8","10")){
    
    if(y >= 2017){
      out <- readRDS("./Objects/cleaned_streamnet_2017.rds")
    } else {
      out <- readRDS(paste0("./Objects/cleaned_streamnet_",y,".rds"))
    }

# Add the stream data
# merge downstream pixel with full data set ------------------------------------------------------------------------------------
# merge info with streamnet pnts file for seasonal wrt
streamnet_pnts <- readRDS("./Data/Traveltime/flacc3000_streamnet_pnts_cleaned_05.02.2024.rds")
# remove old flags
streamnet_pnts[, c("flag_source", "flag_conf",
                   "flag_end", "flag_mouth") := c(NULL, NULL, NULL, NULL)]

# Filter and merge
streamnet_pnts <- streamnet_pnts %>% filter(Year == y & Month == as.numeric(mon))

# make inlets main channels
out[flag_lake.in == 1, flag_main.chan := 1]
out[flag_lake.out == 1, flag_main.chan := 1]

# Get lake WRT estimates
hydrolakes <- readxl::read_excel("./Data/Traveltime/hydrolakes_in_laromaine.xlsx") %>% setDT()
hydrolakes[, Hylak_id := as.character(Hylak_id)]
# Convert to m2
hydrolakes[, wrt_s := Res_time * 86400]
#lakes[, wrt_m := wrt_d * 1440]

# extract the main channel of lakes and assign their WRT estimates to those pixels
main.chan <- out[buf.lake_id != "0", .(main.chan.n = unique(main.chan.n)), by = .(buf.lake_id)]
main.chan <- main.chan[!is.na(main.chan.n),]
main.chan[hydrolakes, wrt_s := i.wrt_s, on = c("buf.lake_id==Hylak_id")]

# calculate wrt in pixel by second and get velocity
main.chan[, wrt.px18_s := wrt_s / main.chan.n]
main.chan[, veloc.px18_ms := 18 / wrt.px18_s]

# calculate residence time in minutes
main.chan[, wrt.px18_min := wrt.px18_s / 60]
main.chan[, flag_main.chan := 1]

# merge into stream network
out <- out[main.chan, c("veloc_px_ms", "wrt_px_min") := list(veloc.px18_ms, i.wrt.px18_min), on = .(buf.lake_id, flag_main.chan)]

# Add reservoir data
res.wrt <- read.csv("./Output/reservoir_wrt_seasonal_perpixel_2024-02-01.csv", sep = ",", stringsAsFactors = F) %>% setDT()

# will just try one year, one month first
temp <- res.wrt %>% dplyr::select(Reservoir:month,wrt_d)
temp[, wrt_s := wrt_d * 86400]

main.chan <- out[grep("RO", buf.lake_id), .(main.chan.n = unique(main.chan.n)), by = .(buf.lake_id)]
main.chan <- main.chan[!is.na(main.chan.n),]

main.chan[temp %>% filter(year == y & month == as.numeric(mon)), wrt_s := i.wrt_s, on = c("buf.lake_id==Reservoir")]

# calculate wrt in pixel by second and get velocity
main.chan[, wrt.px18_s := wrt_s / main.chan.n]
main.chan[, veloc.px18_ms := 18 / wrt.px18_s]

# calculate residence time in minutes
main.chan[, wrt.px18_min := wrt.px18_s / 60]
main.chan[, flag_main.chan := 1]

out <- out[main.chan, c("veloc_px_ms", "wrt_px_min") := list(veloc.px18_ms, i.wrt.px18_min), on = .(buf.lake_id, flag_main.chan)]

# clean out file
out[is.na(flag_lake.in), flag_lake.in := 0]
out[is.na(flag_lake.trib), flag_lake.trib := 0]

out[!is.na(main.chan.n) & buf.lake_id != "0", flag_main.chan := 1]
out[is.na(flag_main.chan), flag_main.chan := 0]

# do all standing water bodies' main channel have a wrt assgined?
out[buf.lake_id != "0" & flag_main.chan == 1 & is.na(wrt_px_min),] # yes

# add stream and river residence times
temp <- streamnet_pnts[out, , on = .(pointid)]

# define reservoirs and lakes
temp[grep("RO",buf.lake_id), water.body := "Reservoir"]
temp[buf.lake_id != "0" & is.na(water.body), water.body := "Lake"]
temp[is.na(water.body) & buf.lake_id == "0", water.body := "Fluvial"]
temp[is.na(water.body),]

# re-arrange
temp <- temp %>% dplyr::select(Year:flow.condition, 
                               pointid, next.pixel, reach_id,
                               water.body, 
                               buf.lake_id,
                               flacc_km2, fllength_dec, strah_ord,
                               velocity_ms, veloc_px_ms,
                               dis_m3s, wrt_min, wrt_px_min,
                               flag_source, flag_conf, flag_mouth,
                               flag_lake.in, flag_lake.trib, flag_lake.out, flag_main.chan,
                               coord_x, coord_y, reach_id)

# correct a few
temp[pointid %in% c(220267, 430489, 430490,453982,454004,480653,480692), 
     water.body := "Fluvial"]

# merge WRT columns
temp[water.body == "Fluvial", wrt_min.wb := wrt_min]
temp[water.body == "Lake" | water.body == "Reservoir", wrt_min.wb := wrt_px_min]

# merge velocity columns
temp[water.body == "Fluvial", vel_ms.wb := velocity_ms]
temp[water.body == "Lake" | water.body == "Reservoir", vel_ms.wb := veloc_px_ms]

temp[reach_id == 2460, c("flag_main.chan", "wrt_min.wb") := list(0,NA)]
# Reservoir confluences that don't have a WRT assigned, fill with fluvial WRT
temp[flag_lake.trib == 1 & is.na(wrt_min.wb), wrt_min.wb := wrt_min]

# sanity check, all important pixels and main channel should have WRT estimates
temp[flag_main.chan == 1 & is.na(wrt_min.wb),] # none, good!
temp[flag_lake.trib == 1 & is.na(wrt_min.wb),]
temp[flag_lake.in == 1 & is.na(wrt_min.wb),]
temp[flag_lake.out == 1 & is.na(wrt_min.wb),]

# one source indication was missing
temp[pointid == 414906, flag_source := 1]

# Look at distribution
# ggplot(test, aes(x = water.body, y = log(wrt_min.wb))) +
#   geom_boxplot()

# re-arrange and save

final <- temp %>% dplyr::select(Year:dis_m3s, fluv.wrt_min = wrt_min,
                                lent.wrt_min = wrt_px_min, vel_ms.wb, wrt_min.wb,
                                flag_source:coord_y)

saveRDS(final, paste0("./Data/Traveltime/flacc3000_streamnet_wWRT_",y,"-",mon,"_", Sys.Date(),".rds"))
# write.table(final, paste0("./Data/Traveltime/flacc3000_streamnet_wWRT_",y,"-",mon,"_", Sys.Date(),".csv"),
#           sep = ",", dec = ".", row.names = F)
rm(test, final, out, streamnet_pnts, temp, res.wrt, main.chan, hydrolakes)
  }}

# # Extract water residence time for FT samples -----------------------------------------------------------------------
# main <- read.csv("./Data/Summary/mainsamples_flacc3000_snap.txt", sep = ";", dec = ",",
#                  stringsAsFactors = F) %>% setDT()
# rest <- read.csv("./Data/Summary/restsamples_flacc3000_snap.txt", sep = ";", dec = ",",
#                  stringsAsFactors = F) %>% setDT()
# 
# wrt <- read.csv("./Data/Traveltime/streamnet_pnts_complete.wrt.csv", 
#                 sep =",", dec = ".", stringsAsFactors = F) %>% setDT()
# 
# # subset by month and year for GIS
# for(i in 2015:2016){
#   for(j in c(6,8)){
#     out <- wrt[Year == i & Month == j,]
#     write.table(out, paste0("./Data/Traveltime/streamnet_pnts_wrt_", i,"-0",j,".csv"),
#                 sep = ",", dec = ".", row.names = F)
#   }
# }
# 
# # ggplot(test, aes(x = water.body, y = wrt_min.wb)) +
# #   geom_boxplot() +
# #   facet_grid(Year~Month)
# 
# wrt[water.body == "RO1",]
# # we extract the water residence time for each sample site
# # some sample sites had corrupt coordinates
# 
# samples <- rbind(main, rest) %>% setDT()
# samples[, coord := paste(NEAR_X, NEAR_Y, sep = " ")]
# wrt[, coord := paste(coord_x, coord_y, sep = " ")]
# 
# samples[wrt, c("velocity_ms", "wrt_min") := list(i.velocity_ms.wb, i.wrt_min.wb), on = .(coord)]
# 
# # raster is in EPSG:4617, which is NAD83 (CSRS)
# # +proj=longlat +ellps=GRS80 +no_defs
# 
# # we will reproject the rater layer into EPSG:2947, which is NAD83 (CSRS) / MTM zone 5
# # Proj4: +proj=tmerc +lat_0=0 +lon_0=-64.5 +k=0.9999 +x_0=304800 +y_0=0 +ellps=GRS80 +units=m +no_defs
# 
# # projection method depends on raster values:
# # - categorical = "ngb", short for nearest neighbour method
# # - numerical / integer = "bilinear"
# 
# library(raster)
# epsg4617 <- "+proj=longlat +ellps=GRS80 +no_defs"
# #mtm5 <- "+proj=tmerc +lat_0=0 +lon_0=-64.5 +k=0.9999 +x_0=304800 +y_0=0 +ellps=GRS80 +units=m +no_defs"
# 
# raster.wrt <- rasterFromXYZ(wrt %>% dplyr::select(x = coord_x, y = coord_y, z = wrt_min.wb),
#                             crs = epsg5617)
# 
# #cdem_mtm5 <- projectRaster(cdem, crs = mtm5, method = "bilinear")
# #crs(cdem_mtm5)
# 
# # 0.00020833333
# 
# #writeRaster(cdem_mtm5, "./GIS/cdem_fill_romain_MTM5.tif", options = c("TFW=YES"))
# # include options = c("TFW=YES") to create auxiliary files needed to read the raster in ArcMap
# 
# # Done!

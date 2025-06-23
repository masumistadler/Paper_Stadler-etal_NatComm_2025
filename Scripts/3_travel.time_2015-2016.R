# Calculate travel time for 2015-06 ---------------------------------------------------------------------------

# R set-up ----------------------------------------------------------------------------------------------------
### Packages -------------------------------------------------------------------------------
.libPaths( c( .libPaths(), "/home/mstadler/projects/def-pauldel/R/x86_64-pc-linux-gnu-library/4.2") )
pckgs <- list("data.table", "tidyverse", # wrangling
              "plyr",
              "doParallel", "doMC","foreach", "doSNOW")
### Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)
# Many packages have to be installed through Bioconductor, please refer to the package websites

### Load
invisible(lapply(pckgs, library, character.only = T))
rm(pckgs, miss.pckgs)

for(y in 2015:2016) {
  for (mon in c("6", "8")){
    # Read in data
    
    # Read in old stream network
    # old <- readRDS(paste0("./Data/Traveltime/flacc3000_streamnet_wWRT_",y,"-0",mon,".rds"))
    # streamnet <-
    #    readRDS(paste0("./Data/Traveltime/flacc3000_streamnet_wWRT_", y, "-", mon, "_2024-02-06.rds"))
    
    old <- readRDS(paste0("./Objects/flacc3000_streamnet_wWRT_",y,"-0",mon,".rds"))
    streamnet <-
      readRDS(paste0("./Objects/flacc3000_streamnet_wWRT_", y, "-", mon, "_2024-02-09.rds"))
    
    streamnet <- merge(old %>% dplyr::select(Year:reach_id, buf.lake_id:strah_ord, flag_source:coord_y),
    streamnet %>% dplyr::select(pointid, water.body, velocity_ms:wrt_min.wb), by = "pointid")

    m <- streamnet
    m <- m[!duplicated(m$pointid), ]
    m <- m[water.body == "Fluvial", flag_main.chan := 1]
    #old <-old[water.body == "Stream" | water.body == "River", flag_main.chan := 1]
    
    #cores <- detectCores()
    cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
    #cl <- makeCluster(cores[1]-1) #not to overload your computer # outfile=""
    #registerDoSNOW(cl)
    registerDoMC(cores = cores)
    #m <- streamnet[streamnet$Year == 2015 & streamnet$Month == 6,] %>% setDT()
    #m <- read.csv("./Data/Traveltime/cumtime_6-2015.csv", sep = ",", stringsAsFactors = F)
    
    #zero <- m[m$travel.time == 0,]
    
    rm(old, streamnet)
    
    # Recode Yves' loop in JMP to R -----------------------------------------------------------------------------------------
    
    # Our data is structured as follows
    # next.pixel = indicates which pointid is the next downstream pixel
    # flag_source = indicates whether current pixel is source
    # flag_conf = indicates whether pixel is a confluence
    # flag_end = indicates whether pixel is just before a confluence
    # flag_mouth = indicates whether pixel is the mouth of the river
    # flag_lake.in = indicates whether pixel is the main inlet of a lake/reservoir
    # flag_lake.trib = indicates whether pixel is entering a lake/reservoir
    # flag_lake.out = indicates whether pixel is a outlet of a lake/reservoir
    # velocity_ms.wb = velocity in metres/sec, considering water bodies (reservoir, riverine lakes)
    # wrt_min.wb = water residence time in pixel, considering water bodies (reservoir, riverine lakes)
    
    # for(yr in 1:length(unique(streamnet$Year))){ #
    #   for(mon in 1:length(unique(streamnet$Month))){ #
    
    #m <- streamnet[streamnet$Month ==  unique(streamnet$Month)[mon] & streamnet$Year == unique(streamnet$Year)[yr],]
    #View(m[pointid %in% dups,])
    # somehow duplicates got inside
    
    
    # First loop ---------------------------------------------------------------------------------
    # Track first stream orders down until a confluence
    # create empty data frame to fill in
    
    # cumtime <- data.frame(pointid = sort(m$pointid), travel.time = rep(0, nrow(m)))
    cumtime <- data.frame()
    #out<- data.frame()
    # extract all pointids that are sources
    sources <- m[m$flag_source == 1 & flag_main.chan == 1, ]$pointid
    
    #pb <- txtProgressBar(max = length(sources), style = 3)
    #progress <- function(n)
    #  setTxtProgressBar(pb, n)
    #opts <- list(progress = progress)
    #sources <- c(1, 1006,1423)
    #359, 1180, 1180
    #m <- setDF(m)
    temp <- foreach(i = 1:length(sources), .combine = rbind) %dopar% {
      #for(i in 1:3){
      if (m[m$pointid == sources[i], "flag_source"] == 1) {
        s <- 1
        k <- sources[i]
        cumtime <- data.frame()
        cumtime[s, "pointid"] <- k
        cumtime[s, "travel.time"] <- m[m$pointid == k, "wrt_min.wb"]
        nx.px <- m[m$pointid == k, ]$next.pixel

        # calculate travel time within first order stream
        while (m[m$pointid == nx.px, "flag_conf"] != 1) {
          # while next pixel is not a confluence, do...
          s <- s + 1
          # calculate cumulative travel time by...
          cumtime[s, "pointid"] <- m[m$pointid == k, ]$next.pixel
          cumtime[s, "travel.time"] <-
            cumtime[cumtime$pointid == k, ]$travel.time + # Time A
            (m[m$pointid == m[m$pointid == k, ]$next.pixel, ]$wrt_min.wb * # Time B
               # adding wrt of next pixel to the weighted wrt of the current pixel
               (1 - (m[m$pointid == m[m$pointid == k, ]$next.pixel, ]$dis_m3s - m[m$pointid == k, ]$dis_m3s) /
                  m[m$pointid == m[m$pointid == k, ]$next.pixel, ]$dis_m3s))
          # Time A + Time B * (1 - (disB - disA)/disB)
          # wrt of next pixel is weighted by discharge of next pixel minus discharge of current pixel
          # divided by discharge of next pixel; equation:
          # m[m$pointid == m[m$pointid == k, ]$next.pixel, ]$wrt_min.wb +
          #   (cumtime[cumtime$pointid == k, ]$travel.time *
          #      # adding wrt of next pixel to the weighted wrt of the current pixel
          #      (1 - (m[m$pointid == m[m$pointid == k, ]$next.pixel, ]$dis_m3s - m[m$pointid == k, ]$dis_m3s) /
          #         m[m$pointid == m[m$pointid == k, ]$next.pixel, ]$dis_m3s))
          # Time B + Time A * (1-(dis_B - dis_A) / dis_B)
          
          # once this is done, overwrite k with the next pixel
          k <-
            m[m$pointid == k, ]$next.pixel # so the loop jumps to the next pixel
          nx.px <- m[m$pointid == k,]$next.pixel
          #print(paste("I'm at pointid", k))
          #if(m[m$pointid == k, "flag_conf"] == 1){
          #cat(paste0("\rReached confluence of source stream. ", length(sources) - i, " streams left."))
          #}
        }
        cumtime <- cumtime[!is.na(cumtime$travel.time), ]
        # out<- rbind(out, cumtime)
        return(cumtime)
      }
    }
    #close(pb)
    #stopCluster(cl)
    
    temp <- temp %>% distinct()
    cumtime <- temp
    #saveRDS(cumtime, "./Data/Traveltime/temp_cumtime_2016-06.rds")
    saveRDS(cumtime,
            paste0("./Objects/temp_strahord1_", y, "-", mon, ".rds"))
    print("Finished stream orders 1")
    #cumtime <- readRDS("./Objects/temp_strahord1_2016-06.rds")
   
    #cumtime <- readRDS( "./Objects/temp_strahord1_2015-06.rds")
    #cumtime <- readRDS( "./Data/Traveltime/temp_strahord1_2015-06.rds")
    # merge with m
    temp <- m %>% dplyr::select(pointid)
    cumtime <- merge(temp, cumtime, all = TRUE)
    # remove all fills of confluences
    # cumtime[pointid %in% m[flag_conf == 1, ]$pointid, travel.time := NA]
    # cumtime <- cumtime %>% distinct()
    #t <- merge(m, cumtime, all = TRUE)
    #ggplot( t %>% filter(strah_ord == 1), aes(x = fllength_dec, y = travel.time)) +
    #  geom_point()
    
    # 2nd loop #
    # takes care of confluences and calculates the time in reaches after confluences
    max.order <- max(m$strah_ord)
   #postponed <- c()
    
    for (s in 2:max.order) {
      print(paste0("Starting loop at stream order ", s))
      #:max.order
      # get all the confluences of the given order
      # order by flow length
      m <- m[order(fllength_dec),]
      conf <- m[m$strah_ord == s & m$flag_conf == 1 & m$flag_main.chan == 1, ]$pointid
      # prep progress bar
      #pb <- txtProgressBar(max = length(conf), style = 3)
      #progress <- function(n) setTxtProgressBar(pb, n)
      #opts <- list(progress = progress)
      
      # extract those confluences where every confluence has a LOWER strahler order
      # = first level confluences
      next.df <- m[next.pixel %in% conf,]
      # for each confluence, pick only the confluences that have merging streams of a lower strahler order
      next.df[, max.ord := max(strah_ord), by = .(next.pixel)]
      # first level confluences = merges only strahler orders smaller
      first.level <- unique(next.df[max.ord < s,]$next.pixel)
      # second level confluences = merges strahler orders smaller and same
      sec.level <- unique(next.df[max.ord == s,]$next.pixel)
      
      print(paste0("Starting parallel loop for stream order ", s," and first level streams"))
      temp <-
        foreach(i = 1:length(first.level), .combine = rbind) %dopar% {
          #.options.snow = opts
          #for(i in 1:length(first.level)){
          if (m[m$pointid == first.level[i], ]$flag_conf == 1 &
              m[m$pointid == first.level[i], ]$strah_ord == s) {
            # if pointid is a confluence and matches the given stream order
            z <- 1
            k <- first.level[i] # get pointid
            temp.df <- data.frame()
            temp.df[z, "pointid"] <- k
            
            # get pointid of pixels flowing into confluence
            conf.df <- m[m$flag_mouth == 0 & m$next.pixel == k, ]
            conf.df <- conf.df[!is.na(wrt_min.wb),]
              
              # if the confluences are empty (= flooded by lake),
              # then track them down until next confluence
              
              # weigh travel time by discharge for pixels flowing into the confluence
              # and take the sum of all weighted travel time of pixels flowing into confluence
              left.points <-
                cumtime[cumtime$pointid %in% conf.df$pointid, ]$pointid
              # keep only confluences that have travel time
              left.points <-
                left.points[!(is.na(cumtime[cumtime$pointid %in% left.points, ]$travel.time))]
              
              # if WRT is available
              if(length(left.points) > 0L){
                ctime <-
                  sum(cumtime[cumtime$pointid %in% left.points, ][order(pointid),]$travel.time *
                       m[m$pointid %in% left.points, ][order(pointid),]$dis_m3s, na.rm = T)
                # take the sum of discharge, of all pixels flowing into confluence
                cda <- sum(m[m$pointid %in% left.points,]$dis_m3s)
                
                # assign new value to cumulative time in confluence pixel based on:
                # ctime = sum of discharge weighted travel time of pixels flowing into the confluence
                # cda = sum of discharge of pixels flowing into the confluence
                # divided by the residence time within the confluence pixel
                # nx.trav <- ctime / cda + m[m$pointid == k, ]$wrt_min.wb
                
                temp.df[z, "travel.time"] <- ctime / cda + m[m$pointid == k, ]$wrt_min.wb
                # if(!is.na(cumtime[cumtime$pointid == k, ]$travel.time)){
                #   # if there is already a value assigned, pick the bigger value
                #   travs <- c(cumtime[cumtime$pointid == k, ]$travel.time, nx.trav)
                #   cumtime[cumtime$pointid == k, ]$travel.time <- travs[which.max(travs)]
                # } else {
                #   cumtime[cumtime$pointid == k, ]$travel.time <-
                #     nx.trav
                # }
                # Once we have corrected the confluence travel time by the discharge of the merging streams,
                # we can trail the stream down until the next confluence
                nx.px <- m[m$pointid == k,]$next.pixel
                
                while (m[m$pointid == nx.px, ]$flag_conf != 1) {
                  # while next pixel is not a confluence, do...
                  # calculate cumulative travel time by...
                  if (m[m$pointid == k, ]$flag_mouth != 1) {
                    #print(paste0("I'm at...", k))
                    z <- z + 1
                    # calculate cumulative travel time by...
                    temp.df[z, "pointid"] <- nx.px
                    temp.df[z, "travel.time"] <- 
                      temp.df[temp.df$pointid == k, ]$travel.time + # Time A
                      (m[m$pointid == nx.px, ]$wrt_min.wb * # Time B
                         # adding wrt of next pixel to the weighted wrt of the current pixel
                         (1 - (m[m$pointid == nx.px, ]$dis_m3s - m[m$pointid == k, ]$dis_m3s) /
                            m[m$pointid == nx.px, ]$dis_m3s))
                    
                    #cumtime[cumtime$pointid == nx.px, ]$travel.time <-
                    
                    # wrt of next pixel is weighted by discharge of next pixel minus discharge of current pixel
                    # divided by discharge of next pixel; equation:
                    # Time B + Time A * (1-(dis_B - dis_A) / dis_B)
                    
                    # once this is done, overwrite k with the next pixel
                    k <-
                      m[m$pointid == k, ]$next.pixel # so the loop jumps to the next pixel
                    nx.px <- m[m$pointid == k,]$next.pixel
                  } else {
                    #if it hits the mouth, stop
                    break
                  }
                  # if(m[m$pointid == k, "flag_conf"] == 1){
                  #   cat(paste0("\rReached confluence of stream of order ", s,". ", length(conf) - i, " streams within order left."))
                  # }
                  
                }
                
                temp.df <- temp.df[!is.na(temp.df$travel.time), ]
                # out<- rbind(out, cumtime)
                return(temp.df)
                # path <- unique(path)
                # return(cumtime[cumtime$pointid %in% path, ])
              
                } else {
            # if no stream flowing into confluence has WRT, skip
            break
          } 
        }
      }
      # save results
      cumtime <-
        cumtime[temp, travel.time := i.travel.time, on = .(pointid)]
      
      # t <- merge(m, cumtime, all = TRUE)
      # ggplot(t %>% filter(strah_ord <= 2), aes(x = fllength_dec, y = log10(travel.time), colour = strah_ord)) +
      #   geom_point()
      
      print("Starting while loop to fill in all other secondary level rivers")  
      while(length(sec.level) > 0L){
      # get next.pixels that have travel time already
      prior <- m[m$next.pixel %in% sec.level,]
      prior <- prior[cumtime, travel.time := i.travel.time, on = .(pointid)][!is.na(wrt_min.wb),]
      prior <- prior[flag_main.chan == 1,]
      # calculate how many streams go in each  confluence
      prior[, n.conf := .N, by = .(next.pixel)]
      #any(prior$strah_ord > s)
      # calculate how many of those have already a travel time assigned
      prior[!is.na(travel.time), n.filled := .N, by = .(next.pixel)]
      prior[order(next.pixel),]
      prior <- prior[!is.na(n.filled),][, .(n.conf = unique(n.conf),
                                        n.filled = unique(n.filled)), by = .(next.pixel)]
      prior <- prior[n.filled == n.conf,]$next.pixel
      first <- sec.level[sec.level %in% prior]
      
      print(paste0("Starting parallel loop for stream orders ", s, " to fill in streams of second level"))
      temp <-
        foreach(i = 1:length(first), .combine = rbind) %dopar% {
          #.options.snow = opts
          #for(i in 1:length(first)){
          if (m[m$pointid == first[i], ]$flag_conf == 1 &
              m[m$pointid == first[i], ]$strah_ord == s) {
            # if pointid is a confluence and matches the given stream order
            z <- 1
            k <- first[i] # get pointid
            temp.df <- data.frame()
            temp.df[z, "pointid"] <- k
            
            # get pointid of pixels flowing into confluence
            conf.df <- m[m$flag_mouth == 0 & m$next.pixel == k, ]
            conf.df <- conf.df[!is.na(wrt_min.wb),]
            
            # if the confluences are empty (= flooded by lake),
            # then track them down until next confluence
            
            # weigh travel time by discharge for pixels flowing into the confluence
            # and take the sum of all weighted travel time of pixels flowing into confluence
            left.points <-
              cumtime[cumtime$pointid %in% conf.df$pointid, ]$pointid
            # keep only confluences that have travel time
            left.points <-
              left.points[!(is.na(cumtime[cumtime$pointid %in% left.points, ]$travel.time))]
            
            # if WRT is available
            if(length(left.points) > 0L){
              ctime <-
                sum(cumtime[cumtime$pointid %in% left.points, ][order(pointid),]$travel.time *
                      m[m$pointid %in% left.points, ][order(pointid),]$dis_m3s, na.rm = T)
              # take the sum of discharge, of all pixels flowing into confluence
              cda <- sum(m[m$pointid %in% left.points,]$dis_m3s)
              
              # assign new value to cumulative time in confluence pixel based on:
              # ctime = sum of discharge weighted travel time of pixels flowing into the confluence
              # cda = sum of discharge of pixels flowing into the confluence
              # divided by the residence time within the confluence pixel
              # nx.trav <- ctime / cda + m[m$pointid == k, ]$wrt_min.wb
              
              temp.df[z, "travel.time"] <- ctime / cda + m[m$pointid == k, ]$wrt_min.wb
              # if(!is.na(cumtime[cumtime$pointid == k, ]$travel.time)){
              #   # if there is already a value assigned, pick the bigger value
              #   travs <- c(cumtime[cumtime$pointid == k, ]$travel.time, nx.trav)
              #   cumtime[cumtime$pointid == k, ]$travel.time <- travs[which.max(travs)]
              # } else {
              #   cumtime[cumtime$pointid == k, ]$travel.time <-
              #     nx.trav
              # }
              # Once we have corrected the confluence travel time by the discharge of the merging streams,
              # we can trail the stream down until the next confluence
              nx.px <- m[m$pointid == k,]$next.pixel
              
              while (m[m$pointid == nx.px, ]$flag_conf != 1) {
                # while next pixel is not a confluence, do...
                # calculate cumulative travel time by...
                # if (m[m$pointid == k, ]$flag_mouth != 1) {
                  #print(paste0("I'm at...", k))
                  z <- z + 1
                  # calculate cumulative travel time by...
                  temp.df[z, "pointid"] <- nx.px
                  temp.df[z, "travel.time"] <- 
                    temp.df[temp.df$pointid == k, ]$travel.time + # Time A
                    (m[m$pointid == nx.px, ]$wrt_min.wb * # Time B
                       # adding wrt of next pixel to the weighted wrt of the current pixel
                       (1 - (m[m$pointid == nx.px, ]$dis_m3s - m[m$pointid == k, ]$dis_m3s) /
                          m[m$pointid == nx.px, ]$dis_m3s))
                  
                  #cumtime[cumtime$pointid == nx.px, ]$travel.time <-
                  
                  # wrt of next pixel is weighted by discharge of next pixel minus discharge of current pixel
                  # divided by discharge of next pixel; equation:
                  # Time B + Time A * (1-(dis_B - dis_A) / dis_B)
                  
                  # once this is done, overwrite k with the next pixel
                  k <-
                    m[m$pointid == k, ]$next.pixel # so the loop jumps to the next pixel
                  nx.px <- m[m$pointid == k,]$next.pixel
                  
                  if(is.na(nx.px)){
                    break
                  }
                # } else {
                #   #if it hits the mouth, stop
                #   break
                # }
                # if(m[m$pointid == k, "flag_conf"] == 1){
                #   cat(paste0("\rReached confluence of stream of order ", s,". ", length(conf) - i, " streams within order left."))
                # }
                
              }
              
              temp.df <- temp.df[!is.na(temp.df$travel.time), ]
              # out<- rbind(out, cumtime)
              return(temp.df)
              # path <- unique(path)
              # return(cumtime[cumtime$pointid %in% path, ])
              
            } else {
              # if no stream flowing into confluence has WRT, skip
              break
            }
          }
        }
      # save results
      cumtime <-
        cumtime[temp, travel.time := i.travel.time, on = .(pointid)]
      
      # t <- merge(m, cumtime, all = TRUE)
      # ggplot(t %>% filter(strah_ord <= 2), aes(x = fllength_dec, y = log10(travel.time), colour = strah_ord)) +
      #   geom_point()
      sec.level <- sec.level[!(sec.level %in% first)]
      }
      
      
      print(paste("Finished stream orders...", s))
      # save before moving on to next order
      
    }
    # save final output
    saveRDS(cumtime, paste0("./Objects/cumtime_", y, "-", mon, "_", Sys.Date(), ".rds"))
    rm(cumtime, streamnet, path, m, temp, sources)
  }}

# # Read in outputs
cumtime <- readRDS("/run/media/mstadler/LAROMAINE/PhD/Analyses/DOM//Data/Traveltime/cumtime_2015-06.rds") %>% setDT()
streamnet <-
  readRDS("/run/media/mstadler/LAROMAINE/PhD/Analyses/DOM//Data/Traveltime/flacc3000_streamnet_wWRT_2015-06.rds")
m <- streamnet
m <- m[!duplicated(m$pointid), ]

 m <- m[water.body == "Stream" | water.body == "River", flag_main.chan := 1]

temp <- merge(m, cumtime, by = "pointid")
temp$water.body <- factor(temp$water.body, 
                          levels = c("Stream","River","Reservoir","Lake"))
ggplot(temp, aes(x = log10(wrt_min.wb), 
                 y = log10(travel.time), colour = water.body)) +
  theme_bw() +
  geom_point() +
  scale_colour_viridis_d(option = "mako", direction = -1) +
  labs(x = "Log WRT (min)", y = "Log FWWA (min)")

ggsave("./Figures/Methods/WRT_vs_FWWA_2015.06.png", last_plot(),
      width = 10, height = 7, units = "cm", dpi = 300)
# 
# temp[pointid == 547207,] # 25
# temp[buf.lake_id == "7808" & flag_lake.out == 1,] #242364
# temp[pointid == 518954,]
# 
# ok <- c(292, 2460, 6207)
# t <- temp[flag_main.chan == 1 & is.na(travel.time) & !(reach_id %in% ok),]
# 
# temp[pointid == 519998,]
# temp[next.pixel == 530898,]
# View(temp[reach_id == 7215,])
# 
# t <- temp %>% filter(water.body == "Lake" & (flag_lake.in == 1 | flag_lake.out == 1)) %>% 
#   dplyr::select(buf.lake_id, travel.time, flag_lake.out, flag_lake.in)
# t[flag_lake.out == 1, flag := "lake.out"]
# t[flag_lake.in == 1, flag := "lake.in"]
# 
# t <- t %>% dplyr::select(buf.lake_id, flag, travel.time) %>%
#   pivot_wider(values_from = travel.time, names_from = flag) %>% setDT()
# 
# t[, lake.diff := lake.out - lake.in]
# 
# View(t[lake.diff < 0,])
# 
# temp[buf.lake_id == "90580" & flag_lake.out == 1,]
# 
# library(RColorBrewer)
# ggplot(temp %>% filter(strah_ord <= 7), aes(x = fllength_dec, y = log10(travel.time), colour = as.factor(strah_ord))) +
#   theme_bw() + 
#   scale_colour_manual(name = "Strahler order", values = brewer.pal(9, "Blues")[-c(1,2)]) +
#   geom_point()
# 
# plot_ly(temp %>% filter(strah_ord <= 2), x = ~fllength_dec, y = ~travel.time, colour = ~strah_ord, name = ~pointid,
#         type = "scatter", mode = "markers") 
# 
# p <- ggplot(temp %>% filter(strah_ord <= 5), aes(x = coord_x, y = coord_y, colour = log10(travel.time))) +
#   geom_point() +
#   theme(axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())
# 
# ggplotly(p)

# reach_id's: 292, 490, 2460, 6207, 7476 (checked, ok)
# reach_id's: 3810, 6387, 6373, 6627, 6948, 7477-7486 (checked, wrong corrected in script 2)

# Clean travel time output

# R set-up ----------------------------------------------------------------------------------------------------
### Packages -------------------------------------------------------------------------------
pckgs <- list("plyr", "data.table", "tidyverse") # wrangling
              
### Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)
# Many packages have to be installed through Bioconductor, please refer to the package websites

### Load
invisible(lapply(pckgs, library, character.only = T))
rm(pckgs, miss.pckgs)

# Read in calculated travel time --------------------------------------------------------------------------------------
# 2015-2016 -----------------------------------------------------------------------------------------------------------

y <- c(2015:2016)
mon <- c("6","8")

# create empty data frames to fill
out <- data.frame()

# Read in files and merge them together
for(i in 1:length(y)){
  for(j in 1:length(mon)){
  # Read in outputs
  cumtime <- readRDS(paste0("./Data/Traveltime/cumtime_",y[i], "-", mon[j],"_2024-02-09.rds")) %>% setDT()
  
  if(mon[j] %in% c("6","8")){
    old <- readRDS(paste0("./Data/Traveltime/flacc3000_streamnet_wWRT_", y[i], "-0", mon[j],".rds"))
  } else {
    old <- readRDS(paste0("./Data/Traveltime/flacc3000_streamnet_wWRT_", y[i], "-", mon[j],".rds"))
  }
  
  streamnet <-
    readRDS(paste0("./Data/Traveltime/flacc3000_streamnet_wWRT_", y[i], "-", mon[j], "_2024-02-09.rds"))
  
  streamnet <- merge(old %>% dplyr::select(Year:reach_id, buf.lake_id:strah_ord, flag_source:coord_y),
                     streamnet %>% dplyr::select(pointid, water.body, velocity_ms:wrt_min.wb), by = "pointid")
  
  m <- streamnet
  m <- m[!duplicated(m$pointid), ]
  m <- m[water.body == "Fluvial", flag_main.chan := 1]
  
  temp <- merge(m, cumtime, by = "pointid")
  
  temp[, water.body := factor(water.body, levels = c("Fluvial", "Reservoir", "Lake"))]
  temp[, travel.time_d := travel.time * 0.00069444]
  
  temp[, coord_x := as.numeric(coord_x)]
  temp[, coord_y := as.numeric(coord_y)]
  
  out <- rbind(out, temp)
  
  # t <- temp[buf.lake_id %in% c("937893"),]
  # out <- rbind(out, t)
  #write.table(temp[flacc_km2 > 4000,], paste0("./Data/Traveltime/toplot_tt_",camps[i],".csv"), sep = ";", row.names = F)
  }}


# 2017-2018 --------------------------------------------------------------------------------------------------------
y <- c(2017:2018)
mon <- c("6","8","10")

# Read in files and merge them together
for(i in 1:length(y)){
  for(j in 1:length(mon)){
    # Read in outputs
    cumtime <- readRDS(paste0("./Data/Traveltime/cumtime_",y[i], "-", mon[j],"_2024-02-09.rds")) %>% setDT()
    
    if(mon[j] %in% c("6","8")){
      old <- readRDS(paste0("./Data/Traveltime/flacc3000_streamnet_wWRT_", y[i], "-0", mon[j],".rds"))
    } else {
      old <- readRDS(paste0("./Data/Traveltime/flacc3000_streamnet_wWRT_", y[i], "-", mon[j],".rds"))
    }
    
    streamnet <-
      readRDS(paste0("./Data/Traveltime/flacc3000_streamnet_wWRT_", y[i], "-", mon[j], "_2024-02-09.rds"))
    
    streamnet <- merge(old %>% dplyr::select(Year:reach_id, buf.lake_id:strah_ord, flag_source:coord_y),
                       streamnet %>% dplyr::select(pointid, water.body, velocity_ms:wrt_min.wb), by = "pointid")
    
    m <- streamnet
    m <- m[!duplicated(m$pointid), ]
    m <- m[water.body == "Fluvial", flag_main.chan := 1]
    
    temp <- merge(m, cumtime, by = "pointid")
    
    temp[, water.body := factor(water.body, levels = c("Fluvial", "Reservoir", "Lake"))]
    temp[, travel.time_d := travel.time * 0.00069444]
    
    temp[, coord_x := as.numeric(coord_x)]
    temp[, coord_y := as.numeric(coord_y)]
    
    out <- rbind(out, temp)
    
    # t <- temp[buf.lake_id %in% c("937893"),]
    # out <- rbind(out, t)
    #write.table(temp[flacc_km2 > 4000,], paste0("./Data/Traveltime/toplot_tt_",camps[i],".csv"), sep = ";", row.names = F)
  }}

max.trav <- out[flag_mouth == 1,]

max.trav[, mouth.trav.time_d := travel.time * 0.00069444]
write.table(max.trav, "./Data/Traveltime/travel.time_byyear_atmouth.csv", sep = ",", dec = ".", row.names = F)

# Assign values to samples ---------------------------------------------------------------------
# read in samples and their pointids
rest <- read.csv("./Data/Traveltime/restsamples_snapped_mainchan_2018_pointid.txt", sep = ";", dec = ",") %>%
  dplyr::select(sample.name = sample_nam, pointid) %>% setDT()

main <- read.csv("./Data/Traveltime/mainsamples_snapped_mainchan_2018_pointid.txt", sep = ";", dec = ",") %>%
  dplyr::select(sample.name = sample_nam, pointid) %>% setDT()

# Read in samples from 2017, 2018 and their snapped pointids
s1718 <- read.csv("./Data/Traveltime/20172018_samples_DOMMic_near3000_pointid.txt", sep = ";", dec = ",") %>%
  dplyr::select(sample.name = sample_nam, pointid = flacc3000_) %>% setDT()

# gather the rest of meta data
mother.main <- read.csv("./Data/Summary/mother_mainsamples.csv", sep = ",", dec = ".") %>% setDT()
mother.rest <- read.csv("./Data/Summary/mother_restsamples.csv", sep = ",", dec = ".") %>% setDT()
mother.meta <- rbind(mother.main, mother.rest)

# combine the samples on main channel and the rest snapped to tributaries
meta <- rbind(main, rest)

# add 2017/18 data into meta
meta <- meta[s1718, pointid := i.pointid, on = .(sample.name)]

# correct a few manually
meta[pointid == 525799, pointid := 526769]
meta[pointid == 457095, pointid := 452871]

# merge the two
meta <- mother.meta[meta, , on = .(sample.name)]
# rm(main, rest, mother.main, mother.rest, mother.meta)
meta <- meta[!is.na(year),]
meta[, month := as.numeric(as.character(factor(campaign, levels = c(1,2,3), labels = c(6,8,10))))]

meta[sample.name %in% meta[sample.type.year == "Marine",]$sample.name, 
     pointid := 547615] # overwrite marine samples with mouth

nrow(meta) #1761

meta[, Year := as.factor(year)]
meta[, Month := as.factor(month)]
out[, Year := as.factor(Year)]
out[, Month := as.factor(Month)]


meta[out, c("flacc_km2", "strahler.order", "dis_m3s", "vel_ms",
                                    "wrt_min", "travel.time_min", "travel.time_days") :=
                         list(i.flacc_km2, i.strah_ord, i.dis_m3s, i.vel_ms.wb, 
                              i.wrt_min.wb, i.travel.time, i.travel.time_d), on = c("Year", "Month", "pointid")]


meta$sample.type.year <- factor(meta$sample.type.year, levels = c("Soil","Sediment",
                                                                  "Soilwater","Hyporheic",
                                                                  "Wellwater","Stream", "Tributary",
                                                                  "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
                                                                  "Upriver", "RO3", "RO2", "RO1","Deep", "Downriver",
                                                                  "Marine"),
                                labels = c("Soil","Sediment",
                                           "Soilwater","Soilwater",
                                           "Groundwater","Stream", "Tributary",
                                           "Riverine Lakes", "Headwater Ponds", "Lake", "Lake",
                                           "Upriver",# "RO3", "RO2", "RO1","Deep",
                                           "Reservoirs","Reservoirs", "Reservoirs","Reservoirs",
                                           "Downriver",
                                           "Estuary"))

# export factors for colouring
sample.factors <- levels(meta$sample.type.year)

# more colour blind friendly
colvec <- c("#FCFDBFFF", #"#FEC589FF", #Soil
            "#FEC589FF", #"#FDEBACFF", #Sediment
            "#F9795DFF", #"#F9795DFF", #Soilwater
            "#DE4968FF", #"#DE4968FF", #Groundwater,
            "skyblue", #Stream
            "#AD347CFF",# Tributary,
            "palegreen", #Riverine Lakes,
            "#7AD151FF", #Headwater Ponds,
            "#FDE725FF",# Lake,
            "#1F9F88FF", # Upriver,
            "orchid", #"#471063FF", #Reservoir,
            #'salmon',
            #"pink",
            #'navy',
            "#375A8CFF", #Downriver,
            "gray40") #Estuary)
names(colvec) <- as.character(sample.factors) # assign names for later easy selection

meta[, Season := factor(campaign, levels = c(1,2,3), labels = c("Spring", "Summer", "Autumn"))]
meta[, Year := factor(year, levels = c(2015:2021))]

write.table(meta, paste0("./Data/Summary/meta_file_traveltime_", Sys.Date(), ".csv"), sep = ",", dec = ".", row.names = F)

# wb.col <- brewer.pal(4, "Set2")
# wb.p <- ggplot(temp %>% filter(!is.na(travel.time)), 
#                aes(x = flacc_km2, y = travel.time * 0.00069444, colour = water.body)) + #as.factor(strah_ord)
#   theme_bw() +
#   labs(x = "Flow length (decimals)", y = "Travel time (days)", title = camps[i]) +
#   scale_colour_manual(name = "Water body type", values = c(wb.col[1],wb.col[3],wb.col[2],wb.col[4])) +
#   geom_point(alpha = 0.7) +
#   geom_rect(aes(xmin = 1.5, xmax = max(temp$fllength_dec), ymin = 0, ymax = 750), colour = "gray20", linetype = "dotted", alpha = 0)
# 
# in.p <- ggplot(temp %>% filter(!is.na(travel.time)), 
#                aes(x = fllength_dec, y = travel.time * 0.00069444, colour = water.body)) + #as.factor(strah_ord)
#   theme_pubr(base_size = 8, border = TRUE) +
#   theme(legend.position = "none", axis.title = element_blank()) +
#   labs(x = "Flow length (decimals)", y = "Travel time (days)") +
#   lims(y = c(0, 750), x = c(1.5, max(temp$fllength_dec))) +
#   scale_colour_manual(name = "Water body type", values = c(wb.col[1],wb.col[3],wb.col[2],wb.col[4])) +
#   geom_point(alpha = 0.7)
# 
# fin.p <- wb.p + annotation_custom(ggplotGrob(in.p), xmin = 1.5, xmax = 6, ymin = 1000, ymax = 4000)
# ggsave(paste0("./Figures/Analysis/LR_trav.time_", camps[i],".png"), fin.p,
#        width = 20, height = 15, units = "cm", dpi = 300)

meta <- read.csv("./Data/Summary/meta_file_traveltime_2024-02-10.csv") %>% setDT()

meta[, d.ex := d2H.h2o - (8 * dO18.h2o)]
meta[, Season := factor(Season, levels = c("Spring","Summer","Autumn"))]

meta[sample.type.year %in% c("Stream","Upriver","Downriver","Tributary") & !is.na(strahler.order), 
     cont := paste("FL order", strahler.order)]
meta[sample.type.year %in% c("Riverine Lakes","Lake","Deep"), cont := "Lake"]
meta[sample.type.year %in% c("Headwater Ponds"), cont := "Pond"]

meta[sample.type.year %in% c("Soil","Soilwater"), cont := "Terrestrial"]
meta[is.na(cont),  cont := sample.type.year]

meta[cont %in% c("FL order 1","FL order 2"),
     cont := "FL order 1-2"]

meta[cont %in% c("FL order 3","FL order 4"),
     cont := "FL order 3-4"]

meta[cont %in% c("FL order 5","FL order 6", "FL order 7"),
     cont := "FL order 5-7"]

ggplot(meta, aes(x = dO18.h2o, y = d2H.h2o)) +
  geom_point(aes(colour = Season)) +
  geom_abline(yintercept = 10, slope = 8)

meta[, cont := factor(cont, levels = c("Terrestrial","FL order 1-2","FL order 3-4", "FL order 5-7","Reservoirs",
                                       "Pond","Lake"))]
colvec <- c("#49C1ADFF", "#357BA2FF", "#3E356BFF",'#B58FC2','#BEBC48') 

(p <- ggplot(meta %>% filter(travel.time_days < 2000 & d.ex > 0 & d.ex < 18 & cont != "Terrestrial" & cont != "Pond"), 
       aes(x = d.ex, y = travel.time_days)) +
  theme_bw() +
  facet_wrap(~Season) +
  geom_point(aes(colour = cont)) +
  #geom_smooth(method = "lm") +
  scale_colour_manual(values = colvec, name = "Habitat") +
  labs(x = "d-excess", y = "Flow-weighted\nwater age (d)") +
  theme(plot.title = element_text(hjust = 0.5,colour = "grey20"),
        axis.text = element_text(colour = "grey50"),
        axis.title = element_text(colour = "grey20"),
        axis.ticks = element_line(colour = "grey50"),
        legend.title = element_text(colour = "grey20"),
        legend.text = element_text(colour = "grey30"),
        strip.background = element_rect(fill = "grey20"),
        strip.text = element_text(colour = "white")))

ggsave("./Figures/Suppl/travel.time_d.ex.png", p,
       dpi = 300, height =5, width =14, units = 'cm')


# entire dataset
cor.res <- cor.test(meta$d.ex, meta$travel.time_days,method = "spearman") # significant
cor.res$p.value
cor.res$estimate

#rho = -0.4322352
#p = 1.50928e-33

ddply(meta %>% filter(travel.time_days < 2000 & d.ex > 0 & d.ex < 18 & cont != "Terrestrial" & cont != "Pond"), 
      .(Season), function(x){
  cor.res <- cor.test(x$d.ex, x$travel.time_days, method = "spearman")
  data.frame(p = cor.res$p.value,
             rho = cor.res$estimate)
})

# Season            p         rho
# 1 Spring 8.644816e-11 -0.40306031
# 2 Summer 1.373443e-08 -0.33395908
# 3 Autumn 6.844951e-01 -0.03685094
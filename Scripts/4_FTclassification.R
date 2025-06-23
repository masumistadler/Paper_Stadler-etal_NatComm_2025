#-- Script for the publication:
#-- Title
#-- Responsible author: Masumi Stadler

#setwd("/home/bioinf/data/Molecular/Masumi/DOM-main")
# R set-up ----------------------------------------------------------------------------------------------------
### Packages -------------------------------------------------------------------------------
pckgs <- list(
  "data.table", "stringr", "tidyverse", # wrangling
  "plyr",
  "rstatix", "mgcv", #stats 
  "doMC", "foreach", # parallel
  'gridExtra', 'ggpubr','plotly') # plotting
### Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)
# Many packages have to be installed through Bioconductor, please refer to the package websites

### Load
invisible(lapply(pckgs, library, character.only = T))
rm(pckgs, miss.pckgs)

### Custom functions -----------------------------------------------------------------------------
source("./Functions/custom_fun.R")

### Parallel set-up -----------------------------------------------------------------------
# register number of cores for parallel computing with 'apply' family
detectCores() # 24 (private PC: 4), we do not have 24 cores. It's 12 cores, 24 threads.
registerDoMC()
getDoParWorkers() # 12 (private PC: 2)

# Read data ---------------------------------------------------------------------------------------------------
cross <- read.csv("./Data/FT/crosstable2015_cor.csv", sep =",", stringsAsFactors = F)
cross16 <- read.csv("./Data/FT/crosstable2016_cor.csv", sep =",", stringsAsFactors = F)

# read in matched mother files
meta15 <- read.csv("./Data/Summary/2015_metadata_upflacc.csv", sep = ",", stringsAsFactors = F)
meta16 <- read.csv("./Data/Summary/2016_metadata_upflacc.csv", sep = ",", stringsAsFactors = F)

meta.all <- read.csv("./Data/Summary/meta_file_traveltime_2024-02-10.csv", sep = ",", stringsAsFactors = F) %>% setDT()

# hist(meta.all[year == 2015,]$travel.time_d)
# hist(meta.all[year == 2016,]$travel.time_d)

# extract only a few samples that we have FT for
meta <- meta.all[sample.name %in% c(meta15$sample.name, meta16$sample.name),] %>% unique()

# join
meta <- left_join(meta,
                  rbind(meta15, meta16) %>% dplyr::select(sample.name, colnames(meta15)[!(colnames(meta15) %in% colnames(meta.all))]),
                  by = "sample.name")

meta[, (n = .N), by = .(year)]

# calculate some more FT variables
cross.tabs <- list(cross, cross16)

# Add FT variables -----------------------------------------------------------------------------------------

ios <- function(x){
  setDT(x)
  x[, molweight := (12.0107*C) + (1.00784*H) + (14.0067*N) + (15.999*O) + (32.065*S)]
  x[, IOS := F]
  x[HC <= 1.3 & HC >= 1.04 & OC <= 0.62 & OC >= 0.42 & molweight <= 388 & molweight >= 332, IOS := T]
  x[HC <= 1.3 & HC >= 1.04 & OC <= 0.62 & OC >= 0.42 & molweight <= 548 & molweight >= 446, IOS := T]
  return(x)
}
cross.tabs <- lapply(cross.tabs, ios)


# Merge and clean --------------------------------------------------------------------------------------
# melt normal samples into long format
cross.tabs[[1]] <- melt.data.table(cross.tabs[[1]], id.vars = c(1:12,15:18,21:26,122:123),
                                   measure.vars = 33:121,
                                   variable.name = "sample.name", value.name = "peak.int")
cross.tabs[[2]] <- melt.data.table(cross.tabs[[2]], id.vars = c(1:12,15:18,21:26,126:127),
                                   measure.vars = 33:125,
                                   variable.name = "sample.name", value.name = "peak.int")

# # merge together
all.df <- bind_rows(cross.tabs)

## Sample renaming ---------------------------------------------------------------------------------------
# some cleaning necessary
# change sample names from . to -
all.df[str_detect(sample.name, "^RO"), sample.name := str_replace(sample.name, "[.]", "-")]
meta[str_detect(sample.name, "^PR"), sample.name := str_replace(sample.name, "[-]", ".")]

all.df[str_detect(sample.name, "^L32"), sample.name := str_replace(sample.name, "[.]", "_")]

# RO1-10
meta[sample.name == "RO1-10D", sample.name := "RO1-10"]
all.df[sample.name == "RO1-10P", sample.name := "RO1-10"]

# RO2-22.r2
meta[sample.name == "RO2-22", sample.name := c("RO2-22","RO2-22.r2")]

# omit estuaries
est <- meta[sample.type.year == "Estuary",]$sample.name
all.df <- all.df[!(sample.name %in% est),]
meta <- meta[sample.type.year != "Estuary",]

#write.table(meta, paste0("./Data/Summary/meta_FTdataset_", Sys.Date(), ".csv"), sep = ",", dec = ".", row.names = F)

# ggplot(meta[sample.type.year != "Lake" & sample.type.year != "Riverine Lakes",], aes(x = sample.type.year, y =wrt_min)) +
#   geom_boxplot()

# merge with meta
all.df <- all.df[meta, , on = .(sample.name)] 

setDT(all.df)

## Define habitats -----------------------------------------------------------------------------------------
all.df$sample.type.year <- factor(all.df$sample.type.year, levels = c(#"Soil","Sediment",
  "Groundwater","Soilwater",#"Hyporheic", "Wellwater",
  "Stream", "Tributary",
  #"IslandLake",
  "Upriver", "Reservoirs",#"Deep",
  "Downriver",
 #"Estuary",
  "Headwater Ponds", "Riverine Lakes", "Lake"))

# save coordinates
out <- all.df %>% dplyr::select(sample.name, Season, year, sample.type.year, lat, long)
out <- out %>% distinct()
#write.table(out, "./Output/coordinates_FT.csv", sep = ",", row.names = F)

# df <- all.df %>% dplyr::select(sample.name, sample.type.year,
#                         DOC, lability, bact.abundance, TP, TN, temp.water.ysi, do.cor.perc.ysi) %>% distinct()
# 
# df <- melt(df, id.vars = c("sample.name","sample.type.year"))
# 
# ggplot(df) +
#   theme_bw()+
#   geom_boxplot(aes(y = value, x = sample.type.year)) +
#   facet_wrap(.~variable, scale = "free")+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# setDT(df)
# options(scipen =9999)
# df[, .(mean = mean(value, na.rm = T)), by = .(variable,  sample.type.year)]


## Remove outliers -------------------------------------------------------------------------------------------
# During plotting (later step in this script), a few outliers were detected
# We apply a function to remove an outlier observation per group
# We identify a group by: year and molecular formula
# If a observation was identified as an outlier, it will be removed
# However, if there are several points of outliers they will be kept
all.df[, ID := paste(year, molecular.formula, sep = "_")]

# detect extreme group outliers, except zeros
all.df[peak.int > 0, ext_outlier := is_extreme(peak.int), by = .(Year, Season, molecular.formula)]

# overwrite all zeros with FALSE
all.df[peak.int == 0, ext_outlier := FALSE]
# check how many outliers
nrow(all.df[ext_outlier == TRUE,]) * 100 / nrow(all.df) # 1.05% outliers
# remove extreme outliers
clean.df <- all.df[ext_outlier == TRUE, peak.int := 0]

## Separate absent and present MF ---------------------------------------------------------------------------
# Define absence
clean.df[, ID := paste(Year, Season, molecular.formula, sep = "_")]
clean.df[, n := .N, by = .(ID)] # number of samples per category
clean.df[, n.obs := nrow(.SD[peak.int > 0,]), by = .(ID)] # number of samples with actual observation

#initiate absence column
clean.df[n.obs == 0, PA := "Absent", by = .(ID)]

levels(factor(clean.df$sample.type.year))

#how many molecular formulae that are absent?
length(levels(factor(clean.df[PA == "Absent",]$molecular.formula))) # 2283
length(levels(factor(clean.df$molecular.formula))) # 13876
# percentage lost
length(levels(factor(clean.df[PA == "Absent",]$molecular.formula))) * 100 / length(levels(factor(clean.df$molecular.formula)))
# 16.5% absent

# Overwrite the rest with present
clean.df[is.na(PA), PA := "Present"]

levels(factor(clean.df[PA == "Present",]$sample.type.year))
levels(factor(clean.df[PA == "Absent",]$sample.type.year))
absent <- clean.df[PA == "Absent",]

# extract only df with present molecular formulae
present <- clean.df[PA == "Present",]
levels(factor(present$sample.type.year))

## Define travel time for non estimated habitats ----------------------------------------------------------
#present <- readRDS(select_newest("./Objects", "presentFT_long"))

#ggplot(present, aes(x = travel.time_d, y = n.obs, colour = sample.type.year)) + geom_point()

present[, travel.time_d := travel.time_days]
present[sample.type.year == "Soilwater", travel.time_d := -50] #
present[sample.type.year == "Groundwater", travel.time_d := -50]
#present[sample.type.year == "Estuary", travel.time_d := 365]

## Transform variables and remove zero observations --------------------------------------------------------
#present[, log.travel.time_d := log10(travel.time_d)]
present[peak.int == 0, peak.int := NA]

ggplot(present[molecular.formula == "C10H10O10N0S0",], aes(x =travel.time_d, y = peak.int, colour = sample.type.year)) +
  geom_point()

#remove samples with no travel time or season
#present <- present[!is.na(log.travel.time_d),]
#present <- present[!is.na(Season),]

# remove some wrongly snapped samples
#present <- present[flacc_npx != -9999,]

# arrange so that the data frame follows travel time
# ID is year_campaign_molecular.formula
present <- present %>% group_by(ID) %>% dplyr::arrange(travel.time_d) %>% setDT()

# save
# saveRDS(clean.df, paste0("./Objects/allFT_long_", Sys.Date(),".rds"))
# saveRDS(present, paste0("./Objects/presentFT_long_", Sys.Date(), ".rds"))
# saveRDS(absent, paste0("./Objects/absentFT_long_", Sys.Date(), ".rds"))

## Z-Scale -------------------------------------------------------------------------------------------------
temp <- dcast(present, sample.name ~ molecular.formula, value.var = "peak.int")
# randomize
#temp <- temp[, (colnames(temp)[-1]) := lapply(.SD, replace_na, 0), .SDcols = colnames(temp)[-1]]
scaled <- cbind(temp[,1], scale(temp[,-1], center = F, scale = T))

# merge back to present
temp <- melt(scaled, id.vars = "sample.name", variable.name = "molecular.formula", value.name = "z.peak.int")
present <- present[temp, , on = .(sample.name, molecular.formula)][!is.na(peak.int),]

# remove Autumn -----------------------------------------------------------------------------------------
present <- present[Season != "Autumn",]

## Filter too rare MF -------------------------------------------------------------------------------------
# Only keep those with more than 8 observations
length(unique(levels(factor(present$ID)))) #42484
hist(present$n.obs, breaks = seq(0, max(present$n.obs) + 10, 1))
# histogram is hard, three major peaks
present <- present[n.obs >= 8,]
length(unique(levels(factor(present$ID)))) #36193

#saveRDS(present, "./Objects/FT_present.rds")

## Randomization approach --------------------------------------------------------------------------------
# For each ID, shuffle the peak intensity along the x-axis and run regression 999 times

# this is actually a mini version of the full function.
# run script 4.1_FT_randomize.R on Compute Canada

# slope.random <- ddply(x, .(ID), function(x) {
# 
#   out <- foreach(i = 1:999, .combine = rbind) %dopar% {
#   # set random iteration seed
#   set.seed(i)
#   # extract data
#   df <- data.frame(x = x$travel.time_d, y = x$z.peak.int)
#   # shuffle y over x
#   df$y <- sample(df$y, replace = F, size = nrow(df))
#   # do linear regression
#   lin <- lm(df$y ~ df$x)
#   # extract run and slope
#   out <- data.frame(run = i, slope = lin$coefficients[2])
#   return(out)
#   }
#   
#   return(out)
# }, .parallel = T)

present <- readRDS("./Objects/FT_present_2024-02-10.rds") %>% setDT()
# read in results
out <- readRDS("./Objects/FT_randomization_output_2024-02-10.rds") %>% setDT()
ci <- readRDS("./Objects/FT_randomization_CI_2024-02-10.rds") %>% setDT()
ci <- ci[order(ID),]

# First filter, calculate the number of reshuffling that worked for a specific ID
out <- out %>% distinct() %>% setDT()
# count the number of runs that did not work
out[is.na(slope), n.na := .N, by = .(ID)]
out[is.na(n.na), n.na := 0]
# count the number of runs that did work
out[!is.na(slope), n.ok := .N, by = .(ID)]
out[is.na(n.ok), n.ok := 0]
# summarize
check <- out[, .(n.na = unique(n.na),
                 n.ok = unique(n.ok)), by = .(ID)]
nrow(check) == length(unique(out$ID))

# IDs that did not work
not.working <- unique(check[n.ok == 0,]$ID)
length(not.working) # 0

# constant but above 0
const <- unique(ci[sd == 0 & mean > 0,]$ID)
length(const) # 0

set.seed(3)
exmpl <- sample(unique(present$ID), 20)

for(i in 1:length(exmpl)){
  # correlation plot
  cor <- ggplot(present[ID == exmpl[i],], aes(x = travel.time_d, y = z.peak.int)) +
    theme_bw() +
    geom_smooth(method = "lm") +
    geom_point()
  # histogram of p-values
  his <- ggplot(out[ID == exmpl[i],], aes(x = slope)) +
    theme_bw() +
    geom_histogram(fill = "white", colour = "gray50") +
    geom_vline(xintercept = ci[ID == exmpl[i] & vars == "slope",]$mean,
               colour = "gray20", linetype = "solid") +
    geom_vline(xintercept = ci[ID == exmpl[i] & vars == "slope",]$upper.bound,
               colour = "tomato", linetype = "dashed") +
    geom_vline(xintercept = ci[ID == exmpl[i] & vars == "slope",]$lower.bound,
               colour = "royalblue", linetype = "dashed")
  # boxplot with notches indicating 95% Confidence interval
  # notch <- ggplot(out[ID == exmpl[i],], aes(y = slope, x = 1)) +
  #   theme_bw() +
  #   geom_boxplot(fill = "white", colour = "gray50", notch = T, width = 0.5) +
  #   geom_hline(yintercept = ci[ID == exmpl[i] & vars == "slope",]$upper.bound,
  #              colour = "tomato", line.type = "dashed") +
  #   geom_hline(yintercept = ci[ID == exmpl[i] & vars == "slope",]$lower.bound,
  #              colour = "royalblue", line.type = "dashed")
  # histogram of R2
  r2 <- ggplot(out[ID == exmpl[i],], aes(x = r2)) +
    theme_bw() +
    geom_histogram(fill = "white", colour = "gray50")
  # histogram of p-values
  pval <- ggplot(out[ID == exmpl[i],], aes(x = p.val)) +
    theme_bw() +
    geom_histogram(fill = "white", colour = "gray50")
  
  # arranged
  arr <- ggarrange(ggarrange(cor, his, align = "v", nrow = 2), ggarrange(r2, pval, align = "hv"), nrow = 2, 
                   heights = c(1,0.5))
  ggsave(paste0("./Figures/Randomization/", exmpl[i],".png"), arr,
         width = 10, height = 15, units = "cm")
}

# For model application refer to script 4.2
# Classification ------------------------------------------------------------------------------------
# read objects
model.ls <- readRDS("./Objects/FT_model.ls_binned_2024-02-11.rds")
model.df <- readRDS("./Objects/FT_model.df_binned_2024-02-11.rds") %>% setDT()
pks.df <- readRDS("./Objects/FT_pks.df_binned_2024-02-11.rds") %>% setDT()

# Plot examples ------------------------------------------------------------------------------------
# remove NA models (not enough observations = none)

# # plot tests
# plots <- lapply(model.ls[1:10], function(x){
#   newd <- data.frame(x =seq(min(x$raw$x),max(x$raw$x), by = 1))
#   newd$y.pred <- predict(x$model, newd)
#   p <- ggplot()+
#     theme_bw() +
#     geom_point(aes(x = x$raw$x, y = x$raw$y), colour = "gray50", alpha = 0.5) +
#     geom_point(aes(x = x$binned$x, y = x$binned$y,
#                    fill = factor(x$binned$n, levels = seq(min(x$binned$n),
#                                                             max(x$binned$n), 1))),
#                size = 5, shape = 21) +
#     scale_fill_brewer(palette = "RdBu", direction = -1, name = "Sample size") +
#     geom_line(aes(x = newd$x, newd$y.pred), colour = "purple") +
#     labs(title = paste("Model = ", x$model.df$model, "\nDirection = ", x$model.df$direction),
#          x = "Travel time (d)",
#          y = "Peak intensity")
# 
#   if(!is.null(x$pks.df)){
#     if(x$pks.df$no.peak[1] > 0){
#     p <- p +
#       geom_point(aes(x = x$pks.df$peak.x, y = x$pks.df$peak.y), colour = "purple", size = 4, shape = 21) +
#       geom_point(aes(x = x$pks.df$infp.x, y = x$pks.df$infp.y, colour = x$pks.df$closest), size = 2) +
#       scale_colour_manual(values = c("tomato","skyblue"), name = "Inflection point")
#   }}
#   return(p)
# })
# 
# # save plots
# for(i in 1:10){
#   ggsave(paste0("./Figures/FT_lm_fit/", names(plots)[i],".png"), plots[[i]])
# }

# Extract raw data
# raw.df <- ldply(model.ls, function(x){
#   x$raw
# }) %>% setDT()

# Clean --------------------------------------------------------------------------------------------
model.df[is.na(p.val), AU := "n.s."]

nrow(model.df[is.na(p.val),]) * 100 / length(model.ls) # we loose 0.58%

model.df[model == "y ~ x", model.type := "LM"]
model.df[model == "y ~ poly(x, 3)" | model == "y ~ poly(x, 2)", model.type := "POLY"]
model.df[model == "y ~ s(x)", model.type := "GAM"]
model.df[is.na(model.type),]

# count numbers of model by type
model.df[, .(n = .N), by = .(model.type)]

model.df <- model.df %>% separate(col = "ID", into = c("year","season","MF"), sep = "[_]", remove = F) %>% setDT()

# Add randomization --------------------------------------------------------------------------------
model.df[ci[vars == "slope",], c("rand.ci.up", "rand.ci.low") := list(i.upper.bound, i.lower.bound), on = .(ID)]

length(unique(model.df$ID)) # 36193 models

# direction is defined 1 == positive, 2 == negative
# see how many pass the randomization filter
model.df[direction == 1 & (lm.slope > rand.ci.up), r.test := TRUE]
model.df[direction == 2 & (lm.slope < rand.ci.low), r.test := TRUE]
model.df[is.na(r.test), r.test := FALSE]

model.df[, .(n = .N), by = .(r.test)] # 35312 pass

# see how many pass p-value filter
model.df[p.val < 0.05, .(n = .N)] # 16073
model.df[p.val < 0.01, .(n = .N)] # 9477

# see how many pass both
model.df[p.val < 0.05 & r.test == TRUE, .(n = .N)] # 15905
model.df[p.val < 0.01 & r.test == TRUE, .(n = .N)] # 9384

# assign flags
model.df[, loose.filt := ifelse(r.test == TRUE, TRUE, FALSE)]
model.df[, cons.filt := ifelse(r.test == TRUE & p.val < 0.05, TRUE, FALSE)]

# Define a few truly flat models ------------------------------------------------------------------------
model.df[cons.filt == TRUE & sd.y == 0, ] # none
model.df[sd.y == 0, AU := "flat"]

# Classification ------------------------------------------------------------------------------------
# classify not significant
model.df[is.na(AU) & is.na(p.val), ns.s := "n.s."]

## Add ranges -----------------------------------------------------------------------------------------------
# Molecules that occur along the whole continuum
limits <- present[peak.int != 0, .(min.x = min(travel.time_d, na.rm = T),
                                max.x = max(travel.time_d, na.rm = T)), by = .(ID)]
limits <- limits[ID %in% names(model.ls),]

model.df <- model.df[limits, , on = .(ID)] #c(".id==ID")

# model.df[min.x <= -15 & max.x >= 5, range := "entire"]
# # Molecules that only occur in freshwater
# model.df[min.x >= -15, range := "freshwater"]
# # Molecules that do not occur in estuary
# model.df[min.x <= -15 & max.x < 6, range := "soil-fresh"]

## Classify models ----------------------------------------------------------------------------------------
# Classify LMs
# classification based on their slope direction
model.df[is.na(AU) & model.type == "LM" & direction == 1, AU := "increase"] #& p.val < 0.05
model.df[is.na(AU) & model.type == "LM" & direction == 2, AU := "decrease"] #& p.val < 0.05
model.df[is.na(AU) & model.type == "LM" & lm.slope == 0, AU := "flat"]

#sanity check
model.df[model.type == "LM" & is.na(AU),]

# One peakers -------------------------------------------------------------------------------------------
# Add a few variables for classification of patterns
model.df[, num.range := max.x - min.x]
model.df[, range.cntr := (min.x + max.x) / 2]
model.df[, buffer.cntr := num.range / 6]
model.df[, buffer.min := range.cntr - buffer.cntr]
model.df[, buffer.max := range.cntr + buffer.cntr]

# Merge model vars to peaks data frame
df <- pks.df[model.df, , on = .(ID)]
df <- df[is.na(AU), ]

# Apply classification settings
df[is.na(AU) & no.peak == 1 &
     peak.x >= buffer.min & peak.x <= buffer.max, AU := "unimodal"]
df[is.na(AU) & no.peak == 1 & peak.x < range.cntr & model.type != "LM", AU := "non-linear decrease"]
df[is.na(AU) & no.peak == 1 & peak.x > range.cntr & model.type != "LM", AU := "non-linear increase"]
# sanity check, one peakers are classified?
df[is.na(AU) & no.peak == 1,]
# which were identified?
unique(df$AU)

# merge into model.df
model.df[ID %in% df[AU == "unimodal",]$ID, AU := "unimodal"]
model.df[ID %in% df[AU == "non-linear decrease",]$ID, AU := "non-linear decrease"]
model.df[ID %in% df[AU == "non-linear increase",]$ID, AU := "non-linear increase"]

# Check how many for each AU
model.df[, .(n = .N), by = .(AU)]

# Bi-peakers ---------------------------------------------------------------------------------------------------
df <- pks.df[model.df, , on = .(ID)]
df <- df[is.na(AU), ]

# Classify two-peakers ------------------------------------------------------------------------------------------
## Quadratic ----------------------------------------------------------------------------------------------------
temp <- df[is.na(AU) & no.peak == 2 & peak == 1 & peak.x <= (min.x + buffer.cntr), ]$ID
# 1st peak is at the beginning of x-axis
temp <- df[(ID %in% temp) & peak == 2 & peak.x >= (max.x - buffer.cntr),]$ID
# 2nd peak is at the end of x-axis
pks.diff <- dcast(df[ID %in% temp,], ID ~ peak, value.var = "peak.y", subset = .(closest == "above"))
pks.diff[, diff := `1`-`2`]
# when the difference is < 0, then the first peak is smaller than second = increasing
# when the difference is > 0, then the first peak is bigger than second = decreasing
model.df[ID %in% pks.diff[diff > 0,]$ID, AU := "decreasing quadratic"]
model.df[ID %in% pks.diff[diff < 0,]$ID, AU := "increasing quadratic"]
df[ID %in% pks.diff[diff > 0,]$ID, AU := "decreasing quadratic"]
df[ID %in% pks.diff[diff < 0,]$ID, AU := "increasing quadratic"]

# Bimodal ------------------------------------------------------------------------------------------------------
pks.diff <- dcast(df[is.na(AU) & no.peak == 2, ], ID ~ peak, value.var = "peak.y", subset = .(closest == "above"))
pks.diff[, diff := `1` - `2`]

model.df[ID %in% pks.diff[diff > 0,]$ID, AU := "decreasing bimodal"]
model.df[ID %in% pks.diff[diff < 0,]$ID, AU := "increasing bimodal"]
df[ID %in% pks.diff[diff > 0,]$ID, AU := "decreasing bimodal"]
df[ID %in% pks.diff[diff < 0,]$ID, AU := "increasing bimodal"]

# sanity check, all two peakers are classified?
df[is.na(AU) & no.peak == 2,]

# how many peaks max?
max(pks.df$no.peak) # 5

# what models were selected?
unique(model.df$AU)

# Classify multi-peakers ------------------------------------------------------------------------------------------
pks.diff <- dcast(df[is.na(AU), ], ID ~ peak, value.var = "peak.y", subset = .(closest == "above"))
pks.diff <- melt(pks.diff, id.vars = "ID", 
                 measure.vars = c("1","2","3","4"), 
                 variable.name = "peak", value.name = "y" ) %>% arrange(ID)
pks.diff <- pks.diff[!is.na(y),]

# calculate the difference between the first and last peak
opposite.diff <- ddply(pks.diff, .(ID), function(x){
  data.frame(max.diff = x[x$peak == 1,]$y - x[which.max(x$peak),]$y)
}) %>% setDT()

opposite.diff[max.diff < 0, AU := "non-linear increase"]
opposite.diff[max.diff > 0, AU := "non-linear decrease"]
opposite.diff[is.na(AU),] # none

model.df[ID %in% unique(opposite.diff[AU == "non-linear increase",]$ID), AU := "non-linear increase"]
model.df[ID %in% unique(opposite.diff[AU == "non-linear decrease",]$ID), AU := "non-linear decrease"]

#sanity check
model.df[is.na(AU),] # none, all classified

# add info whether non-linear models have humps
cast.pks <- dcast(pks.df[no.peak != 0,], ID + peak ~ closest, value.var = "infp.slope")
model.df[ID %in% cast.pks[!is.na(above) & !is.na(below),]$ID, hump := TRUE]
model.df[is.na(hump), hump := FALSE]

# identify those that are significant
model.df[!is.na(ns.s), c.ns.s := "n.s."]
model.df[is.na(ns.s) & loose.filt == F, ns.s := "n.s."]
model.df[is.na(c.ns.s) & cons.filt == F, c.ns.s := "n.s."]
model.df[is.na(ns.s) & loose.filt == T, ns.s := "sig."]
model.df[is.na(c.ns.s) & cons.filt == T, c.ns.s := "sig."]
model.df[AU == "no.p" | AU == "stable", ns.s := "n.s."]
model.df[AU == "stable" | AU == "flat",]
# quality filter
model.df[is.nan(p.val), c("ns.s", "c.ns.s") := list("n.s.", "n.s.")]
# sanity check
model.df[is.na(ns.s),] #none
model.df[is.na(c.ns.s),] #none

model.df[model.type == "none", AU := "flat"]


# Save classification -----------------------------------------------------------------------
model.df[is.na(p.val), AU := "n.s."]
unique(model.df$AU)
model.df[, AU := factor(AU, levels = c("n.s.", "flat",
                                       "increase","non-linear increase", "decreasing quadratic","decreasing bimodal",
                                       "unimodal",
                                       'decrease','non-linear decrease', 'increasing quadratic','increasing bimodal'))]

model.df <- model.df %>% dplyr::select(ID:p.val, min.x:buffer.max, rand.ci.up:cons.filt, model.type, hump, AU, ns.s, c.ns.s)

write.table(model.df,paste0("./Data/Summary/FT_patterns.classified_", Sys.Date(), ".csv"), sep = ",", row.names = F)


# Extra code -------------------------------------------------------------------------------------------------
# # Get a grasp on proportion of each AU -----------------------------------------------------------------
# n.au <- model.df[order(AU), .(n = .N), by = .(AU)]
# 
# # get groups of AU
# groups <- n.au$AU
# 
# folders<-list.dirs("./Figures/FT_lm_fit/", full.names = F)[-1]
# groups<-as.vector(groups[!(groups %in% folders)])
# 
# # multiple models in one plot -----------------------------------------------------------------------
# set.seed(3)
# for(i in 1:length(groups)){
#   dir.create(paste0("./Figures/FT_lm_fit/", groups[i]))
#   print(paste0("Plotting models of group: ", groups[i]))
#   # get all IDs from the same group
#   temp <- model.df[AU == groups[i],]$ID
#   
#   if(length(temp) > 160){
#     temp <- sample(temp, 160)
#   }
#   if(groups[i] == "flat" | groups[i] == "n.s."){
#     plots <- lapply(model.ls[names(model.ls) %in% temp], function(x){
#       p <- ggplot()+
#         theme_bw() +
#         theme(legend.position = "none") +
#         geom_point(aes(x = x$raw$x, y = x$raw$y), colour = "gray50", alpha = 0.5) +
#         geom_point(aes(x = x$binned$x, y = x$binned$y,
#                        fill = factor(x$binned$n, levels = seq(min(x$binned$n),
#                                                               max(x$binned$n), 1))),
#                    size = 5, shape = 21) +
#         scale_fill_brewer(palette = "RdBu", direction = -1, name = "Sample size") +
#         #geom_line(aes(x = newd$x, newd$y.pred), colour = "purple") +
#         labs(subtitle = paste("Model = ", x$model.df$model),
#              x = "Travel time (d)",
#              y = "Peak intensity (z-scaled)")
#       
#       if(!is.null(x$pks.df)){
#         if(x$pks.df$no.peak[1] > 0){
#           p <- p +
#             geom_point(aes(x = x$pks.df$peak.x, y = x$pks.df$peak.y), colour = "purple", size = 4, shape = 21) +
#             geom_point(aes(x = x$pks.df$infp.x, y = x$pks.df$infp.y, colour = x$pks.df$closest), size = 2) +
#             scale_colour_manual(values = c("tomato","skyblue"), name = "Inflection point")
#         }}
#       return(p)
#     })
#     
#   } else {
#     
#     plots <- lapply(model.ls[names(model.ls) %in% temp], function(x){
#       newd <- data.frame(x =seq(min(x$raw$x, na.rm = T),max(x$raw$x, na.rm = T), by = 1))
#       newd$y.pred <- predict(x$model, newd)
#       p <- ggplot()+
#         theme_bw() +
#         theme(legend.position = "none") +
#         geom_point(aes(x = x$raw$x, y = x$raw$y), colour = "gray50", alpha = 0.5) +
#         geom_point(aes(x = x$binned$x, y = x$binned$y,
#                        fill = factor(x$binned$n, levels = seq(min(x$binned$n),
#                                                               max(x$binned$n), 1))),
#                    size = 5, shape = 21) +
#         scale_fill_brewer(palette = "RdBu", direction = -1, name = "Sample size") +
#         geom_line(aes(x = newd$x, newd$y.pred), colour = "purple") +
#         labs(subtitle = paste("Model = ", x$model.df$model),
#              x = "Travel time (d)",
#              y = "Peak intensity (z-scaled)")
#       
#       if(!is.null(x$pks.df)){
#         if(x$pks.df$no.peak[1] > 0){
#           p <- p +
#             geom_point(aes(x = x$pks.df$peak.x, y = x$pks.df$peak.y), colour = "purple", size = 4, shape = 21) +
#             geom_point(aes(x = x$pks.df$infp.x, y = x$pks.df$infp.y, colour = x$pks.df$closest), size = 2) +
#             scale_colour_manual(values = c("tomato","skyblue"), name = "Inflection point")
#         }}
#       return(p)
#     })
#   }
#   
#   # Add ID to plot title
#   named <- mapply(function(x, y){
#     p <- x + labs(title=y)
#     return(p)
#   }, x = plots, y = names(plots), SIMPLIFY = F)
#   
#   plot.vec <- seq(from = 1, to = length(named), by = 16)
#   
#   # plot 9 models in one plot
#   for(j in 1:length(plot.vec)){
#     if(j != length(plot.vec)){
#       p <- do.call(grid.arrange, named[plot.vec[j]:(plot.vec[j+1]-1)])
#       ggsave(paste0("./Figures/FT_lm_fit/", groups[i],"/models_", plot.vec[j], "-", (plot.vec[j+1]-1),".pdf"), p,
#              height = 30, width = 35, units = "cm")
#     } else {
#       p <- do.call(grid.arrange, named[plot.vec[j]:length(named)])
#       ggsave(paste0("./Figures/FT_lm_fit/", groups[i],"/models_", plot.vec[j], "-", length(named),".pdf"), p,
#              height = 30, width = 35, units = "cm")
#     }
#   }
# }

# GAM: Multi-peakers -------------------------------------------------------------------------------------------
# df <- pks.df[model.df, , on = .(ID)]
# df <- df[is.na(AU), ]
# unique(df$no.peak) # how many peaks 3-5
# 
# # Classify generally declining and increasing multi-peakers
# df[is.na(AU) & no.peak >= 3 & direction == 2, AU := "declining_multimodal"]
# model.df[ID %in% unique(df[AU == "declining_multimodal",]$ID), AU := "declining_multimodal"]
# df[is.na(AU) & no.peak >= 3 & direction == 1, AU := "increasing_multimodal"]
# model.df[ID %in% unique(df[AU == "increasing_multimodal",]$ID), AU := "increasing_multimodal"]
# 
# # sanity check
# df[is.na(AU),] # all classified
# model.df[is.na(AU),] # all classified
# 
# # Check number of models classified into each category
# model.df[, .(n = .N), by = .(AU)]
# 
# # Plot them all ----------------------------------------------------------------------------------------------------------
# 
# groups<- c('declining_quadratic', "increasing_quadratic",
#            'declining_bimodal', "increasing_bimodal",
#            'declining_multimodal', "increasing_multimodal")
# 
# # multiple models in one plot -----------------------------------------------------------------------
# for(i in 1:length(groups)){
#   dir.create(paste0("./Figures/MO_lm_fit/", groups[i]))
#   print(paste0("Plotting models of group: ", groups[i]))
#   # get all IDs from the same group
#   temp <- model.df[AU == groups[i],]$ID
#   if(length(temp) > 500){
#     temp <- sample(temp, 500)
#   }
#   
#   # plot
#   plots <- lapply(model.ls[names(model.ls) %in% temp], function(x){
#     newd <- data.frame(x =seq(min(x$raw$x),max(x$raw$x), by = 1))
#     newd$y.pred <- predict(x$model, newd)
#     p <- ggplot()+
#       theme_bw() +
#       theme(legend.position = "none") +
#       geom_line(aes(x = newd$x, newd$y.pred), colour = "purple") +
#       geom_point(aes(x = x$raw$x, y = x$raw$y), colour = "gray50", alpha = 0.5, size = 3) +
#       labs(subtitle = paste("Model = ", x$model.df$model),
#            x = "Log Travel time (d)",
#            y = "Reads")
#     
#     if(!is.null(x$pks.df)){
#       if(x$pks.df$no.peak[1] > 0){
#         p <- p +
#           geom_point(aes(x = x$pks.df$peak.x, y = x$pks.df$peak.y), colour = "purple", shape = 21,size = 4) +
#           geom_point(aes(x = x$pks.df$infp.x, y = x$pks.df$infp.y, colour = x$pks.df$closest), size = 2, alpha = 0.5)
#       }}
#     return(p)
#   })
#   
#   # Add ID to plot title
#   named <- mapply(function(x, y){
#     p <- x + labs(title=y)
#     return(p)
#   }, x = plots, y = names(plots), SIMPLIFY = F)
#   
#   plot.vec <- seq(from = 1, to = length(named), by = 16)
#   
#   # plot 16 models in one plot
#   for(j in 1:length(plot.vec)){
#     if(j != length(plot.vec)){
#       p <- do.call(grid.arrange, named[plot.vec[j]:(plot.vec[j+1]-1)])
#       ggsave(paste0("./Figures/MO_lm_fit/", groups[i],"/models_", plot.vec[j], "-", (plot.vec[j+1]-1),".pdf"), p,
#              height = 30, width = 35, units = "cm")
#     } else {
#       p <- do.call(grid.arrange, named[plot.vec[j]:length(named)])
#       ggsave(paste0("./Figures/MO_lm_fit/", groups[i],"/models_", plot.vec[j], "-", length(named),".pdf"), p,
#              height = 30, width = 35, units = "cm")
#     }
#   }
# }

# classify where is the biggest peak among multi-peakers
# df[is.na(AU), diff.x := abs(shift(peak.x,2) - peak.x), by = .(.id)]
# df[is.na(AU), diff.y := abs(shift(peak.y,2) - peak.y), by = .(.id)]
# 
# df[diff.y]
# df[.id %in% df[is.na(AU) & no.peak == 2 & diff.y > mean.y & diff.x >= 15 & peak.y <= 0 & direction == 2,]$.id, AU := "2p-decay"]
# model.df[.id %in% df[AU == "2p-decay",]$.id, AU := "2p-decay"]
# 
# df[AU == "2p-unimodal" & diff.y > mean.y & direction == 2,]

#df[.id %in% df[is.na(AU) & no.peak == 2 & diff.x >= 15,]$.id, AU := "quadratic"]

#df[.id %in% df[is.na(AU) & no.peak == 2 & diff.x >= 10,]$.id, AU := "bimodal"]
#model.df[.id %in% df[AU == "quadratic",]$.id, AU := "quadratic"]
#model.df[.id %in% df[AU == "bimodal",]$.id, AU := "bimodal"]

# groups<- c('declining_bimodal', "increasing_bimodal",
#            'declining_multimodal', "increasing_multimodal")
# 
# for(i in 1:length(groups)){
#   dir.create(paste0("./Figures/MO_lm_fit/", groups[i]))
#   #sub <- model.df[AU == groups[i] & direction == 2,]$.id
#   sub <- sample(model.df[AU == groups[i],]$.id, 10)
#   #sub <- model.df[AU == groups[i],]$.id
#   plots <- lapply(model.ls[names(model.ls) %in% sub], function(x){
#     newd <- data.frame(x =seq(min(x$raw$x),max(x$raw$x), by = 1))
#     newd$y.pred <- predict(x$model, newd)
#     p <- ggplot()+
#       theme_bw() +
#       geom_point(aes(x = x$raw$x, y = x$raw$y), colour = "gray50", alpha = 0.5) +
#       geom_line(aes(x = newd$x, newd$y.pred), colour = "purple") +
#       labs(title = paste("Model = ", x$model.df$model, "\nDirection = ", x$model.df$direction),
#            x = "Log Travel time (d)",
#            y = "Reads")
#     
#     if(!is.null(x$pks.df)){
#       p <- p +
#         geom_point(aes(x = x$pks.df$peak.x, y = x$pks.df$peak.y), colour = "purple", shape = 21,size = 4) +
#         geom_point(aes(x = x$pks.df$infp.x, y = x$pks.df$infp.y, colour = x$pks.df$closest), size = 2, alpha = 0.5)
#     }
#     return(p)
#   })
#   
#   # save plots
#   for(j in 1:length(plots)){
#     ggsave(paste0("./Figures/MO_lm_fit/", groups[i],"/", "GAM_", names(plots)[j],".png"), plots[[j]],
#            height = 10, width = 10, units = "cm")
#   }
#   
# }
# 



# # Explore patterns of AU --------------------------------------------------------------------------------
# # split seasons
# model.df <- model.df %>% separate(col = "ID", into = c("Year","Season","MF"), sep = "[_]", remove = F)
# setDT(model.df)
# # overall percentage
# model.df[, n := .N]
# model.df[, n.au.all := .N, by = .(AU)]
# model.df[, perc.all := (n.au.all * 100) / n]
# 
# # sanity check
# model.df %>% dplyr::select(AU, perc.all) %>% distinct() %>% summarize(sum = sum(perc.all))
# 
# # overall percentage
# model.df[, n.camp := .N, by = .(Season)]
# model.df[, n.au.camp := .N, by = .(AU, Season)]
# model.df[, perc.camp := (n.au.camp * 100) / n.camp]
# 
# # sanity check
# t <- model.df %>% dplyr::select(AU, Season, perc.camp) %>% distinct() %>% arrange(Season, AU)
# t[, .(sum = sum(perc.camp)), by = .(Season)]
# 
# # export table
# out <- model.df %>% dplyr::select(Season, 
#                                   AU, perc.all, n.au.all, perc.camp) %>% distinct()
# 
# setDT(out)
# 
# out[, Season := factor(Season, levels = c("Spring", "Summer"))]
# # save to merge with FT data
# write.table(out, "./Data/Summary/AU_perc_overall_ft.csv", sep = ",", row.names = F)
# 
# # make a bar plot
# (p <- ggplot(out) +
#     theme_bw() +
#     geom_col(aes(x = Season, y = perc.camp, fill = AU)) +
#     #geom_text(aes(y = -1.5, label= unique(out$n.by.range), x = unique(out$range)), colour = "gray40", size = 2) +
#     scale_fill_viridis_d(name = 'Spatial pattern', option = "inferno", direction = -1)+
#     labs(x = "Season", y = "Percentage of spatial patterns (%)") +
#     theme(panel.grid = element_blank(),
#           axis.text.x = element_text(angle = 45, hjust = 1)))
# ggsave("./Figures/mic_perc.range.pattern.season_DNARNA.png", p, width = 14, height = 10, units = "cm")
# out.all <- model.df %>% dplyr::select(AU, perc.all) %>% distinct()
# # continuum


# add growth and decay rate -----------------------------------------------------------------------------
# rate.travtime <- function(x, y, min.x, max.x, break.n){
#   # remove zeros
#   temp <- data.frame(x = x, y = y)
#   #temp <- temp[!is.na(temp$y),]
#   #temp <- temp[temp$y > 0,]
#   avg <- mean(temp$x, na.rm = T)
#   # aggregate
#   bins <- tapply(temp$y, cut(round(temp$x,0), seq(min.x-1,
#                                                   max.x+1, by = break.n)), mean, na.rm = TRUE)
#   
#   out <- data.frame(y = bins)
#   setDT(out, keep.rownames = "bins")
#   out[, bins := str_remove(bins, "[(]")]
#   out[, bins := str_remove(bins, "[]]")]
#   out <- separate(out, col = "bins", into = c("start", "end"), sep = ",") %>% setDT()
#   out[, c('start', 'end') := list(as.numeric(start), as.numeric(end))]
#   out[, x := rowMeans(.SD), .SDcols = c("start", "end")]
#   out[is.na(y), y := 0]
#   
#   #y <- (out$y /avg)
#   y <- out$y
#   x <- out$x
#   
#   #apply linear model
#   lin <- try(lm(y ~ x), silent = T)
#   
#   lin.df <- data.frame(lin.intercept = coef(lin)[[1]],
#                        lin.slope = coef(lin)[[2]],
#                        lin.p.val = summary(lin)$coefficients[2,4])
#   return(lin.df)
# }
# 
# # apply on all MF
# 
# lin.ls <- dlply(present, .(ID), function(x) {
#   min.x <- minmax.trav[Season == unique(x$Season) & Year == unique(x$Year),]$min.trav
#   max.x <-  minmax.trav[Season == unique(x$Season) & Year == unique(x$Year),]$max.trav
#   
#   model.ls <- rate.travtime(x$log.travel.time_d, x$peak.int, min.x, max.x, break.n = 2)
#   
#   return(model.ls)
#   
# }, .parallel = T)
# 
# sub.ls <- lin.ls[!is.na(lin.ls)]
# length(lin.ls) - length(sub.ls) # we loose 21

## Range filter ------------------------------------------------------------------------------------------
# # remove those with extremely narrow ranges
# present[peak.int > 0, min.range := min(.SD$travel.time_d), by = .(ID)]
# present[peak.int > 0, max.range := max(.SD$travel.time_d), by = .(ID)]
# present[, range.span := max.range - min.range]
# hist(present$range.span)
# # we want at least 10
# hist(present$range.span)
# # set based on histogram, there are not that many below 10
# #present <- present[ID %in% present[range.span >= 4,]$ID,]
# #length(unique(levels(factor(present$ID)))) #17207
# 
# # include a filter for distribution of sample sizes across range 
# 
# check.bin.n <- function(x){
#   steps <- unique(x$range.span) / 3
#   bins <- c(unique(x$min.range), unique(x$min.range) + (steps), unique(x$min.range) + ((steps)*2), 
#             unique(x$max.range))
#   n.bins <- c(sum(x$travel.time_d >= bins[1] & x$travel.time_d < bins[2]),
#     sum(x$travel.time_d >= bins[2] & x$travel.time_d < bins[3]),
#     sum(x$travel.time_d >= bins[3] & x$travel.time_d <= bins[4]))
#   
#   return(any((n.bins >= 2) == F)) # if there are any bins with n < 2, return TRUE
#   }
# 
# present <- present[, bin.check := check.bin.n(.SD), by = .(ID)]
# length(unique(present$ID)) #36193
# # apply filter
# present <- present[bin.check == FALSE,]
# length(unique(present$ID)) #31810
# 
# rm(temp, t, scaled)
# 
# # Slope filter --------------------------------------------------------------------------------------------
# # Run lm on all models and classify those that are slope == 0 as stable even before running model decision tree
# slope.df <- ddply(present, .(ID), function(x) {
#   lin <- lm(x$z.peak.int ~ x$log.travel.time_d)
#   df <- data.frame(slope = lin$coefficients[2])
#   return(df)
# }, .parallel = T)
# 
# # identify those that have zero slope
# setDT(slope.df)
# stable <- slope.df[abs(slope) == 0,]$ID
# length(stable) # 0
# #present[ID %in% stable, AU := "stable"]
# 
# hist(present$log.travel.time_d) # better distribution than MO, a lot of points in water
# # go with binned version

# Add rates based on non z-scaled data --------------------------------------------------------------------------------
# slope.log2 <- ddply(present, .(ID), function(x) {
#   lin <- lm(x$peak.int ~ log2(x$travel.time_d))
#   df <- data.frame(slope.log2 = lin$coefficients[2])
#   return(df)
# }, .parallel = T)
# 
# slope.z.log2 <- ddply(present, .(ID), function(x) {
#   lin <- lm(x$z.peak.int ~ log2(x$travel.time_d))
#   df <- data.frame(slope.log2.z = lin$coefficients[2])
#   return(df)
# }, .parallel = T)
# 
# stat.peaks <- present[, .(mean.peak.int = mean(peak.int, na.rm = T),
#                           var.peak.int = var(peak.int, na.rm = T),
#                           median.peak.int = median(peak.int, na.rm = T)), by = .(ID)]
# 
# stat.peaks <- stat.peaks[slope.log2, , on = .(ID)]
# stat.peaks <- stat.peaks[slope.z.log2, , on = .(ID)]


# Testing and defining function ---------------------------------------------------------------------------

# temp <- present %>% sample_n(1)
# test <- present[ID %in% temp$ID,]
# 
# # set vectors
# x <- test$log.travel.time_d
# y <- test$z.peak.int
# 
# # try different models
# loess.mod <- loess(y ~ x, span = 0.7)
# lin <- lm(y ~ x)
# poly.lin <- lm(y ~ poly(x, degree = 2))
# poly.lin3 <- lm(y ~ poly(x, degree = 3))
# gam <- mgcv::gam(y ~ x)
# gam.poly <- mgcv::gam(y ~ s(x))
# 
# # predict in shorter intervals
# frq.x <- seq(min(x),max(x), by = 0.5)
# newd <- data.frame(x = frq.x)
# newd$y.loess <- as.numeric(predict(loess.mod, newdata = newd))
# newd$y.lm <- as.numeric(predict(lin, newdata = newd))
# newd$y.lmpoly <- as.numeric(predict(poly.lin, newdata = newd))
# newd$y.lmpoly3 <- as.numeric(predict(poly.lin3, newdata = newd))
# newd$y.gam <- as.numeric(predict(gam, newdata = newd))
# newd$y.gampoly <- as.numeric(predict(gam.poly, newdata = newd))

# #plot 
# ggplot()+
#   theme_bw() +
#   geom_point(data = test, aes(x = x, y = y, colour = sample.type.year), alpha = 0.5) +
#   geom_line(aes(x = newd$x, newd$y.loess), colour = "gray50") +
#   geom_line(aes(x = newd$x, y = newd$y.lm), colour = "tomato") +
#   geom_line(aes(x = newd$x, y = newd$y.lmpoly), colour = "gold") +
#   geom_line(aes(x = newd$x, y = newd$y.lmpoly3), colour = "navy") +
#   geom_line(aes(x = newd$x, y = newd$y.gam), colour = "pink") +
#   geom_line(aes(x = newd$x, y = newd$y.gampoly), colour = "purple") +
#   geom_point(aes(x = newd$x[localMaxima(newd$y.gampoly)],
#                  y = newd$y.gampoly[localMaxima(newd$y.gampoly)]),colour = "purple", size = 4)

# GAM seems best in terms of smoothness
# remove unncesseary objects
# rm(newd, lin, loess.mod, poly.lin, poly.lin3, temp, test, frq.x, x, y, gam , gam.poly)
# 
# lin.df <- ldply(lin.ls, function(x){
#   x
# })
# setDT(lin.df)
# # merge to general model data frame
# model.df <- model.df[lin.df, , on = c(".id==ID")]
# cross <- present %>% dplyr::select(molecular.formula, mz:molweight) %>% distinct() %>% setDT()
# 
# # extract a few columns
# 
# out <- model.df %>% dplyr::select(.id, Year, Season, MF, AU, model.type, mean.y, var.y, dev.expl, p.val, range, 
#                            lin.intercept:lin.p.val, n:perc.camp) %>% setDT()
# out <- out[cross, , on = c("MF==molecular.formula")]
# 
# write.table(out, "./Data/Summary/cross_AU_FT.csv", sep = ",", dec=  ".", row.names = F)

# P-value correction --------------------------------------------------------------------------------
# # Check p-value distribution
# hist(model.df[p.val < 0.1,]$p.val)
# 
# # number of models that worked
# (n <- nrow(model.df[is.na(AU),])) # 31708
# nrow(model.df) #31810
# 
# # create p-value bins
# p.bins <-  c(min(model.df$p.val, na.rm = T), 0.00001, 0.0001, 0.001, 0.01, 0.05, max(model.df$p.val, na.rm = T))
# 
# # theoretical false-positives
# fp <- round(p.bins * n, 0 )
# 
# # number of models positive within each bin
# ranges <- paste(head(p.bins,-1), p.bins[-1], sep=" - ")
# freq <- hist(model.df[is.na(AU),]$p.val, breaks = p.bins, include.lowest = T, plot = F)
# p.freq <- data.frame(range = ranges, frequency = freq$counts, freq.fp = fp[-1]) %>% setDT()
# p.freq[, freq.diff := frequency - freq.fp]
# p.freq[, p.bins := c(1:5,NA)]
# 
# # sanity check
# sum(p.freq$frequency) == n # should be TRUE
# 
# sum(p.freq[-nrow(p.freq),]$freq.diff)
# # we keep 6968 models this way
# 
# # Identify models to keep
# model.df <- model.df[order(rank(p.val)),]
# # separate models that are not included in this
# out.df <- model.df[!is.na(AU),]
# 
# model.df <- model.df[is.na(AU),]
# model.df[p.val <= 0.00001, p.bins := 1]
# model.df[p.val > 0.00001 & p.val <= 0.0001, p.bins := 2]
# model.df[p.val > 0.0001 & p.val <= 0.001, p.bins := 3]
# model.df[p.val > 0.001 & p.val <= 0.01, p.bins := 4]
# model.df[p.val > 0.01 & p.val <= 0.05, p.bins := 5]
# 
# p.freq <- p.freq[!is.na(p.bins),]
# 
# # overwrite p.val of models that didn't pass
# for(i in 1:nrow(p.freq)){
#   temp <- model.df[p.bins == p.freq$p.bins[i], head(.SD,p.freq$freq.diff[i])]$ID
#   model.df[p.bins == p.freq$p.bins[i] & !(ID %in% temp), c("p.cor") :=
#              list("remove")]
# }
# 
# #sanity check
# nrow(model.df[p.cor == "remove",]) == sum(p.freq$freq.fp) # should be TRUE
# 
# # bind the good and "bad" models back
# #model.df<-rbind(model.df %>% dplyr::select(-p.bins), out.df)
# model.df<-rbind(model.df %>% dplyr::select(-c(p.bins)), out.df %>% mutate(p.cor = "NA"))
# 
# # save corrected model.df
# saveRDS(model.df, paste0("./Objects/FT/model.df_FT_scaled.binned_p.corr_", Sys.Date(),".rds"))
# 
# rm(p.freq, ranges, freq, out.df)

# Check
# library(plotly)
# plot_ly(x = model.df[!is.na(AU) & p.val < 0.05 & (AU == "increase" | AU == "decrease") & model == "y ~ x",]$lm.slope, type = 'histogram')
# # maybe slope's less than abs(<0.006)
# # smaller than 0.02
# 
# temp <- model.df[!is.na(AU) & p.val < 0.05 & (AU == "increase" | AU == "decrease") & model == 'y ~ x',][abs(lm.slope) < 0.008,]$ID
# temp <- model.df[!is.na(AU) & p.val < 0.05 & (AU == "increase" | AU == "decrease"),][model == 'y ~ x' & var.y < abs(lm.slope),]$ID
# temp <- model.df[!is.na(AU) & p.val < 0.05 & (AU == "increase" | AU == "decrease"),][model == 'y ~ x' & abs(lm.slope) >= 0.01,]$ID
# temp <- model.df[!is.na(AU) & p.val < 0.05,][dev.expl <= 0.3,]$ID
# 
# 
# plot_ly(alpha = 0.9) %>%
#   add_histogram(x = ~lm.slope,
#                 data = model.df[!is.na(AU) & p.val < 0.05 & (AU == "increase" | AU == "decrease") & model == 'y ~ x',],
#                 name = "All linear models") %>%
#   add_histogram(x = ~lm.slope,
#                 data = model.df[ID %in% temp,], name = "slope > var") %>%
#   layout(yaxis = list(title = "Frequency")
#          )
#   
# plot_ly(alpha = 0.9) %>%
#   add_histogram(x = ~dev.expl,
#                 data = model.df[!is.na(AU) & p.val < 0.05 & (AU == "increase" | AU == "decrease") & model == 'y ~ x',],
#                 name = "All linear models") %>%
#   add_histogram(x = ~dev.expl,
#                 data = model.df[ID %in% temp,], name = "slope > var") %>%
#   layout(yaxis = list(title = "Frequency")
#   )
# 
# plot_ly(data = model.df[p.val < 0.05,], x = ~dev.expl, type = "histogram")
# 
# plot_ly(alpha = 0.9) %>%
#   add_histogram(x = ~dev.expl,
#                 data = model.df[p.val < 0.05,],
#                 name = "Dev. explained", xaxis = "x1") %>%
#   add_histogram(x = ~mean.y,
#                 data = model.df[p.val <0.05,],
#                 name = "Mean", xaxis = "x2") %>%
#   layout(yaxis = list(title = "Frequency"),
#          xaxis = list(side = "bottom", title = "Deviance explained"),
#          xaxis2 = list(side = 'top', overlaying = "x", title = "Mean")
#   )
# 
# plot_ly(alpha = 0.9) %>%
#   add_histogram(x = ~dev.expl,
#                 data = model.df[p.val < 0.05 & model == "y ~ x",],
#                 name = "Linear", xaxis = "x1") %>%
#   add_histogram(x = ~dev.expl,
#                 data = model.df[p.val < 0.05 & model == "y ~ s(x)",],
#                 name = "GAM", xaxis = "x1") %>%
#   add_histogram(x = ~dev.expl,
#                 data = model.df[p.val < 0.05 & model == "y ~ poly(x, 2)",],
#                 name = "Poly2", xaxis = "x1") %>%
#   add_histogram(x = ~dev.expl,
#                 data = model.df[p.val < 0.05 & model == "y ~ poly(x, 3)",],
#                 name = "Poly3", xaxis = "x1") %>%
#   layout(yaxis = list(title = "Frequency"),
#          xaxis = list(side = "bottom", title = "Deviance explained")
#   )
# 
# temp <- model.df[!is.na(AU) & p.val < 0.05 & (AU == "increase" | AU == "decrease"),][model == 'y ~ x' & dev.expl < 0.01,]$ID
# temp <- model.df[!is.na(AU) & p.val < 0.05 & (AU == "increase" | AU == "decrease"),][model == 'y ~ x' & dev.expl < 0.15 & abs(lm.slope) <= 0.01,]$ID
# temp <- model.df[!is.na(AU) & p.val < 0.05 & (AU == "increase" | AU == "decrease"),][model == 'y ~ x' & dev.expl < 0.185 & abs(lm.slope) <= 0.01,]$ID
# 
# #good?
# temp <- model.df[!is.na(AU) & p.val < 0.05 & (AU == "increase" | AU == "decrease"),][model == 'y ~ x' & dev.expl >= 0.385 & abs(lm.slope) > 0.01,]$ID
# 


# Visually check classification
# Looks good, some models seem spurious - update manual correction later on with all other model corrections

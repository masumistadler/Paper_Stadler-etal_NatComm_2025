#-- Script for the publication:
#-- Title
#-- Responsible author: Masumi Stadler

# This script is the third of a series of scripts that were used to analyse the data
# used in the publication.

#setwd("/home/bioinf/data/Molecular/Masumi/DOM-main")
# 1. R set-up ------------------------------------------------------------------------------
### Packages -------------------------------------------------------------------------------
pckgs <- list("phyloseq", "plyr", "tidyverse", "data.table", # wrangling
              "ggpubr",  "cowplot", "plotly", "gridExtra", # plotting, "ggnewscale","ggpmisc", 
              "kableExtra",# "xlsx", # making tables for manuscript (LaTeX) and export as excel (for ISME)
              "doMC", # parallel computing
              "vegan", "ape", "ade4", "rstatix", "mgcv") # statistics

### Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)
# Many packages have to be installed through Bioconductor, please refer to the package websites

### Load
invisible(lapply(pckgs, library, character.only = T))
rm(pckgs, miss.pckgs)

# Note: For package "xlsx" Java has to be installed
# For Linux users install dependencies: default-jre, default-jdk

### Functions -----------------------------------------------------------------------------
source("./Functions/custom_fun.R")

### Parallel set-up -----------------------------------------------------------------------
# register number of cores for parallel computing with 'apply' family
detectCores() # 24 (private PC: 4), we do not have 24 cores. It's 12 cores, 24 threads.
registerDoMC()
getDoParWorkers() # 12 (private PC: 2)

# prepare for parallel processing (with 'parallel') for 'vegan' functions
numCores <- detectCores()
cl <- makeCluster(numCores, type = "FORK") # using forking

# Set seed for session and reproducibility of permutations
# (NMDS were done, but not part of main workflow anymore)
set.seed(3)

# 2. Read and prepare data ---------------------------------------------------------------
# Output of first paper. Workflow in https://github.com/CarBBAS/Paper_Stadler-delGiorgio_ISMEJ_2021
# Output of script: 2_analysis.R
# do we have several files per object? -> take newest version
# ASV CSS transformed table
otu.tab <- select_newest("./Data/Microbial", "201520162017_fin_css_otu99_phantomcor_no.dups_")
otu.tab <- as.matrix(read.csv(
  otu.tab,
  sep = ";",
  dec = ".",
  row.names = 1,
  stringsAsFactors = F
))
#otu.tab <- otu.tab[, 1:1000]

# row orders need to match between tax.tab and otu.tab
otu.tab <- otu.tab[order(row.names(otu.tab)),]

# Output of first paper. Workflow in https://github.com/CarBBAS/Paper_Stadler-delGiorgio_ISMEJ_2021
# Output of script: 1_correctByBin_css.R
# Taxonomy table
tax.tab <- select_newest("./Data/Microbial", "201520162017_tax_otu99_table_no.dups_")
tax.tab <-
  as.matrix(read.csv(
    tax.tab,
    sep = ";",
    dec = ".",
    row.names = 1,
    stringsAsFactors = F
  ))
# orders need to match between tax.tab and otu.tab
tax.tab <- tax.tab[order(row.names(tax.tab)),]
#tax.tab <- tax.tab[,1:1000]

# Output of first paper. Workflow in https://github.com/CarBBAS/Paper_Stadler-delGiorgio_ISMEJ_2021
# Output of script: 1_correctByBin_css.R
# Meta data
met.df <-
  select_newest(path = "./Data/Microbial", file.pattern = "201520162017_meta_otu99_no.dups_")
met.df <-
  read.csv(
    met.df,
    sep = ";",
    dec = ".",
    row.names = 1,
    stringsAsFactors = F
  ) %>% setDT(keep.rownames = "row.names")


## Sample renaming ---------------------------------------------------------------------------------------
# add travel time
meta.all <- read.csv("./Data/Summary/meta_file_traveltime_2024-02-10.csv", sep = ",", stringsAsFactors = F) %>% setDT()
#ft.df <- read.csv("./Data/Summary/meta_FTdataset_2022-05-31.csv", sep = ",", dec = ".", stringsAsFactors = F) %>% setDT()

# omit sediment samples
met.df <- met.df[sample.type.year != "Sediment",]
met.df <- met.df[year <= 2016,]
met.df <- met.df[campaign <= 2,]

# omit estuary
met.df <- met.df[sample.type.year != "Marine",]

meta.all[, mic.name := sample.name]

meta.all[str_detect(mic.name, "^RO"), mic.name := str_replace(mic.name, "[-]", "")]
meta.all[str_detect(mic.name, "^RO"), mic.name := str_replace(mic.name, "[_]", ".")]
meta.all[str_detect(mic.name, "hypo"), mic.name := str_replace(mic.name, "hypo", ".Hypo")]
meta.all[str_detect(mic.name, "^T49"), mic.name := str_replace(mic.name, "T49", "TR49")]
meta.all[str_detect(mic.name, "^T5"), mic.name := str_replace(mic.name, "T5", "TR5")]
meta.all[str_detect(sample.name, "^L32"), mic.name := str_replace(mic.name, "[_]", "")]
meta.all[str_detect(sample.name, "^L33"), mic.name := str_replace(mic.name, "[.]", "")]
meta.all[sample.name == "RO2-111hypo", mic.name := "RO2111.90m"]
meta.all[sample.name == "T49_C2", mic.name := "TR49"]
meta.all[sample.name == "L330_C1", mic.name := "L330.C1"]
meta.all[sample.name == "L330_C2", mic.name := "L330"]
meta.all[sample.name == "LR09.1", mic.name := "LR09"]
meta.all[sample.name == "SWPR02.1", mic.name := "SWPR02"]

met.df$dr_match_name[!(met.df$dr_match_name %in% meta.all$mic.name)]

setDT(met.df)
# merge travel time to mic meta data
met.df[meta.all, c("Season","pointid", "strahler.order","flacc_km2",
                   "vel_ms","dis_m3s","wrt_min", "travel.time_d") := # "velocity_ms", 
         list(i.Season, i.pointid, i.strahler.order, i.flacc_km2,
              i.vel_ms, i.dis_m3s,  i.wrt_min, i.travel.time_days), on = c("dr_match_name==mic.name")] #i.velocity_ms,



#temp[,mic.samples := str_replace(sample.name, "[-]","")]
#temp[,ft.samples := str_replace(sample.name, "[-]",".")]

# setDT(met.df)
# met.df <- met.df[temp, c("flacc_npx") := list(i.flacc_npx), on = .(sample.name)]
setDF(met.df, rownames = met.df$row.names) # convert back to data frame
#row.names(met.df) <- samples # add row names for phyloseq

## Define habitats -----------------------------------------------------------------------------------------
# merge some sample types
met.df$sample.type.year <- factor(met.df$sample.type.year, levels = c("Soil","Sediment",
                                                                      "Soilwater","Hyporheic",
                                                                      "Wellwater","Stream", "Tributary",
                                                                      "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
                                                                      "Upriver", "RO3", "RO2", "RO1","Deep", "Downriver",
                                                                      "Marine","Unknown"),
                                  labels = c("Soil","Sediment",
                                             "Soilwater","Soilwater",
                                             "Groundwater","Stream", "Tributary",
                                             "Riverine \nLakes", "Headwater \nPonds", "Lake", "Lake",
                                             "Upriver",
                                             "Reservoirs","Reservoirs", "Reservoirs","Reservoirs", "Downriver",
                                             "Estuary","Unknown"))

met.df$Season <- factor(met.df$Season, levels = c("Spring","Summer","Autumn"),
                        labels = c("Spring","Summer","Autumn"))

met.df$dna_type <- factor(met.df$dna_type, levels = c("DNA","cDNA"),
                          labels = c("DNA","RNA"))

# save coordinates
out <- met.df %>% dplyr::select(seq_name, dr_match_name, dna_type, season, year, sample.type.year, lat, long)
write.table(out, "./Output/coordinates_MO.csv", sep = ",", row.names = F)

## Converge all data and prepare for processing -------------------------------------------------------
# Construct phyloseq object
pb <- phyloseq(otu_table(otu.tab[row.names(otu.tab) %in% row.names(met.df),], taxa_are_rows = F),
               sample_data(met.df),
               tax_table(tax.tab))

pb <- prune_taxa(!taxa_sums(pb) == 0, pb)
pb <- prune_samples(!sample_sums(pb) == 0, pb)

 #saveRDS(pb, './Objects/cleaned.mic_phyloseq.rds')
 # write.table(met.df, './Data/Microbial/cleaned.mic_meta.df.csv', sep = ",", dec = '.', row.names = F)

# format into long
# melt community matrix for tidy data set
commat <- melt.data.table(
  setDT(as.data.frame(otu_mat(pb)), keep.rownames = "Sample"), #otu_mat(pb)
  #scale(otu_mat(pb), center = F, scale = T))
  id.vars = "Sample",
  measure.vars = patterns("^OTU_"),
  variable.name = "OTU",
  value.name = "reads"
)

# Some stats
minhead(otu_mat(pb))
nrow(otu_mat(pb)) # 410 samples
ncol(otu_mat(pb)) # 10589 OTUs

# merge back with meta data
mic <- commat[setDT(met.df %>% dplyr::select(seq_name, dr_match_name, replicate, dna_type, seq_depth, year,
                                             Season, sample.type.year, sample.name, flacc_km2, dis_m3s,
                                             wrt_min,
                                             travel.time_d), keep.rownames = "Sample"),
              , on = .(Sample)]

# filter out some irrelevant samples
#mic <- mic[year <= 2016 & seq_depth == "Shallow" & dna_type != "RNA" & Season != "Autumn",]
# remove some wrongly snapped samples
#mic <- mic[flacc_npx != -9999,]

# update flacc_npx (= number of pixel) to catchment area
# multiply by pixel area and convert in km (pixel in metres)
#mic[,catch.area_km2 := (flacc_npx *(14*14))* 0.000001]
# we have to separate a few samples that are not a direct part of the continuum
# create a big gap in the x-axis for clear seperation of strong terrestrially influenced and estuary sites
# mic[sample.type.year == "Groundwater", catch.area_km2 := -5000]
# mic[sample.type.year == "Soilwater", catch.area_km2 := -5000]
# mic[sample.type.year == "Hyporheic", catch.area_km2 := -2000]
# mic[sample.type.year == "Estuary", catch.area_km2 := 20000]

## Remove outliers -------------------------------------------------------------------------------------------
# detect extreme group outliers, except zeros
mic[reads > 0, ext_outlier := is_extreme(reads), by = .(year, Season, dna_type, OTU)]

# sanity check
any(mic[is.na(ext_outlier),]$reads > 0)
# overwrite all zeros with FALSE
mic[reads == 0, ext_outlier := FALSE]

# check how many outliers
nrow(mic[ext_outlier == TRUE,]) * 100 / nrow(mic) # 0.53% outliers

# remove extreme outliers
clean.df <- mic[ext_outlier == FALSE,]

## Separate absent and present MF ---------------------------------------------------------------------------
# Define absence
clean.df[, ID := paste(year, Season, dna_type, OTU, sep = "_")] # year,
clean.df[, n := .N, by = .(ID)] # number of samples per category
clean.df[, n.obs := nrow(.SD[reads > 0,]), by = .(ID)] # number of samples with actual observation

#initiate absence column
clean.df[n.obs == 0, PA := "Absent", by = .(ID)]

#how many microbes that are absent?
length(levels(factor(clean.df[PA == "Absent",]$OTU))) # 94258
length(levels(factor(clean.df$OTU))) # 10589
# percentage lost
length(levels(factor(clean.df[PA == "Absent",]$OTU))) * 100 / length(levels(factor(clean.df$OTU)))
# 87.43% absent
# this doesn't mean that all OTUs never appear in this dataset. This means that only 11% appear in all years, campaigns
# and DNA vs RNA

# Overwrite the rest with present
clean.df[is.na(PA), PA := "Present"]
# extract only df with present OTUs
present <- clean.df[PA == "Present",]
levels(factor(present$sample.type.year))

## Define travel time for non estimated habitats ----------------------------------------------------------
present[sample.type.year == "Soil", travel.time_d := -50]
present[sample.type.year == "Soilwater", travel.time_d := -50]
present[sample.type.year == "Groundwater", travel.time_d := -50] #0.0000001
#present[sample.type.year == "Estuary", travel.time_d := 600] # 365

## Transform variables and remove zero observations --------------------------------------------------------
#present[, log.travel.time_d := log(travel.time_d)]
present[reads== 0, reads := NA]
# remove weird sample types
present <- present[sample.type.year != "Unknown", ]
present[, ID := paste(year, Season, dna_type, OTU, sep = ".")] #year,
# arrange so that the data frame follows travel time
present <- present %>% group_by(ID) %>% dplyr::arrange(travel.time_d)
#remove samples with no travel time or season
setDT(present)
present <- present[!is.na(travel.time_d),]
present <- present[!is.na(Season),]

absent <- clean.df[PA == "Absent",]

# save
# saveRDS(clean.df, paste0("./Objects/allMO_long_",Sys.Date(),".rds"))
# saveRDS(present, paste0("./Objects/presentMO_long_", Sys.Date(), ".rds"))
# saveRDS(absent, paste0("./Objects/absentMO_long_", Sys.Date(),".rds"))

# Check WRT - Travel time relationship ------------------------------------------------------------------
#' plot.df <- present %>% dplyr::select(Sample, sample.type.year,
#'                                      travel.time_d, wrt_min) %>% unique() %>% filter(sample.type.year != "Unknown") %>%
#'   filter(wrt_min > 0.002224513)
#' # define habitat types
#' plot.df$sample.type.year <- factor(plot.df$sample.type.year, levels = c("Soil",
#'                                                                         "Soilwater",
#'                                                                         "Groundwater","Stream", "Tributary",
#'                                                                         "Riverine \nLakes", "Headwater \nPonds", "Lake",
#'                                                                         "Upriver",
#'                                                                         "Reservoirs", "Downriver",
#'                                                                         "Estuary"))
#' # export factors for colouring
#' sample.factors <- levels(plot.df$sample.type.year)
#'
#' # more colour blind friendly
#' colvec <- c("#FCFDBFFF", #"#FEC589FF", #Soil
#'             #"#FEC589FF", #"#FDEBACFF", #Sediment
#'             "#F9795DFF", #"#F9795DFF", #Soilwater
#'             "#DE4968FF", #"#DE4968FF", #Groundwater,
#'             "skyblue", #Stream
#'             "#AD347CFF",# Tributary,
#'             "palegreen", #Riverine Lakes,
#'             "#7AD151FF", #Headwater Ponds,
#'             "#FDE725FF",# Lake,
#'             "#1F9F88FF", # Upriver,
#'             "orchid", #"#471063FF", #Reservoir,
#'             #'salmon',
#'             #"pink",
#'             #'navy',
#'             "#375A8CFF", #Downriver,
#'             "gray40") #Estuary)
#' names(colvec) <- as.character(sample.factors) # assign names for later easy selection
#'
#' # set wake WRT for terrestrial and estuary
#' plot.df[sample.type.year == "Soil", wrt_min := 0.006]
#' plot.df[sample.type.year == "Soilwater", wrt_min := 0.008]
#' plot.df[sample.type.year == "Groundwater", wrt_min := 0.01]
#' plot.df[sample.type.year == "Estuary", wrt_min := 5000]
#'
#' (p <- ggplot(plot.df, aes(x = log(travel.time_d), y = log(wrt_min))) +
#'   theme_bw() +
#'   geom_point(aes(colour = sample.type.year), size = 3) +
#'   scale_colour_manual(values = colvec[names(colvec) %in% as.character(levels(plot.df[!is.na(travel.time_d),]$sample.type.year))],
#'                       name = "Habitat types") +
#'     labs(x = "Travel time (days; log-scale)",
#'          y = "Water residence time (min; log-scale)") +
#'   theme(panel.grid = element_blank(),
#'         strip.background = element_rect(fill = "gray20"),
#'         axis.text = element_text(size = 14),
#'         axis.title = element_text(size = 16)))
#'
#' ggsave("./Figures/travel.time_vs_wrt.png", p,
#'        height = 12, width = 15, units = "cm")



## Z-Scale -------------------------------------------------------------------------------------------------
temp <- dcast(present, Sample ~ OTU, value.var = "reads")
#temp <- temp[, (colnames(temp)[-1]) := lapply(.SD, replace_na, 0), .SDcols = colnames(temp)[-1]]
scaled <- cbind(temp[,1], scale(temp[,-1], center = F, scale = T))

# merge back to present
temp <- melt(scaled, id.vars = "Sample", variable.name = "OTU", value.name = "z.reads")
present <- present[temp, , on = .(Sample, OTU)][!is.na(reads),]

# # Create randomized dataset ---------------------------------------------------------------------------------
# random.ls <- list()
# begin <- c(1, 501)
# end <- c(500, 1000)
# for(j in 1:length(end)) {
#   random.mat <- temp[, 1:2]
#   cols <- colnames(scaled)[-1]
#
#   print(paste("Round", j, sep = " "))
#   pb <- txtProgressBar(min = 0, max = 500, style = 3)
#   for (i in begin[j]:end[j]) {
#     setTxtProgressBar(pb, i)
#     if(j == 2){
#       set.seed(i+50)
#     } else {
#       set.seed(i)
#     }
#
#     r <- copy(scaled)
#     r <-
#       scaled[, (cols) := lapply(.SD, sample, size = nrow(scaled), replace = F), .SDcols = cols]
#     r <-
#       melt.data.table(
#         r,
#         id.vars = "Sample",
#         variable.name = "OTU",
#         value.name = paste0("random_", i)
#       )
#     #random.mat <- random.mat[r, , on = .(Sample, OTU)]
#     random.ls[[i]] <- r
#     Sys.sleep(time = 1)
#   }
#   close(pb)
#   saveRDS(random.ls, paste0("./Output/random.matrix_mic_run", j,".rds"))
#   rm(random.mat, r)
# }

## Filter too rare MOs -------------------------------------------------------------------------------------
# only keep those that have more or equal then 7 observations
length(unique(levels(factor(present$ID)))) #28814
hist(present$n.obs, breaks = seq(0, max(present$n.obs) + 10, 1))
# decided based on histogram, n.obs reaches plateau at 7
present <- present[n.obs >= 7,]
length(unique(levels(factor(present$ID)))) #18471

saveRDS(present, "./Objects/MO_present.rds")

## Randomization approach --------------------------------------------------------------------------------
# For each ID, shuffle the peak intensity along the x-axis and run regression 999 times

# this is actually a mini version of the full function.
# run script 5.1_MO_randomize.R on Compute Canada

# x <- present[ID %in% levels(as.factor(present$ID))[320],]
# slope.random <- ddply(x, .(ID), function(x) {
#   
#   out <- foreach(i = 1:999, .combine = rbind) %dopar% {
#     # set random iteration seed
#     set.seed(i)
#     # extract data
#     df <- data.frame(x = x$travel.time_d, y = x$z.reads)
#     # shuffle y over x
#     df$y <- sample(df$y, replace = F, size = nrow(df))
#     # do linear regression
#     lin <- try(lm(df$y ~ df$x))
#    
#     if (!inherits(lin, "try-error") & !is.na(lin$coefficients[[2]])) {
#     # extract run and slope
#     out <- data.frame(run = i, slope = lin$coefficients[2],
#                       r2 = summary(lin)$r.squared,
#                       p.val = summary(lin)$coefficients[2,4])
#     } else {
#       out <- data.frame(run = i, slope = NA)
#     }
#     return(out)
#   }
#   
#   return(out)
# }, .parallel = T)


present <- readRDS("./Objects/MO_present_2024-02-10.rds") %>% setDT()
# read in results
out <- readRDS("./Objects/MO_randomization_output_2024-02-10.rds") %>% setDT()
ci <- readRDS("./Objects/MO_randomization_CI_2024-02-10.rds") %>% setDT()
ci <- ci[order(ID),]

# # try another way
# t.test(out[ID == "2015.Spring.DNA.OTU_1",]$slope)
# ci[ID == "2015.Spring.DNA.OTU_1" & vars == "slope",] # two methods are the same

# re-calcuate CI, the one we calculate is based on t-score, which only works with sample sizes smaller than 30
# calculate z-score instead

# ci[, alpha := 0.01]
# ci[, t.score := qt(p=alpha/2, df = df, lower.tail = F)]
# ci[, margin.error := t.score * se]
# ci[, c("upper.bound", "lower.bound") := list(mean + margin.error,
#                                              mean - margin.error)]

# try another way


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
length(not.working) # 532

# constant but above 0
const <- unique(ci[sd == 0 & mean > 0,]$ID)
length(const) # 1387

# remove those IDs that did not work
present <- present[!(ID %in% not.working),]

set.seed(3)
exmpl <- sample(unique(present$ID), 20)

for(i in 1:length(exmpl)){
  # correlation plot
  cor <- ggplot(present[ID == exmpl[i],], aes(x = travel.time_d, y = z.reads)) +
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
  ggsave(paste0("./Figures/Randomization/MO_", exmpl[i],".png"), arr,
         width = 10, height = 15, units = "cm")
}

# Apply to all MOs -------------------------------------------------------------------------------------
# for applying function see script 5.2

model.ls <- readRDS("./Objects/MO_model.ls_binned_2024-02-10.rds")
model.df <- readRDS("./Objects/MO_model.df_binned_2024-02-10.rds") %>% setDT()
pks.df <- readRDS("./Objects/MO_pks.df_binned_2024-02-10.rds") %>% setDT()

#raw.df <- ldply(model.ls, function(x){
#  x$raw
#}) %>% setDT()

# Clean --------------------------------------------------------------------------------------------
model.df[is.na(p.val), AU := "n.s."]

nrow(model.df[is.na(p.val),]) * 100 / length(model.ls) # we loose 20.1%

model.df[model == "y ~ x", model.type := "LM"]
model.df[model == "y ~ poly(x, 3)" | model == "y ~ poly(x, 2)", model.type := "POLY"]
model.df[model == "y ~ s(x)", model.type := "GAM"]
model.df[is.na(model.type),]

# count numbers of model by type
model.df[, .(n = .N), by = .(model.type)] # 2237 NA models

model.df <- model.df %>% separate(col = "ID", into = c("year","season","dna.type","OTU"), sep = "[.]", remove = F) %>% setDT()

# Add randomization --------------------------------------------------------------------------------
model.df[ci[vars == "slope",], c("rand.ci.up", "rand.ci.low") := list(i.upper.bound, i.lower.bound), on = .(ID)]

length(unique(model.df$ID)) # 17939 models

# direction is defined 1 == positive, 2 == negative
# see how many pass the randomization filter
model.df[direction == 1 & (lm.slope > rand.ci.up), r.test := TRUE]
model.df[direction == 2 & (lm.slope < rand.ci.low), r.test := TRUE]
model.df[is.na(r.test), r.test := FALSE]

model.df[, .(n = .N), by = .(r.test, dna.type)] # 11308 DNA and 3968 RNA pass the test

# see how many pass p-value filter
model.df[p.val < 0.05, .(n = .N)] # 1949
model.df[p.val < 0.01, .(n = .N)] # 672

# see how many pass both
model.df[p.val < 0.05 & r.test == TRUE, .(n = .N), by = .(dna.type)] # D 1483, R 451
model.df[p.val < 0.01 & r.test == TRUE, .(n = .N), by = .(dna.type)] # D 535, R 131

# assign flags
model.df[, loose.filt := ifelse(r.test == TRUE, TRUE, FALSE)]
model.df[, cons.filt := ifelse(r.test == TRUE & p.val < 0.05, TRUE, FALSE)]

# Define a few truly flat models ------------------------------------------------------------------------
model.df[cons.filt == TRUE & sd.y == 0, ] # none
model.df[sd.y == 0, AU := "flat"]

# classify not significant
model.df[is.na(p.val), ns.s := "n.s."]
#model.df[p.cor == "remove", ns.s := "n.s."]

## Add ranges -----------------------------------------------------------------------------------------------
# Molecules that occur along the whole continuum
limits <- present[reads != 0, .(min.x = min(travel.time_d, na.rm = T),
            max.x = max(travel.time_d, na.rm = T)), by = .(ID)]
limits <- limits[ID %in% names(model.ls),]

model.df <- model.df[limits, , on = .(ID)] #c(".id==ID")

# model.df[min.x <= -15 & max.x >= 5, range := "entire"]
# # Molecules that only occur in freshwater
# model.df[min.x >= -15, range := "freshwater"]
# # Molecules that do not occur in estuary
# model.df[min.x <= -15 & max.x < 6, range := "soil-fresh"]


## Define model.type ------------------------------------------------------------------------
# model.df[model == "y ~ x", model.type := "LM"]
# model.df[model == "y ~ poly(x, 2)", model.type := "poly2"]
# model.df[model == "y ~ poly(x, 3)", model.type := "poly3"]
# model.df[model == "y ~ s(x)", model.type := "GAM"]
model.df[is.na(model.type), model.type := "none"]
# all defined?
model.df[is.na(model.type),] # yes

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
                 measure.vars = c("1","2","3","4","5"), 
                 variable.name = "peak", value.name = "y" ) %>% arrange(ID)
pks.diff <- pks.diff[!is.na(y),]

# calcualte the difference between the first and last peak
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
# sanity check
model.df[is.na(ns.s),] #none
model.df[is.na(c.ns.s),] #none

model.df[model.type == "none", AU := "flat"]

# Save results ----------------------------------------------------------------------------------------------------
unique(model.df$AU)
model.df[, AU := factor(AU, levels = c("n.s.", "flat",
                                       "increase","non-linear increase", "decreasing quadratic","decreasing bimodal",
                                       "unimodal",
                                       'decrease','non-linear decrease', 'increasing quadratic','increasing bimodal'))]

model.df <- model.df %>% dplyr::select(ID:p.val, min.x:buffer.max, rand.ci.up:cons.filt, model.type, hump, AU, ns.s, c.ns.s)

write.table(model.df,paste0("./Data/Summary/MO_patterns.classified_", Sys.Date(), ".csv"), sep = ",", row.names = F)


# Extra code ---------------------------------------------------------------------------------------------
# # Get a grasp on proportion of each AU -----------------------------------------------------------------
# n.au <- model.df[order(AU), .(n = .N), by = .(AU)]
# 
# # get groups of AU
# groups <- as.vector(n.au$AU)#[-1]
# 
# folders<-list.dirs("./Figures/MO_lm_fit/", full.names = F)[-1]
# groups<-as.vector(groups[!(groups %in% folders)])
# 
# # multiple models in one plot -----------------------------------------------------------------------
# set.seed(3)
# for(i in 1:length(groups)){
#   dir.create(paste0("./Figures/MO_lm_fit/", groups[i]))
#   print(paste0("Plotting models of group: ", groups[i]))
#   # get all IDs from the same group
#   temp <- model.df[AU == groups[i],]$ID
# 
#   if(length(temp) > 160){
#     temp <- sample(temp, 160)
#   }
#   if(groups[i] == "flat"){
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
#              y = "Reads (z-scaled)")
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
#       newd <- data.frame(x =seq(min(x$raw$x),max(x$raw$x), by = 1))
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
#              y = "Reads (z-scaled)")
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
#       ggsave(paste0("./Figures/MO_lm_fit/", groups[i],"/models_", plot.vec[j], "-", (plot.vec[j+1]-1),".pdf"), p,
#              height = 30, width = 35, units = "cm")
#     } else {
#       p <- do.call(grid.arrange, named[plot.vec[j]:length(named)])
#       ggsave(paste0("./Figures/MO_lm_fit/", groups[i],"/models_", plot.vec[j], "-", length(named),".pdf"), p,
#              height = 30, width = 35, units = "cm")
#     }
#   }
# }


# range(df[!is.na(diff.x),]$diff.x)
#
# df[no.peak >= 2, max.peak := .SD[which.max(peak.y),]$peak, by = .(.id)]
# df[is.na(AU) & peak == max.peak, diff.inf := infp.slope - shift(infp.slope, 1), by = .(.id)]
# flat <- df[is.na(AU) & peak == max.peak & abs(diff.inf) < 0.01, ]$.id
#
# decay <- df[is.na(AU) & peak == max.peak & closest == "below",][peak.x < -8 & direction == 2 & infp.slope < 0,]$.id
# increase <- df[is.na(AU) & peak == max.peak & closest == "before",][peak.x > -2 & direction == 1,]$.id
# hump <- df[is.na(AU) & peak == max.peak & closest == "below",][peak.x <= -2 & peak.x >= -8,]$.id

# Classify peakers
# one peakers have decaying and increasing fractions, classify by peak location
# gam.decay <- df[closest == "above" & is.na(infp.x),][model.type == "GAM" & no.peak == 1 & peak.x < -10,]$.id
# gam.increase <- df[closest == "below" & is.na(infp.x),][model.type == "GAM" & no.peak == 1 & peak.x > 0,]$.id
#
# df[.id %in% gam.decay & p.val < 0.05, AU := "decay"]
# df[.id %in% gam.increase & p.val < 0.05, AU := "increase"]
# df[.id %in% c(gam.decay, gam.increase) & p.val >= 0.05, AU := "stable"]
#
# # Classify the rest as dynamic (= peaking at)
# df[model.type == "GAM" & no.peak == 1 & p.val < 0.05 & is.na(AU), AU :=  "dynamic"]
# df[is.na(AU) & no.peak > 1 & p.val < 0.05 & model.type == "GAM", AU := "dynamic"]
# df[model.type == "GAM" & is.na(AU) & p.val >= 0.05, AU := "stable"]
# df[model.type == "GAM" & is.na(AU) & is.nan(p.val), AU := "stable"]
#
# # sanity check
# df[is.na(AU),]

# split seasons
#model.df <- model.df %>% separate(col = ".id", into = c("Year","Season","dna_type","OTU"), sep = "[.]", remove = F)
# setDT(model.df)
# # overall percentage
# model.df[, n := .N]
# model.df[, n.au.all := .N, by = .(AU)]
# model.df[, perc.all := (n.au.all * 100) / n]
# 
# # overall percentage
# model.df[, n.camp := .N, by = .(season,dna.type)]
# model.df[, n.au.camp := .N, by = .(AU, season, dna.type)]
# model.df[, perc.camp := (n.au.camp * 100) / n.camp]
# 
# # export table
# out <- model.df %>% dplyr::select(season, dna.type,
#                             AU, perc.all, n.au.all, perc.camp) %>% distinct()
# setDT(out)
# # out[,AU := factor(AU, levels = c("stable", "unimodal",
# #                                  "increase", "increasing_bimodal", "increasing_multimodal",
# #                                  "decay", "declining_bimodal", "declining_multimodal"),
# #                   labels = c("stable", "unimodal",
# #                              "increase", "increase", "increase",
# #                              "decay", "decay", "decay"))]
# out[, season := factor(season, levels = c("Spring", "Summer"))]
# # save to merge with FT data
# write.table(out, "./Data/Summary/AU_perc_overall_mic.csv", sep = ",", row.names = F)

# sub <- sample(df[AU == "dynamic",]$.id, 50)
# ggplot(present[ID %in% sub,], aes(x = log.travel.time_d, y = reads)) +
#   geom_point() +
#   geom_smooth(method = "gam") +
#   facet_wrap(~ID, scale = "free_y") +
#   guides(colour = "none")

# add growth and decay rate -----------------------------------------------------------------------------
# rate.travtime <- function(x, y, min.x, max.x, break.n){
#   # remove zeros
#   temp <- data.frame(x = x, y = y)
#   #temp <- temp[!is.na(temp$y),]
#   #temp <- temp[temp$y > 0,]
# 
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
# # apply on all OTUs
# 
# lin.ls <- dlply(present, .(ID), function(x) {
#   min.x <- minmax.trav[Season == unique(x$Season) & year == unique(x$year),]$min.trav
#   max.x <-  minmax.trav[Season == unique(x$Season) & year == unique(x$year),]$max.trav
# 
#   model.ls <- rate.travtime(x$log.travel.time_d, x$reads, min.x, max.x, break.n = 3)
# 
#   return(model.ls)
# 
# }, .parallel = T)
# 
# sub.ls <- lin.ls[!is.na(lin.ls)]
# length(lin.ls) - length(sub.ls) # we loose 21
# 
# lin.df <- ldply(lin.ls, function(x){
#   x
# })
# setDT(lin.df)
# # merge to general model data frame
# model.df <- model.df[lin.df, , on = c(".id==ID")]
# tax <- setDT(as.data.frame(tax.tab), keep.rownames = "OTU")
# tax[!is.na(species),] # species is empty
# tax[, species := NULL]
# 
# # extract a few columns
# 
# out <- model.df %>% dplyr::select(ID, year, season, dna.type, OTU, AU, model.type, mean.y, var.y, dev.expl, p.val, range,
#                                   lin.intercept:lin.p.val, n:perc.camp) %>% setDT()
# out <- out[tax, , on = .(OTU)]
# 
# write.table(out, "./Data/Summary/cross_AU_MO.csv", sep = ",", dec=  ".", row.names = F)

# Done ----------------------------------------------------------------------------------------------------------------

# Rest are figures for GRIL 2022 poster

# make a bar plot
# (p <- ggplot(out) +
#     theme_bw() +
#     facet_grid(.~Season) +
#     geom_col(aes(x = dna_type, y = perc.camp, fill = AU)) +
#     #geom_text(aes(y = -1.5, label= unique(out$n.by.range), x = unique(out$range)), colour = "gray40", size = 2) +
#     scale_fill_viridis_d(name = 'Spatial pattern', option = "inferno", direction = -1)+
#     labs(x = "Nucleic acid type", y = "Percentage of spatial patterns (%)") +
#     theme(panel.grid = element_blank(),
#           axis.text.x = element_text(angle = 45, hjust = 1)))
# ggsave("./Figures/mic_perc.range.pattern.season_DNARNA.png", p, width = 14, height = 10, units = "cm")
# out.all <- df %>% dplyr::select(AU, perc.all) %>% distinct()
# # continuum
# 
# 
# # merge with meta data
# mic <- readRDS("./Objects/mic_long.rds")
# 
# # arrange that row order follows catchment area
# # ID is year_campaign_molecular.formula
# mic <- mic %>% group_by(ID) %>% dplyr::arrange(catch.area_km2)
# setDT(mic)
# 
# # pick example spatial behaviours
# temp <- sample_n_by(df[range == "entire",], AU, size = 10)
# setDT(temp)
# temp.ls <- model.ls[names(model.ls) %in% temp$.id]
# 
# raw.data <- ldply(temp.ls, function(x){x$raw})
# frq.data <- ldply(temp.ls, function(x){x$frq.data})
# setDT(raw.data); setDT(frq.data)
# 
# raw.data[df, "AU" := i.AU, on = .(.id)]
# frq.data[df, "AU" := i.AU, on = .(.id)]
# 
# sel <- c("Spring.OTU_14060", # stable
#          "Spring.OTU_1696", # increase
#          "Summer.OTU_1904", #decay
#          "Spring.OTU_15561") #dynamic
# 
# 
# temp.ls <- model.ls[names(model.ls) %in% sel]
# raw.data <- ldply(temp.ls, function(x){x$raw})
# setDT(raw.data)
# 
# raw.data[df, c("AU") := list(i.AU), on = .(.id)]
# raw.data[,AU := factor(AU, levels = c("stable","increase","decay",
#                                       "dynamic"),
#                        labels = c("Stable","Increase","Decay",
#                                   "Dynamic"))]
# 
# # plot tests
# (p <- ggplot(raw.data, aes(x = x, y = y))+
#     theme_bw() +
#     facet_wrap(AU~.) +
#     #geom_point(aes(x = x, y = y), colour = "gray50", alpha = 0.5) +
#     geom_point(alpha = .5) +
#     geom_smooth(aes(colour = AU), formula = y ~ s(x), method = "gam", se = F) +
#     labs(
#       x = expression(paste("Catchment area [km"^2,"]")),
#       y = "OTU reads (z-scaled)") +
#     scale_colour_viridis_d(name = "Spatial\nbehvaiour", option = "inferno", direction = -1) +
#     theme(panel.grid = element_blank(),
#           axis.text.x = element_text(angle = 45, hjust = 1)) +
#     guides(colour = "none") +
#     theme(axis.title = element_text(size = 14),
#           strip.background = element_rect(fill = "gray20"),
#           strip.text = element_text(colour = "white", size = 14)))
# ggsave("./Figures/spatbehav_ex_mic.png", width = 13, height = 14, units = "cm")
# 
# 
# # Link both microbes and FT overall percentage
# mic <- read.csv("./Data/Summary/AU_perc_overall_mic.csv", sep = ",")
# ft <- read.csv("./Data/Summary/AU_perc_overall_ft.csv", sep = ",")
# setDT(mic); setDT(ft)
# ft[, season := factor(campaign, levels = c(1,2),
#                       labels = c("Spring", "Summer"))]
# ft[,AU := factor(AU, levels = c("stable","increase","decay",
#                                 "dynamic"),
#                  labels = c("Stable","Increase","Decay",
#                             "Dynamic"))]
# 
# # add category
# ft[, cat := "DOM"]
# mic[, cat := "DNA"]
# mic[ft, , on = .(AU, season)]
# 
# comb <- rbind(mic %>% select(cat, AU, perc.all) %>% distinct(), ft %>% select(cat, AU, perc.all) %>% distinct())
# setDT(comb)
# comb[,AU := factor(AU, levels =  c("Stable","Increase","Decay",
#                                    "Dynamic"))]
# 
# (p <- ggplot() +
#     theme_bw() +
#     geom_col(aes(x = cat, y = perc.all, fill = AU), colour = "gray20", data = comb) +
#     #geom_text(aes(y = -1.5, label= unique(out$n.by.range), x = unique(out$range)), colour = "gray40", size = 2) +
#     scale_fill_viridis_d(name = 'Spatial pattern', option = "inferno", direction = -1)+
#     labs(x = "Assemblage", y = "Percentage of patterns (%)") +
#     theme(panel.grid = element_blank(),
#           axis.text.x = element_text(size = 14),
#           strip.background = element_rect(fill = "gray20"),
#           strip.text = element_text(colour = "white", size = 14)))
# ggsave("./Figures/overall_perc_DOMDNA.png", p, height = 10, width = 12, units = 'cm')
# 
# 
# # seasonal
# comb <- rbind(mic %>% select(cat,season, AU, perc.camp) %>% distinct(),
#               ft %>% select(cat, season, AU, perc.camp = perc.by.camp) %>% distinct())
# setDT(comb)
# comb[,AU := factor(AU, levels =  c("Stable","Increase","Decay",
#                                    "Dynamic"))]
# 
# (p <- ggplot() +
#     theme_bw(base_size = 24) +
#     facet_wrap(season~.)+
#     geom_col(aes(x = cat, y = perc.camp, fill = AU), colour = "gray20", data = comb) +
#     #geom_text(aes(y = -1.5, label= unique(out$n.by.range), x = unique(out$range)), colour = "gray40", size = 2) +
#     scale_fill_viridis_d(name = 'Spatial pattern', option = "inferno", direction = -1)+
#     labs(x = "Pool", y = "Percentage of patterns (%)") +
#     theme(panel.grid = element_blank(),
#           strip.background = element_rect(fill = "gray20"),
#           strip.text = element_text(colour = "white"),
#           legend.position = "none"))
# ggsave("./Figures/camp_perc_DOMDNA_poster.png", p, height = 12, width = 15, units = 'cm')
# 
# 

# make groups factors
# pca.df <- df
# pca.df <- pca.df %>% select(.id, no.peak, model.type, direction, intercept, lm.slope, dev.expl, adj.r.sq, p.val, range)
# pca.df[, model.type := as.numeric(as.factor(model.type))]
# pca.df[, range := as.numeric(as.factor(range))]
#
# for(j in 1:ncol(pca.df)) set(pca.df, which(is.infinite(pca.df[[j]])), j, NA)
# for(j in 1:ncol(pca.df)) set(pca.df, which(is.na(pca.df[[j]])), j, 0)
#

# pca.df <- pca.df %>% unique()
# ordiplot(prcomp(pca.df[,-1]))
# d <- dist(pca.df[,-1])
# hcl <-hclust(d)
# plot(hcl)
# abline(h = 50, col = "red")

# k-means
# set.seed(3)
# ks <- 1:10
# kmeans(pca.df[,-1], centers = 5, nstart = 10)
# tot_within_ss <- sapply(ks, function(k) {
#   cl <- kmeans(pca.df[,-1], k, nstart = 10)
#   cl$tot.withinss
# })
# plot(ks, tot_within_ss, type = "b")
#
# k5 <- kmeans(pca.df[,-1], centers = 5, nstart = 10)
# k8 <- kmeans(pca.df[,-1], centers = 8, nstart = 10)
#
# plot(pca.df[,-1]$lm.slope, pca.df$adj.r.sq, col = k5$cluster)
# plot(pca.df[,-1]$lm.slope, pca.df$adj.r.sq, col = k8$cluster)
#
# pca.df$k5 <- k5$cluster
# pca.df$k8 <- k8$cluster

#df[pca.df, c('k5','k8') := list(i.k5, i.k8), on = .(.id)]

# merge peak info to model.df
#df <- pks.df[model.df, , on = .(.id)]

# calculate peak differences to find GAMs that are flat
#df[, diff.y.pks := peak.y - shift(peak.y, 1), by = .(.id, closest)]

# update lm peak number
# pca.df <- model.df
# pca.df <- pca.df %>% select(.id, model.type:lm.slope, mean.y:max.x)
#
# # make model choice as factor
# pca.df[, model.type := as.numeric(factor(model.type))] #1 = GAM
# pca.df <- na.omit(pca.df)
# pca.df <- pca.df[, lapply(.SD, function(x) replace(x, is.infinite(x), 10))]
#
# # k-means
# set.seed(3)
# ks <- 1:10
# tot_within_ss <- sapply(ks, function(k) {
#   cl <- kmeans(pca.df[,-1], k, nstart = 10)
#   cl$tot.withinss
# })
# plot(ks, tot_within_ss, type = "b") #4
#
# k4 <- kmeans(pca.df[,-1], centers = 4, nstart = 200)
# df <- data.frame(cluster = k4$cluster, pca.df)
# setDT(df)
#
# ggplot(df, aes(x = lm.slope, y = mean.y, colour = as.character(cluster))) +
#   geom_point(shape = 19)
# # df[model.type == "LM", no.peak := 1]
# # df[model.type == "LM", peak := 1]
# sub <- sample(df[cluster == 3,]$.id, 30)
# ggplot(present[ID %in% sub,], aes(x = log.travel.time_d, y = reads)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   facet_wrap(~ID, scale = "free_y") +
#   guides(colour = "none")
#
# ggplot(df, aes(y = var.y, x = mean.y)) +
#   geom_point()

# # remove those with extremely narrow ranges
# present[reads > 0, min.range := min(.SD$log.travel.time_d), by = .(ID)]
# present[reads > 0, max.range := max(.SD$log.travel.time_d), by = .(ID)]
# present[, range.span := max.range - min.range]
# hist(present$range.span)
# # we want at least 10
# hist(present$range.span)
# # set based on histogram, there are not that many below 10
# #present <- present[ID %in% present[range.span >= 10,]$ID,]
# #length(unique(levels(factor(present$ID)))) #17207
# 
# # include a filter for distribution of sample sizes across range
# 
# check.bin.n <- function(x){
#   steps <- unique(x$range.span) / 3
#   bins <- c(unique(x$min.range), unique(x$min.range) + (steps), unique(x$min.range) + ((steps)*2),
#             unique(x$max.range))
#   n.bins <- c(sum(x$log.travel.time_d >= bins[1] & x$log.travel.time_d < bins[2]),
#               sum(x$log.travel.time_d >= bins[2] & x$log.travel.time_d < bins[3]),
#               sum(x$log.travel.time_d >= bins[3] & x$log.travel.time_d <= bins[4]))
# 
#   return(any((n.bins >= 2) == F)) # if there are any bins with n < 2, return TRUE
# }
# 
# present <- present[, bin.check := check.bin.n(.SD), by = .(ID)]
# length(unique(present$ID)) #18126
# # apply filter
# present <- present[bin.check == FALSE,]
# length(unique(present$ID)) #10064
# 
# rm(temp, scaled)
# 
# # Slope filter --------------------------------------------------------------------------------------------
# # Run lm on all models and classify those that are slope == 0 as stable even before running model decision tree
# slope.df <- ddply(present, .(ID), function(x) {
#   lin <- lm(x$z.reads ~ x$log.travel.time_d)
#   df <- data.frame(slope = lin$coefficients[2])
#   return(df)
# }, .parallel = T)
# 
# # identify those that have zero slope
# setDT(slope.df)
# stable <- slope.df[abs(slope) == 0,]$ID
# length(stable) # 99
# present[ID %in% stable, AU := "stable"]
# 
# hist(present$log.travel.time_d) # very skewed and sparsely distributed
# 
# # Add rates based on non z-scaled data --------------------------------------------------------------------------------
# slope.log2 <- ddply(present, .(ID), function(x) {
#   lin <- lm(x$reads ~ log2(x$travel.time_d))
#   df <- data.frame(slope.log2 = lin$coefficients[2])
#   return(df)
# }, .parallel = T)
# 
# slope.z.log2 <- ddply(present, .(ID), function(x) {
#   lin <- lm(x$z.reads ~ log2(x$travel.time_d))
#   df <- data.frame(slope.log2.z = lin$coefficients[2])
#   return(df)
# }, .parallel = T)
# 
# stat.reads <- present[, .(mean.reads = mean(reads, na.rm = T),
#             var.reads = var(reads, na.rm = T),
#             median.reads = median(reads, na.rm = T)), by = .(ID)]
# 
# stat.reads <- stat.reads[slope.log2, , on = .(ID)]
# stat.reads <- stat.reads[slope.z.log2, , on = .(ID)]

# y ~ x
# response ~ predictor

# log of 2
# slope corresponds to an increase in y, when x is increased by two-fold
# doubling rate

# natural log interpretation
# y = a + b*ln(x) + e
# A 1% change in X is associated with a change in Y of 0.01*b
# Taken from Introduction to Econometrics from Stock and Watson, 2003, p. 215

# OR
# The intercept is the expected value of the response when the predictor is 1.
# The slope measures the expected change in the response when the predictor increases by a fixed percentage
# Reads = 4906 + 298.6 * log(1)
# 4906 + 298.6 * log(1) = 4906
# 4906 reads/d

# Other semi-log and log-log transformations:
# ln(Y)=B0 + B1*X + u ~ A change in X by one unit (∆X=1) is associated with a (exp(B1) - 1)*100 % change in Y
# ln(Y)=B0 + B1*ln(X) + u ~ A 1% change in X is associated with a B1% change in Y, so B1 is the elasticity of Y with respect to X.

# Testing and defining function ---------------------------------------------------------------------------
# Function to evaluate whether a molecule has a linear, or polynomial curve
# select a hundred per compound category
# dirty plotting
set.seed(3)

# present <- readRDS("./Objects/presentMO_long.rds")
#
# temp <- present %>% sample_n(1)
# test <- present[ID %in% temp$ID,]
#
# # set vectors
# x <- test$log.travel.time_d
# y <- test$z.reads
#
#
# #saveRDS(mic, "./Objects/mic_long.rds")
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
#
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
#
# # GAM seems best in terms of smoothness
# # remove unncesseary objects
# rm(newd, lin, loess.mod, poly.lin, poly.lin3, temp, test, frq.x, x, y, gam , gam.poly)


# test with a 100 OTUs
# temp <- unique(present[is.na(AU),]$ID) %>% sample(100)
# temp <- present[ID %in% temp,]
#
# #trouble-shooting
# # x <- present[ID %in% levels(factor(present$ID))[4685],]
# #
# x <- present[ID %in% levels(factor(present$ID))[3],]$log.travel.time_d
# y <- present[ID %in% levels(factor(present$ID))[3],]$reads
#
# x <- present[ID %in% sub[2],]$log.travel.time_d
# y <- present[ID %in% sub[2],]$reads

# x <- temp[ID %in% "2015_Spring_C30H29O19N1S0",]$log.travel.time_d
# y <- temp[ID %in% "2015_Spring_C30H29O19N1S0",]$reads

# # retrieve max and min of that season/year
# minmax.trav <- present[, .(min.trav = round(min(log.travel.time_d, na.rm = T),0),
#                            max.trav = round(max(log.travel.time_d, na.rm = T),0)), by = .(year, Season)]
#
# # retrieve best model fit
# model.ls <- dlply(temp, .(ID), function(x) {
#   min.x <- minmax.trav[Season == unique(x$Season) & year == unique(x$year),]$min.trav
#   max.x <-  minmax.trav[Season == unique(x$Season) & year == unique(x$year),]$max.trav
#
#   #z.mean <- z.scaling[year == x$year,]$overall.mean
#
#   model.ls <- lm.or.gam(x$log.travel.time_d, x$z.reads, min.x, max.x) #, break.n = 3
#
#   return(model.ls)
#
# }, .parallel = T)
#
# # # remove NA models (not enough observations = only 2)
# model.ls <- model.ls[!is.na(model.ls)]
#
# # plot tests
# plots <- lapply(model.ls, function(x){
#   newd <- data.frame(x =seq(min(x$raw$x),max(x$raw$x), by = 1))
#   newd$y.pred <- predict(x$model, newd)
#   p <- ggplot()+
#     theme_bw() +
#     geom_point(aes(x = x$raw$x, y = x$raw$y), colour = "gray50", alpha = 0.5) +
#     geom_line(aes(x = newd$x, newd$y.pred), colour = "purple") +
#     labs(title = paste("Model = ", x$model.df$model, "\nDirection = ", x$model.df$direction),
#          x = "Log Travel time (d)",
#          y = "Reads")
#
#   if(!is.null(x$pks.df)){
#     if(x$pks.df$no.peak[1] > 0){
#       p <- p +
#         geom_point(aes(x = x$pks.df$peak.x, y = x$pks.df$peak.y), colour = "purple", size = 4, alpha = 0.5) +
#         geom_point(aes(x = x$pks.df$infp.x, y = x$pks.df$infp.y, colour = x$pks.df$closest), size = 2, alpha = 0.5)
#     }
#   }
#   return(p)
# })
#
# # save plots
# for(i in 1:20){
#   ggsave(paste0("./Figures/MO_lm_fit/", names(plots)[i],".png"), plots[[i]])
# }

# looks good

# P-value correction --------------------------------------------------------------------------------
# # Check p-value distribution
# hist(model.df[p.val < 0.1,]$p.val)
# 
# # number of models that worked
# (n <- nrow(model.df[is.na(AU),])) # 9364
# nrow(model.df) # out of binned 10064
# 
# # out of which, how many are RNA?
# nrow(model.df[grep("RNA", ID),]) # 2202
# nrow(model.df[grep("DNA", ID),]) # 7862
# 
# rna <- model.df[grep("RNA", ID),]
# dna <- model.df[grep("DNA", ID),]
# (n <- nrow(dna[is.na(AU),]))  # 7211
# 
# # create p-value bins
# p.bins <-  c(min(dna$p.val, na.rm = T), 0.00001, 0.0001, 0.001, 0.01, 0.05, max(dna$p.val, na.rm = T))
# 
# # theoretical false-positives
# fp <- round(p.bins * n, 0 )
# 
# # number of models positive within each bin
# ranges <- paste(head(p.bins,-1), p.bins[-1], sep=" - ")
# freq <- hist(dna[is.na(AU),]$p.val, breaks = p.bins, include.lowest = T, plot = F)
# p.freq <- data.frame(range = ranges, frequency = freq$counts, freq.fp = fp[-1]) %>% setDT()
# p.freq[, freq.diff := frequency - freq.fp]
# p.freq[, p.bins := c(1:5,NA)]
# 
# # sanity check
# sum(p.freq$frequency) == n # should be TRUE
# 
# sum(p.freq[-nrow(p.freq),]$freq.diff)
# # we keep 524 models this way
# 
# # Identify models to keep
# dna <- dna[order(rank(p.val)),]
# # separate models that are not included in this
# out.df <- dna[!is.na(AU),]
# 
# dna <- dna[is.na(AU),]
# dna[p.val <= 0.00001, p.bins := 1]
# dna[p.val > 0.00001 & p.val <= 0.0001, p.bins := 2]
# dna[p.val > 0.0001 & p.val <= 0.001, p.bins := 3]
# dna[p.val > 0.001 & p.val <= 0.01, p.bins := 4]
# dna[p.val > 0.01 & p.val <= 0.05, p.bins := 5]
# 
# # also add FDR correction
# dna[, q.val := p.adjust(p.val, method = "fdr")]
# nrow(dna[q.val < 0.05,]) # worse than our method, 2 models left
# 
# p.freq <- p.freq[!is.na(p.bins),]
# 
# # overwrite p.val of models that didn't pass
# for(i in 1:nrow(p.freq)){
#   temp <- dna[p.bins == p.freq$p.bins[i], head(.SD,p.freq$freq.diff[i])]$ID
#   dna[p.bins == p.freq$p.bins[i] & !(ID %in% temp), c("p.cor") :=
#                                           list("remove")]
# }
# 
# #sanity check
# nrow(dna[p.cor == "remove",]) == sum(p.freq$freq.fp) # should be TRUE
# 
# # ok we have all models that we want to keep
# 
# # bind dna and rna together
# model.df <- bind_rows(dna %>% dplyr::select(-c(p.bins, q.val)), rna)
# 
# # bind the good and "bad" models back
# model.df<-rbind(model.df, out.df %>% mutate(p.cor = "NA"))
# 
# # save corrected model.df
# saveRDS(model.df, paste0("./Objects/MO/model.df_mic_scaled.binned_p.corr_", Sys.Date(),".rds"))

# Classification ------------------------------------------------------------------------------------
# read objects
# model.ls <- readRDS(select_newest("./Objects/", "model.ls_mic_scaled"))
# pks.df <- readRDS(select_newest("./Objects/", "pks.df_mic_scaled")) %>% setDT()

#model.df <- readRDS(select_newest("./Objects/MO", "model.df_mic_scaled.binned_p.corr_")) %>% setDT()

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


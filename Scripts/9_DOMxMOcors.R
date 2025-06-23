## -------------------------------------------------------------------------
##
## Script name: 9_DOMxMOcors.R
##
## Purpose of script: 
##
## Author: Masumi Stadler
##
## Date Finalized: 
##
## Copyright (c) Masumi Stadler, 2025
## Email: m.stadler.jp.at@gmail.com
##
## -------------------------------------------------------------------------
##
## Notes:
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
pckgs <- list("plyr", "tidyverse", "data.table", # wrangling
              "RColorBrewer", "ggpubr" # plotting
) # statistics

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
# Get final melted correlation matrix 
cordf <- readRDS("./Objects/all_cors_2024-03-04.rds")

#\ Stats ------------------------------------------------------------------
# some stats
# how many correlations overall?
nrow(cordf) #101,280,418
# how many didn't return robust results?
nrow(cordf[is.na(pval),]) # 12,180,924

# how many MF and OTU?
length(unique(cordf$MF)) #13236
length(unique(cordf$OTU)) #6019

# how many cors by year and season?
stats.cor <- cordf[!is.na(pval), .(n = .N), by = .(Year, Season)]
# Year Season        n
# 2015 Spring  4537606
# 2015 Summer 23694454
# 2016 Spring 20473655
# 2016 Summer 40393779

# Format data ------------------------------------------------------------
# remove any NAs in pvalues, and any rho == 0
cordf <- cordf[!is.na(pval),]
cordf <- cordf[rho != 0,]
# how many left?
nrow(cordf) #54472074

# check distribution of correlation coefficients
hist(cordf$rho)
hist(cordf[pval < 0.05,]$rho)
nrow(cordf[abs(rho) == 1,]) # 3,104,575

cordf <- cordf[abs(rho) != 1,]
# apply FDR correction for multiple comparison
#cordf[, p.cor := p.adjust(pval, method = "fdr")]
#hist(cordf[p.cor < 0.05,]$rho)

# add categories do identify significant and ns. correlations according raw p-value and adjusted p-val
cordf[, sig.group := ifelse(pval > 0.05, "n.s.", "sig.")]
#cordf[, fdr.group := ifelse(p.cor > 0.05, "n.s.", "sig.")]

# number of correlations by year and season
cordf[, n.all := .N, by = .(Year, Season)]

#calculate percentage of significance group by raw p and adjusted p
stats.cor <- cordf[, .(n = .N,
                       n.all = unique(n.all)), by = .(Year, Season, sig.group)]
# stats.cor.fdr <- cordf[, .(n = .N,
#                            n.all = unique(n.all)), by = .(Year, Season, fdr.group)]
stats.cor[, perc := n * 100/ n.all]

# plot
ggplot(stats.cor, aes(x = paste(Year, Season), y = perc)) +
  geom_col(aes(fill = sig.group))


# reduce dataset to all correlations that are significant
cordf <- cordf[pval < 0.05, ]

rm(stats.cor, stats.cor.fdr)

#\ Merge with spatial pattern data -----------------------------------------------------------
# get spatial patterns categorization
mo <- read.csv("./Data/Summary/MO_patterns.classified_cleaned_2024-02-12.csv", sep = ",", stringsAsFactors = F) %>% setDT()
#ft <- read.csv("./Data/Summary/FT_patterns.classified_cleaned_2024-02-12.csv", sep = ",", stringsAsFactors = F) %>% setDT()

ft <- readRDS("./Objects/ft_cleaned_withclusters.rds")
#mo <- readRDS("./Objects/mo_cleaned_withphy.rds")

# Add taxonomy
phy <- readRDS("./Data/phyloseq_otu.tax.samp.phylo_SINA_chapter2.rds")
tax <- as.data.frame(phyloseq::tax_table(phy))[,1:5] %>% setDT(keep.rownames = "OTU")
mo[, domain := NULL]
mo <- mo[tax,, on = .(OTU)]

# Unify column names and factors
colnames(cordf)[5:6] <- c("year","season")
cordf[, c("year", "season") := list(factor(year, levels = c("2015","2016")),
                                    factor(season, levels = c("Spring","Summer")))]
mo[, c("year", "season") := list(factor(year, levels = c("2015","2016")),
                                 factor(season, levels = c("Spring","Summer")))]
ft[, c("year", "season") := list(factor(year, levels = c("2015","2016")),
                                 factor(season, levels = c("Spring","Summer")))]

# unify mo and dom sAU categories
mo[, sAU := factor(sAU, levels = c("increase", "non-linear increase", "multimodal increase",
                                   "unimodal",
                                   "multimodal decrease", "non-linear decrease", "decrease"),
                   labels =  c("increase", "non-linear increase", "non-linear increase",
                               "unimodal",
                               "non-linear decrease", "non-linear decrease", "decrease"))]
ft[, sAU := factor(sAU, levels = c("increase", "non-linear increase", "multimodal increase",
                                   "unimodal",
                                   "multimodal decrease", "non-linear decrease", "decrease"),
                   labels =  c("increase", "non-linear increase", "non-linear increase",
                               "unimodal",
                               "non-linear decrease", "non-linear decrease", "decrease"))]

# remove OTUs that dropped out during modelling exercise
mo <- mo[!is.na(season),]

# cordf[mo[dna.type == "DNA"], c("dna.nss", "dna.au", "dna.sau", "dna.slope", "dna.range", "mo.bin", "domain",
#                                "phylum", "class", "rna.rho", "rna.pval") := list(i.c.ns.s, i.minAU, i.sAU, i.lm.slope, i.tt.range,
#                                                                                  i.bin.cat,
#                                                                                  i.domain, i.phylum, i.class, i.dr.rho, i.dr.pval), on = .(season, year, OTU)]
# cordf[ft, c("dom.nss", "dom.au", "dom.sau", "dom.slope", "dom.range","dom.bin", "dom.clust", "dom.order","dom.cat",
#             "nosc", "dbeo", "mz","aimod","hc","oc") := list(i.c.ns.s, i.minAU, i.sAU, i.lm.slope, i.tt.range,
#                                                             i.bin.cat, i.cluster, i.order,
#                                                             i.categories, i.NOSC, i.DBE_O, i.mz, i.AImod, i.HC, i.OC), 
#       on = .(season, year, MF)]

cordf[mo[dna.type == "DNA"], c("dna.nss", "dna.au", "dna.sau", "dna.slope", "dna.range", "mo.bin", "domain",
                               "phylum", "class") := list(i.c.ns.s, i.minAU, i.sAU, i.lm.slope, i.tt.range,
                                                                                 i.bin.cat,
                                                                                 i.domain, i.phylum, i.class), on = .(season, year, OTU)]
cordf[ft, c("dom.nss", "dom.au", "dom.sau", "dom.slope", "dom.range","dom.bin", "dom.clust", "dom.order","dom.cat",
            "nosc", "dbeo", "mz","aimod","hc","oc") := list(i.c.ns.s, i.minAU, i.sAU, i.lm.slope, i.tt.range,
                                                            i.bin.cat, i.cluster, i.order,
                                                            i.categories, i.NOSC, i.DBE_O, i.mz, i.AImod, i.HC, i.OC), 
      on = .(season, year, MF)]

cordf[, dna.au := factor(dna.au,
                         levels = c('increase','unimodal','decrease'))]
cordf[, dna.end.eco := factor(dna.range, levels = c("soil-mouth","stream-mouth","soil-lake","stream-lake"),
                              labels = c("Mouth","Mouth","Lake","Lake"))]
cordf[, dna.origin := factor(dna.range, levels = c("soil-mouth","soil-lake","stream-lake","stream-mouth"),
                             labels = c("Terrestrial","Terrestrial","Aquatic","Aquatic"))]

cordf[, dom.au := factor(dom.au,
                         levels = c('increase','unimodal','decrease'))]
cordf[, dom.end.eco := factor(dom.range, levels = c("soil-mouth","stream-mouth","soil-lake","stream-lake"),
                              labels = c("Mouth","Mouth","Lake","Lake"))]
cordf[, dom.origin := factor(dom.range, levels = c("soil-mouth","soil-lake","stream-lake","stream-mouth"),
                             labels = c("Terrestrial","Terrestrial","Aquatic","Aquatic"))]

# Update MO with no phyla entry
cordf[is.na(phylum), phylum := "Unclassified"]

# remove some OTU and MFs that never had any significant correlations
cordf <- cordf[!is.na(dna.au) & !is.na(dom.au),]

# Sanity check: Do the correlations signs and direction of spatial patterns match? --------------
nrow(cordf) # 2,119,512, significant correlations

sig.cor <- copy(cordf)
all.cor <- copy(cordf)

# With detailed spatial pattern categories --------------------------------
#\ Significant dataset - aka reactive -------------------------------------
# Make a subset of those correlations of OTUs and MFs that have a sign. spatial pattern
sig.cor <- sig.cor[dna.nss == "sig." & dom.nss == "sig.",]
nrow(sig.cor) # 541,875

# add a category that separates positive and negative correlations
sig.cor[, rho.group := ifelse(rho > 0, "pos", "neg")]
sig.cor[, rho.group := factor(rho.group, levels = c("pos","neg"))]

# calculate percentage of positive vs negative correlations in each relation group
sig.cor[, all.n := .N, by = .(dom.sau, dna.sau)]
rho.df <- sig.cor[, .(n = .N,
                      all.n = unique(all.n)), by = .(dom.sau, dna.sau, rho.group)]
rho.df[, perc := n * 100/all.n]

# sanity check
rho.df[, .(sum = sum(perc)), by = .(dom.sau, dna.sau)]
sum(rho.df$n) == nrow(sig.cor)

# add dataset identifier
rho.df[, dataset := "sig"]

#\ With bulk dataset - aka bulk ---------------------------------------------
all.cor <- copy(cordf)
nrow(all.cor) # 2,368,307

all.cor[, rho.group := ifelse(rho > 0, "pos", "neg")]
all.cor[, rho.group := factor(rho.group, levels = c("pos","neg"))]

# calculate percentage
all.cor[, all.n := .N, by = .(dom.sau, dna.sau)]
df <- all.cor[, .(n = .N,
                  all.n = unique(all.n)), by = .(dom.sau, dna.sau, rho.group)]
df[, perc := n * 100/all.n]

# sanity check
df[, .(sum = sum(perc)), by = .(dom.sau, dna.sau)]
sum(df$n) == nrow(all.cor)

df[, rho.group := factor(rho.group, levels = c("pos","neg"))]
df[, dataset := "all"]

rho.df <- rbind(rho.df, df)
rho.df[, dataset := factor(dataset, levels = c("all","sig"), labels = c("Bulk", "Reactive"))]

(p <- ggplot(rho.df, aes(x = dataset, y = perc)) +
    theme_custom() +
    theme(legend.position = "bottom") +
    facet_grid(dom.sau~dna.sau) +
    geom_col(aes(fill = rho.group), width = .8) +
    geom_text(data = rho.df %>% dplyr::select(dom.sau, dna.sau, dataset, all.n),
              aes(x = dataset, y = 50, label = all.n), colour = "white", size = 3) +
    scale_fill_manual(values = c("tomato","skyblue"), name = "Correlation\ncoefficient", labels = c(expression(paste(rho, " > 0")),
                                                                                                    expression(paste(rho, " < 0")))) +
    labs(x = "Dataset", y = "Percent"))
p <- annotate_figure(p, top = "Microbial spatial pattern", right = "Molecular spatial pattern")

ggsave("./Figures/Analysis/DOMxMO_correlation_stacked_rho_sAU_byMF.png", p, height = 20, width = 20, units = "cm")



# Is there a difference in contribution to correlation categories by cluster? ----------------
head(cordf)
hist(cordf$rho)

cordf[, rho.group := ifelse(rho > 0, "pos", "neg")]
cordf[, rho.group := factor(rho.group, levels = c("pos","neg"))]

by.all <- copy(cordf)
by.all <- by.all[dna.nss == "sig." & dom.nss == "sig.",]
by.all <- by.all[dom.au != "unimodal" & dna.au != "unimodal",]
by.all[, id := paste(dom.au, dna.au, rho.group)]
by.all <- by.all[id %in% c('increase increase pos', 'decrease decrease pos',
                           'increase decrease neg', 'decrease increase neg'),]

by.all[, n := .N, by = .(rho.group, dna.au, season, year)]
rho.df <- by.all[, .(n.all = unique(n),
                     n = .N), by = .(dom.clust, rho.group, dna.au, season, year)]
rho.df[, perc := n * 100 / n.all]
# sanity check
rho.df[, .(sum = sum(perc)), by = .(rho.group, dna.au, season, year)]
rho.df[, variable := paste(dna.au, rho.group, sep = ".")]
rho.df[, variable := factor(variable, levels = c("increase.pos","increase.neg", "decrease.neg", "decrease.pos"),
                            labels = c("Increase - Positive", "Increase - Negative",
                                       "Decrease - Negative", "Decrease - Positive"))]
total.perc <- copy(rho.df)
# order.vec <- rho.df[rho.group == "pos" & dna.au == "increase" & dom.clust == "1",]
# order.vec <- order.vec[order(-perc),]$phylum
# 
# rho.df[, phylum := factor(phylum, levels = unique(order.vec))]
# rho.df[, rho.group := factor(rho.group, levels = c("pos", "neg"),
#                              labels = c(expression(paste(rho, " > 0")),
#                                         expression(paste(rho, " < 0"))))]

cl.col <- viridisLite::viridis(6, option = "magma", direction = -1)[2:6]
(p <- ggplot(rho.df, aes(x = perc, y = variable)) +
    theme_custom() +
    #facet_wrap(.~season+year, scales = "free_y", ncol = 2,  strip.position = "left") +
    facet_grid(season~year, scales = "free_y") +
    geom_col(aes(fill = dom.clust)) +
    geom_vline(xintercept = 50, linetype = "dashed", colour = "white") +
    scale_y_discrete(limits = rev) +
    scale_x_continuous(breaks = c(0,50, 100)) +
    scale_fill_manual(values = cl.col, 
                      name = "Cluster") +
    labs(x = "Relative contribution to \n total number of significant correlations", y = "Correlation category") +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 10), #angle = 90, vjust = 0.5
          axis.title = element_text(size = 12),
          legend.justification = c(0,0.5),
          #plot.title = element_text(size = 20, hjust = 0.5),
          strip.text = element_text(colour = "white", size = 14)))

ggsave("./Figures/Analysis/all_rho.group_proportion_clust_seasonal.pdf", p, height = 13, width = 25, units = "cm")
# For publication this figure's y-axis was edited in Inkscape


# Proportion of DOM cluster correlations by phyla -----------------------------------------------------
by.phy <- copy(cordf)

by.phy <- by.phy[dna.nss == "sig." & dom.nss == "sig.",]
by.phy <- by.phy[dom.au != "unimodal" & dna.au != "unimodal",]
by.phy[, id := paste(dom.au, dna.au, rho.group)]
by.phy <- by.phy[id %in% c('increase increase pos', 'decrease decrease pos',
                           'increase decrease neg', 'decrease increase neg'),]

# by.phy <- by.phy[!(phylum %in% c("CSP1-3","FCPU426","Gennatimonadota","Cloacimonadota", "Latescibacterota",
#                                  "Fibrobacterota", "Methylomirabilota", "Nitrospirota","Deinococcota", "Dependentiae",
#                                  "Armatimonadota","Desulfobacterota_B", "Firmicutes_C", "Eremiobacterota")),]
nrow(by.phy) # 258,982

by.phy[, n := .N, by = .(phylum, rho.group, dna.au, season, year)]

by.phy %>% group_by(phylum, rho.group, dna.au, season, year) %>%
  summarise(n = unique(n)) %>% ungroup() %>% summarise(n = sum(n))
# 258,982

rho.df <- by.phy[, .(n.all = unique(n),
                     n = .N), by = .(phylum, dom.clust, rho.group, dna.au, season, year)]
rho.df[, perc := n * 100 / n.all]
# sanity check
rho.df[, .(sum = sum(perc)), by = .(phylum, rho.group, dna.au, season, year)]
rho.df[, perc := perc / 100]
# order.vec <- rho.df[rho.group == "pos" & dna.au == "increase" & dom.clust == "1",]
# order.vec <- order.vec[order(-perc),]$phylum

rho.df %>% group_by(phylum, rho.group, dna.au, season, year) %>%
  summarise(n = unique(n.all)) %>% ungroup() %>% summarise(n = sum(n))
# 258,982

rho.df[, variable := paste(dna.au, rho.group, sep = ".")]
rho.df[, variable := factor(variable, levels = c("increase.pos","increase.neg", "decrease.neg", "decrease.pos"),
                            labels = c("Increase - Positive", "Increase - Negative",
                                       "Decrease - Negative", "Decrease - Positive"))]
# Usual bootstrap of the ratio of means using the city data
setDF(rho.df)
mean.boot <- function(data, i){
  d <- data[i,]
  return(mean(d))
}
set.seed(3)
boot.df <- ddply(rho.df, .(dom.clust), function(x){
  bo <- boot::boot(x[, "perc", drop = FALSE], mean.boot, R = 1000, stype = "i")
  bo.conf <- boot::boot.ci(bo, conf=0.95, type="bca")
  
  data.frame(mean = bo$t0,
             ci.lower = bo.conf$bca[4],
             ci.upper = bo.conf$bca[5])
})

setDT(rho.df)
rho.df <- rho.df[boot.df, , on = .(dom.clust)]

rho.df[perc > mean + ci.upper, cat := "above"]
rho.df[perc < mean - ci.lower, cat := "below"]
rho.df[is.na(cat), cat := "below"]
rho.df[cat == "above", cat := season]
rho.df[, cat := factor(cat, levels = c("Spring","Summer", "below"))]

rho.df %>% group_by(phylum, rho.group, dna.au, season, year) %>%
  summarise(n = unique(n.all)) %>% ungroup() %>% summarise(n = sum(n))
# 258,982

rho.df %>% filter(cat != "below") %>% group_by(phylum, rho.group, dna.au, season, year) %>%
  summarise(n = unique(n.all)) %>% ungroup() %>% summarise(n = sum(n))
# 79,115

# only keep phyla that ever have any big shift
rho.df <- rho.df[phylum %in% unique(rho.df[cat %in% c("Spring","Summer"),]$phylum),]

rho.df %>% filter(cat != "below") %>% group_by(phylum, rho.group, dna.au, season, year) %>%
  summarise(n = unique(n.all)) %>% ungroup() %>% summarise(n = sum(n))
# 54,865

# simplify phyla categories
rho.df[, simp.phy := str_split_i(phylum, pattern = "_", i = 1)]

rho.df[, simp.phy := factor(simp.phy, levels = c("Chloroflexota","Dormibacterota","Armatimonadota", "Eremiobacterota",
                                                 "Firmicutes",
                                                 "Actinobacteriota",
                                                 "Cyanobacteria",
                                                 # until here "Terrabacteria"
                                                 "Patescibacteria",
                                                 "Deinococcota",
                                                 # from here "Hydrobacteria"
                                                 "Elusimicrobiota", "Omnitrophota",
                                                 "Dependentiae",
                                                 "Planctomycetota",
                                                 "Verrucomicrobiota",
                                                 "Latescibacterota",
                                                 "Fibrobacterota",
                                                 "Bacteroidota",
                                                 "Proteobacteria", "Acidobacteriota",
                                                 "Bdellovibrionota", "Desulfobacterota",
                                                 "Myxococcota","Methylomirabilota",
                                                 "Nitrospirota",
                                                 "Unclassified"))]
rho.df <- rho.df[simp.phy != "Unclassified"]

phy.cons <- rho.df[(variable == "Increase - Positive" & dom.clust %in% c(5)) & perc > mean + ci.upper,.(n = .N), by = .(simp.phy, season, dom.clust)] %>% arrange(n)
phy.cons$n <- NULL
phy.cons <- phy.cons %>% distinct()

(p <- ggplot() +
    theme_custom() +
    geom_rect(data = boot.df,
              aes(xmin = mean - ci.lower, xmax = mean + ci.upper, ymin = 0 , ymax = length(unique(rho.df$simp.phy)) + 1),
              fill = "#DDDDDD") +
    #geom_hline(yintercept = c(13.5,15.5), linetype = "dashed", colour = "grey50", size = 0.25) +
    geom_vline(data = boot.df, aes(xintercept = mean), linetype = "dashed", colour = "grey20") +
    facet_grid(variable~dom.clust, scales = "free_x") +
    geom_point(data = rho.df[cat == "below",], aes(y = simp.phy, x = perc), colour = "grey70") +
    geom_point(data = rho.df[cat != "below",], aes(y = simp.phy, x = perc, colour = cat), size = 2.5) +
    #geom_point(aes(y = simp.phy, x = perc), size = 2.5) +
    scale_y_discrete(limits=rev) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c(0, 25, 50, 75, 100)) +
    labs(x = "Relative contribution to total number \n of correlations by campaign",
         y = "Phylum") +
    scale_colour_manual(values = c("#44AA99","#DDCC77"), name = "Season") +
    theme(axis.text.x = element_text(size = 8), #angle = 90, vjust = 0.5
          axis.text.y = element_text(size = 7), #angle = 90, vjust = 0.5
          axis.title = element_text(size = 12),
          legend.justification = c(0,0.5),
          #plot.title = element_text(size = 20, hjust = 0.5),
          strip.text = element_text(colour = "white", size = 14)))

(p <- annotate_figure(p, top = text_grob("DOM Cluster", size = 14), 
                      right = text_grob("Correlation category", size = 14, rot = -90)))

ggsave("./Figures/Analysis/phyla_dom.clust_percent_season_new.pdf", p, height = 25, width = 25, units = "cm")


# calculate number of above average links by phyla and cat
n.phy <- rho.df[cat != "below", .(n = .N) , by = .(phylum, variable)]

n.phy <- n.phy %>% arrange(variable, -n)

sum.n <- n.phy %>% group_by(phylum) %>% summarise(sum = sum(n))
nrow(sum.n) # 26 phyla

length(unique(cordf$phylum)) # 36 phyla

# calculate number of above average links by dom cluster
n.dom <- rho.df[cat != "below", .(n = .N) , by = .(dom.clust, variable)]

# calculate how many by season & category
rho.df[perc > mean, side := "above"]
rho.df[perc < mean, side := "below"]
temp <- rho.df[cat != "below", .(n = .N), by = .(cat, dom.clust, variable, side)]

tt <- 
  temp[, .(sum = sum(n)), by = .(cat, variable, side)] %>% arrange(variable, cat, side)
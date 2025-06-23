# Bring FT and MO datasets together and clean

# 1. R set-up ------------------------------------------------------------------------------
### Packages -------------------------------------------------------------------------------
pckgs <- list("phyloseq", "plyr", "tidyverse", "data.table", # wrangling
              "ape") # statistics

### Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)

# remotes::install_github("hrbrmstr/waffle")
# install.packages("waffle", repos = "https://cinc.rud.is")
# Many packages have to be installed through Bioconductor, please refer to the package websites

### Load
invisible(lapply(pckgs, library, character.only = T))
rm(pckgs, miss.pckgs)

### Functions -----------------------------------------------------------------------------
source("./Functions/custom_fun.R")

# 2. Read in data -------------------------------------------------------------------------
# output from model classification
mo <- read.csv("./Data/Summary/MO_patterns.classified_2024-02-12.csv", sep = ",", stringsAsFactors = F) %>% setDT()
ft <- read.csv("./Data/Summary/FT_patterns.classified_2024-02-12.csv", sep = ",", stringsAsFactors = F) %>% setDT()

pks.df <- readRDS(select_newest("./Objects", "MO_pks.df_binned_")) %>% setDT()
mo[pks.df, no.peak := i.no.peak, on = .(ID)]
pks.df <- readRDS(select_newest("./Objects", "FT_pks.df_binned_")) %>% setDT()
ft[pks.df, no.peak := i.no.peak, on = .(ID)]

# # read in database tree
# pb <- readRDS("./Data/phyloseq_otu.tax.samp.phylo_2015-18_SINA.rds")
# otus <- t(otu_table(pb))
# # subset OTUs
# sub.otus <- subset(otus, rownames(otus) %in% unique(sapply(str_split(mo$ID, "[.]"), "[[",4)))
# # subset samples
# cleaned.meta <- read.csv('./Data/Microbial/cleaned.mic_meta.df.csv', stringsAsFactors = F)
# samples <- unique(cleaned.meta$seq_name)
# meta <- sample_data(pb)
# sub.samples <- subset(meta, rownames(meta) %in% samples)
# # extract final phyloseq object
# 
# # specify root of tree
# # Function from: john-quensen.com/r/unifrac-and-tree-roots/
# root.outgroup <- function(unrooted.tree){
#   treeDT <- cbind(data.table(unrooted.tree$edge),
#                   data.table(length = unrooted.tree$edge.length))[1:ape::Ntip(unrooted.tree)] %>%
#     cbind(data.table(id = unrooted.tree$tip.label))
#   # Take the longest terminal branch as outgroup
#   new.outgroup <- treeDT[which.max(length)]$id
#   return(new.outgroup)
# }
# 
# outgroup <- root.outgroup(phy_tree(pb))
# rooted.tree <- ape::root(phy_tree(pb), outgroup = outgroup, resolve.root = T)
# 
# physeq <- merge_phyloseq(sub.otus, tax_table(pb), sub.samples, rooted.tree)
# saveRDS(physeq, "./Data/phyloseq_otu.tax.samp.phylo_SINA_chapter2.rds")
# ape::write.tree(phy_tree(physeq), file = "./Output/chapter2_rooted_SINA.tree")

# if tree figrue is used in publication cite
#Letunic I and Bork P (2021) Nucleic Acids Res doi: 10.1093/nar/gkab301 Interactive Tree Of Life (iTOL) v5: an online tool for phylogenetic tree display and annotation

# all models as list
# ft.ls <- readRDS(select_newest("./Objects","FT_model.ls"))
# mo.ls <- readRDS(select_newest("./Objects","MO_model.ls"))

# 3. Clean data ------------------------------------------------------------------------------------------------------
## Add factor identifiers ------------------------------------------------------------------------------------------
ft[, dataset := "DOM"]

## split ID info
#ft <- ft %>% separate(col = "ID", into = c("year", "season", "MF"), sep = "[_]", remove = F) %>% setDT()
#mo <- mo %>% separate(col = "ID", into = c("year", "season", "dna.type", "OTU"), sep = "[.]", remove = F) %>% setDT()
mo[dna.type == "DNA", dataset := "DNA"]
mo[dna.type == "RNA", dataset := "RNA"]

## Add travel time ranges -----------------------------------------------------------------------------------------
# extract mouth travel time
tt.mouth <- read.csv("./Data/Traveltime/travel.time_byyear_atmouth.csv", sep = ",", dec = ".", stringsAsFactors = F) %>% setDT()
tt.mouth[, season := factor(Month, levels = c("6","8","10"), labels = c("Spring","Summer","Autumn"))]
tt.mouth[, c("year","Year") := list(Year, NULL)]
mo[tt.mouth, tt.at.mouth := i.mouth.trav.time_d, on = .(year, season)]
ft[tt.mouth, tt.at.mouth := i.mouth.trav.time_d, on = .(year, season)]

## Add travel time range categories
ft[x.start == -50 & x.end <= tt.at.mouth, tt.range := "soil-mouth"]
ft[x.start >= 0 & x.end <= tt.at.mouth, tt.range := "stream-mouth"]
ft[x.start == -50 & x.end > tt.at.mouth, tt.range := "soil-lake"]
ft[x.start >= 0 & x.end > tt.at.mouth, tt.range := "stream-lake"]

mo[x.start == -50 & x.end <= tt.at.mouth, tt.range := "soil-mouth"]
mo[x.start >= 0 & x.end <= tt.at.mouth, tt.range := "stream-mouth"]
mo[x.start == -50 & x.end > tt.at.mouth, tt.range := "soil-lake"]
mo[x.start >= 0 & x.end > tt.at.mouth, tt.range := "stream-lake"]

# identify one-peaking non-linear models
mo[no.peak > 1 & AU == "non-linear increase", AU := "multimodal increase"]
ft[no.peak > 1 & AU == "non-linear increase",  AU := "multimodal increase"]

mo[no.peak > 1 & AU == "non-linear decrease", AU := "multimodal increase"]
ft[no.peak > 1 & AU == "non-linear decrease", AU := "multimodal increase"]

# make factors
# ft[, year := factor(year, levels = c("2015","2016"))]
# ft[, season := factor(season, levels = c("Spring","Summer"))]
# mo[, year := factor(year, levels = c("2015","2016"))]
# mo[, season := factor(season, levels = c("Spring","Summer"))]
#mo <- mo[!(dna.type == "RNA"),]v

## simplify trends category ------------------------------------------------------------------------
# ft[ , sAU := factor(AU, levels = c("n.s.", "increase", "non-linear increase",
#                                    "increasing quadratic", "increasing bimodal",
#                                    "unimodal",
#                                    "decreasing bimodal", "decreasing quadratic",
#                                    "non-linear decrease", "decrease"), 
#                     labels = c("n.s.", "increase", "non-linear increase",
#                                "non-linear increase", "non-linear increase",
#                                "unimodal",
#                                "non-linear decrease", "non-linear decrease",
#                                "non-linear decrease", "decrease"))]
ft[ , sAU := factor(AU, levels = c("n.s.", "increase", "non-linear increase",
                                   "increasing quadratic", "increasing bimodal", "multimodal increase",
                                   "unimodal",
                                   "multimodal decrease",
                                   "decreasing bimodal", "decreasing quadratic",
                                   "non-linear decrease", "decrease"), 
                    labels =c("n.s.", "increase", "non-linear increase",
                              "non-linear increase", "multimodal increase", "multimodal increase",
                              "unimodal",
                              "multimodal decrease",
                              "multimodal decrease", "non-linear decrease",
                              "non-linear decrease", "decrease"))]

# mo[ , sAU := factor(AU, levels = c("n.s.", "flat", "increase", "non-linear increase",
#                                    "increasing quadratic", "increasing bimodal",
#                                    "unimodal",
#                                    "decreasing bimodal", "decreasing quadratic",
#                                    "non-linear decrease", "decrease"), 
#                     labels = c("n.s.", "n.s.", "increase", "non-linear increase",
#                                "non-linear increase", "non-linear increase",
#                                "unimodal",
#                                "non-linear decrease", "non-linear decrease",
#                                "non-linear decrease", "decrease"))]
mo[ , sAU := factor(AU, levels = c("n.s.", "flat","increase", "non-linear increase",
                                   "increasing quadratic", "increasing bimodal", "multimodal increase",
                                   "unimodal",
                                   "multimodal decrease",
                                   "decreasing bimodal", "decreasing quadratic",
                                   "non-linear decrease", "decrease"), 
                    labels =c("n.s.", "flat","increase", "non-linear increase",
                              "non-linear increase", "multimodal increase", "multimodal increase",
                              "unimodal",
                              "multimodal decrease",
                              "multimodal decrease", "non-linear decrease",
                              "non-linear decrease", "decrease"))]

# ft[ , ns.s := factor(ns.s, levels = c("n.s.", "sign"), 
#                      labels = c("n.s.", "significant"))]
# 
# mo[ , ns.s := factor(ns.s, levels = c("n.s.", "sign"), 
#                      labels = c("n.s.", "significant"))]

## Remove a few taxa ------------------------------------------------------------------------------
phy <- readRDS("./Data/phyloseq_otu.tax.samp.phylo_SINA_chapter2.rds")
tax <- tax_mat(phy)
tax <- as.data.frame(tax) %>% setDT(keep.rownames = "OTU")
mo <- merge(mo, tax %>% dplyr::select(OTU:domain), by = "OTU", all.x = T)

# number of Archaea and Unclassified
length(unique(mo[domain == "Archaea",]$OTU)) # 55
length(unique(mo[is.na(domain),]$OTU)) #2019
length(unique(mo$OTU)) # out of 6103

# How many of significant models were either one of these categories
nrow(mo[is.na(domain) | domain == "Archaea",][c.ns.s == "sig.",]) # 486 significant models were either archaea or Unclassified
nrow(mo[c.ns.s == "sig.",]) # 1934
(468*100) / 1711 # that's 27%

# Remove
mo <- mo[!(is.na(domain) | domain == "Archaea"),]

## How many of the models actually didn't work? -------------------------------------------------------
nrow(mo[sAU == "flat",]) #1782
nrow(mo[sAU == "n.s.",]) #608

mo[(n = .N), by = .(dna.type)]
nrow(mo[(sAU == "flat" | sAU == "n.s.") & dna.type == "DNA",]) #1924
nrow(mo[(sAU == "flat" | sAU == "n.s.") & dna.type == "RNA",]) #466

nrow(ft[sAU == "n.s.",]) #211

nrow(mo); nrow(ft)

# Just remove
mo <- mo[sAU != "n.s.",][sAU != "flat",]
ft <- ft[sAU != "n.s.",]

## Add a minimal spatial pattern column -------------------------------------------------------------------------------
mo[, minAU := factor(sAU,
                     levels = c('increase', "non-linear increase", "multimodal increase",
                                'unimodal',
                                'multimodal decrease', 'non-linear decrease', 'decrease'),
                     labels = c('increase','increase','increase','unimodal',
                                'decrease', 'decrease', 'decrease'))]
ft[, minAU := factor(sAU,
                     levels = c('increase', "non-linear increase", "multimodal increase",
                                'unimodal',
                                'multimodal decrease', 'non-linear decrease', 'decrease'),
                     labels = c('increase','increase','increase','unimodal',
                                'decrease', 'decrease', 'decrease'))]

## Identify 'true' reactives -------------------------------------------------------------------------
# histogram, calculate quantiles
ft.quar.pos <- quantile(ft[sAU %in% c('increase', 'non-linear increase') & c.ns.s == 'sig.' & lm.slope > 0, ]$lm.slope,
                        probs = c(0,0.33333333, 0.6666666, 1))
ft.quar.neg <- quantile(ft[sAU %in% c('decrease', 'non-linear decrease') & c.ns.s == 'sig.' & lm.slope < 0, ]$lm.slope,
                        probs = c(0,0.33333333, 0.6666666, 1))

(ft.hist <- ggplot(ft[sAU %in% c('increase', 'non-linear increase', 
                                 'non-linear decrease', 'decrease') & c.ns.s == "sig.",],
                   aes(x = lm.slope)) +
    theme_bw() +
    geom_histogram(bins = 100, fill = "white", colour = "gray20") +
    geom_vline(xintercept = ft.quar.pos, colour = "tomato", linetype = 'dashed') +
    geom_vline(xintercept = ft.quar.neg, colour = "royalblue", linetype = 'dashed') +
    labs(x = "z-scaled slope of sign. linear models", y = "Frequency", title = "DOM"))

mo.quar.pos <- quantile(mo[sAU %in% c('increase', 'non-linear increase') & c.ns.s == 'sig.' & 
                             lm.slope > 0 & dna.type == "DNA", ]$lm.slope,
                        probs = c(0,0.33333333, 0.6666666, 1))
mo.quar.neg <- quantile(mo[sAU %in% c('decrease', 'non-linear decrease') & c.ns.s == 'sig.' & lm.slope < 0 &
                             dna.type == "DNA", ]$lm.slope,
                        probs = c(0,0.33333333, 0.6666666, 1))


(mo.hist <- ggplot(mo[sAU %in% c('increase', 'non-linear increase', 
                                 'non-linear decrease', 'decrease') & c.ns.s == "sig.",],
                   aes(x = lm.slope)) + 
    theme_bw() +
    geom_histogram(bins = 100, fill = "white", colour = "gray20") +
    geom_vline(xintercept = mo.quar.pos, colour = "tomato", linetype = 'dashed') +
    geom_vline(xintercept = mo.quar.neg, colour = "royalblue", linetype = 'dashed') +
    labs(x = "z-scaled slope of sign. linear models", y = "Frequency", title = "MO"))


ggsave("./Figures/Analysis/slope_histogram_lm_only_detailed.png",ggarrange(mo.hist, ft.hist, nrow = 2),
       height = 12, width = 15, units = "cm")
#(pat <- pat + geom_smooth(aes(x = x, y = y.pred), method="custom.smooth", se = F, colour = "black"))

#ggsave("./Figures/Analysis/spatial.patterns_modom_poster.png", pat, dpi = 300, height = 15, width = 28, unit = 'cm')
#ggsave("./Figures/Analysis/spatial.patterns_modom_nolabs.png", pat, dpi = 500, height = 10, width = 23, unit = 'cm')
#ggsave("./Figures/Analysis/spatial.patterns_modom_nolabs.pdf", pat, dpi = 500, height = 10, width = 23, unit = 'cm')

# add bins to main data frame
mo[lm.slope <= mo.quar.neg[2], bin.cat := -3]
mo[lm.slope > mo.quar.neg[2] & lm.slope <= mo.quar.neg[3], bin.cat := -2]
mo[lm.slope > mo.quar.neg[3] & lm.slope <= 0, bin.cat := -1]
mo[lm.slope >= 0 & lm.slope <= mo.quar.pos[2], bin.cat := 1]
mo[lm.slope > mo.quar.pos[2] & lm.slope <= mo.quar.pos[3], bin.cat := 2]
mo[lm.slope > mo.quar.pos[3], bin.cat := 3]
# sanity check
mo[is.na(bin.cat) & !is.na(lm.slope) & dna.type == "DNA",]

ft[lm.slope <= ft.quar.neg[2], bin.cat := -3]
ft[lm.slope > ft.quar.neg[2] & lm.slope <= ft.quar.neg[3], bin.cat := -2]
ft[lm.slope > ft.quar.neg[3] & lm.slope <= 0, bin.cat := -1]
ft[lm.slope >= 0 & lm.slope <= ft.quar.pos[2], bin.cat := 1]
ft[lm.slope > ft.quar.pos[2] & lm.slope <= ft.quar.pos[3], bin.cat := 2]
ft[lm.slope > ft.quar.pos[3] , bin.cat := 3]
# sanity check
ft[is.na(bin.cat) & !is.na(lm.slope),]

mo[, bin.cat := factor(bin.cat, levels = c("-3","-2","-1", "1", "2", "3"),
                       labels = c("Q1","Q2", "Q3","Q4","Q5","Q6"))]
ft[, bin.cat := factor(bin.cat, levels = c("-3","-2","-1", "1", "2", "3"),
                       labels = c("Q1","Q2", "Q3","Q4","Q5","Q6"))]

## Reorder -----------------------------------------------------------------------------------
mo <- mo %>% dplyr::select(dataset, ID, OTU, domain, year:dna.type, tt.at.mouth, tt.range, 
                           model:hump, no.peak, ns.s:c.ns.s, AU, sAU, minAU, bin.cat)
ft <- ft %>% dplyr::select(dataset, ID, MF, year:season, tt.at.mouth, tt.range, 
                           model:hump, no.peak, ns.s:c.ns.s, AU, sAU, minAU, bin.cat)

# 4. Save clean version ---------------------------------------------------------------------------------
write.table(mo, paste0("./Data/Summary/MO_patterns.classified_cleaned_",Sys.Date(), ".csv"), sep = ",",
            row.names = F)
write.table(ft, paste0("./Data/Summary/FT_patterns.classified_cleaned_",Sys.Date(), ".csv"), sep = ",",
            row.names = F)

## -------------------------------------------------------------------------
##
## Script name: 8_phylo.analysis.R
##
## Purpose of script: Phylogenetic analyses
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
pckgs <- list("phyloseq", "plyr", "tidyverse", "data.table", # wrangling
              "ggpubr", "plotly", "ggforce","waffle", "gridExtra", "cowplot", "RColorBrewer", 'patchwork', # plotting,
              "kableExtra",# "xlsx", # making tables for manuscript (LaTeX) and export as excel (for ISME)
              "doMC", # parallel computing
              "vegan", "ape", "ade4", "rstatix", "mgcv", "picante") # statistics

## Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)

## Load
invisible(lapply(pckgs, library, character.only = T))
rm(pckgs, miss.pckgs)

## Load custom functions --------------------------------------------------
funs <- list.files("./Functions", full.names = T)
invisible(lapply(funs, source))

# Colour palettes ---------------------------------------------------------
medcont.col <- c( '#EECC66','#EE99AA', '#6699CC', '#997700', '#994455',  '#004488')
light.col <- c('#77AADD', '#EE8866', '#EEDD88', '#FFAABB', '#99DDFF', '#44BB99', '#BBCC33', '#AAAA00', '#DDDDDD')
col.vec <- brewer.pal(7, "RdBu")

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
# output from model classification
mo <- read.csv("./Data/Summary/MO_patterns.classified_cleaned_2024-02-12.csv", sep = ",", stringsAsFactors = F) %>% setDT()
ft <- read.csv("./Data/Summary/FT_patterns.classified_cleaned_2024-02-12.csv", sep = ",", stringsAsFactors = F) %>% setDT()

# Get some statistics ----------------------------------------------------------------------------
  ## How many MF/OTU overall ------------------------------------------------------------------------
(core.mf <- length(unique(ft$MF))) # how many FT molecular formulae in the whole dataset?
# 13251
(core.mo <- length(unique(mo[!(dna.type == "RNA"),]$OTU)))
# 3333

## How many MF/OTU in each season / year? ----------------------------------------------------------
camp.n.ft <- ft[, .(n.camp = .N), by = .(year, season)]
camp.n.mo <- mo[!(dna.type == "RNA"), .(n.camp = .N), by = .(year, season)]

## How many are unique within each season / year ---------------------------------------------------
# Get all unique MFs
mf.ls <- dlply(ft, .(year, season), function(x){
  unique(x$MF)
})

## Count how many MF/OTU are unique within each campaign --------------------------------------------
ft.uni <- data.frame()
names(mf.ls) # order
for(i in 1:length(mf.ls)){
  ft.uni[i,"ID"]<- names(mf.ls)[i]
  ft.uni[i, "uni.n"] <- length(mf.ls[[i]][!(mf.ls[[i]] %in% unique(Reduce(c, mf.ls[-i])))])
}

ft.uni<- ft.uni %>% separate(col = "ID", into = c("year","season")) %>% setDT()

# same for microbes
# Get all unique MFs
mf.ls <- dlply(mo[!(dna.type == "RNA"),], .(year, season), function(x){
  unique(x$OTU)
})

# Count how many are unique within each campaign
mo.uni <- data.frame()
names(mf.ls) # order
for(i in 1:length(mf.ls)){
  mo.uni[i,"ID"]<- names(mf.ls)[i]
  mo.uni[i, "uni.n"] <- length(mf.ls[[i]][!(mf.ls[[i]] %in% unique(Reduce(c, mf.ls[-i])))])
}

mo.uni<- mo.uni %>% separate(col = "ID", into = c("year","season")) %>% setDT()

# check how many of non linear models do not pass the randomization filter (since it's based on linear modelling)
length(unique(ft[model.type != "LM" & r.test == FALSE & p.val < 0.05,]$ID)) # 168
length(unique(mo[model.type != "LM" & r.test == FALSE & p.val < 0.05,]$ID)) # 12

# Microbial ------------------------------------------------------------
#\ Phylogenetic difference in spatial patterns -------------------------

#\\ UniFrac -------------------------------------------------------------
set.seed(3)
# read in chapter 2 phyloseq object
phy <- readRDS("./Data/phyloseq_otu.tax.samp.phylo_SINA_chapter2.rds")
# extract OTU matrix
otu.mat <- otu_mat(phy)

# cast microbial data frame by spatial pattern x year x month
uni.df <- dcast(mo[dna.type == "DNA",],
                OTU ~ year + season + sAU + c.ns.s, value.var = "lm.slope")

all.df <- dcast(mo[dna.type == "DNA",],
                OTU ~ year + season + sAU, value.var = "lm.slope")

# Make sure we have all 0 and 1s (= presence absence)
uni.df[,(colnames(uni.df)[-1]) := lapply(.SD, function(x) ifelse(is.na(x), 0, 1)), .SDcols = colnames(uni.df)[-1]]
uni.df <- setDF(uni.df, rownames = uni.df$OTU)

all.df[,(colnames(all.df)[-1]) := lapply(.SD, function(x) ifelse(is.na(x), 0, 1)), .SDcols = colnames(all.df)[-1]]
all.df <- setDF(all.df, rownames = all.df$OTU)

# specify root of tree
# Function from: john-quensen.com/r/unifrac-and-tree-roots/
root.outgroup <- function(unrooted.tree){
  treeDT <- cbind(data.table(unrooted.tree$edge),
                  data.table(length = unrooted.tree$edge.length))[1:ape::Ntip(unrooted.tree)] %>%
    cbind(data.table(id = unrooted.tree$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup)
}

outgroup <- root.outgroup(phy@phy_tree)
rooted.tree <- ape::root(phy@phy_tree, outgroup = outgroup, resolve.root = T)
# overwrite unrooted tree
phy_tree(phy) <- rooted.tree
# general hypothesis: the groups are the same
# make a phyloseq object for UniFrac distances
uniphy <- phyloseq(otu_table(uni.df[-1], taxa_are_rows = T), rooted.tree)
unifr <- phyloseq::distance(uniphy, method = "unifrac")

uniphy.all <- phyloseq(otu_table(all.df[-1], taxa_are_rows = T), rooted.tree)
unifr.all <- phyloseq::distance(uniphy.all, method = "unifrac")
# Ordinate: PCoA
pcoa.uni <- ordinate(uniphy, method = "PCoA", distance = unifr)
dbrda.uni <- ordinate(uniphy, method = "RDA", distance = unifr)

# And NMDS, stress should be below 0.2
#nmds.uni <- ordinate(uniphy, method = "NMDS", k = 5, distance = unifr) # low enough dimensions with acceptable stress
nmds.uni <- ordinate(uniphy, method = "NMDS", k = 2, distance = unifr)
nmds.uni.all <- ordinate(uniphy.all, method = "NMDS", k = 2, distance = unifr.all)

# quick plot
plot_ordination(uniphy, pcoa.uni)
plot_ordination(uniphy.all, nmds.uni.all)
plot_ordination(uniphy, nmds.uni)
plot_ordination(uniphy, dbrda.uni)

# Cleaner plot, colour by spatial patterns and significance
pcoa.df <- data.frame(pcoa.uni$vectors)
pcoa.df$cat <- row.names(pcoa.df)
pcoa.df <- pcoa.df %>% separate(cat, sep = "[_]", into = c("year","season","sAU","c.ns.s"))
pcoa.df$sAU <- factor(pcoa.df$sAU, levels = c("increase","non-linear increase", "multimodal increase",
                                              "unimodal",
                                              'multimodal decrease', "non-linear decrease", "decrease"))

# Find hulls
find_hull <- function(x, axes){x[chull(x[,paste0("Axis.", axes[1])], x[,paste0("Axis.", axes[2])]),]}
hulls <- ddply(pcoa.df, "c.ns.s", find_hull, axes = c(1,2))

(pcoa.p <- ggplot(pcoa.df, aes(x = Axis.1, y = Axis.2)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_hline(yintercept = 0, colour = "gray70") +
    geom_vline(xintercept = 0, colour = "gray70") +
    geom_polygon(data = hulls, alpha = 0, aes(linetype = c.ns.s), fill = "white", colour = "grey30") +
    geom_point(aes(fill = sAU, size = c.ns.s), shape = 21, colour = "gray70") +
    scale_fill_manual(values = col.vec, name = "Spatial pattern") +
    scale_shape_manual(values = c(4,19), name = "Significance") +
    scale_size_manual(values = c(3,5)) +
    scale_linetype_manual(values = c("dashed","solid")) +
    #scale_shape_manual(values = c(4,19), name = "Significance") +
    annotate(x = c(-0.12,0.25), y = c(0.26,0.23), geom = "text", label = c("sign.","n.s."), colour = "grey30") +
    labs(title = "PCoA: Unweighted UniFrac",
         x = paste("PC1 [", round(pcoa.uni$values$Relative_eig[1],3) * 100,"%]"),
         y = paste("PC2 [", round(pcoa.uni$values$Relative_eig[2],3) * 100,"%]")) +
    guides(linetype = "none", size = "none", 
           fill = guide_legend(override.aes = list(size = 3))))


# Same for NMDS
nmds.df <- data.frame(nmds.uni$points)
nmds.df$cat <- row.names(nmds.df)
nmds.df <- nmds.df %>% separate(cat, sep = "[_]", into = c("year","season","sAU","c.ns.s"))
nmds.df$sAU <- factor(nmds.df$sAU, levels = c("increase","non-linear increase", "multimodal increase",
                                              "unimodal",
                                              'multimodal decrease', "non-linear decrease", "decrease"))
# Find hulls
find_hull <- function(x, axes){x[chull(x[,paste0("MDS", axes[1])], x[,paste0("MDS", axes[2])]),]}
hulls <- ddply(nmds.df, "c.ns.s", find_hull, axes = c(1,2))

(nmds.p <- ggplot(nmds.df, aes(x = MDS1, y = MDS2)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_hline(yintercept = 0, colour = "gray70") +
    geom_vline(xintercept = 0, colour = "gray70") +
    geom_polygon(data = hulls, alpha = 0, aes(linetype = c.ns.s), fill = "white", colour = "grey30") +
    geom_point(aes(fill = sAU, size = c.ns.s), shape = 21, colour = "gray70") +
    scale_fill_manual(values = col.vec, name = "Spatial pattern") +
    scale_size_manual(values = c(3,5)) +
    scale_linetype_manual(values = c("dashed","solid")) +
    #scale_shape_manual(values = c(4,19), name = "Significance") +
    annotate(x = -.25, y = 0.3, geom = "text", label = paste("k = 2, stress = ", round(nmds.uni$stress, 3)),
             colour = "gray20") +
    annotate(x = c(-0.25, 0.5), y = c(0.15,0.15), geom = "text", label = c("n.s.","sign."), colour = "grey30") +
    labs(title = "NMDS: Unweighted UniFrac") +
    guides(linetype = "none", size = "none", 
           fill = guide_legend(override.aes = list(size = 3))))

ggsave("./Figures/Analysis/UniFrac_NMDS_DNA.png", nmds.p, height = 15, width = 20, units = 'cm')

nmds.df.all <- data.frame(nmds.uni.all$points)
nmds.df.all$cat <- row.names(nmds.df.all)
nmds.df.all <- nmds.df.all %>% separate(cat, sep = "[_]", into = c("year","season","sAU"))
nmds.df.all$sAU <- factor(nmds.df.all$sAU, levels = c("increase","non-linear increase", "multimodal increase",
                                                      "unimodal",
                                                      'multimodal decrease', "non-linear decrease", "decrease"))
# Find hulls
# find_hull <- function(x, axes){x[chull(x[,paste0("MDS", axes[1])], x[,paste0("MDS", axes[2])]),]}
# hulls <- ddply(nmds.df.all, "c.ns.s", find_hull, axes = c(1,2))

(nmds.p.all <- ggplot(nmds.df.all, aes(x = MDS1, y = MDS2)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_hline(yintercept = 0, colour = "gray70") +
    geom_vline(xintercept = 0, colour = "gray70") +
    #geom_polygon(data = hulls, alpha = 0, aes(linetype = c.ns.s), fill = "white", colour = "grey30") +
    geom_point(aes(fill = sAU), shape = 21, colour = "gray70", size = 5) +
    scale_fill_manual(values = col.vec, name = "Spatial pattern") +
    scale_linetype_manual(values = c("dashed","solid")) +
    #scale_shape_manual(values = c(4,19), name = "Significance") +
    labs(title = "NMDS: Unweighted UniFrac with all models") +
    guides(linetype = "none", size = "none", 
           fill = guide_legend(override.aes = list(size = 3))))

ggsave("./Figures/Analysis/UniFrac_NMDS_DNA_all.png", nmds.p.all, height = 15, width = 20, units = 'cm')
# save plot
ggsave("./Figures/Analysis/UniFrac_NMDSPCOA_DNA_24-02.png", 
       ggarrange(nmds.p.all, nmds.p, pcoa.p, ncol = 3, common.legend = T, legend = "right"),
       height = 13, width = 45, units = 'cm')

#\\\ PERMANOVA ----------------------------------------------------------------
# add bias.adjust since it's not a balanced sampling design
set.seed(3)
adonis2(formula = unifr ~ sAU + c.ns.s + year + season, data = nmds.df, permutations = 9999, parallel = 10)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = unifr ~ sAU + c.ns.s + year + season, data = nmds.df, permutations = 9999, parallel = 10)
# Df SumOfSqs      R2      F Pr(>F)    
# sAU       6   2.7588 0.16567 1.7510 0.0001 ***
#   c.ns.s    1   1.4522 0.08721 5.5302 0.0001 ***
#   year      1   0.2183 0.01311 0.8315 0.7885    
# season    1   0.1438 0.00863 0.5476 1.0000    
# Residual 46  12.0794 0.72538                  
# Total    55  16.6526 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis2(formula = unifr.all ~ sAU + year + season, data = nmds.df.all, permutations = 9999, parallel = 10)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = unifr.all ~ sAU + year + season, data = nmds.df.all, permutations = 9999, parallel = 10)
# Df SumOfSqs      R2      F Pr(>F)    
# sAU       6   2.6099 0.37426 2.0297 0.0001 ***
#   year      1   0.1739 0.02494 0.8115 0.8639    
# season    1   0.1178 0.01689 0.5497 0.9999    
# Residual 19   4.0719 0.58390                  
# Total    27   6.9735 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# calculate multivariate dispersions
mod <- betadisper(unifr,
                  group = as.character(nmds.df$sAU), bias.adjust = T)

mod2 <- betadisper(unifr,
                   group = as.character(nmds.df$c.ns.s), bias.adjust = T)
mod3 <- betadisper(unifr,
                   group = as.character(nmds.df$year), bias.adjust = T)
mod4 <- betadisper(unifr,
                   group = as.character(nmds.df$season), bias.adjust = T)

# Perform test
perm1 <- anova(mod, permutations = 9999, pariwise = T, parallel = 10)
perm2 <- anova(mod2, permutations = 9999, pariwise = T, parallel = 10)
perm3 <- anova(mod3, permutations = 9999, pariwise = T, parallel = 10)
perm4 <- anova(mod4, permutations = 9999, pariwise = T, parallel = 10)
# dispersion are all equal
# = PERMANOVA assumption fulfilled

# Calculate phylogenetic diversity by spatial pattern
# Faith's diversity is highly correlated to species richness, so it's not particularly helpful in our case.
# faith <- pd(t(otu_mat(uniphy)), tree = phy_tree(uniphy), include.root = T)
# set.seed(3)
# cor.faith <- ses.pd(t(otu_mat(uniphy)), tree = phy_tree(uniphy), include.root = T)
# 
# faith <- setDT(faith, keep.rownames = "ID") %>% separate(col = "ID", into = c("year",'season','sAU','c.ns.s'), sep = "_") %>% setDT()
# faith <- melt(faith, id.vars = c('year','season','sAU','c.ns.s'), variable.name = "diversity", value.name = "value")
# faith[, sAU := factor(sAU, 
#                       levels = c("increase", "non-linear increase",
#                                  "multimodal increase",
#                                  "unimodal",
#                                  "multimodal decrease",
#                                  "non-linear decrease", "decrease"))]
# 
# ggplot(faith, aes(x = sAU, y = log10(value))) + 
#   theme_bw() +
#   geom_boxplot(aes(colour = c.ns.s)) +
#   facet_grid(diversity~., scale ="free_y") +
#   theme(axis.text.x = element_text(angle = 90, vjust = .5))
#unifr <- picante::unifrac(uni.df[,-1], phy@phy_tree) # unweighted, presence absence

rm(nmds.uni.all, pcoa.df, pcoa.p, pcoa.uni, perm1, perm2, perm3, perm4, mod, mod2, mod3, mod4, test.df, uni.df,
   copy.ls, dbrda.uni, hulls, match.df,
   nmds.df, nmds.df.all, nmds.p, nmds.p.all, nmds.uni, p, df)

#\\ NTI vs NRI -------------------------------------------------------------------------------
# Calculate Nearest taxon index (NTI) & Net relatednss index (NRI)
# These are standardised values, standardised to null model.
# mpd = mean pairwise distance -> mean phylogenetic distance (i.e. branch length) among all pairs of species within a community
# is thought to reflect phylogenetic structuring across the whole phylogeny
# mntd = mean nearest taxon distance -> mean distance between each species within a community and its closest relative
# is thought to reflect phylogenetic structure closer to the tips

# NRI = -1 * (MPDobs - mean MPDnull / sd MPDnull)
# NTI = -1 * (MNTDobs - mean MNTDnull / sd MNTDnull)

# in picante, rather than calculating NRI and NTI, standardised effect sizes (SES) are reported.
# These values are equivalent to -1 x NRI or NTI

# cast microbial data frame by spatial pattern x year x month
uni.df <- dcast(mo[dna.type == "DNA",],
                OTU ~ year + season + minAU + c.ns.s, value.var = "lm.slope")

all.df <- dcast(mo[dna.type == "DNA",],
                OTU ~ year + season + minAU, value.var = "lm.slope")

# Make sure we have all 0 and 1s (= presence absence)
uni.df[,(colnames(uni.df)[-1]) := lapply(.SD, function(x) ifelse(is.na(x), 0, 1)), .SDcols = colnames(uni.df)[-1]]
uni.df <- setDF(uni.df, rownames = uni.df$OTU)

all.df[,(colnames(all.df)[-1]) := lapply(.SD, function(x) ifelse(is.na(x), 0, 1)), .SDcols = colnames(all.df)[-1]]
all.df <- setDF(all.df, rownames = all.df$OTU)

# make a phyloseq object for UniFrac distances
uniphy <- phyloseq(otu_table(uni.df[-1], taxa_are_rows = T), rooted.tree)
uniphy.all <- phyloseq(otu_table(all.df[-1], taxa_are_rows = T), rooted.tree)

# calculate phylogenetic distance matrix
phy.dist <- cophenetic(phy_tree(uniphy))
phy.dist.all <- cophenetic(phy_tree(uniphy.all))

# calculate NRI and NTI
set.seed(3)
mpd <- ses.mpd(t(otu_mat(uniphy)), phy.dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 100)
mpd.all <- ses.mpd(t(otu_mat(uniphy.all)), phy.dist.all, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 100)
# mpd.obs.z = standardizied MPD (equivalent to -NRI)
# mpd.obs.p = p-value of observed MPD vs null communities
# positive values of mpd.obs.z indicate phylogenetic evenness 
# (i.e. species within the community are more distantly related than expected by chance)
# negative values of mpd.obs.z indicate phylogenetic clustering
# (i.e. species within the community are more closely related than expected by chance)

# Calculate NRI
nri <- as.matrix(-1 * ((mpd[,2] - mpd[,3]) / mpd[,4]))
rownames(nri) <- row.names(mpd)
colnames(nri) <- "NRI"

nri.all <- as.matrix(-1 * ((mpd.all[,2] - mpd.all[,3]) / mpd.all[,4]))
rownames(nri.all) <- row.names(mpd.all)
colnames(nri.all) <- "NRI"
# negative NRI = overdispersion
# positive NRI = clustering

# Calculate Nearest Taxon Index (NTI)
set.seed(3)
mntd <- ses.mntd(t(otu_mat(uniphy)), phy.dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 100)
mntd.all <- ses.mntd(t(otu_mat(uniphy.all)), phy.dist.all, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 100)

# Calculate NTI
nti <- as.matrix(-1 * ((mntd[,2] - mntd[,3]) / mntd[,4]))
rownames(nti) <- row.names(mntd)
colnames(nti) <- "NTI"
nti.all <- as.matrix(-1 * ((mntd.all[,2] - mntd.all[,3]) / mntd.all[,4]))
rownames(nti.all) <- row.names(mntd.all)
colnames(nti.all) <- "NTI"
# negative NRI = overdispersion
# positive NRI = clustering

# NTI is more sensitive to the clustering or overdispersion happening closer to the tips of a tree
# NRI represents an overall pattern of structuring across the whole phylogenetic tree
nti <- setDT(data.frame(nti), keep.rownames = "ID")
nri <- setDT(data.frame(nri), keep.rownames = "ID")
nti.all <- setDT(data.frame(nti.all), keep.rownames = "ID")
nri.all <- setDT(data.frame(nri.all), keep.rownames = "ID")

phy <- merge(nri, nti, by = "ID")
phy.all <- merge(nri.all, nti.all, by = "ID")

# Let's plot this
phy <- phy %>% separate(ID, into = c("year","season","minAU","c.ns.s"), sep = "_") %>% setDT()
phy.all <- phy.all %>% separate(ID, into = c("year","season","minAU"), sep = "_") %>% setDT()
phy <- melt(phy, id.vars = 1:4, measure.vars = 5:6, value.name = "phyl.rel", variable.name = "index")
phy.all <- melt(phy.all, id.vars = 1:3, measure.vars = 4:5, value.name = "phyl.rel", variable.name = "index")
phy.all$c.ns.s <- "all"

phy <- bind_rows(phy.all, phy)
# phy[, minAU := factor(sAU, levels = c("increase", "non-linear increase",
#                                       "multimodal increase",
#                                       "unimodal",
#                                       "multimodal decrease",
#                                       "non-linear decrease", "decrease"),
#                       labels = c("increase", "increase",
#                                  "increase",
#                                  "unimodal",
#                                  "decrease",
#                                  "decrease", "decrease"))]
# 
# phy[, sAU := factor(sAU,
#                       levels = c("increase", "non-linear increase",
#                                  "multimodal increase",
#                                  "unimodal",
#                                  "multimodal decrease",
#                                  "non-linear decrease", "decrease"))]

phy[, minAU := factor(minAU, levels = c('increase', "unimodal", 'decrease'))]

phy[, c.ns.s := factor(c.ns.s, levels = c('sig.', 'n.s.', 'all'))]

(nri.p <- ggplot(phy[c.ns.s != "n.s." & index == "NRI",]) +
    theme_custom() +
    geom_hline(yintercept = 0, colour = "grey50", linetype = 'dashed') +
    geom_boxplot(aes(x = minAU, y = phyl.rel, fill = c.ns.s), outlier.size = 1, outlier.colour = "grey50") +
    geom_bracket(xmin = c(0.81,1.81), xmax = c(2.81,2.81), label = c("*","*"), y.position = c(-1.1,1.3),
                 tip.length = 0, colour = "darkorange") +
    #facet_grid(index~., scale ="free_y") +
    scale_fill_manual(values = c("orange", "white"), labels = c('Reactive', 'Bulk'), name = "Dataset") +
    scale_x_discrete() +
    scale_y_continuous(limits = c(-1.8,1.8)) +
    annotate(geom = "segment", y = c(0.1,-0.1), yend = c(1.5,-1.5), x = c(0.55, 0.55), xend = c(0.55,0.55),
             arrow = arrow(type = "open", length = unit(0.02, 'npc')), colour = "grey50") +
    annotate(geom = "text", y = c(-1.7,1.7), x = c(0.5,0.5), label = c("Overdispersion","Clustering"),
             colour = "grey50", hjust = 0) +
    labs(y = "NRI", x = "Spatial pattern",
         title = "Microbial phylogeny"))
#ggsave("./Figures/Analysis/NRI_minAU_dunnTest_24-06.png", nti.p, width = 10, height = 10, units = 'cm')

#t <- ggarrange(nri.p, b, ncol = 2, common.legend = T, legend = "right", align = "hv")
#ggsave("./Figures/Analysis/NTI_16Scopynum_minAU.png", t, width = 17, height = 10, units = 'cm')

#\\\ Statistical Testing ----------------------------------------------------------------------------
test.df <- phy[c.ns.s != "n.s." & index == "NTI",]
t.test(test.df[minAU == "increase" & index == "NTI",]$phyl.rel ~ test.df[minAU == "increase" & index == "NTI",]$c.ns.s)
# t.test(test.df[minAU == "non-linear increase" & index == "NTI",]$phyl.rel ~ test.df[minAU == "non-linear increase" & index == "NTI",]$c.ns.s)
# t.test(test.df[minAU == "multimodal increase" & index == "NTI",]$phyl.rel ~ test.df[minAU == "multimodal increase" & index == "NTI",]$c.ns.s)
t.test(test.df[minAU == "unimodal" & index == "NTI",]$phyl.rel ~ test.df[minAU == "unimodal" & index == "NTI",]$c.ns.s)
# t.test(test.df[minAU == "multimodal decrease" & index == "NTI",]$phyl.rel ~ test.df[minAU == "multimodal decrease" & index == "NTI",]$c.ns.s)
# t.test(test.df[minAU == "non-linear decrease" & index == "NTI",]$phyl.rel ~ test.df[minAU == "non-linear decrease" & index == "NTI",]$c.ns.s)
t.test(test.df[minAU == "decrease" & index == "NTI",]$phyl.rel ~ test.df[minAU == "decrease" & index == "NTI",]$c.ns.s)
kruskal.test(phyl.rel ~ c.ns.s, data = test.df) # no

test.df <- phy[c.ns.s == "sig." & index == "NRI",]
#kruskal.test(phyl.rel ~ sAU, data = test.df) # no
kruskal.test(phyl.rel ~ minAU, data = test.df) # yes
dunn.res <- dunn_test(data = test.df, phyl.rel ~ minAU, p.adjust.method = "BH") # increase vs. decrease, unimodal vs. decrease

test.df <- phy[c.ns.s == "sig." & index == "NTI",]
kruskal.test(phyl.rel ~ sAU, data = test.df) # yes
kruskal.test(phyl.rel ~ minAU, data = test.df) # no
dunn.res <- dunn_test(data = test.df, phyl.rel ~ minAU, p.adjust.method = "BH") # none
#dunn.res <- dunn_test(data = test.df, phyl.rel ~ minAU, p.adjust.method = "BH")

test.df <- phy[c.ns.s == "all" & index == "NTI",]
#kruskal.test(phyl.rel ~ sAU, data = test.df) # no
kruskal.test(phyl.rel ~ minAU, data = test.df) # no

test.df <- phy[c.ns.s == "all" & index == "NRI",]
#kruskal.test(phyl.rel ~ sAU, data = test.df) # no
kruskal.test(phyl.rel ~ minAU, data = test.df) # no

# DOM ---------------------------------------------------------------------------
# Read in FT-ICR MS cross tables
cross <- read.csv("./Data/FT/crosstable2015_cor.csv", sep =",", stringsAsFactors = F) %>% dplyr::select(molecular.formula:Class,AI:fl_sugar)
cross16 <- read.csv("./Data/FT/crosstable2016_cor.csv", sep =",", stringsAsFactors = F) %>% dplyr::select(molecular.formula:Class,AI:fl_sugar)

cross <- bind_rows(cross, cross16) %>% setDT()
cross <- cross %>% distinct()

#\ Hierarchical clustering of FT -------------------------------------------------------------------
# following uc-r.github.io/hc_clustering tutorial
library(cluster); library(factoextra)
mf.df <- setDF(cross %>% dplyr::select(molecular.formula:mz, C:S, AImod, NOSC), rownames = cross$molecular.formula)
mf.df$molecular.formula <- NULL
mf.df <- scale(mf.df)

d <- dist(mf.df, method = 'euclidean')
hc1 <- hclust(d, method = 'complete')
dend <- as.dendrogram(hc1)
ft.tree <-  as.phylo(hc1)

#write.tree(ft.tree,"./Output/FT_hiercl.tree")

# Plot tree
# no.label <- function(x) {
#   if (stats::is.leaf(x)){
#     attr(x, "label") <- NULL
#   }
#   return(x)
# }
#test <- dendrapply(dend, no.label)

#clusters <- cutree(hc1, 5)
#clusters <- cutree(hc1, 4)
#plot(test, cex = 0.5, col = "grey50")
#rect.hclust(tree = hc1, k = 5, which = 1:5, cluster =  clusters, border = c("#999933", "#882255","gray","#88CCEE","#AA4499"))


# make it circular and pretty
#library(circlize); library(dendextend)
library(viridisLite)
cl.col <- viridis(6, option = "magma", direction = -1)[2:6]
cl.col.unord <-c(cl.col[1],cl.col[5],cl.col[3],cl.col[4],cl.col[2])
#Cluster order: 1, 5, 3, 4, 2
# small, high, medium, high, small mass
# fav, fav, fav, unfav, unfav
# arom, arom, low arom, low arom, low arom


# Check if the clusters are significantly different
ct <- cutree(hc1, 5) # tried with 5, however, two clusters not statistically significantly different
cl.df <- data.frame(cluster = ct, molecular.formula = names(ct)) %>% setDT()
# save order as well
order.df <- data.frame(order = hc1$order, molecular.formula = hc1$labels) %>% setDT()
cl.df <- cl.df[order.df, , on = .(molecular.formula)]
cl.df[, cluster := factor(cluster, levels = c(1,5,3,4,2), labels = c(1:5))]

# Plot boxplots per clusters
mf.df <- setDF(cross %>% dplyr::select(molecular.formula, mz, C:S, AImod, NOSC))
mf.df <- merge(mf.df, cl.df, by = "molecular.formula") %>% setDT()
mf.df[, CN := C / N]
mf.df[, OC := O / C]
mf.df[, HC := H / C]
mf.df[is.infinite(CN), CN := C]
# melt
box.df <- melt(mf.df[,order :=NULL], id.vars = c("molecular.formula", "cluster"))
box.df <- box.df[, variable := factor(variable, levels = c('mz','C','H','O','N','S',
                                                           'HC','OC','CN', 'AImod', 'NOSC'),
                                      labels = c('Mass*(mz)','C','H','O','N','S', 
                                                 'H/C','O/C','C/N', "'AI'['mod']", 'NOSC'))]

# Do tests
# Variances are not equal. Go for Kruskal-Wallis
data <- box.df
variable.list <- split(data, f = data$variable)
krusk <- lapply(variable.list, function(x) kruskal.test(value ~ cluster, data = x))
dunn <- lapply(variable.list, function(x) rstatix::dunn_test(data = x, value ~ cluster, p.adjust.method = "BH"))

# Extract letters by significance
dunn.mult <- function(x){
  dif <- x$p.adj
  names(dif) <- paste(x$group1, x$group2, sep = "-")
  out <- as.vector(multcompView::multcompLetters(dif))
  out <- data.frame(groups = names(out$Letters), sig.let = out$Letters)
  return(out)
}

# Extract significance of kruskal wallis test
krusk.sig <- function(x){
  return(x$p.value)
}

# get p value of tests
krusk.res <- lapply(krusk, krusk.sig)
krusk.res <- data.frame(variable = names(krusk.res), krusk.p.val = unlist(krusk.res))
any(krusk.res[krusk.res$krusk.p.val >= 0.05,]) # FALSE, all variables significantly different

# get matching letters from post-hoc test
dunn.let <- lapply(dunn, dunn.mult)
group.n <- nrow(dunn.let[[1]])
dunn.let <- data.frame(variable = rep(names(dunn.let), each = group.n),
                       bind_rows(dunn.let)) %>% setDT()

# Merge to plot df
sum.box <- box.df[,.(max.val = max(value),
                     cluster = unique(cluster)), by = .(variable)]
sum.box[, add := max.val *1.1]
sum.box <- merge(sum.box, krusk.res, by = "variable")
sum.box[, cluster := as.factor(cluster)]; dunn.let[, groups := as.factor(groups)]
sum.box <- merge(sum.box, dunn.let, by.x = c("variable","cluster"), by.y = c("variable","groups"))
sum.box <- sum.box[, variable := factor(variable, levels =c('Mass*(mz)','C','H','O','N','S', 'P',
                                                            'H/C','O/C','C/N', "'AI'['mod']", 'NOSC'))]

box.df <- box.df %>% filter(variable != "S")
sum.box <- sum.box %>% filter(variable != "S")
# Plot

(bh <- ggplot() +
    theme_custom() +
    geom_violin(data = box.df, aes(x = as.factor(cluster), y = value,
                                   fill = as.factor(cluster)),
                alpha = 0.5, colour = "white") +
    geom_boxplot(data = box.df, aes(x = as.factor(cluster), y = value,
                                    fill = as.factor(cluster)),
                 alpha = 0.8, width = 0.5, outlier.color = "gray70", outlier.alpha = 0.5, outlier.size = 1) +
    geom_text(data = sum.box, aes(x = as.factor(cluster), y = add, label = sig.let), size = 2) +
    facet_wrap(.~variable, scales = "free", labeller = label_parsed, ncol = 5) +
    scale_fill_manual(values = cl.col, name = "") +
    guides(fill = "none") +
    labs(x = "Cluster", y = ""))
#ggsave("./Figures/Suppl/FT_clusters_boxpot_char_dunn.png", bh, height = 10, width = 20, units = 'cm')

#write.table(cl.df, file = "./Output/dom_hier.cluster_res.csv", sep = ",", row.names = F)

rm(order.df, krusk, dunn, dunn.let, bh, box.df, sum.box,
   variable.list)

cross <- cross[cl.df, , on = .(molecular.formula)]
#cross[, cluster := factor(cluster, levels = c(1,5,3,4,2), labels = c(1,2,3,4,5))]
cl.col <- viridis(6, option = "magma", direction = -1)[2:6]
#cl.col <- viridis(5, option = "magma", direction = -1)

(van1 <- ggplot(cross) +
    theme_custom()+
    geom_point(aes(y = HC, x = OC, colour = as.factor(cluster)), alpha = 0.7) +
    scale_colour_manual(values = cl.col, name = "Cluster")+
    guides(colour = guide_legend(position = "left", override.aes = list(size = 3))) +
    theme(axis.title  = element_text(size = 14),
          legend.title = element_text(size = 12), legend.text = element_text(size = 10),
          legend.key.size = unit(0.5, 'cm')) +
    labs(x = "O/C", y = "H/C"))

cross[, categories := factor(categories,
                             levels = c("aliphatic","high O unsaturated", "low O unsaturated", "polyphenols","condensed aromatics"),
                             labels = c("aliphatic","high O\nunsaturated", "low O\nunsaturated", "polyphenols","condensed\naromatics"))]

(van2 <- ggplot(cross) +
    theme_custom()+
    geom_point(aes(y = HC, x = OC, colour = categories)) +
    scale_colour_brewer(palette = "YlOrBr", name = "Compound\ncategory",) +
    guides(colour = guide_legend(position = "right", override.aes = list(size = 3))) + #nrow = 2, byrow = T, 
    theme(axis.title  = element_text(size = 14),
          legend.title = element_text(size = 12), legend.text = element_text(size = 10),
          legend.key.size = unit(0.5, 'cm'), legend.key.spacing.y = unit(5, 'pt')) +
    labs(x = "O/C", y = "H/C"))
#(cl.van <- ggarrange(van1, van2, ncol = 2, labels = "auto",align = "hv"))

(cl.van <- (van1 + van2 + plot_layout(guides = 'collect', axes = 'collect') +
              plot_annotation(tag_levels = "a") & 
              theme(plot.tag = element_text(size = 14, face = "bold"))))

ggsave("./Figures/Analysis/FT_hier.cl_van.krev_24-06.png", cl.van, width = 23, height = 10, units = 'cm')

# Chemical similarity ----------------------------------------------------------------------

# merge results to FT dataset
clusters <- cutree(hc1, 5)
clusters <- data.frame(MF = names(clusters), clusters = clusters)
ft <- merge(ft, clusters, by = "MF")

# cast FT data frame by spatial pattern x year x month
uni.df <- dcast(ft[c.ns.s == "sig.",],
                MF ~ year + season + minAU + c.ns.s, value.var = "lm.slope")

all.df <- dcast(ft,
                MF ~ year + season + minAU, value.var = "lm.slope")



# Make sure we have all 0 and 1s (= presence absence)
uni.df[,(colnames(uni.df)[-1]) := lapply(.SD, function(x) ifelse(is.na(x), 0, 1)), .SDcols = colnames(uni.df)[-1]]
uni.df <- setDF(uni.df, rownames = uni.df$MF)

all.df[,(colnames(all.df)[-1]) := lapply(.SD, function(x) ifelse(is.na(x), 0, 1)), .SDcols = colnames(all.df)[-1]]
all.df <- setDF(all.df, rownames = all.df$MF)

uniphy <- phyloseq(otu_table(uni.df[,-1], taxa_are_rows = T), ft.tree)
allphy <- phyloseq(otu_table(all.df[,-1], taxa_are_rows = T), ft.tree)

# calculate NRI and NTI
set.seed(3)
phy.dist <- cophenetic(phy_tree(uniphy))
phy.dist.all <- cophenetic(phy_tree(allphy))

# calculate NRI and NTI
set.seed(3)
mpd <- ses.mpd(t(otu_mat(uniphy)), phy.dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 100)
mpd.all <- ses.mpd(t(otu_mat(allphy)), phy.dist.all, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 100)
# mpd.obs.z = standardizied MPD (equivalent to -NRI)
# mpd.obs.p = p-value of observed MPD vs null communities
# positive values of mpd.obs.z indicate phylogenetic evenness 
# (i.e. species within the community are more distantly related than expected by chance)
# negative values of mpd.obs.z indicate phylogenetic clustering
# (i.e. species within the community are more closely related than expected by chance)

# Calculate NRI
nri <- as.matrix(-1 * ((mpd[,2] - mpd[,3]) / mpd[,4]))
rownames(nri) <- row.names(mpd)
colnames(nri) <- "NRI"

nri.all <- as.matrix(-1 * ((mpd.all[,2] - mpd.all[,3]) / mpd.all[,4]))
rownames(nri.all) <- row.names(mpd.all)
colnames(nri.all) <- "NRI"
# negative NRI = overdispersion
# positive NRI = clustering

# Calculate Nearest Taxon Index (NTI)
set.seed(3)
mntd <- ses.mntd(t(otu_mat(uniphy)), phy.dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 100)
mntd.all <- ses.mntd(t(otu_mat(allphy)), phy.dist.all, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 100)

# Calculate NTI
nti <- as.matrix(-1 * ((mntd[,2] - mntd[,3]) / mntd[,4]))
rownames(nti) <- row.names(mntd)
colnames(nti) <- "NTI"
nti.all <- as.matrix(-1 * ((mntd.all[,2] - mntd.all[,3]) / mntd.all[,4]))
rownames(nti.all) <- row.names(mntd.all)
colnames(nti.all) <- "NTI"
# negative NRI = overdispersion
# positive NRI = clustering

# NTI is more sensitive to the clustering or overdispersion happening closer to the tips of a tree
# NRI represents an overall pattern of structuring across the whole phylogenetic tree
nti <- setDT(data.frame(nti), keep.rownames = "ID")
nri <- setDT(data.frame(nri), keep.rownames = "ID")
nti.all <- setDT(data.frame(nti.all), keep.rownames = "ID")
nri.all <- setDT(data.frame(nri.all), keep.rownames = "ID")

phy <- merge(nri, nti, by = "ID")
phy.all <- merge(nri.all, nti.all, by = "ID")

# Let's plot this
phy <- phy %>% separate(ID, into = c("year","season","sAU","c.ns.s"), sep = "_") %>% setDT()
phy.all <- phy.all %>% separate(ID, into = c("year","season","sAU"), sep = "_") %>% setDT()
phy <- melt(phy, id.vars = 1:4, measure.vars = 5:6, value.name = "phyl.rel", variable.name = "index")
phy.all <- melt(phy.all, id.vars = 1:3, measure.vars = 4:5, value.name = "phyl.rel", variable.name = "index")
phy.all$c.ns.s <- "all"

phy <- bind_rows(phy.all, phy)
setDT(phy)
phy[, minAU := factor(sAU, levels = c("increase", "non-linear increase",
                                      "multimodal increase",
                                      "unimodal",
                                      "multimodal decrease",
                                      "non-linear decrease", "decrease"),
                      labels = c("increase", "increase",
                                 "increase",
                                 "unimodal",
                                 "decrease",
                                 "decrease", "decrease"))]

phy[, sAU := factor(sAU,
                    levels = c("increase", "non-linear increase",
                               "multimodal increase",
                               "unimodal",
                               "multimodal decrease",
                               "non-linear decrease", "decrease"))]

phy[, c.ns.s := factor(c.ns.s, levels = c('sig.', 'n.s.', 'all'))]

phy <- phy[!is.na(minAU),]

(nti.p.ft <- ggplot(phy[c.ns.s != "n.s." & index == "NTI",]) +
    theme_custom() +
    geom_hline(yintercept = 0, colour = "grey50", linetype = 'dashed') +
    geom_boxplot(aes(x = minAU, y = phyl.rel, fill = c.ns.s), outlier.size = 1, outlier.colour = "grey50") +
    geom_bracket(xmin = 0.81, xmax = 1.81, label = "*", y.position = 5, tip.length = 0, colour = "darkorange") +
    #facet_grid(season~year, scale ="free_y") +
    scale_fill_manual(values = c("orange", "white"), labels = c('Reactive', 'Bulk'), name = "Dataset") +
    scale_x_discrete() +
    #scale_y_continuous(limits = c(-3,3)) +
    annotate(geom = "segment", y = c(1,-1), yend = c(48,-3), x = c(0.55, 0.55), xend = c(0.55,0.55),
             arrow = arrow(type = "open", length = unit(0.02, 'npc')), colour = "grey50") +
    annotate(geom = "text", y = c(-5,50), x = c(0.5,0.5), label = c("Overdispersion","Clustering"),
             colour = "grey50", hjust = 0) +
    labs(y = "NTI", x = "Spatial pattern",
         title = "Chemical similarity"))
ggsave("./Figures/Analysis/FT_NTI_minAU_dunnTest.png", nti.p, width = 10, height = 10, units = 'cm')

# (nri.p <- ggplot(phy[c.ns.s != "n.s." & index == "NRI",]) +
#     theme_custom() +
#     geom_hline(yintercept = 0, colour = "grey50", linetype = 'dashed') +
#     geom_boxplot(aes(x = minAU, y = phyl.rel, fill = c.ns.s), outlier.size = 1, outlier.colour = "grey50") +
#     geom_bracket(xmin = 2.2, xmax = 3.2, label = "**", y.position = 68, tip.length = 0) +
#     #facet_grid(index~., scale ="free_y") +
#     scale_fill_manual(values = c("orange", "white"), labels = c('Reactive', 'Bulk'), name = "Dataset") +
#     scale_x_discrete() +
#     #scale_y_continuous(limits = c(-3,3)) +
#     annotate(geom = "segment", y = c(2,-2), yend = c(70.5,-5), x = c(0.55, 0.55), xend = c(0.55,0.55),
#              arrow = arrow(type = "open", length = unit(0.02, 'npc')), colour = "grey50") +
#     annotate(geom = "text", y = c(-9,75), x = c(0.5,0.5), label = c("Overdispersion","Clustering"),
#               colour = "grey50", hjust = 0) +
#     labs(y = "Nearest Relative Index (NRI)", x = "Spatial pattern",
#          title = "Chemical similarity"))
ggsave("./Figures/Analysis/FT_NRI_minAU_dunnTest.png", nri.p, width = 10, height = 10, units = 'cm')

ft.t <- ggarrange(nti.p, b, ncol = 2, common.legend = T, legend = "right", align = "hv")
#ggsave("./Figures/Analysis/NTI_16Scopynum_minAU.png", t, width = 17, height = 10, units = 'cm')

(t <- (nri.p + nti.p.ft + plot_layout(guides = 'collect', axes = 'collect') +
    plot_annotation(tag_levels = "a") & 
    theme(plot.tag = element_text(size = 14, face = "bold"))))

ggsave("./Figures/Analysis/MO-NRI_FT-NTI.png", t, width = 17, height = 9.5, units = 'cm')

# test between bulk vs reactive
test.df <- phy[c.ns.s != "n.s." & index == "NTI",]
t.test(test.df[sAU == "increase" & index == "NTI",]$phyl.rel ~ test.df[sAU == "increase" & index == "NTI",]$c.ns.s)
# t.test(test.df[sAU == "non-linear increase" & index == "NTI",]$phyl.rel ~ test.df[sAU == "non-linear increase" & index == "NTI",]$c.ns.s)
# t.test(test.df[sAU == "multimodal increase" & index == "NTI",]$phyl.rel ~ test.df[sAU == "multimodal increase" & index == "NTI",]$c.ns.s)
t.test(test.df[sAU == "unimodal" & index == "NTI",]$phyl.rel ~ test.df[sAU == "unimodal" & index == "NTI",]$c.ns.s)
# t.test(test.df[sAU == "multimodal decrease" & index == "NTI",]$phyl.rel ~ test.df[sAU == "multimodal decrease" & index == "NTI",]$c.ns.s)
# t.test(test.df[sAU == "non-linear decrease" & index == "NTI",]$phyl.rel ~ test.df[sAU == "non-linear decrease" & index == "NTI",]$c.ns.s)
t.test(test.df[sAU == "decrease" & index == "NTI",]$phyl.rel ~ test.df[sAU == "decrease" & index == "NTI",]$c.ns.s)

# bulk NRI test
test.df <- phy[c.ns.s == "all" & index == "NRI",]
kruskal.test(phyl.rel ~ minAU, data = test.df) # no
#dunn.res <- dunn_test(data = test.df, phyl.rel ~ sAU, p.adjust.method = "BH") # none
#dunn.res <- dunn_test(data = test.df, phyl.rel ~ minAU, p.adjust.method = "BH") # unimodal vs. decrease p = 0.016

# bulk NTI test
test.df <- phy[c.ns.s == "all" & index == "NTI",]
kruskal.test(phyl.rel ~ minAU, data = test.df) # no

# reactive NRI test
test.df <- phy[c.ns.s == "sig." & index == "NRI",]
kruskal.test(phyl.rel ~ minAU, data = test.df) # no
# reactive NTI test
test.df <- phy[c.ns.s == "sig." & index == "NTI",]
kruskal.test(phyl.rel ~ minAU, data = test.df) # yes
dunn.res <- dunn_test(data = test.df, phyl.rel ~ minAU, p.adjust.method = "BH") # increase vs. unimodal p = 0.016

## Functional similarity in DOM -----------------------------------------------------------------------------------

ft[mf.df, c("AImod","NOSC") := list(i.AImod, i.NOSC), on = c("MF==molecular.formula")]
reac <- ft[c.ns.s == "sig.",]
reac[, bulk := "Reactive"]
ft[, bulk := "Bulk"]

all.ft <- bind_rows(reac, ft)

all.ft[, bulk := factor(bulk, levels = c("Reactive","Bulk"))]
all.ft[, minAU := factor(minAU, levels = c("increase","unimodal","decrease"))]

all.ft <- all.ft[!is.na(minAU),]
(nosc <- ggplot(all.ft, aes(x = minAU, y = NOSC)) +
    theme_custom() +
    geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
    #facet_grid(season~year) +
    geom_boxplot(aes(fill = bulk), outlier.size = 1, outlier.colour = "grey50") +
    geom_bracket(xmin = c(1.2), xmax = c(3.2), label = c("***"), y.position = c(2.2), tip.length = 0, colour = "grey50",
                 label.size = 2.5) +
    geom_bracket(xmin = c(1.2, 2.2), xmax = c(2.15,3.2), label = c("***","***"), y.position = c(2,2), tip.length = 0,
                 colour = "grey50", label.size = 2.5) +
    geom_bracket(xmin = c(.8, 1.85), xmax = c(1.8,2.8), label = c("***","***"), y.position = c(1.6,1.6), tip.length = 0,
                 colour= 'darkorange', label.size = 2.5) +
    geom_bracket(xmin = c(.8), xmax = c(2.8), label = c("***"), y.position = c(1.8), tip.length = 0,
                 colour= 'darkorange', label.size = 2.5) +
    scale_fill_manual(values = c("orange","white"), name = "Dataset") +
    labs(x = "Spatial pattern", y = "NOSC"))

test.df <- all.ft[bulk == "Bulk",]
kruskal.test(NOSC ~ minAU, data = test.df) # yes, p < 0.016
dunn.res <- dunn_test(data = test.df, NOSC ~ minAU, p.adjust.method = "BH") # unimodal vs. decrease, unimodal vs. increase

test.df <- all.ft[bulk == "Reactive",]
kruskal.test(NOSC ~ minAU, data = test.df) # yes, p < 0.016
dunn.res <- dunn_test(data = test.df, NOSC ~ minAU, p.adjust.method = "BH") # all

# calculate means and sd
all.ft[, .(mean = mean(NOSC), sd = sd(NOSC)), by = .(minAU)]


(ai <- ggplot(all.ft, aes(x = minAU, y = AImod)) +
    theme_custom() +
    geom_boxplot(aes(fill = bulk), outlier.size = 1, outlier.colour = "grey50") +
    geom_bracket(xmin = c(1.2), xmax = c(3.2), label = c("***"), y.position = c(1.2), tip.length = 0, colour = "grey50",
                 label.size = 2.5) +
    geom_bracket(xmin = c(1.2, 2.2), xmax = c(2.15,3.2), label = c("***","***"), y.position = c(1.1), tip.length = 0,
                 colour = "grey50", label.size = 2.5) +
    geom_bracket(xmin = c(.8, 1.85), xmax = c(1.8,2.8), label = c("***","***"), y.position = c(0.9,0.9), tip.length = 0,
                 colour= 'darkorange') +
    geom_bracket(xmin = c(.8), xmax = c(2.8), label = c("***"), y.position = c(1), tip.length = 0,
                 colour= 'darkorange') +
    scale_fill_manual(values = c("orange","white"), name = "Dataset") +
    labs(x = "Spatial pattern", y = expression(paste("AI"["mod"]))))

test.df <- all.ft[bulk == "Bulk",]
kruskal.test(AImod ~ minAU, data = test.df) # yes, p < 0.016
dunn.res <- dunn_test(data = test.df, AImod ~ minAU, p.adjust.method = "BH") # all

test.df <- all.ft[bulk == "Reactive",]
kruskal.test(AImod ~ minAU, data = test.df) # yes, p < 0.016
dunn.res <- dunn_test(data = test.df, AImod ~ minAU, p.adjust.method = "BH") # all

(func.p <- (nosc + ai) + plot_layout(axes = "collect", guides = "collect") +
    plot_annotation(tag_levels = "a") & 
    theme(plot.tag = element_text(size = 14, face = "bold")))

ggsave("./Figures/Analysis/FT-NOSC_AImod.png", func.p, width = 17, height = 9.5, units = 'cm')


# Number of models per travel time range --------------------------------------------------------------------
setDT(mo); setDT(ft)
mo[, n.camp := .N, by = .(year, season, dataset)]
ft[, n.camp := .N, by = .(year, season)]
range.df <- bind_rows(mo[, .(n.range = .N), by = .(tt.range, year, season, dataset, n.camp)],
                      ft[, .(n.range = .N), by = .(tt.range, year, season, dataset, n.camp)])
# calculate percentage
range.df[, perc := round(n.range * 100 / n.camp, 1)]

range.df[, tt.range := factor(tt.range, levels = c("soil-lake", "soil-mouth", "stream-lake", "stream-mouth"))]
range.df[, dataset := factor(dataset, levels = c("DOM", "DNA", "RNA"))]

(p <- ggplot(range.df, aes(x = season, y = perc)) +
    theme_custom() +
    facet_grid(dataset~year) +
    geom_col(aes(fill = tt.range), width = 0.5) +
    geom_text(aes(label = n.camp), y = 104, size = 2, colour = 'gray50') +
    scale_fill_manual(values = brewer.pal(4, "BrBG")[c(1:2,4,3)], name = "Ecosystem\nrange") +
    scale_y_continuous(breaks = c(0,50,100), limits = c(0,104)) +
    labs(x = "", y = "Proportion of models\nby travel time range") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)))

ggsave('./Figures/Analysis/travel.time_range_bymodel.png', height = 10, width = 10, units = 'cm', dpi = 300)
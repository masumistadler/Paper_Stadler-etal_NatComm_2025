## -------------------------------------------------------------------------
##
## Script name: 0_FTcleaning.R
##
## Purpose of script: Clean DOM data derived from FT-ICR MS.
##
## Author: Masumi Stadler
## Based on recommendations of: Francois Guillemette, Trista Vick-Majors
##
## Date Finalized: 05.06.2025
##
## Copyright (c) Masumi Stadler, 2025
## Email: m.stadler.jp.at@gmail.com
##
## -------------------------------------------------------------------------
##
## Notes:
##
## Cleaning steps ##
## Following steps were deemed necessary by Trista Vick Majors and Francois Guillemette:
##
## 1. Remove molecular formulae deemed "rare"
##    A formula that occurs less than 9 times
##    Which corresponds to the minimum number of sites within a given sample type
##    We have 91 samples in 2015, and any formulae should be present in > ~10% of samples
##
## 2. Merging FT categories
##   - "Black Carbon" will be removed
##   - Merge "High" and "Low" OC categories
##   - Merge "highly unsaturated" and "phenolics"
##   - Merge "with" and "wihout" sugars
## Final category list should be:
## Peptides, Polyphenolics, Sugars, Condensed Aromatics, Highly Unsaturated and Phenolics, and Aliphatics
##
## 3. Remove contaminants (= outliers)
## 2015: LRH08 and LRH15 (LR08 is a bit odd but we keep it)
## 2016: SW34, SW36, RO110, LR70, RO240, PR32, LR63, LR67, LR59, LR55, LR57, M18, HW42, T48, RO251, LR58, LR60, LR62
## Outliers were removed based on ordination plots
##
## Cleaning was done partially in excel, partially in R
## To unify/track the cleaning process, the above mentioned steps were implemented in this script
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
              "data.table", "stringr", "plyr", "tidyverse", # wrangling
              "vegan", # multivariate
              "plotly"
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

# 2015 -----------------------------------------------------------------------------------------------------------
# Read in raw data ------------------------------------------------------------------------------------------------
cross <- read_excel("./Data/FT/Raw/Cross-Table - FT-MS - Quebec 2015.xlsx")
compclass <- read_excel("./Data/FT/Raw/Compound Classes Table - FT-MS - Quebec 2015.xlsx", sheet = 1)
leg <- read_excel("./Data/FT/Raw/Compound Classes Table - FT-MS - Quebec 2015.xlsx", sheet = 2)

# Clean column names -----------------------------------------------------------------------------------------
colnames(cross)[1:3] <- c("rID", "molecular.formula","mz")
colnames(compclass)[1] <- "orig"

first.pos <-which(str_detect(colnames(cross),"ESI") == T)[1]
# what is column "3"? = found in compclass
out<-which(str_detect(colnames(cross),"ESI") == F)
colnames(cross)[out[out > first.pos]] <- "ESI.Neg_Quebec_PR03_sum100.42.csv.RelAbun"

orig <- colnames(cross)[first.pos:ncol(cross)]
new.col <- sapply(strsplit(colnames(cross)[first.pos:ncol(cross)], split = "_"), "[", 3)
# one sample not separated by "."
new.col[new.col == "RO2"] <- "RO2.22t7"

# # read in Tristy's data to check
# tristy <- read.csv("./Data/FT/Tristy/CrossTable2015.csv", sep = ",", stringsAsFactors = F)
# tristy.nc <- read.csv("./Data/FT/Tristy/CrossTable2016_Final_TVM_noRare_noContam.csv", sep = ",", stringsAsFactors = F)
# colnames(tristy)[!(colnames(tristy) %in% new.col)]

clean.col <- data.frame(col = new.col)
clean.col <- clean.col %>% separate(col = "col", into = c("a","b"),sep = "[.]")

# unify sample naming strategy to match with master file
setDT(clean.col)
clean.col[, typ := sapply(str_extract_all(clean.col$a, "[A-Z]+"),"[",1)] # extract capital letters
clean.col[!is.na(b), typ := a] # overwrite RO2, as number is needed
clean.col[typ == "RO", typ := a] # one assign mistake
clean.col[is.na(b), no := parse_number(a)] # get numbers
clean.col[is.na(no), no := as.numeric(b)] # add RO2 sample number
clean.col[!str_detect(typ, pattern = "WELL"), cor.no := sprintf("%02d", no)] # add 0 to single digit numbers
# except WELL, meta data is without extra zero
clean.col[str_detect(typ, pattern = "WELL"), cor.no := no]
clean.col[typ != "RO2", cor.col := paste0(typ, cor.no)] # paste type and number together
clean.col[typ == "RO2", cor.col := paste(typ, cor.no, sep = ".")] # add a dot for RO samples

# clean one sample, RO2_22t7 = duplicate with RO2.22?
clean.col[cor.col == "RO2.NA", cor.col := "RO2.22-r2"]

# compare
col.rename <- data.frame(orig = str_remove(orig,".csv.RelAbun"), old = new.col, new = clean.col$cor.col)
# # all good
# colnames(tristy)[!(colnames(tristy) %in% col.rename$new)]

# Overwrite sample names ---------------------------------------------------------------------------------------
colnames(cross)[first.pos:ncol(cross)] <- col.rename$new

# correct compclass sample names
setDT(compclass); setDT(col.rename)
compclass[col.rename, samples := i.new, on = .(orig)]

setcolorder(compclass, neworder = "samples") # put to first
compclass[, orig := NULL]

# 2016 --------------------------------------------------------------------------------------------------------
# Read in raw data ------------------------------------------------------------------------------------------------
cross16 <- read_excel("./Data/FT/Raw/FT-ICRMS-2016-CrossTable_Final.xlsx", sheet = 1)
compclass16 <- read_excel("./Data/FT/Raw/FT-ICRMS-2016-SummaryTable_Final.xlsx", sheet = 1)
sample.ref <- read_excel("./Data/FT/Raw/FT_ICRMS_Samples2016.xlsx", sheet = "coord")

colnames(cross16)[1:3] <- c("rID", "molecular.formula","mz")
colnames(compclass16)[1] <- "orig"

# separate numbers from column names
colnames(cross16)

first.pos <-which(str_detect(colnames(cross16),"ESI") == T)[1]
orig <- colnames(cross16)[first.pos:ncol(cross16)]

new.col <- sapply(strsplit(orig, split = "_"), "[", 2)

# merge with reference table of sample names corresponding to vial numbers
clean.col <- data.frame(vial = as.numeric(new.col), orig = str_remove(orig,".csv.RelAbun"))
setDT(clean.col); setDT(sample.ref)
clean.col <- sample.ref[clean.col, , on = .(vial)]

# some have a weird "_" at the end
# unify ROx- to ROx.
clean.col[, new := site]
clean.col[, new := str_replace(new, "-", ".")]
clean.col[, new := str_replace(new, "_$", "")]
# sanity check
any(duplicated(clean.col$new) == T) # FALSE

# Overwrite sample names ---------------------------------------------------------------------------------------
colnames(cross16)[first.pos:ncol(cross16)] <- clean.col$new

# correct compclass sample names
setDT(compclass16)
compclass16[clean.col, samples := i.new, on = .(orig)]

setcolorder(compclass16, neworder = "samples") # put to first
compclass16[, orig := NULL]

# 1. Remove rare molecular formulae ------------------------------------------------------------------------------
setDT(cross); setDT(cross16)
sub.cross <- cross[Count >= 9,]
nrow(cross) - nrow(sub.cross) # removed 5303 MF

sub.cross16 <- cross16[Count >= 9,]
nrow(cross16) - nrow(sub.cross16) # removed 7306 MF

# overwrite main cross table
cross <- sub.cross
cross16 <- sub.cross16

# 2. Merge FT categories ----------------------------------------------------------------------------------------
colnames(compclass) == colnames(compclass16) # not the same among years
#   - "Black Carbon" will be removed
#   - Merge "High" and "Low" OC categories
#   - Merge "highly unsaturated" and "phenolics"
#   - Merge "with" and "wihout" sugars
# Final category list should be:
# Peptides, Polyphenolics, Sugars, Condensed Aromatics, Highly Unsaturated and Phenolics, and Aliphatics

# 2015 ----------------------------------------------------------------------------------------------------------
# melt
melt.comp <- melt(compclass, id.vars = "samples", variable.name = "Comp")
levels(factor(melt.comp$Comp))
# remove Black carbon
melt.comp <- melt.comp[!str_detect(Comp, "BC.X"),]
# merge high and low OC categories
melt.comp[, gen.comp := str_remove(Comp, pattern = "_HighOC|_LowOC")]
# merge with and without sugars
melt.comp[, gen.comp := str_remove(gen.comp, pattern = "_wo_X|_w_X")]
# remove "like"
melt.comp[, gen.comp := str_remove(gen.comp, pattern = "_like")]
# correct some inconsistencies
melt.comp[gen.comp == "Condensed_Aromatics_(%RA)", gen.comp := "Condensed_aromatics_(%RA)"]
melt.comp[gen.comp == "Aliphatic_R=(%RA)", gen.comp := "Aliphatic_(%RA)"]

# How should we summarise the now merged categories?
# % and # count = sum
# average = average

averages <- melt.comp[str_detect(melt.comp$gen.comp, pattern = "Ave"),]
sums <- melt.comp[!str_detect(melt.comp$gen.comp, pattern = "Ave"),]

# sanity check
levels(factor(averages$gen.comp))
levels(factor(sums$gen.comp))

# Summarise
sums <- sums[, .(sum = sum(value)), by = .(samples, gen.comp)]
averages <- averages[, .(mean = mean(value)), by = .(samples, gen.comp)]

# cast back
sums <- dcast(sums, samples ~ gen.comp, value.var = "sum")
averages <-  dcast(averages, samples ~ gen.comp, value.var = "mean")

# merge
n.compclass <- sums[averages, , on = .(samples)]

# Apply new categories to cross table
# remove Black carbon
cross <- cross[!str_detect(Compound.category, "BC.X"),]
# merge high and low OC categories
cross[, Compound.category := str_remove(Compound.category, pattern = "_HighOC|_LowOC")]
# merge with and without sugars
cross[, Compound.category := str_remove(Compound.category, pattern = "_wo_X|_w_X")]
# remove "like" in peptides
cross[, Compound.category := str_remove(Compound.category, pattern = "_like")]

# 2016 --------------------------------------------------------------------------------------------------------
# melt
melt.comp <- melt(compclass16, id.vars = "samples", variable.name = "Comp")
levels(factor(melt.comp$Comp))
# no black carbon in 2016
levels(factor(melt.comp$gen.comp))
# merge high and low OC categories
melt.comp[, gen.comp := str_remove(Comp, pattern = "_HighOC|_LowOC")]
# merge with and without sugars
melt.comp[, gen.comp := str_remove(gen.comp, pattern = "_wo_X|_w_X")]
# remove "like"
melt.comp[, gen.comp := str_remove(gen.comp, pattern = "_like")]
# correct some inconsistencies
melt.comp[gen.comp == "Condensed_Aromatics_(%RA)", gen.comp := "Condensed_aromatics_(%RA)"]
melt.comp[, gen.comp := str_remove(gen.comp, pattern = ".X")]

# How should we summarise the now merged categories?
# % and # count = sum
# average = average

averages <- melt.comp[str_detect(melt.comp$gen.comp, pattern = "Ave"),]
sums <- melt.comp[!str_detect(melt.comp$gen.comp, pattern = "Ave"),]

# sanity check
levels(factor(averages$gen.comp))
levels(factor(sums$gen.comp))

# Summarise
sums <- sums[, .(sum = sum(value)), by = .(samples, gen.comp)]
averages <- averages[, .(mean = mean(value)), by = .(samples, gen.comp)]

# cast back
sums <- dcast(sums, samples ~ gen.comp, value.var = "sum")
averages <-  dcast(averages, samples ~ gen.comp, value.var = "mean")

# merge
n.compclass16 <- sums[averages, , on = .(samples)]

# Apply new categories to cross table
# merge high and low OC categories
cross16[, Compound.category := str_remove(Compound.category, pattern = "_HighOC|_LowOC")]
# merge with and without sugars
cross16[, Compound.category := str_remove(Compound.category, pattern = "_wo_X|_w_X")]
# remove "like" in peptides
cross16[, Compound.category := str_remove(Compound.category, pattern = "_like")]
cross16[, Compound.category := str_remove(Compound.category, pattern = ".X")]

levels(factor(cross16$Compound.category)) == levels(factor(cross$Compound.category)) # TRUE

# Final tables at this point are:
head(n.compclass); head(cross)
head(n.compclass16); head(cross16)

# Calculate weighted indices ----------------------------------------------------------------------------------

head(n.compclass)
sample.mat <- cross[,19:ncol(cross)]
cross <- cross[,1:18]

sample.mat16 <- cross16[,19:ncol(cross16)]
cross16 <- cross16[,1:18]

# Categorization used in Peterborough lab
cross[, aliphatic := ifelse(HC >= 1.5, 1, 0)]
cross[, low.O.unsaturated := ifelse(HC < 1.5 & AImod < 0.5 & OC < 0.5, 1, 0)]
cross[, high.O.unsaturated := ifelse(HC < 1.5 & AImod < 0.5 & OC >= 0.5, 1, 0)]
cross[, polyphenols := ifelse(AImod > 0.5 & AImod < 0.67, 1, 0)]
cross[, condensed.aromatics := ifelse(AImod >= 0.67, 1, 0)]
cross16[, aliphatic := ifelse(HC >= 1.5, 1, 0)]
cross16[, low.O.unsaturated := ifelse(HC < 1.5 & AImod < 0.5 & OC < 0.5, 1, 0)]
cross16[, high.O.unsaturated := ifelse(HC < 1.5 & AImod < 0.5 & OC >= 0.5, 1, 0)]
cross16[, polyphenols := ifelse(AImod > 0.5 & AImod < 0.67, 1, 0)]
cross16[, condensed.aromatics := ifelse(AImod >= 0.67, 1, 0)]

# Categorization used in Florida lab
cross[, fl_condensed.aromatics := ifelse(AImod > 0.66, 1, 0)]
cross[, fl_polyphenolic := ifelse(AImod > 0.5 & AImod <= 0.66 &
                                    fl_condensed.aromatics == 0, 1, 0)]
cross[, fl_highly.unsat.phenolic := ifelse(AImod <= 0.5 & HC < 1.5 &
                                             fl_condensed.aromatics == 0 &
                                             fl_polyphenolic == 0, 1, 0)]
cross[, fl_aliphatic := ifelse(HC <= 2 & HC >= 1.5 & N == 0 &
                                 fl_condensed.aromatics == 0 &
                                 fl_polyphenolic == 0 & fl_highly.unsat.phenolic == 0, 1, 0)]
cross[, fl_peptide := ifelse(HC <= 2 & HC >= 1.5 & N > 0 &
                               fl_condensed.aromatics == 0 &
                               fl_polyphenolic == 0 & fl_highly.unsat.phenolic == 0 &
                               fl_aliphatic == 0, 1, 0)]
cross[, fl_sugar := ifelse(OC > 0.9 &
                             fl_condensed.aromatics == 0 &
                             fl_polyphenolic == 0 & fl_highly.unsat.phenolic == 0 &
                             fl_aliphatic == 0 & fl_peptide == 0, 1, 0)]
cross16[, fl_condensed.aromatics := ifelse(AImod > 0.66, 1, 0)]
cross16[, fl_polyphenolic := ifelse(AImod > 0.5 & AImod <= 0.66 &
                                    fl_condensed.aromatics == 0, 1, 0)]
cross16[, fl_highly.unsat.phenolic := ifelse(AImod <= 0.5 & HC < 1.5 &
                                             fl_condensed.aromatics == 0 &
                                             fl_polyphenolic == 0, 1, 0)]
cross16[, fl_aliphatic := ifelse(HC <= 2 & HC >= 1.5 & N == 0 &
                                 fl_condensed.aromatics == 0 &
                                 fl_polyphenolic == 0 & fl_highly.unsat.phenolic == 0, 1, 0)]
cross16[, fl_peptide := ifelse(HC <= 2 & HC >= 1.5 & N > 0 &
                               fl_condensed.aromatics == 0 &
                               fl_polyphenolic == 0 & fl_highly.unsat.phenolic == 0 &
                               fl_aliphatic == 0, 1, 0)]
cross16[, fl_sugar := ifelse(OC > 0.9 &
                             fl_condensed.aromatics == 0 &
                             fl_polyphenolic == 0 & fl_highly.unsat.phenolic == 0 &
                             fl_aliphatic == 0 & fl_peptide == 0, 1, 0)]

# add as a category
cross[AImod > 0.66, fl.categories := "condensed.aromatics"]
cross[AImod > 0.5 & AImod <= 0.66 &
        fl_condensed.aromatics == 0, fl.categories := "polyphenolic"]
cross[AImod <= 0.5 & HC < 1.5 &
        fl_condensed.aromatics == 0 &
        fl_polyphenolic == 0, fl.categories := "highly.unsat.phenolic"]
cross[HC <= 2 & HC >= 1.5 & N == 0 &
        fl_condensed.aromatics == 0 &
        fl_polyphenolic == 0 & fl_highly.unsat.phenolic == 0, fl.categories := "aliphatic"]
cross[HC <= 2 & HC >= 1.5 & N > 0 &
        fl_condensed.aromatics == 0 &
        fl_polyphenolic == 0 & fl_highly.unsat.phenolic == 0 &
        fl_aliphatic == 0, fl.categories := "peptide"]
cross[OC > 0.9 &
        fl_condensed.aromatics == 0 &
        fl_polyphenolic == 0 & fl_highly.unsat.phenolic == 0 &
        fl_aliphatic == 0 & fl_peptide == 0, fl.categories := "sugar"]
cross16[AImod > 0.66, fl.categories := "condensed.aromatics"]
cross16[AImod > 0.5 & AImod <= 0.66 &
        fl_condensed.aromatics == 0, fl.categories := "polyphenolic"]
cross16[AImod <= 0.5 & HC < 1.5 &
        fl_condensed.aromatics == 0 &
        fl_polyphenolic == 0, fl.categories := "highly.unsat.phenolic"]
cross16[HC <= 2 & HC >= 1.5 & N == 0 &
        fl_condensed.aromatics == 0 &
        fl_polyphenolic == 0 & fl_highly.unsat.phenolic == 0, fl.categories := "aliphatic"]
cross16[HC <= 2 & HC >= 1.5 & N > 0 &
        fl_condensed.aromatics == 0 &
        fl_polyphenolic == 0 & fl_highly.unsat.phenolic == 0 &
        fl_aliphatic == 0, fl.categories := "peptide"]
cross16[OC > 0.9 &
        fl_condensed.aromatics == 0 &
        fl_polyphenolic == 0 & fl_highly.unsat.phenolic == 0 &
        fl_aliphatic == 0 & fl_peptide == 0, fl.categories := "sugar"]

# Rossel, P. E.;Stubbins,A.; Hach,P.F.; Dittmar, T. 
# Bioavailability and molecular composition of dissolved organic matter from a diffuse hydrothermal system. Mar. Chem. 2015, 177, 257−266
# Kellerman, A. M.; Kothawala, D. N.; Dittmar, T.; Tranvik, L. J. 
# Persistence of dissolved organic matter in lakes related to its molecular characteristics. Nat. Geosci. 2015, 8 (6), 454−457.

# add as a category
cross[HC >= 1.5, categories := "aliphatic"]
cross[HC < 1.5 & AImod < 0.5 & OC < 0.5, categories := "low O unsaturated"]
cross[HC < 1.5 & AImod < 0.5 & OC >= 0.5, categories := "high O unsaturated"]
cross[AImod >= 0.5 & AImod < 0.67, categories := "polyphenols"]
cross[AImod >= 0.67, categories := "condensed aromatics"]
cross16[HC >= 1.5, categories := "aliphatic"]
cross16[HC < 1.5 & AImod < 0.5 & OC < 0.5, categories := "low O unsaturated"]
cross16[HC < 1.5 & AImod < 0.5 & OC >= 0.5, categories := "high O unsaturated"]
cross16[AImod >= 0.5 & AImod < 0.67, categories := "polyphenols"]
cross16[AImod >= 0.67, categories := "condensed aromatics"]

# add a few indices
cross[, DBE_O := DBE - O]
cross[, GFE := 60.3 - 28.5 * NOSC]
cross16[, DBE_O := DBE - O]
cross16[, GFE := 60.3 - 28.5 * NOSC]

# Use Ryan's code -->
nsample <- ncol(sample.mat)
nformula <- nrow(sample.mat)
sample.mat <- setDF(sample.mat)

weighted<-as.data.frame(matrix(nrow = nsample))
for(i in 1:nsample){
  weighted$sample[i]<-colnames(sample.mat[i])
  weighted$mz[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE)*cross$mz),na.rm = TRUE)
  weighted$HC[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE))*cross$HC,na.rm = TRUE)
  weighted$OC[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE))*cross$OC,na.rm = TRUE)
  weighted$DBE[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE))*cross$DBE,na.rm = TRUE)
  weighted$DBE_O[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE))*cross$DBE_O,na.rm = TRUE)
  weighted$AImod[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE))*cross$AImod,na.rm = TRUE)
  weighted$NOSC[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE))*cross$NOSC,na.rm = TRUE)
  weighted$GFE[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE))*cross$GFE,na.rm = TRUE)
  # 1st compound classes
  weighted$aliphatic[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE))*cross$aliphatic,na.rm = TRUE)*100
  weighted$low.O.unsaturated[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE))*(cross$low.O.unsaturated),na.rm = TRUE)*100
  weighted$high.O.unsaturated[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE))*(cross$high.O.unsaturated),na.rm = TRUE)*100
  weighted$polyphenols[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE))*(cross$polyphenols),na.rm = TRUE)*100
  weighted$condensed.aromatics[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE))*(cross$condensed.aromatics),na.rm = TRUE)*100
  #2nd compound classes
  weighted$fl_condensed.aromatics[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE))*cross$fl_condensed.aromatics,na.rm = TRUE)*100
  weighted$fl_polyphenolic[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE))*(cross$fl_polyphenolic),na.rm = TRUE)*100
  weighted$fl_highly.unsat.phenolic[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE))*(cross$fl_highly.unsat.phenolic),na.rm = TRUE)*100
  weighted$fl_aliphatic[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE))*(cross$fl_aliphatic),na.rm = TRUE)*100
  weighted$fl_peptide[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE))*(cross$fl_peptide),na.rm = TRUE)*100
  weighted$fl_sugar[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE))*(cross$fl_sugar),na.rm = TRUE)*100
  #others
  weighted$H[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE)*cross$H),na.rm = TRUE)
  weighted$O[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE)*cross$O),na.rm = TRUE)
  weighted$C[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE)*cross$C),na.rm = TRUE)
  weighted$N[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE))*(cross$N),na.rm = TRUE)
  weighted$S[i]<-sum((sample.mat[,i]/sum(as.numeric(sample.mat[,i]),na.rm = TRUE))*(cross$S),na.rm = TRUE)
  weighted$sum[i]<-sum(as.numeric(sample.mat[,i]),na.rm = TRUE)*1e-9
  weighted$nformula[i]<-nformula-sum(is.na(sample.mat[,i]))
}

nsample <- ncol(sample.mat16)
nformula <- nrow(sample.mat16)
sample.mat16 <- setDF(sample.mat16)

weighted16<-as.data.frame(matrix(nrow = nsample))
for(i in 1:nsample){
  weighted16$sample[i]<-colnames(sample.mat16[i])
  weighted16$mz[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE)*cross$mz),na.rm = TRUE)
  weighted16$HC[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE))*cross$HC,na.rm = TRUE)
  weighted16$OC[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE))*cross$OC,na.rm = TRUE)
  weighted16$DBE[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE))*cross$DBE,na.rm = TRUE)
  weighted16$DBE_O[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE))*cross$DBE_O,na.rm = TRUE)
  weighted16$AImod[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE))*cross$AImod,na.rm = TRUE)
  weighted16$NOSC[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE))*cross$NOSC,na.rm = TRUE)
  weighted16$GFE[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE))*cross$GFE,na.rm = TRUE)
  # 1st compound classes
  weighted16$aliphatic[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE))*cross$aliphatic,na.rm = TRUE)*100
  weighted16$low.O.unsaturated[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE))*(cross$low.O.unsaturated),na.rm = TRUE)*100
  weighted16$high.O.unsaturated[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE))*(cross$high.O.unsaturated),na.rm = TRUE)*100
  weighted16$polyphenols[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE))*(cross$polyphenols),na.rm = TRUE)*100
  weighted16$condensed.aromatics[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE))*(cross$condensed.aromatics),na.rm = TRUE)*100
  #2nd compound classes
  weighted16$fl_condensed.aromatics[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE))*cross$fl_condensed.aromatics,na.rm = TRUE)*100
  weighted16$fl_polyphenolic[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE))*(cross$fl_polyphenolic),na.rm = TRUE)*100
  weighted16$fl_highly.unsat.phenolic[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE))*(cross$fl_highly.unsat.phenolic),na.rm = TRUE)*100
  weighted16$fl_aliphatic[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE))*(cross$fl_aliphatic),na.rm = TRUE)*100
  weighted16$fl_peptide[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE))*(cross$fl_peptide),na.rm = TRUE)*100
  weighted16$fl_sugar[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE))*(cross$fl_sugar),na.rm = TRUE)*100
  #others
  weighted16$H[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE)*cross$H),na.rm = TRUE)
  weighted16$O[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE)*cross$O),na.rm = TRUE)
  weighted16$C[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE)*cross$C),na.rm = TRUE)
  weighted16$N[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE))*(cross$N),na.rm = TRUE)
  weighted16$S[i]<-sum((sample.mat16[,i]/sum(as.numeric(sample.mat16[,i]),na.rm = TRUE))*(cross$S),na.rm = TRUE)
  weighted16$sum[i]<-sum(as.numeric(sample.mat16[,i]),na.rm = TRUE)*1e-9
  weighted16$nformula[i]<-nformula-sum(is.na(sample.mat16[,i]))
}
# --> Ryan code ends

# 3. Remove outliers -------------------------------------------------------------------------------------------

# 2015 --------------------------------------------------------------------------------------------------------
# extract "community matrix"
set.seed(3)
# dim 1 = row = sample, dim 2 = col = species
setDF(cross)
com.mat <- data.frame(molecular.formula = cross[,"molecular.formula"], sample.mat)
row.names(com.mat) <- com.mat$molecular.formula
com.mat <- com.mat[,-1]
com.mat <- t(as.matrix(com.mat))
com.mat <- as.data.frame(com.mat)

# Remove GW
# Add meta data
meta <- read.csv("./Data/LaRomaine_MasterFile_withDist_2021-02-22.csv", sep = ",", stringsAsFactors = F)
setDT(meta); setDT(com.mat, keep.rownames = "sample.name")
meta[,ft.samples := str_replace(sample.name, "[-]",".")]

# duplicate "RO2.22"
meta <- rbind(meta, meta[sample.name == "RO2-22",])
meta[nrow(meta), ft.samples := "RO2.22.r2"]

nrow(com.mat[!(sample.name %in% meta$ft.samples),]) # all

# Merge
com.mat <- com.mat[meta[ft.samples %in% com.mat$sample.name], sample.type.year := i.sample.type.year, on = c("sample.name" = "ft.samples")]
# Remove groundwater
com.mat <- com.mat[sample.type.year != "Wellwater",]

# Format back to matrix
com.mat[, sample.type.year := NULL]
com.mat <- setDF(com.mat[,-1], rownames = com.mat$sample.name)

# Run NMDS
nmds <- metaMDS(com.mat, k = 2) # default is bray curtis
stressplot(nmds)
ordiplot(nmds, display = 'sites')

nmds.df <- as.data.frame(nmds$points)
nmds.df$samples <- row.names(nmds.df)

setDT(nmds.df)
# merge
nmds.df <- nmds.df[meta[ft.samples %in% nmds.df$samples], , on = c("samples" = "ft.samples")]

plot_ly(nmds.df, mode = "markers", type = "scatter", x = ~MDS1, y = ~MDS2, text = ~samples, color = ~sample.type.year)
# 2015: LRH08 and LRH15 (LR08 is a bit odd but we keep it)

nmds.df <- nmds.df[!(samples %in% c("LRH08","LRH15")),]
nmds.df <- nmds.df[!(samples %in% c("RO2.22-r2", "PR19.1","PR23.1"))]
#nmds.df[, Season := factor(campaign, levels = c(1,2), labels = c("Spring","Summer"))]
nmds.df[, Season := factor(campaign, levels = c(1,2), labels = c("June","August"))]
nmds.df[, sample.type.year := factor(sample.type.year, levels = c("Wellwater", "Stream",  "Tributary", "HeadwaterLakes",
                                                                  "PRLake" ,"IslandLake",
                                                                  "Upriver", "RO2", "Downriver", "Marine" ),
                                     labels = c("Groundwater","Stream","Tributary", "Riverine Lakes",
                                                "Headwater Ponds", "Lake", "Upriver","RO2", "Downriver", "Estuary"))]
# nmds.df[, sample.type.year := factor(sample.type.year, levels = c("Wellwater", "Stream",  "Tributary", "HeadwaterLakes",
#                                                                   "PRLake" ,"IslandLake",
#                                                                   "Upriver", "RO2", "Downriver", "Marine"),
#                                      labels = c("Groundwater","Petite Romaine","Tributary", "Lake",
#                                                 "Lake", "Lake", "Upstream","RO2", "Downstream", "Estuary"))]

colvec <- c(#"#F9795DFF", #"#F9795DFF", #Wellwater,
                      #"#FCFDBFFF", #"#FEC589FF", #Soilwater
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

# Plot hulls
# Make hulls around DNA and RNA
find_hull <- function(x, axes){x[chull(x[,paste0("MDS",axes[1])], x[,paste0("MDS",axes[2])]),]}
hulls <- plyr::ddply(nmds.df, .(sample.type.year), find_hull, axes = c(1,2))

# Calculate centroids
cent.df <- nmds.df[, .(MDS1 = mean(MDS1),
            MDS2 = mean(MDS2)), by = .(sample.type.year, Season)]

(p <- ggplot(cent.df, aes(x = MDS1, y = MDS2)) +
  theme_bw() +
  geom_hline(yintercept = 0, colour = "gray50") +
  geom_vline(xintercept = 0, colour = "gray50") +
    geom_polygon(data = hulls, alpha = 0.3, aes(colour = sample.type.year, fill = sample.type.year)) +
    geom_point(aes(fill = sample.type.year, shape = Season), colour = "gray50", size = 3) +
    scale_fill_manual(values = colvec, name = "Sample type") +
    scale_colour_manual(values = colvec, name = "Sample type") +
    scale_shape_manual(values = c(21, 23, 24)) +
    annotate(geom = "text", x = min(nmds.df$MDS1) + .1,
             y = max(nmds.df$MDS2), label = paste0("Stress = ", round(nmds$stress, 2))) +
    guides(fill = guide_legend(override.aes = list(shape = 21))))

#ggsave("./Figures/NMDS_FTall_2015.png", p)

nmds.df <- nmds.df %>% dplyr::select(-MDS1, -MDS2)
# save as 2015 meta file
write.table(nmds.df, "./Data/Summary/2015_metadata.csv", sep = ",", row.names = F)

# 2016 --------------------------------------------------------------------------------------------------------
# extract "community matrix"
set.seed(3)
# dim 1 = row = sample, dim 2 = col = species
setDF(cross16)
com.mat <- data.frame(molecular.formula = cross16$molecular.formula, sample.mat16)
row.names(com.mat) <- com.mat$molecular.formula
com.mat <- com.mat[,-1]
com.mat <- t(as.matrix(com.mat))
com.mat <- as.data.frame(com.mat)

# Remove GW
# Add meta data
meta <- read.csv("./Data/LaRomaine_MasterFile_withDist_2021-02-22.csv", sep = ",", stringsAsFactors = F)
setDT(meta); setDT(com.mat, keep.rownames = "sample.name")
meta[,ft.samples := str_replace(sample.name, "[-]",".")]

# some sample names do not match, correct in meta file
meta[grep("L328", ft.samples), ft.samples := str_replace(ft.samples, "_", ".")]
meta[grep("L329", ft.samples), ft.samples := str_replace(ft.samples, "_", ".")]
# assuming that RO1-10P is a hypolimnion sample
meta[ft.samples == "RO1.10D", ft.samples := "RO1.10P"]

com.mat[!(sample.name %in% meta$ft.samples),] # all

# Merge
com.mat <- com.mat[meta[ft.samples %in% com.mat$sample.name], c("sample.type.year",
                                                                "campaign") := list(i.sample.type.year, i.campaign), on = c("sample.name" = "ft.samples")]
temp <- meta[ft.samples %in% com.mat$sample.name,]

median(temp[!(sample.type %in% c("Soilwater","Marine","Tributary")),]$DOC) # 6
sd(temp[!(sample.type %in% c("Soilwater","Marine","Tributary")),]$DOC)
median(temp[(sample.type %in% c("Marine")),]$DOC) # 1.5
median(temp[(sample.type %in% c("Tributary")),]$DOC) # 9.4
median(temp[(sample.type %in% c("Soilwater")),]$DOC) # 9.4

ggplot(temp[sample.type != "Soilwater",], aes(x = sample.type, y = DOC)) +
  geom_boxplot()

# Remove groundwater
com.mat <- com.mat[sample.type.year != "Wellwater" & sample.type.year != "Soilwater" & sample.type.year != "Hyporheic" & sample.type.year != "Deep",]
com.mat <- com.mat[campaign != 3,]

# Format back to matrix
com.mat[, sample.type.year := NULL]
com.mat <- setDF(com.mat[,-1], rownames = com.mat$sample.name)

nmds <- metaMDS(com.mat, k = 2) # default is bray curtis
stressplot(nmds)
ordiplot(nmds, display = 'sites')

nmds.df <- as.data.frame(nmds$points)
nmds.df$samples <- row.names(nmds.df)

# Add meta data
setDT(meta); setDT(nmds.df)

# merge
nmds.df <- nmds.df[meta[ft.samples %in% nmds.df$samples], , on = c("samples" = "ft.samples")]

plot_ly(nmds.df, mode = "markers", type = "scatter", x = ~MDS1, y = ~MDS2, text = ~samples, color = ~sample.type.year)
# 2016:
# Obviously, SW34 and SW36 are extreme outliers
# HW42, M18, T48, PR32
# RO1.10, RO2.40, RO2.51
# LR55, LR57, LR58, LR59, LR60, LR63, LR67, LR70
# Outlier decision by Francois Guillemette and Tristy Vick-Majors

# T47 listed as potentially contaminated, but not removed? Seems a bit odd in NMDS.
nmds.df <- nmds.df[!(samples %in% c("SW34","SW36","HW42",
                                    "M18","T48","PR32",
                                    "RO1.10", "RO2.40", "RO2.51",
                                    "LR55","LR57","LR58","LR59", "LR60", "LR63","LR67","LR70")),]

nmds.df[, Season := factor(campaign, levels = c(1,2,3), labels = c("Spring","Summer","Autumn"))]
nmds.df[, sample.type.year := factor(sample.type.year, levels = c("Soilwater", "Hyporheic", 
                                                                  "Stream",  "Tributary", "HeadwaterLakes",
                                                                  "Lake",
                                                                  "Upriver", "RO2", "RO1", "Deep", "Downriver", "Marine" ),
                                     labels = c("Soilwater", "Hyporheic", 
                                                "Stream",  "Tributary", "Riverine Lakes",
                                                "Lake", "Upriver","RO1", "RO2", "Hypolimnion",
                                                "Downriver", "Estuary"))]

colvec <- c(
            #"#FCFDBFFF", #"#FEC589FF", #Soilwater,
            #"#F9795DFF", #"#F9795DFF", #Wellwater/Hyporheic
            "skyblue", #Stream
            "#AD347CFF",# Tributary, 
            "palegreen", #Riverine Lakes, 
            "#FDE725FF",# Lake, 
            "#1F9F88FF", # Upriver, 
            "salmon",
            "orchid", #"#471063FF", #Reservoir,
            #"navy", # deep
            #"pink",
            #'navy',
            "#375A8CFF", #Downriver,
            "gray40") #Estuary)

# Plot hulls
find_hull <- function(x, axes){x[chull(x[,paste0("MDS",axes[1])], x[,paste0("MDS",axes[2])]),]}
hulls <- ddply(nmds.df, .(sample.type.year), find_hull, axes = c(1,2))

# Calculate centroids
cent.df <- nmds.df[, .(MDS1 = mean(MDS1),
                       MDS2 = mean(MDS2)), by = .(sample.type.year, Season)]

(p <- ggplot(cent.df, aes(x = MDS1, y = MDS2)) +
    theme_bw() +
    geom_hline(yintercept = 0, colour = "gray50") +
    geom_vline(xintercept = 0, colour = "gray50") +
    geom_polygon(data = hulls, alpha = 0, aes(colour = sample.type.year), fill = "white") +
    geom_point(aes(fill = sample.type.year, shape = Season), colour = "gray50", size = 3) +
    scale_fill_manual(values = colvec, name = "Sample type") +
    scale_colour_manual(values = colvec, name = "Sample type") +
    scale_shape_manual(values = c(21, 23, 24)) +
    annotate(geom = "text", x = min(nmds.df$MDS1) + .1,
             y = max(nmds.df$MDS2), label = paste0("Stress = ", round(nmds$stress, 2))) +
    guides(fill = guide_legend(override.aes = list(shape = 21))))
#ggsave("./Figures/NMDS_FTall_2016.png", p)

nmds.df <- nmds.df %>% dplyr::select(-MDS1, -MDS2)

# save as 2016 meta file
write.table(nmds.df, "./Data/Summary/2016_metadata.csv", sep = ",", row.names = F)

# Remove outliers from cross table
# 2015
sample.mat[,c("LRH08","LRH15")] <- NULL
# removed two samples, now
ncol(sample.mat)
# 89 samples

# 2016
sample.mat16[,c("SW34","SW36","HW42",
           "M18","T48","PR32",
           "RO1.10", "RO2.40", "RO2.51",
           "LR55","LR57","LR58","LR59", "LR60", "LR62", "LR63","LR67","LR70")] <- NULL
# removed 17 samples, now
ncol(sample.mat16)
# 93 samples

# merge samples and cross table
cross <- cross %>% dplyr::select(molecular.formula:S, Class, Count, AI:NOSC, DBE_O:GFE, 
                        Compound.category, fl.categories, categories, aliphatic:fl_sugar)
cross16 <- cross16 %>% dplyr::select(molecular.formula:S, Class, Count, AI:NOSC, DBE_O:GFE, 
                                 Compound.category, fl.categories, categories, aliphatic:fl_sugar)

cross <- cbind(cross, sample.mat)
cross16 <- cbind(cross16, sample.mat16)

# Save output --------------------------------------------------------------------------------------------
write.table(cross, "./Data/FT/crosstable2015_cor.csv", sep =",", dec = ".", row.names = F)
write.table(cross16, "./Data/FT/crosstable2016_cor.csv", sep =",", dec = ".", row.names = F)

write.table(weighted[,-1], "./Data/FT/weighted2015.csv", sep =",", dec = ".", row.names = F)
write.table(weighted16[,-1], "./Data/FT/weighted2016.csv", sep =",", dec = ".", row.names = F)

write.table(n.compclass, "./Data/FT/compoundtable2015_cor.csv", sep =",", dec = ".", row.names = F)
write.table(n.compclass16, "./Data/FT/compoundtable2016_cor.csv", sep =",", dec = ".", row.names = F)






# # calculate sum of rel. abundance to check if we have similar cross tables
# setDT(tristy.nc); setDT(cross16)
# melt.t <- melt.data.table(tristy.nc, id.vars = "Molecular.Formula", 
#                 measure.vars = 19:ncol(tristy.nc),variable.name = "sample", value.name = "rel.abun")
# melt.c <- melt.data.table(cross16, id.vars = "molecular.formula", 
#                           measure.vars = 19:ncol(cross16),variable.name = "sample", value.name = "rel.abun")
# 
# c<-melt.c[,.(sumc = sum(rel.abun)), by = .(sample)]
# t<-melt.t[,.(sumt = sum(rel.abun)), by = .(sample)]
# 
# check <- cbind(c,t)
# check[sumc != sumt,]
# identical!

sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.2 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
# 
# locale:
#   [1] LC_CTYPE=en_CA.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_CA.UTF-8        LC_COLLATE=en_CA.UTF-8    
# [5] LC_MONETARY=en_CA.UTF-8    LC_MESSAGES=en_CA.UTF-8   
# [7] LC_PAPER=en_CA.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_CA.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] plyr_1.8.6        mgcv_1.8-33       nlme_3.1-149      plotly_4.9.2.1   
# [5] vegan_2.5-6       lattice_0.20-41   permute_0.9-5     forcats_0.5.0    
# [9] dplyr_1.0.2       purrr_0.3.4       readr_1.4.0       tidyr_1.1.2      
# [13] tibble_3.0.4      ggplot2_3.3.2     tidyverse_1.3.0   stringr_1.4.0    
# [17] data.table_1.14.0 readxl_1.3.1     
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.6         lubridate_1.7.9.2  assertthat_0.2.1   digest_0.6.27     
# [5] R6_2.5.0           cellranger_1.1.0   backports_1.2.0    reprex_0.3.0      
# [9] httr_1.4.2         pillar_1.4.7       rlang_0.4.9        lazyeval_0.2.2    
# [13] rstudioapi_0.13    Matrix_1.2-18      labeling_0.4.2     splines_4.0.3     
# [17] htmlwidgets_1.5.2  munsell_0.5.0      broom_0.7.2        compiler_4.0.3    
# [21] modelr_0.1.8       pkgconfig_2.0.3    htmltools_0.5.0    tidyselect_1.1.0  
# [25] fansi_0.4.1        viridisLite_0.3.0  crayon_1.3.4       dbplyr_2.0.0      
# [29] withr_2.3.0        MASS_7.3-53.1      grid_4.0.3         jsonlite_1.7.1    
# [33] gtable_0.3.0       lifecycle_0.2.0    DBI_1.1.0          magrittr_2.0.1    
# [37] scales_1.1.1       cli_2.2.0          stringi_1.5.3      farver_2.0.3      
# [41] fs_1.5.0           xml2_1.3.2         ellipsis_0.3.1     generics_0.1.0    
# [45] vctrs_0.3.5        RColorBrewer_1.1-2 tools_4.0.3        glue_1.4.2        
# [49] crosstalk_1.1.0.1  hms_0.5.3          yaml_2.2.1         parallel_4.0.3    
# [53] colorspace_2.0-0   cluster_2.1.0      rvest_0.3.6        haven_2.3.1   
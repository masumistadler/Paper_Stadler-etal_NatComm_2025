# Randomize data and calculate CI for slopes ------------------------------------------------------------------

# R set-up ----------------------------------------------------------------------------------------------------
### Packages -------------------------------------------------------------------------------
.libPaths( c( .libPaths(), "/home/mstadler/projects/def-pauldel/R/x86_64-pc-linux-gnu-library/4.2") )
pckgs <- list("data.table", "tidyverse", # wrangling
              "plyr",
              "doMC","foreach")
### Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)
# Many packages have to be installed through Bioconductor, please refer to the package websites

### Load
invisible(lapply(pckgs, library, character.only = T))
rm(pckgs, miss.pckgs)

# Set up parallel environment ------------------------------------------------------------------------------
cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
registerDoMC(cores = cores)

# Read in data ---------------------------------------------------------------------------------------------
present <- readRDS("./Objects/MO_present.rds") %>% setDT()

# Run function ----------------------------------------------------------------------------------------------

random <- ddply(present, .(ID), function(x) {
  
  out <- foreach(i = 1:999, .combine = rbind) %dopar% {
    # set random iteration seed
    set.seed(i)
    # extract data
    df <- data.frame(x = x$travel.time_d, y = x$z.reads)
    # shuffle y over x, needs to be with replacement for bootstrapping
    df$y <- sample(df$y, replace = T, size = nrow(df))
    # do linear regression
    lin <- try(lm(df$y ~ df$x), silent = T)
    # extract run and slope
    if (!inherits(lin, "try-error") & !is.na(lin$coefficients[[2]])) {
    out <- data.frame(run = i, slope = lin$coefficients[2], r2 = summary(lin)$r.squared,
                      p.val = summary(lin)$coefficients[2,4])
    } else {
    out <- data.frame(run = 1, slope = NA, r2 = NA,
                      p.val = NA)
    }
    return(out)
  }
  
  return(out)
}, .parallel = T)

# Save intermediate output -----------------------------------------------------------------------------
saveRDS(random,"./Objects/MO_randomization_output.rds")
write.table(random,"./Objects/MO_randomization_output.csv", sep = ',', dec = ".", row.names = F)

# Calculate confidence interval ------------------------------------------------------------------------
setDT(random)
random <- melt(random, id.vars = 1:2, measure.vars = 3:5, variable.name = "vars", value.name = "value")

df <- random[, .(mean = mean(value, na.rm =T),
                 sd = sd(value, na.rm = T),
                 n = .N), by = .(ID, vars)]

df[, c("se","alpha", "df") := list(sd / sqrt(n),
                                   0.05,
                                   n - 1)]
df[, t.score := qt(p=alpha/2, df = df, lower.tail = F)]
df[, margin.error := t.score * se]
df[, c("upper.bound", "lower.bound") := list(mean + margin.error,
                                             mean - margin.error)]

# Save -----------------------------------------------------------------------------------------------
saveRDS(df,"./Objects/MO_randomization_CI.rds")
write.table(df,"./Objects/MO_randomization_CI.csv", sep = ',', dec = ".", row.names = F)

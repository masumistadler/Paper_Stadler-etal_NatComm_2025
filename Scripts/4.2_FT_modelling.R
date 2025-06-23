# Run modelling decision tree ---------------------------------------------------------------------------------

# R set-up ----------------------------------------------------------------------------------------------------
### Packages -------------------------------------------------------------------------------
.libPaths( c( .libPaths(), "/home/mstadler/projects/def-pauldel/R/x86_64-pc-linux-gnu-library/4.2") )
pckgs <- list("data.table", "tidyverse", # wrangling
              "plyr", "mgcv",
              "doMC","foreach")
### Check if packages are installed, output packages not installed:
(miss.pckgs <- unlist(pckgs)[!(unlist(pckgs) %in% installed.packages()[,"Package"])])
#if(length(miss.pckgs) > 0) install.packages(miss.pckgs)
# Many packages have to be installed through Bioconductor, please refer to the package websites

### Load
invisible(lapply(pckgs, library, character.only = T))
rm(pckgs, miss.pckgs)

### Set up parallel environment ------------------------------------------------------------------------------
cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
registerDoMC(cores = cores)

# Read in data ---------------------------------------------------------------------------------------------
# present <- readRDS("./Objects/FT_present_2024-02-10.rds") %>% setDT()
# # read in randomization results
# out <- readRDS("./Objects/FT_randomization_output_2024-02-10.rds") %>% setDT()


present <- readRDS("./Objects/FT_present.rds") %>% setDT()
# read in randomization results
out <- readRDS("./Objects/FT_randomization_output.rds") %>% setDT()


# Identify IDs for which modelling does not work and remove ------------------------------------------------
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

# remove those IDs that did not work
present <- present[!(ID %in% not.working),]
rm(out, check, not.working)

# Define function ---------------------------------------------------------------------------
localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  # Use .Machine$integer.max instead of Inf for integer
  y <- diff(c(-Inf, x)) > 0L
  # returns logical vector with which numbers in 1st derivative are positive
  rle(y)$lengths
  # returns how often the same value (i.e. TRUE) repeats in a sequence without being interrupted
  # and when there is a switch from T to F
  y <- cumsum(rle(y)$lengths)
  # returns vector with cumulating the sum of TRUE and FALSE observations = returns location of switch from T to F
  y <- y[seq.int(1L, length(y), 2L)]
  # samples the every second location (i.e. switch to TRUE)
  
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

closest <- function(x,value){
  b <- x[x < value]
  if(length(b) == 0){
    below <- NA
  } else {
    below <- b[which.min(abs(b-value))]
  }
  
  a <- x[x > value]
  if(length(a) == 0){
    above <- NA
  } else {
    above <- a[which.min(abs(a-value))]
  }
  
  if(length(c(a,b)) == 0L){
    above <- x; below <- x
  }
  
  return(c(below, above))
}

lm.or.gam <- function(x, y) { #, min.x, max.x
  # bin
  raw <- data.table(x,y)
  brks <- c(seq(-50, max(x, na.rm = T), by = 50), max(x, na.rm = T)) %>% unique()
  bins <- data.table(brks, bin.id = 1:length(brks))
  raw[, bin.id := findInterval(x, bins$brks)]
  raw <- raw[order(x),]
  df <- raw[bins, , on = .(bin.id)]
  df <- df[,.(x = unique(brks),
              y = mean(y, na.rm = T),
              n = .N), by = .(bin.id)]
  
  # remove n = 1 if it is NA
  df[is.na(y), n := 0]
  df <- df[!is.na(y),]
  
  if(nrow(df) >= 3){
    
    #y <- (out$y /avg)
    y <- df$y
    x <- df$x
    n <- df$n
    
    # model decision tree start --------------------------------------------------------
    lin <- try(gam(y ~ x, weights = n/mean(n)), silent = T)
    
    if (!inherits(lin, "try-error") & !is.na(coef(lin)[[2]]) & sd(y) != 0) {
      # try all models
      poly.2 <- try(gam(y ~ poly(x, 2), weights = n/mean(n)), silent = T)
      poly.3 <- try(gam(y ~ poly(x, 3), weights = n/mean(n)), silent = T)
      gam <- try(gam(y ~ s(x), select = T, weights = n/mean(n)), silent = T) #constrain wiggliness
      
      model.ls <-
        list(
          lin = lin,
          poly.2 = poly.2,
          poly.3 = poly.3,
          gam = gam
        )
      # error handling, remove models that didn't work
      errors <-
        names(which(sapply(model.ls, inherits, 'try-error') == T))
      model.ls <- model.ls[!(names(model.ls) %in% errors)]
      
      # which has the smallest AIC?
      aic.df <-
        data.frame(models = as.vector(names(sapply(model.ls, AIC))),
                   AIC = as.vector(sapply(model.ls, AIC))) %>% arrange(AIC) %>% filter(!is.infinite(AIC))
      
      if (nrow(aic.df) != 0L) {
        # if any model goes through filtering, do...
        # if the best model is linear, keep the linear model
        if (aic.df[1,]$models == "lin") {
          #re-fit with lm
          model <- lm(y ~ x)
        } else {
          # if it's not linear, test against linearity
          best.model <-
            model.ls[names(model.ls) %in% aic.df[1,]$models][[1]]
          chisq.test <- try(mgcv::anova.gam(lin,        
                                            best.model, test = "Chisq"), silent = T)
          # aov.test <- anova(lin,
          #                   best.model)
          
          
          if (inherits(chisq.test, "try-error")){
            # non-linearity is not justified
            model <- lm(y ~ x)
          } else if(is.na(chisq.test$`Pr(>Chi)`[2]) |
                    chisq.test$`Pr(>Chi)`[2] > 0.05){
            # non linearity is not justified
            model <- lm(y ~ x)
          } else {
            # non-linearity is justified
            # re-fit with lm if poly
            if (formula(best.model) == "y ~ poly(x, 2)") {
              model <- lm(y ~ poly(x, 2))
            } else if (formula(best.model) == "y ~ poly(x, 3)") {
              model <- lm(y ~ poly(x, 3))
            } else {
              # if gam
              model <- best.model
            }
          }
        }
      } else {
        # if there is no model that doesn't have a AIC value, return NA
        model <- NA
      }
    } else {
      # if even the linear model doesn't work return NA
      model <- NA
    }
    # model decision tree finish -------------------------------------------------------
  } else {model <- NA}
  if (length(model) > 1L) {
    # if a real model was returned do...
    # we have our best model now, extract coefficients and get peak values
    # model peak extraction start -----------------------------------------------------
    if (grepl("s(", deparse(formula(model)), fixed = T) |
        grepl("poly(", deparse(formula(model)), fixed = T)) {
      # if GAM or poly then do ...
      if(max(x) - min(x) < 10){step <- 1} else {step <- 10}
      # make a x-axis with a higher frequency of points to make rollapply work
      frq.x <- c(seq(min(x), max(x), by = step), max(x)) %>% unique()
      
      # predict y-values along high frequency x axis using the best model
      newd <- data.frame(x = frq.x)
      y.pred <- as.numeric(predict(model, newdata = newd))
      
      # apply peak detection function
      i.max <- localMaxima(y.pred)
      
      # error handling code, if no peaks fill final data frames with NA
      if (length(i.max) != 0L) {
        # if there are peaks
        # Find inflection points
        deriv3rd <- diff(diff(diff(y.pred)))
        d <- sign(deriv3rd)
        d <- cumsum(rle(d)$lengths)
        inf.pts <- d[seq.int(1L, length(d), 2L)]
        
        # Take inflection points before and after each peak
        df.ls <- list()
        
        for (i in 1:length(i.max)) {
          clos.inf <- closest(inf.pts, i.max[i])
          
          temp <-
            data.frame(
              no.peak = rep(length(i.max), times = length(clos.inf)),
              peak = rep(i, times = length(clos.inf)),
              peak.x = rep(frq.x[i.max[i]], times = length(clos.inf)),
              peak.y = rep(y.pred[i.max[i]], times = length(clos.inf)),
              closest = c("above", "below"),
              infp.x = frq.x[clos.inf],
              infp.y = y.pred[clos.inf],
              infp.slope = diff(y.pred)[clos.inf],
              stringsAsFactors = F
            )
          df.ls[[i]] <- temp
        }
        
        pk.df <- do.call(rbind, df.ls)
        
        #output
        frq.data <- data.frame(
          x.pred = frq.x,
          y.pred = y.pred,
          stringsAsFactors = F
        )
      }
    }
    # model peak extraction end ------------------------------------------------------
    
    # get model coefficients depending on the model used -----------------------------
    if (grepl("s(", deparse(formula(model)), fixed = T)) {
      # gather final output for GAM
      model.df <- data.frame(
        model = deparse(formula(model)),
        direction = ifelse(lin$coefficients[2] > 0, 1, 2),
        intercept = model$coefficients[1],
        lm.slope = lin$coefficients[2],
        x.start = min(x),
        x.end = max(x),
        mean.y = mean(y, na.rm = T),
        median.y = median(y, na.rm = T),
        sd.y = sd(y, na.rm = T),
        var.y = var(y, na.rm = T),
        dev.expl = summary(model)$dev.expl,
        adj.r.sq = summary(model)$r.sq,
        p.val = summary(model)$s.table[, 4]
      )
    } else {
      # if it's LM or poly get these....
      # same statistics as GAM just linear equivalents
      model.df <- data.frame(
        model = deparse(formula(model)),
        direction = ifelse(lin$coefficients[2] > 0, 1, 2),
        intercept = model$coefficients[1],
        lm.slope = lin$coefficients[2],
        x.start = min(x),
        x.end = max(x),
        mean.y = mean(y, na.rm = T),
        median.y = median(y, na.rm = T),
        sd.y = sd(y, na.rm = T),
        var.y = var(y, na.rm = T),
        dev.expl = summary(model)$r.squared,
        adj.r.sq = summary(model)$adj.r.squared,
        p.val = mean(summary(model)$coefficients[-1, 4]) # calculate the mean p-val across polys
      )
    }
  }
  try(rm(model.ls))
  # export list
  model.ls <- list(
    raw = raw,
    binned = df,
    frq.data = if(exists("frq.data")) {
      frq.data
    } else {
      data.frame(x.pred = NA,
                 y.pred = NA)
    },
    model = model,
    model.df = if(exists("model.df")){
      model.df
    } else {
      data.frame(
        model = NA,
        direction = NA,
        intercept = NA,
        lm.slope = NA,
        x.start = min(x),
        x.end = max(x),
        mean.y = mean(y, na.rm = T),
        median.y = median(y, na.rm = T),
        sd.y = sd(y, na.rm = T),
        var.y = var(y, na.rm = T),
        dev.expl = NA,
        adj.r.sq = NA,
        p.val = NA
      )
    },
    pks.df = if(exists("pk.df")){
      pk.df
    } else {
      data.frame(
        no.peak = 0,
        peak = NA,
        peak.x = NA,
        peak.y = NA,
        closest = NA,
        infp.x = NA,
        infp.y = NA,
        infp.slope = NA,
        auc = NA
      )
    }
  )
  # remove files
  try(rm(pk.df, model.df, frq.data, model, lin, gam, poly2, poly3), silent = T)
  # return everything
  return(model.ls)
}

# d <- present[ID == levels(factor(present$ID))[14953],]
#  x <- d$travel.time_d
#  y <- d$z.peak.int

# Apply function ----------------------------------------------------------------------------
model.ls <- dlply(present, .(ID), function(x) {
  model.ls <- lm.or.gam(x$travel.time_d, x$z.peak.int) 
  return(model.ls)
}, .parallel = T)

# Extract results --------------------------------------------------------------------------
# get model stats data frame
model.df <- ldply(model.ls, function(x){
  x$model.df
}) %>% setDT()

# Peaks data frame
pks.df <- ldply(model.ls, function(x){
  x$pks.df
}) %>% setDT()

# Save ------------------------------------------------------------------------------------
saveRDS(model.ls, paste0("./Objects/FT_model.ls_binned_", Sys.Date(),".rds"))
saveRDS(model.df, paste0("./Objects/FT_model.df_binned_", Sys.Date(),".rds"))
saveRDS(pks.df, paste0("./Objects/FT_pks.df_binned_", Sys.Date(),".rds"))

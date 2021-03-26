# Front-end ----
# . Library load ----
  library(Rcpp)
  library(myCpp) # devtools::install_github("danStich/myCpp")
  library(rstan)
  library(snowfall)
  library(rlecuyer)
  options(scipen=999)
  
# . Set options ----
  options(mc.cores = parallel::detectCores())  
  rstan_options(auto_write = TRUE)
  Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

# Data load ----
# Run the data script if not done yet
  ### Throws warnings about conversion of weight and 
  ### length to numeric from factor, but okay otherwise
  if(!exists('flathead_mort')) source('data/data_prep.R')
  
# Get the data for this species
  # Length, weight, age data
    gdat <- channel_growth
  # Catch curve data
    cdat <- data.frame(
      channel_mort[ ,1],
      stack(channel_mort[ ,2:ncol(channel_mort)])
    )  
    cdat <- cdat[!is.na(cdat$values), ]
    cdat$gear <- as.numeric(as.factor(cdat$ind))
    
# Growth and mortality ----
# . Stan model(s) ----
# Package the data for stan
  vb_data = list(
    length = gdat$Length,
    age = gdat$Age,
    n1 = nrow(gdat),
    log10length = gdat$lLength,
    log10weight = gdat$lWeight,
    N = log(cdat$values),
    age2 = cdat$Age,
    gear = cdat$gear,
    ngear = length(unique(cdat$gear)),
    n2 = nrow(cdat),
    s_vb = 10,
    s_lw = 10,
    s_cc = 10,
    s_cc_2 = 10,
    ages = seq(1, max(gdat$Age), 1),
    n_ages = length(seq(1, max(gdat$Age), 1))
    )

# Fit the model with stan
# This will take a while the first time you
# compile the model. Right now it is commented
# out and results from previous run loaded.
# Uncomment to reproduce analysis.
  # channel <- stan(file = 'models/population_analysis_wcc.stan',
  #             data = vb_data,
  #             chains = 3,
  #             iter = 10000,
  #             warmup = 9000,
  #             control = list(
  #               adapt_delta = 0.80,
  #               max_treedepth = 10
  #             ), save_warmup = FALSE
  # )
  # # Save the fitted model object to file for use below
  # # and in plotting scripts
  # save(channel, file="results/channel_posts.rda")

# Load results from previously fit model
  load("results/channel_posts.rda")  
  
# Print model summary
#  print(channel, digits=3)
  
# Extract posteriors
  pars <- rstan::extract(channel)
  list2env(pars, envir = .GlobalEnv)

# . Derived quantities ----
# Age at length from posteriors
  t330 = -1/K * log(1 - 330.2/Linf) + t0 # Age for min Length
  t711 = -1/K * log(1 - 711.6/Linf) + t0 # Age at 29"
  t889 = -1/K * log(1 - 889.4/Linf) + t0 # Age at 36"
  
  
# W-infinity from Linf and L-W regression
  Winf <- 10 ^ (a + log10(Linf) * b)
  
# Derive omega from VBGF params (Galucci and Quinn 1979)
  w = K*Linf  
  
# Natural Mortality
  # Life-history (indirect) estimators
  # CW_NM and QD_NM are funky, but for different reasons?
    C_NM <- 4.31 * ((t0 - (log(0.05) / K))^-1.01) 
    D_NM <- 0.8598 * (Winf^-0.0302) * (K^0.5280)# k is from length von b, ok?
    J_NM <- 1.50 * K
    L_NM <- 3.00 * (Winf^-0.288)
    CW_NM <- (1 / (17 - 3)) * log(exp(K * 17) - exp(K * t0))/(exp(K * 3) - exp(K * t0)) #set t0 same age range as mort
    QD_NM <- -log(.01)/17
    NM <- (C_NM + D_NM + J_NM + L_NM + CW_NM + QD_NM) / 6 # element-wise average of all estimates
    
  # Catch curve estimate  
    CC_M <- -beta_2
    
  # Put them all together
    morts <- data.frame(
      do.call(cbind, mget(ls()[grep(pattern="_NM", x=ls())])),
      NM,
      CC_M
      )
    
# . Data write ----  
# Build data.frame of mean and 95% HDI for output data from
# Stan models
  out_data = data.frame(cbind(
    Winf, Linf, K, t0, w, # VBGF params
    t330, t711, t889,     # Age at size (mm)
    a, b,                 # L-W regression params
    NM                    # L-H mortality
    ))
    
# Add names
  colnames(out_data) = c(
    "Winf", "Linf", "K",
    "t0", "w", "t330",
    "t711", "t889", "a",
    "b", "Life_hist_M"
    )  
  
# Summarize  
  out_data = tidyr::gather(out_data) %>% group_by(key) %>% 
    summarise(mean = mean(value, na.rm=T), 
              L_ci = quantile(value, probs = 0.025, na.rm=T),
              U_ci = quantile(value, probs = 0.975, na.rm=T)
    )
  
# Add catch-curve mortality (Z) 
  cc_Z <- c('Z', mean(-beta_2), quantile(-beta_2, probs=c(0.025,0.975)))
  names(cc_Z) <- names(out_data)
  out_data <- rbind(out_data, cc_Z)
  
# Add estimated fishing mortality (F)
  FMM <- apply(-beta_2, 2, "-", NM)
  FMM <- c('F', mean(FMM), quantile(FMM, probs=c(0.025,0.975)))
  names(FMM) <- names(out_data)  
  out_data <- rbind(out_data, FMM)
  
# Write the table out to a file  
  write.table(out_data, file="results/channel_summary_stan.csv",
                    quote=FALSE, sep=",", row.names = FALSE)  
  
  
# YPR From age-structured models (Dippold et al. 2016) ----
# . Age and size from posts ----
# Age classes in population
  age <- seq(1, max(gdat$Age))
  
  # Length at age from VBGF
  length_at_age <- matrix(data = NA, nrow = length(K), ncol=length(age))
  for(i in 1:length(K)){
    for(t in 1:length(age)){
      length_at_age[i, t] <- Linf[i] * (1 - exp(-K[i]*(age[t] - t0[i])))    
    }
  }
  
  # Weight at age from L-W regression
  log_length <- apply(length_at_age, c(1, 2), log10)
  Wt <- matrix(data = NA, nrow = length(K), ncol=length(age))
  for(i in 1:length(K)){
    for(t in 1:length(age)){
      Wt[i,t] <- 10^(a[i] + b[i] * log_length[i, t])
    }
  }  
  
# . Parallel settings ----
# Number of cores
  ncpus <- 7
  
# Initialize snowfall
  #cl <- makeCluster(ncpus, type="SOCK")
  sfInit(parallel = TRUE, cpus=ncpus, type="SOCK")
  
# . sim function ----
  sim <- function(x){  
  
  # Fishing mortality 
    fM <- sample(seq(0, 1.00, 0.01), 1)
    
  # Ages corresponding to length limits  
    Tc <- mean(t330)
    
    if(is.na(mean(t711, na.rm=TRUE))){
      Tmax <- max(age)
    } else {
      Tmax <- min(c(max(age), mean(t711, na.rm = TRUE)))
    }
    
  # Mortality rates for each mgmt scenario
    scenarios <- c('none', 'min', 'slot')
    FM_none <- rep(fM, max(age))
    
    FM_min <- c(rep(0, mean(floor(Tc)-1)),
                mean(Tc-floor(Tc))*fM,
                rep(fM, mean((max(age)-(floor(Tc)))))
                )
    
    FM_slot <- c(rep(0, mean(floor(Tc)-1)),
                 mean(Tc-floor(Tc))*fM,
                 rep(fM, mean(Tmax-floor(Tc))),
                 mean(Tmax-floor(Tmax)) * fM,
                 rep(0, max(0,mean(length(age)-ceiling(Tmax))))
                 )
    
    scen <- sample(scenarios, 1)
    FM <- unlist(mget(paste0('FM_', scen)))
  
  # Population projection with demographic error from pop analysis
  # Vector pre-allocation
    N <- matrix(data = NA, nrow=length(K), ncol=max(gdat$Age))  
    N[ , 1] <- 1000
  
    res <- yprC(fm=FM, nm=morts$NM, N=N, Wt=Wt)$ypr  
    
    outs <- 
      data.frame(
        FM = rep(fM, length(K)),
        scen = rep(scen, length(K)),
        ypr = rowSums(res)
    )
    
  }

# . Load libraries on workers -----
  sfLibrary(myCpp)
  sfLibrary(Rcpp)
  sfExport(
      'cdat',
      'gdat',
      'K',
      'age',
      't330',
      't711',
      't889',
      'NM',
      'Wt',
      'morts'
  )
  
# . Distribute to workers -----
# Number of simulations to run
  niterations <- 1000  
  
# . Run simulation in parallel ----
# Store system time
  start <- Sys.time()
  
# Run simulation
  result <- sfLapply(1:niterations, sim) 
  
# Calculate elapsed time  
  Sys.time()-start
  
# Stop snowfalll
  sfStop()    

# . Summary statistics----
# Process output into a data.table  
  resypr <- data.table::rbindlist(result)
  
# Summarize ypr by F and scenario 
  ypr <- plyr::ddply(resypr,
             c('FM', 'scen'),
             summarize,
             avg=mean(ypr),
             upper=quantile(ypr, probs=0.975),
             lower=quantile(ypr, probs=0.025)
             )
  
# Calculate mean and 95% CRI for F_msy by scenario 
  # Filter data by scenario
  none_lim <- filter(ypr, scen=='none')
  min_lim <- filter(ypr, scen=='min')
  slot_lim <- filter(ypr, scen=='slot')
  
  # Calculate mean and 95% CRI for F_msy by scenario 
  none_fmsy <- none_lim$FM[apply(apply(none_lim[ , 3:5], 2, diff), 2, function(x) first(which(x<0)))] 
  min_fmsy <- min_lim$FM[apply(apply(min_lim[ , 3:5], 2, diff), 2, function(x) first(which(x<0)))] 
  slot_fmsy <- slot_lim$FM[apply(apply(slot_lim[ , 3:5], 2, diff), 2, function(x) first(which(x<0)))] 

  # Calculate mean and 95% CRI for MSY by scenario 
  none_msy <- none_lim$avg[apply(apply(none_lim[ , 3:5], 2, diff), 2, function(x) first(which(x<0)))] 
  min_msy <- min_lim$lower[apply(apply(min_lim[ , 3:5], 2, diff), 2, function(x) first(which(x<0)))] 
  slot_msy <- slot_lim$upper[apply(apply(slot_lim[ , 3:5], 2, diff), 2, function(x) first(which(x<0)))] 

  
# . Data write ----  
# Write ypr summary for plotting
  save(ypr, file='results/channel_ypr.rda')
  
# Write F_msy estimates to a file
  channel_F <- do.call(rbind, mget(ls()[grep(pattern='_fmsy', x=ls())]))
  row.names(channel_F) <- gsub(
    pattern="_fmsy", x=row.names(channel_F), replacement = ""
    )
  colnames(channel_F) <- c('mean', 'lwr', 'upr')
  write.csv(channel_F, file='results/channel_fmsy.csv', quote=FALSE)
  
# Write MSY estimates to a file
  channel_msy <- do.call(rbind, mget(ls()[grep(pattern='_msy', x=ls())]))
  row.names(channel_msy) <- gsub(
    pattern="_msy", x=row.names(channel_msy), replacement = ""
    )
  colnames(channel_msy) <- c('mean', 'lwr', 'upr')
  write.csv(channel_msy, file='results/channel_msy.csv', quote=FALSE)
  
  
# . Preliminary YPR plots ----
# Plot the mean and 95% CRI for posterior ypr 
# predictions.
  # No size limit
  plot(ypr$FM[ypr$scen=='none'],
       ypr$avg[ypr$scen=='none'],
       type='l', col='red', lwd=2,
       ylab = "Yield per recruit (g)",
       xlab = "F",
       ylim = c(0, 1e6),
       main='channel catfish'
       )
  polygon(x=c(ypr$FM[ypr$scen=='none'],
              rev(ypr$FM[ypr$scen=='none'])),
          y=c(ypr$lower[ypr$scen=='none'],
              rev(ypr$upper[ypr$scen=='none'])),
          col = rgb(1, 0, 0, 0.5),
          border=NA)
  # Minimum size limit
  lines(ypr$FM[ypr$scen=='min'],
       ypr$avg[ypr$scen=='min'],
       col='blue', lwd=2
  )
  polygon(x=c(ypr$FM[ypr$scen=='min'],
              rev(ypr$FM[ypr$scen=='min'])),
          y=c(ypr$lower[ypr$scen=='min'],
              rev(ypr$upper[ypr$scen=='min'])),
          col = rgb(0, 0, 1, 0.45),
          border=NA)  
  # Slot limit (736 mm)
  lines(ypr$FM[ypr$scen=='slot'],
       ypr$avg[ypr$scen=='slot'],
       col = 'gray90'
  )
  polygon(x=c(ypr$FM[ypr$scen=='slot'],
              rev(ypr$FM[ypr$scen=='slot'])),
          y=c(ypr$lower[ypr$scen=='slot'],
              rev(ypr$upper[ypr$scen=='slot'])),
          col = rgb(.7, .7, 0.7, 0.5),
          border=NA)
  legend(x = 0.70, y=8.5e5, fill = c('red', 'blue', 'gray'),
         legend = c('none', 'min', 'slot'), border = NA,
         box.col = 'white')
  

  
  

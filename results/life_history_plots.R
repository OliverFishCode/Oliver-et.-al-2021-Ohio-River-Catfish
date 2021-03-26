# Plotting code for growth and mort predictions
# We make a 3 x 2 plot of L-H predictions
# with species in rows, growth on left,
# and catch-curve on the right.

# Package load ----
  library(rstan)
  library(plyr)

# Data load ----
# Run the data script if not done yet
  ### Throws warnings about conversion of weight and 
  ### length to numeric from factor, but okay otherwise
  if(!exists('flathead_mort')) source('data/data_prep.R')  
  

# Get a list of files containing ypr preds
  files <-  dir('results')[grep(pattern = '_posts.rda', x=dir('results'))]

# Graphs with polygons -----  
# Open file connection
  png(
    filename = paste0('results/life-history.png'),
    height = 1350,
    width = 1350,
    pointsize = 32
  )
# Graphical parameters
  par(mar=c(4,4,1,1))    
# Set graphical parameters
  par(mfrow=c(3,2), oma=c(5,1,2,0), mar=c(1,5,1,1))
  
# Loop through species and mortality estimators      
  for(i in 1:length(files)){
  
  # Get ith species
    species = c('blue', 'channel', 'flathead')[i]
    
  # Load data file for ith species
    load(file=paste0('results/', files[i]))
    mod <- mget(species)
    pars <- rstan::extract(mod[[1]])
    list2env(pars, envir = .GlobalEnv)
    
  # Load raw data
    rdata <- mget(ls()[grep(pattern=paste0(species,"_combined"),ls())])[[1]]
    cdata <- mget(ls()[grep(pattern=paste0(species,"_mort"),ls())])[[1]]
    cdata <- data.frame(
      cdata[ ,1],
      stack(cdata[ ,2:ncol(cdata)])
    )  
    cdat <- cdata[!is.na(cdata$values), ]
    cdat$gear <- as.numeric(as.factor(cdat$ind))
    gdata <- mget(ls()[grep(pattern=paste0(species,"_growth"),ls())])[[1]]
    
  # Make length-at-age predictions  
    # Von bertalanffy
      ages = seq(0, max(rdata$Age), 0.1)
      preds = matrix(NA, nrow=length(K), ncol=length(ages))
      for(j in 1:length(K)){
        for(t in 1:length(ages)){
          preds[j, t] <- Linf[j] * (1 - exp(-K[j]*(ages[t]-t0[j])))
        }
      }
      
      vbmeans <- apply(preds, 2, mean)
      vblwr <- apply(preds, 2, quantile, 0.025) 
      vbupr <- apply(preds, 2, quantile, 0.975)
    
    # Plot
      plot(ages, vbmeans, 
           type='l', col='black', lwd=1,
           ylab = "",
           xlab = "",
           ylim = c(0, 1100),
           xlim = c(0, 25),
           yaxt = 'n',
           xaxt='n'
           )
      polygon(x=c(ages,rev(ages)),
              y=c(vblwr, rev(vbupr)),
              col = rgb(0.4, 0.4, 0.4, 0.35),
              border=NA)
      points(rdata$Age, rdata$Length, pch=21, col=NA,
             bg= rgb(0.2,0.2,0.2,0.15))
      
      # X-axis ticks and labels
      if(i < 3){
        axis(side = 1, at=seq(0,25,5), tick = T, labels = F)
      }
      if(i == 3){
        axis(side = 1, at=seq(0,25,5), labels = seq(0,25,5))
      }
      
      # Y-axis ticks and labels
      axis(side = 2, at=seq(0,2000,200),
       labels = sprintf("%.0f", seq(0,2000,200)),
       las = 2)
      
      # Y-axis title
      if(i == 2){
        mtext('Total length (mm)', side=2, line=4)  
      }
    
    # Catch-curve mortality
      # Predictions
        ages = seq(0, max(cdat$Age), 1)
        preds = array(NA, c(nrow(beta_2), length(ages), ncol(beta_2)))
        for(p in 1:nrow(beta_2)){
          for(t in 1:length(ages)){
            for(j in 1:ncol(beta_2)){
              
              preds[p, t, j] <- alpha_2[p, j] + ages[t]*beta_2[p, j]
              
            }  
          }
        }
      
      ccmeans <- (apply(preds, c(2), mean))
      cclwr <- (apply(preds, c(2), quantile, 0.025)) 
      ccupr <- (apply(preds, c(2), quantile, 0.975))        
        
      # Plot raw catch by age and gear
      # e15, e60, Hoopnet, Trotline (square, circle, triangle, plus sign)
        plot(x=cdat$Age, y = log(cdat$values), 
             pch = c(0,1,2,3)[as.numeric(as.factor(cdat$ind))],
             col= 'black',
             ylab = "",
             xlim = c(0, 25),
             ylim = c(-0.01, 5.5),
             xlab = "",
             yaxt = 'n',
             xaxt='n'
             ) 
      polygon(x=c(ages+min(cdata$Age), rev(ages+min(cdata$Age))),
              y=c(ccupr, rev(cclwr)),
              col = rgb(0.4, 0.4, 0.4, 0.15),
              border=NA)
      lines(x=ages+min(cdata$Age), y=ccmeans, lty=1, lwd=1, col='black')
      # lines(x=ages+min(cdata$Age), y=cclwr)
      # lines(x=ages+min(cdata$Age), y=ccupr)      
      
    # X-axis ticks and labels
      if(i < 3){
        axis(side = 1, at=seq(0,25,5), tick = T, labels = F)
      }
      if(i==3){
        axis(side = 1, at=seq(0,25,5), labels = seq(0,25,5))
      }

      
    # Y-axis ticks and labels
      axis(side = 2, at=seq(0,5.5,.5),
       labels = sprintf("%.1f", seq(0,5.5,.5)),
       las = 2)
      if(i==2){
        mtext(side=2, text=expression(paste('log'[italic(e)], 'N')), line=3.5)
      }

    # Species labels
      text(
        x=25,
        y=4.5, 
        paste(c("Blue", "Channel", "Flathead")[i], " Catfish"),
        pos=2
        )

  }
  # Add x-axis label  
    mtext("Age (years)", side=1, adj=-.5, line=3)  
  
dev.off()  

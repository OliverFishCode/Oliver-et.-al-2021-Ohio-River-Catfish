# Plotting code for ypr predictions
# We make a 3 x 2 plot of YPR predictions
# with species in rows and mortality estimator
# in columns.

# Get a list of files containing ypr preds
  files <-  dir('results')[grep(pattern = '_ypr.rda', x=dir('results'))]
  sums <- dir('results')[grep(pattern = '_summary_stan.csv', x=dir('results'))]
  
# Graphs with polygons -----  
# Open file connection
  png(
    filename = paste0('results/ypr.png'),
    height = 1350,
    width = 1350,
    pointsize = 32
  )
# Graphical parameters
  par(mar=c(4,4,1,1))    
# Set graphical parameters
  par(mfrow=c(3,2), oma=c(4,5,2,0), mar=c(1,1,1,1))
  
# Loop through species and mortality estimators      
  for(i in 1:length(files)){
  
  # Load data file for ith scenario
    load(file=paste0('results/', files[i]))  
    ypr[,3:5] <- ypr[,3:5]/1000
    sumi <- read.csv(paste0('results/', sums[i]))
    
  # No size limit
    plot(ypr$FM[ypr$scen=='none'],
         ypr$avg[ypr$scen=='none'],
         type='l', col='black', lwd=1,
         ylab = "Yield per recruit (g)",
         xlab = "F",
         ylim = c(0, .5e3),
         yaxt = 'n',
         xaxt='n'
         )
    polygon(x=c(ypr$FM[ypr$scen=='none'],
                rev(ypr$FM[ypr$scen=='none'])),
            y=c(ypr$lower[ypr$scen=='none'],
                rev(ypr$upper[ypr$scen=='none'])),
            col = rgb(0.4, 0.4, 0.4, 0.35),
            border=NA)
    lines(ypr$FM[ypr$scen=='none'],
     ypr$avg[ypr$scen=='none'],
     col='black', lwd=1, lty=1
    )    
    lines(x = c(0, sumi[13, 4]), y = c(10,10), lwd=2)
    if(sumi[13, 'mean'] >= 0){
      points(x = sumi[13, 'mean'], y = 10, pch=19, col='black')
    } else {
      points(x = 0, y = 10, pch=19, col='black')
    }
    # X-axis ticks and labels
    if(i %in% c(1,2)){
      axis(side = 1, at=seq(0,1,.1), tick = T, labels = F)
    } else {
      axis(side = 1, at=seq(0,1,.1),
           labels = sprintf("%.1f", seq(0,1,.1)))
    }
    # Y-axis ticks and labels
    axis(side = 2, at=seq(0,1000,100),
     labels = sprintf("%.0f", seq(0,1000,100)),
     las = 2)
    # Y-axis title
    if(i == 2){
      mtext('Yield per recruit (kg)', side=2, line=4)  
    }
    if(i == 1){
      mtext('No limit', side=3, line=1) 
    }    
      
  # Minimum size limit
    plot(ypr$FM[ypr$scen=='min'],
         ypr$avg[ypr$scen=='min'],
         type='l', col='black', lwd=1,
         ylab = "Yield per recruit (g)",
         xlab = "F",
         ylim = c(0, .5e3),
         yaxt = 'n',
         xaxt='n'
         )
    polygon(x=c(ypr$FM[ypr$scen=='min'],
                rev(ypr$FM[ypr$scen=='min'])),
            y=c(ypr$lower[ypr$scen=='min'],
                rev(ypr$upper[ypr$scen=='min'])),
            col = rgb(0.4, 0.4, 0.4, 0.4),
            border=NA)
    lines(x = c(0, sumi[13, 4]), y = c(10,10), lwd=2)
    if(sumi[13, 'mean'] >= 0){
      points(x = sumi[13, 'mean'], y = 10, pch=19, col='black')
    } else {
      points(x = 0, y = 10, pch=19, col='black')
    }

  # Slot limit (species-specific)
    if(length(grep(pattern='blue', x=files[i])) > 0|
       length(grep(pattern='flathead', x=files[i])) > 0
       ){
      polygon(x=c(ypr$FM[ypr$scen=='slot'],
                  rev(ypr$FM[ypr$scen=='slot'])),
              y=c(ypr$lower[ypr$scen=='slot'],
                  rev(ypr$upper[ypr$scen=='slot'])),
              col = rgb(0.25, 0.25, 0.25, 0.6),
              border=NA)
      # Means    
      lines(ypr$FM[ypr$scen=='min'],
           ypr$avg[ypr$scen=='min'],
           col='black', lwd=1, lty=1
      )  
      lines(ypr$FM[ypr$scen=='slot'],
           ypr$avg[ypr$scen=='slot'],
           col = 'black', lty=2
      )
      lines(x = c(0, sumi[13, 4]), y = c(10,10), lwd=2)
      if(sumi[13, 'mean'] >= 0){
        points(x = sumi[13, 'mean'], y = 10, pch=19, col='black')
      } else {
        points(x = 0, y = 10, pch=19, col='black')
      }
    }

    # X-axis ticks and labels
      if(i %in% c(1,2)){
        axis(side = 1, at=seq(0,1,.1), tick = T, labels = F)
      } else {
        axis(side = 1, at=seq(0,1,.1),
             labels = sprintf("%.1f", seq(0,1,.1)))
      }
      axis(side = 2, at=seq(0,1000,100), tick = T, labels = F)
      
    # Column title
    if(i == 1){
      mtext('Size limit', side=3, line=1)  
    }
    # Species labels
      text(
        x=1,
        y=450, 
        paste(c("Blue", "Channel", "Flathead")[i], " Catfish"),
        pos=2
        )

  }
  # Add x-axis label  
    mtext(expression(italic(F)), side=1, adj=-.10, line=3)  
  
dev.off()  

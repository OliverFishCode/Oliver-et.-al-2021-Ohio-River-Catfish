# YPR predictions from Beverton and Holt (1957) ----
Tc <- t330               # Age at first capture - Not sure what this should be
Tr <- rep(8, length(Tc)) # Age at recruitment - Same here
M <- morts$e15#apply(morts[,c(10, 11)], 1, mean)
FM <- seq(0, 2, 0.05)
Z <- vector(mode='numeric', length=length(M))
ypr <- matrix(data=NA, nrow=length(K), ncol=length(FM))
for(i in 1:length(K)){
  for(t in 1:length(FM)){
    Z[i] <- FM[t] + M[i]
    ypr[i,t] <- FM[t] * exp( -M[i] * (Tc[i] - Tr[i]) ) *
      Winf[i] * ( 1/Z[i] - b[i]*exp(-K[i] * (Tc[i] - t0[i]))/ (Z[i]+K[i]) +
                    b[i]*exp(-K[i] * (Tc[i] - t0[i]))^2 / (Z[i] +2*K[i]) -
                    exp(-K[i] * (Tc[i] - t0[i]))^3 / (Z[i] + b[i]*K[i]))   
  }
}
# Plot the result
par(mar=c(5,4,1,1))
plot(FM, apply(ypr, 2, mean), type='l', ylim=c(0, 2500), 
     xlab = expression(italic('F')), ylab='YPR', yaxt='n')
axis(2, las=2)
polygon(c(FM, rev(FM)),
        c(apply(ypr, 2, quantile, probs=0.025),
          rev(apply(ypr, 2, quantile, probs=0.975))),
        col=rgb(.8, .8, .8, 0.5),
        border=NA
)
box()

### ISSUES IDENTIFIED SO FAR ----
#
#   A COUPLE OF LENGTHS AT AGE IN CHAN AND FLAT DATA
#   THAT DON'T MAKE BIOLOGICAL SENSE
#
#   A COUPLE OF POINTS IN THE L-W REGRESSION THAT ARE
#   WAY TOO HEAVY FOR THE CORRESPONDING LENGTHS AND
#   DO NOT MAKE ANY BIOLOGICAL SENSE



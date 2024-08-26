# An example of how to generate a periodic changepoint structure time series and how to run the SNcirc algorithm

# change to an appropriate directory
# load library
library(changepoint.periodic)

# Create a periodic changepoint time series (with stable periodicity)
period.len=24 # number of observations within one fixed period. 24 => hourly data for daily patterns
per.cpts = list(c(7,22))
params = list()
params[[1]] = list(mean=c(2,10),var=c(1,3)) # if dist=Bernoulli, no need to include the var argument
n = period.len*14 # total number of observations. 14 periods of data, i.e. 2 weeks
set.seed(283)
dat = pgcptSim(period=period.len,dist="Normal",glob.cpts=c(0,n),per.cpts=per.cpts,params=params,reps=1)

# Plot the data
# plot(dat,pch=16) # linear time plot
# withinPeriod_plot(dat,period.len=period.len,test.stats="Normal meanvar",circData=TRUE) # circular time plot

# Run SNcirc
result = sncirc(dat, period.len,dist="Normal mean",max.cpts=5,minseglen=1,pen.val=3*log(length(dat)), useClass = FALSE)
result$sncirc.results$op.cpt.loc # optimal changepoint locations
## Plot the data and the identified periodic changepoints
withinPeriod_plot(dat,period.len=period.len,test.stats="Normal meanvar",circData=TRUE) # circular time plot
abline(v=result$op.cpt.loc,col='blue') # plot the optimal changepoint locations onto the circular data


# Extract optimal number of changepoints
m=result$sncirc.results$op.ncpts # no. of changepoints
result$sncirc.results$cps[,m]
withinPeriod_plot(dat,period.len=period.len,test.stats="Normal meanvar",circData=TRUE)
abline(v=result$sncirc.results$cps[,m],col='blue')






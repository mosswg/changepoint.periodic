

test_that("sncirc.norm mean", {
  period.len=24 # number of observations within one fixed period. 24 => hourly data for daily patterns
  params = list()
  params[[1]] = list(mean=c(2,10),var=c(1,3)) # if dist=Bernoulli, no need to include the var argument
  n = period.len*14 # total number of observations. 14 periods of data, i.e. 2 weeks

  for (i in 1:23) {
    for (j in (i+1):24) {
      set.seed(283040) # occasionally sncirc will find an additional change point with random seed
      changepoints = c(i, j)
      per.cpts = list(changepoints)
      dat = pgcptSim(period=period.len,dist="Normal",glob.cpts=c(0,n),per.cpts=per.cpts,params=params,reps=1)

      # Run SNcirc
      data.input = data_circ(dat,period.len,n)
      out = sncirc.norm(data=data.input,period.len=period.len,max.cpts=5,minseglen=1,pen.val=3*log(length(dat)),dist="Normal mean")
      expect_equal(unname(out$op.cpt.loc), changepoints)
    }
  }
})


test_that("sncirc.norm meanvar", {
  period.len=24 # number of observations within one fixed period. 24 => hourly data for daily patterns
  params = list()
  params[[1]] = list(mean=c(2,10),var=c(1,3)) # if dist=Bernoulli, no need to include the var argument
  n = period.len*14 # total number of observations. 14 periods of data, i.e. 2 weeks

  for (i in 1:23) {
    for (j in (i+1):24) {
      set.seed(283040) # occasionally sncirc will find an additional change point with random seed
      changepoints = c(i, j)
      per.cpts = list(changepoints)
      dat = pgcptSim(period=period.len,dist="Normal",glob.cpts=c(0,n),per.cpts=per.cpts,params=params,reps=1)

      # Run SNcirc
      data.input = data_circ(dat,period.len,n)
      out = sncirc.norm(data=data.input,period.len=period.len,max.cpts=5,minseglen=1,pen.val=3*log(length(dat)),dist="Normal meanvar")
      expect_equal(unname(out$op.cpt.loc), changepoints)
    }
  }
})


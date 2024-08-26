

test_that("param.ests Normal", {
  period.len=24 # number of observations within one fixed period. 24 => hourly data for daily patterns
  changepoints = c(7, 22)
  per.cpts = list(changepoints)
  params = list()
  params[[1]] = list(mean=c(2,10),var=c(1,3)) # if dist=Bernoulli, no need to include the var argument
  n = period.len*14 # total number of observations. 14 periods of data, i.e. 2 weeks
  set.seed(283040) # results change wildly but stay around 2 and 10 with set seed
  data = pgcptSim(period=period.len,dist="Normal",glob.cpts=c(0,n),per.cpts=per.cpts,params=params,reps=1)

  dat = data_circ(data,period.len,n)
  data2plot = dat[,1:2]
  params = param.ests(data2plot, period.len, params2est="Mean", changepoints, circdata=FALSE)
  expect_equal(params, c(2.114654, 9.851897), tolerance = .0001) # this is probably a bad way of doing this test but since i dont know what exactly what this function calculates its hard to write a better one
})


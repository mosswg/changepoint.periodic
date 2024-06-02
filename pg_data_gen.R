## Code to simulate periodic-global changepoint series (if no global changes, then a time series with stable periodicity is generated)

pgcptSim=function(period=12,dist="Normal",glob.cpts=NA,per.cpts=list(),
                  params=list(),reps=1){
  # Function to simulate data with a periodic-global structure for different distributions
  
  if(!is.numeric(period)){stop("period must be numeric")}
  if(period<3){stop("period must be atleast 3 to contain periodic changepoints")}
  
  dists = c("Normal","Bernoulli")
  if(all(dist!=dists)){stop(paste0("Only distributions which are currently implemented are: ",paste0(dists,collapse=',')))}
  
  #if(is.na(glob.cpts)){stop("Global changepoints must contain the length of the data to be simulated as the final element")}
  if(any(!is.numeric(glob.cpts))){stop("Global changepoints must be numeric")}
  if(glob.cpts[length(glob.cpts)]<(2*period)){stop("Length of the data (final glob.cpts element) must be atleast 2*period")}
  
  if(length(per.cpts)==1 && (length(glob.cpts)-1)!=1){ 
    # if the per.cpts are the same for all segments then repeat in a list for ease later on
    tmp=per.cpts
    per.cpts=list()
    for(i in 1:(length(glob.cpts)-1)){
      per.cpts[[i]]=tmp
    }
  }
  if(any(!is.numeric(unlist(per.cpts)))){stop("Periodic changepoints must be numeric")}
  if(max(unlist(per.cpts))>period){stop("All per.cpts vectors must be in 1:period")}
  if(min(unlist(per.cpts)<1)){stop("All per.cpts vectors must be in 1:period")}
  if(any(unlist(per.cpts)>period)){stop("Periodic changepoints can only take indices in 1:period")}
  
  if(!is.numeric(reps)){stop("reps must be numeric")}
  if(reps<1){stop("number of repetitions (reps) must be positive")}
  
  if(length(params)!=(length(glob.cpts)-1)){stop("There must be the same number of segments (length(params)) as there are glob.cpts.")}
  
  if(dist=="Normal"){
    out=mmcptSim.norm(period=period,glob.cpts=glob.cpts,per.cpts=per.cpts,
                  params=params,reps=reps)
  } else if(dist=="Bernoulli"){
    out=mmcptSim.bern(period=period,glob.cpts=glob.cpts,per.cpts=per.cpts,
                      params=params,reps=reps)
  }
  return(out)
}


mmcptSim.norm=function(period,glob.cpts,per.cpts,params,reps){
  # call the circcptSim.norm function to generate data within a global-segment
  # given glob.cpts, generate periodic level data for that segment
  # glob.cpts=vector. glob.cpts starts with 1 and ends with n=length of data.
  # per.cpts=list of vectors.
  # period = length of single period
  # params=list of list of parameters. For now, we are just simulating normal data. So param[[i]][1]=mean, param[[i]][2]=var, i=ith segment.
  # reps = number of replication of time series for simulations, i.e. number of data sets for a number of simulations
  
  oneMMdataset <- function(period=period,glob.cpts=glob.cpts,per.cpts=per.cpts,params=params){
  # out = NULL
  # for(i in 1:reps){
    dat = NULL
    tlen <<- length(glob.cpts)
    for(i in 1:(length(glob.cpts)-1)){
      l=glob.cpts[i+1]-glob.cpts[i] # length of this global-segment
      # need to shift the periodic-cpts and the params here in relation to the global-cpt position
      # do this here and then put the new cpts and params into the circcptSim function 
      if(glob.cpts[i]%%period!=0){
        params[[i]]$mean = c(params[[i]]$mean[which(per.cpts[[i]]>=(glob.cpts[i]%%period))],params[[i]]$mean[which(per.cpts[[i]]<(glob.cpts[i]%%period))]) # shift the position of the means
        params[[i]]$var = c(params[[i]]$var[which(per.cpts[[i]]>=(glob.cpts[i]%%period))],params[[i]]$var[which(per.cpts[[i]]<(glob.cpts[i]%%period))]) # shift the position of the vars
        per.cpts.i = sort((per.cpts[[i]]-glob.cpts[i])%%period) # shift them
      } else{
        per.cpts.i = per.cpts[[i]]
      }
      
      seg.dat = circcptSim.norm(period,length=l,per.cpts=per.cpts.i,params=params[[i]],reps=1)
      dat = c(dat,seg.dat)
    }
    # out <- cbind(out, dat)
    return(dat)
  # }
  }
  
  out=replicate(reps,oneMMdataset(period,glob.cpts,per.cpts,params))
  return(out)
}

circcptSim.norm=function(period,length,per.cpts,params,reps){
  # Function is to be intended to be used within mmcptSim only, not to be called directly
  # as no error checking!
  
  # period=period length
  # length=length of the data for one series, n.
  # per.cpts = location of periodic cpts
  # params=list of parameters. For now, we are just simulating normal data. So param[1]=mean, param[2]=var
  # reps = number of replication of time series for simulations, i.e. number of data sets for a number of simulations
  
  if(!any(names(params)=="mean")){stop("Normal distribution must have a params list called 'mean'")}
  if(!any(names(params)=="var")){stop("Normal distribution must have a params list called 'var'")}
  
  onedataset=function(period,length,per.cpts,params,global.start.index){
    means.per=params$mean
    if(length(means.per)!=length(per.cpts)){means.per=rep(means.per,length(per.cpts))}
    # if the mean doesn't change (i.e. 1 value for all segments) then repeat it to ease calculations later
    vars.per=params$var
    if(length(vars.per)!=length(per.cpts)){vars.per=rep(vars.per,length(per.cpts))}
    # if the var doesn't change (i.e. 1 value for all segments) then repeat it to ease calculations later

    if(per.cpts[length(per.cpts)]!=period){
      per.cpts=c(per.cpts,period)
      means.per=c(means.per,means.per[1])
      vars.per=c(vars.per,vars.per[1])
    }
    # add the period to per.cpts but put the same mean as the start of the period
    # this gets the circular effect but using linear time (as data are independent)
    per.cpts=c(0,per.cpts)
    means.per=rep(means.per,times=diff(per.cpts)) # get the means per timepoint within a period
    means.per=rep(means.per,length.out=length) # repeat the means across periods to the desired length
    vars.per=rep(vars.per,times=diff(per.cpts)) # get the vars per timepoint within a period
    vars.per=rep(vars.per,length.out=length) # repeat the vars across periods to the desired length
    
    data=rnorm(length,mean=means.per,sd=sqrt(vars.per))
    return(data)
  }
  
  out=replicate(reps,onedataset(period,length,per.cpts,params))
  return(out)
}




mmcptSim.bern=function(period,glob.cpts,per.cpts,params,reps){
  # call the circcptSim.norm function to generate data within a global-segment
  # given glob.cpts, generate Periodic level data for that segment
  # glob.cpts=vector. glob.cpts starts with 1 and ends with n=length of data.
  # per.cpts=list of vectors.
  # period = length of single period
  # params=list of list of parameters. For now, we are just simulating normal data. So param[[i]][1]=mean, param[[i]][2]=var, i=ith segment.
  # reps = number of replication of time series for simulations, i.e. number of data sets for a number of simulations
  
  oneMMdataset <- function(period=period,glob.cpts=glob.cpts,per.cpts=per.cpts,params=params){
    # out = NULL
    # for(i in 1:reps){
    dat = NULL
    tlen <<- length(glob.cpts)
    for(i in 1:(length(glob.cpts)-1)){
      l=glob.cpts[i+1]-glob.cpts[i] # length of this global-segment
      # need to shift the Periodic-cpts and the params here in relation to the global-cpt position
      # do this here and then put the new cpts and params into the circcptSim function 
      if(glob.cpts[i]%%period!=0){
        params[[i]]$mean = c(params[[i]]$mean[which(per.cpts[[i]]>=(glob.cpts[i]%%period))],params[[i]]$mean[which(per.cpts[[i]]<(glob.cpts[i]%%period))]) # shift the position of the means
        per.cpts.i = sort((per.cpts[[i]]-glob.cpts[i])%%period) # shift them
      } else{
        per.cpts.i = per.cpts[[i]]
      }
      
      seg.dat = circcptSim.bern(period,length=l,per.cpts=per.cpts.i,params=params[[i]],reps=1)
      dat = c(dat,seg.dat)
    }
    # out <- cbind(out, dat)
    return(dat)
    # }
  }
  
  out=replicate(reps,oneMMdataset(period,glob.cpts,per.cpts,params))
  return(out)
}


circcptSim.bern=function(period,length,per.cpts,params,reps){
  # Function is to be intended to be used within mmcptSim only, not to be called directly
  # as no error checking!
  
  # period=period length
  # length=length of the data for one series, n.
  # per.cpts = location of periodic cpts
  # params=list of parameters. For now, we are just simulating normal data. So param[1]=mean, param[2]=var
  # reps = number of replication of time series for simulations, i.e. number of data sets for a number of simulations
  
  if(!any(names(params)=="mean")){stop("Bernoulli distribution must have a params list called 'mean'")}
  
  onedataset=function(period,length,per.cpts,params){
    means.per=params$mean
    if(length(means.per)!=length(per.cpts)){means.per=rep(means.per,length(per.cpts))}
    # if the mean doesn't change (i.e. 1 value for all segments) then repeat it to ease calculations later
    
    if(per.cpts[length(per.cpts)]!=period){
      per.cpts=c(per.cpts,period)
      means.per=c(means.per,means.per[1])
    }
    # add the period to per.cpts but put the same mean as the start of the period
    # this gets the circular effect but using linear time (as data are independent)
    per.cpts=c(0,per.cpts)
    means.per=rep(means.per,times=diff(per.cpts)) # get the means per timepoint within a period
    means.per=rep(means.per,length.out=length) # repeat the means across periods to the desired length
    
    
    data=rbinom(length,size=1,prob=means.per)
    return(data)
  }
  
  out=replicate(reps,onedataset(period,length,per.cpts,params))
  return(out)
}


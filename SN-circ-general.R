# SNcirc function for user
sncirc <- function(data, period.len=96, dist="Normal meanvar",max.cpts=5,minseglen=1,pen.val=3*log(length(data)),cost="Likelihood",circData=TRUE){
  # Function to run SNcirc on periodic data for different distributions
  
  # assume data comes in two forms, matrix or vector.
  # vector: just the data. matrix/dataframe: two columns, first is data
  if(is.null(dim(data))){
    n = length(data)
  } else if(is.numeric(dim(data))){
    n = nrow(data)
  } else{stop("Data must be a vector or matrix/dataframe with 1st column as data")}
  
  if(!is.numeric(period.len)){stop("Period length must be numeric")}
  if(period.len<3){stop("Period length must be at least 3 to contain changepoints")}
  
  if(!is.numeric(max.cpts)){stop("Maximum number of changepoints must be numeric")}
  if(max.cpts>floor(period.len/minseglen)){stop(paste('M is larger than the maximum number of changepoints: ',floor(period.len/minseglen)))}
  
  if(!is.numeric(minseglen)){stop("Minimum segment length must be numeric")}
  if(minseglen<1){stop("Minimum segment length needs to be at least 1.")}
  
  # circData: does the data need a circular time index column? If TRUE, then add a circular time index column using the code below
  if(isTRUE(circData)){
    data.input = data_circ(data,period.len,n)
  } else{ data.input = data}
  
  dist.list = c("Normal mean","Normal meanvar","Bernoulli") #list of distributions that can be done
  # penalty.list = c("SIC"); pen.vals.list = c(log(n))
  # pen.val = pen.vals.list[which(penalty.list==penalty)]
  
  if(dist=="Normal mean"){
    out = sncirc.norm(data=data.input,period.len=period.len,max.cpts=max.cpts,minseglen=minseglen,pen.val=pen.val,dist="Normal mean")
  } else if(dist=="Normal meanvar"){
    out = sncirc.norm(data=data.input,period.len=period.len,max.cpts=max.cpts,minseglen=minseglen,pen.val=pen.val,dist="Normal meanvar")
  } else if(dist=="Bernoulli"){
    out = sncirc.bern(data=data.input,period.len=period.len,max.cpts=max.cpts,minseglen=minseglen,pen.val=pen.val,dist="Bernoulli")
  }
  
  return(list(sncirc.results=out,period.len=period.len,dist=dist,max.cpts=max.cpts,minseglen=minseglen,pen.val=pen.val,cost=cost))
}

# function to add within-period time index column to data. 
# Output: 1st column=within-period time index. 2nd column=data.
data_circ <- function(data,period.len=96,n.obs){
  if(missing(n.obs)){
    if(is.null(dim(data))){
      n.obs = length(data)
    } else if(is.numeric(dim(data))){
      n.obs = nrow(data)
    } else{stop("Data must be a vector or matrix/dataframe with 1st column as data")}
  }
  time.index = rep(1:period.len, floor(n.obs/period.len))
  if(n.obs%%period.len!=0){time.index=c(time.index, 1:(n.obs%%period.len))}
  out = cbind(time.index, data)
  return(out)
}


# SNcirc distribution functions -------------------------------------------

# function to run SNcirc with normal distribution
sncirc.norm <- function(data,period.len=96,max.cpts=5,minseglen=1,pen.val=0,dist="Normal meanvar",all.seg){
  if(missing(all.seg)){all.seg = sncirc.allseg.norm(data,period.len,dist)}
  N = period.len
  M = max.cpts
  
  like.M <- array(0,c(M,N,N)) #in the format like.M[m,j,k]
  #initialise like.M for m=1
  for(k in 1:N){
    jvals = 1:N
    for(j in jvals[-k]){ #09/06/22: I think k->k should be included in the likelihoods now
      like.M[1,j,k] = all.seg[k,j]
    }
  }
  
  cp = array(NA,c(M,N,N)) # cp[m,j,k], last cpt location prior to j (mth changepoint), given m cpts and starting at k
  for(k in 1:N){
    for(m in 2:M){
      for(j in (k+m*minseglen):(k-1+N)){
        like=NULL
        v=(k+(m-1)*minseglen):(j-minseglen)%%N # potential positions for the (m-1)th cpt prior to j
        v[which(v==0)] = N
        v2=(v+1)%%N # the start of the next segment
        v2[which(v2==0)] = N
        j.circ = j%%N
        j.circ[j.circ==0] = N
        like=like.M[m-1,v,k]+all.seg[v2,j.circ]
        
        like.M[m,j.circ,k] = max(like,na.rm=TRUE) 
        add.cpt.loc = which(like==max(like,na.rm=TRUE))[1]+(k+(m-1)*minseglen-1) # cpt location prior to j
        if(add.cpt.loc>N){add.cpt.loc = add.cpt.loc%%N}
        cp[m,j.circ,k]=add.cpt.loc
      }
    }
  }
  
  #collapse the current like.M dimension from 3 to 2
  like.M.coll = matrix(ncol=N,nrow=M,0)
  op.k = NULL #optimal starting positions for each m
  #future work: speed up the below by using the apply function
  for(m in 1:M){
    likes.k = NULL
    for(k in 1:N){
      wrapAround = (k-1+N)%%N
      if(wrapAround==0){wrapAround = N}
      likes.k = c(likes.k, like.M[m,wrapAround,k])
    }
    k.opt = which(likes.k==max(likes.k,na.rm=TRUE))[1] #find the optimal k for this m
    like.M.coll[m,] = like.M[m,,k.opt]
    op.k = c(op.k, k.opt)
  }
  
  cps.M=matrix(NA,ncol=m,nrow=m) #goes back and finds the optimal cpt locations for each m
  #k=start of a segment. So "first" changepoint is at k-1.
  f.cpts = (op.k[1]-1)%%N
  if(f.cpts==0){f.cpts=N}
  for(m in 2:M){
    f = (op.k[m] - 1)%%N
    if(f==0){f=N}
    f.cpts = c(f.cpts,f) #save the first cpts for each m in a list.
    cps.M[m,1]=cp[m,f,op.k[m]]
    for(i in 1:(m-1)){
      cps.M[m,(i+1)]=cp[(m-i),cps.M[m,i],op.k[m]]
    }
  }
  cps.M = cbind(f.cpts, cps.M[,-ncol(cps.M)])
  
  
  op.ncps=NULL #the optimal number cpts for a given penalty
  h=c(0,2:M) #1cpt=0cpts so we should have a 0 penalty at first.
  
  lv = NULL #likelihood vector
  for(i in 1:M){
    k = op.k[i]
    k.end = (k-1+N)%%N
    if(k.end == 0){k.end = N}
    lv = c(lv, like.M.coll[i,k.end])
  }
  criterion=-2*lv+h*pen.val #likelihood plus the penalty term for each m
  op.ncps<-h[which(criterion==min(criterion,na.rm=T))[1]]
  
  if(op.ncps==(M)){warning('The number of segments identified is M, it is advised to increase M to make sure changepoints have not been missed.')}
  if(op.ncps==0){cpts=N}
  else{
    cpts = c(sort(cps.M[op.ncps,][cps.M[op.ncps,]>0])) #cpt locations for the optimal no. of cpts
  }
  
  if(op.ncps==0){op.like=criterion[1]}else{op.like=criterion[op.ncps]}
  return(list(cps=apply(cps.M,1,sort,na.last=TRUE),op.ncpts=op.ncps,op.cpt.loc=cpts,op.like=op.like, like.M=like.M, period.len=period.len, pen.val=pen.val,
              cp=cp, like.M.coll=like.M.coll, op.k=op.k, op.likes.m=lv, op.likes.m.pen=criterion, max.cpts=max.cpts, all.seg=all.seg))
}

# function to run SNcirc with bernoulli distribution
sncirc.bern <- function(data,period.len=96,max.cpts=5,minseglen=1,pen.val=0,dist="Bernoulli",all.seg){
  if(missing(all.seg)){all.seg = sncirc.allseg.bern(data,period.len)}
  N = period.len
  M = max.cpts
  
  like.M <- array(0,c(M,N,N)) #in the format like.M[m,j,k]
  #initialise like.M for m=1
  for(k in 1:N){
    jvals = 1:N
    for(j in jvals[-k]){ #09/06/22: I think k->k should be included in the likelihoods now
      like.M[1,j,k] = all.seg[k,j]
    }
  }
  
  cp = array(NA,c(M,N,N)) # cp[m,j,k], last cpt location prior to j (mth changepoint), given m cpts and starting at k
  for(k in 1:N){
    for(m in 2:M){
      for(j in (k+m*minseglen):(k-1+N)){
        like=NULL
        v=(k+(m-1)*minseglen):(j-minseglen)%%N # potential positions for the (m-1)th cpt prior to j
        v[which(v==0)] = N
        v2=(v+1)%%N # the start of the next segment
        v2[which(v2==0)] = N
        j.circ = j%%N
        j.circ[j.circ==0] = N
        like=like.M[m-1,v,k]+all.seg[v2,j.circ]
        
        like.M[m,j.circ,k] = max(like,na.rm=TRUE)
        add.cpt.loc = which(like==max(like,na.rm=TRUE))[1]+(k+(m-1)*minseglen-1) # cpt location prior to j
        if(add.cpt.loc>N){add.cpt.loc = add.cpt.loc%%N}
        cp[m,j.circ,k]=add.cpt.loc
      }
    }
  }
  
  #collapse the current like.M dimension from 3 to 2
  like.M.coll = matrix(ncol=N,nrow=M,0)
  op.k = NULL #optimal starting positions for each m
  #future work: speed up the below by using the apply function
  for(m in 1:M){
    likes.k = NULL
    for(k in 1:N){
      wrapAround = (k-1+N)%%N
      if(wrapAround==0){wrapAround = N}
      likes.k = c(likes.k, like.M[m,wrapAround,k])
    }
    k.opt = which(likes.k==max(likes.k,na.rm=TRUE))[1] #find the optimal k for this m
    like.M.coll[m,] = like.M[m,,k.opt]
    op.k = c(op.k, k.opt)
  }
  
  cps.M=matrix(NA,ncol=m,nrow=m) #goes back and finds the optimal cpt locations for each m
  #k=start of a segment. So "first" changepoint is at k-1.
  f.cpts = (op.k[1]-1)%%N
  if(f.cpts==0){f.cpts=N}
  for(m in 2:M){
    f = (op.k[m] - 1)%%N
    if(f==0){f=N}
    f.cpts = c(f.cpts,f) #save the first cpts for each m in a list.
    cps.M[m,1]=cp[m,f,op.k[m]]
    for(i in 1:(m-1)){
      cps.M[m,(i+1)]=cp[(m-i),cps.M[m,i],op.k[m]]
    }
  }
  cps.M = cbind(f.cpts, cps.M[,-ncol(cps.M)])
  
  
  op.ncps=NULL #the optimal number cpts for a given penalty
  h=c(0,2:M) #1cpt=0cpts so we should have a 0 penalty at first.
  
  lv = NULL #likelihood vector
  for(i in 1:M){
    k = op.k[i]
    k.end = (k-1+N)%%N
    if(k.end == 0){k.end = N}
    lv = c(lv, like.M.coll[i,k.end])
  }
  criterion=-2*lv+h*pen.val #likelihood plus the penalty term for each m
  op.ncps<-h[which(criterion==min(criterion,na.rm=T))[1]]
  
  if(op.ncps==(M)){warning('The number of segments identified is M, it is advised to increase M to make sure changepoints have not been missed.')}
  if(op.ncps==0){cpts=N}
  else{
    cpts = c(sort(cps.M[op.ncps,][cps.M[op.ncps,]>0])) #cpt locations for the optimal no. of cpts
  }
  
  if(op.ncps==0){op.like=criterion[1]}else{op.like=criterion[op.ncps]}
  return(list(cps=apply(cps.M,1,sort,na.last=TRUE),op.ncpts=op.ncps,op.cpt.loc=cpts,op.like=op.like, like.M=like.M, period.len=period.len, pen.val=pen.val,
              cp=cp, like.M.coll=like.M.coll, op.k=op.k, op.likes.m=lv, op.likes.m.pen=criterion, max.cpts=max.cpts, all.seg=all.seg))
}



# all.seg calculation functions -------------------------------------------

# function to calculate all.seg with normal distribution
sncirc.allseg.norm <- function(data,period.len=96,dist="Normal meanvar"){
  #do the subset stuff before all.seg to reduce computational time
  N=period.len
  
  subset.mat = matrix(0,ncol=3,nrow=period.len) #col1=len, col2=sumx, col3=ssq
  for(t in 1:N){ #go through for each within-period time index
    dat = subset(data,data[,1]==t)
    len = nrow(dat); sumx = sum(dat[,2]); ssq = sum(dat[,2]^2)
    subset.mat[t,] = c(len,sumx,ssq)
  }
  
  all.seg=matrix(0,ncol=N,nrow=N) #ncol=N instead of N+(N-1)
  for(i in 1:N){
    ssq=0; sumx=0; jend = N+i-1; len = 0
    for(j in i:jend){
      j.circ=j%%N; j.circ[j.circ==0]=N;
      len = len + subset.mat[j.circ,1]
      sumx = sumx + subset.mat[j.circ,2]
      ssq = ssq + subset.mat[j.circ,3]
      if(dist=="Normal meanvar"){
        sigsq=ssq-(sumx^2)/len
        if(sigsq<=0){sigsq=0.0000000000000001}
        
        all.seg[i,j.circ]=-0.5*len*(log(2*pi)+log(sigsq)-log(len)+1)
      } else if(dist=="Normal mean"){
        all.seg[i,j.circ]=-0.5*(ssq-(sumx^2)/len)
      }
    }
  }
  
  return(all.seg)
}

# function to calculate all.seg with bernoulli distribution
sncirc.allseg.bern <- function(data,period.len=96){
  #do the subset stuff before all.seg to reduce computational time
  N=period.len
  
  subset.mat = matrix(0,ncol=2,nrow=N) #col1=len, col2=sumx, col3=ssq
  for(t in 1:N){ #go through for each within-period time index
    dat = subset(data,data[,1]==t)
    len = nrow(dat); sumx = sum(dat[,2]);
    subset.mat[t,] = c(len,sumx)
  }
  
  all.seg=matrix(0,ncol=N,nrow=N) #ncol=N instead of N+(N-1)
  for(i in 1:N){
    ssq=0; sumx=0; jend = N+i-1; len = 0
    for(j in i:jend){
      j.circ=j%%N; j.circ[j.circ==0]=N;
      len = len + subset.mat[j.circ,1]
      sumx = sumx + subset.mat[j.circ,2]
      if(sumx==len){sl = 1; reverse.sl = 1e-323
      } else if(sumx==0){sl = 1e-323; reverse.sl = 1
      } else{sl = sumx/len; reverse.sl = 1-sl}
      all.seg[i,j.circ]= sumx*log(sl)+(len-sumx)*log(reverse.sl)
    }
  }
  
  return(all.seg)
}



# Plot function -----------------------------------------------------------

withinPeriod_plot <- function(data, period.len=96,test.stats='Normal meanvar',title='',param.lines=FALSE,cpts=c(period.len),params2plot="Mean",circData=TRUE,param.est.col='red',defaultAxes=TRUE){
  # assume data comes in two forms, matrix or vector.
  # vector: just the data. matrix/dataframe: two columns, first is data
  if(is.null(dim(data))){
    n = length(data)
  } else if(is.numeric(dim(data))){
    n = nrow(data)
  } else{stop("Data must be a vector or matrix/dataframe with 1st column as data")}
  
  
  if(isTRUE(circData)){
    # TRUE = need to put the data into circular time
    dat = data_circ(data,period.len,n)
  } else{dat = data}
  #assumes first column is within period time index, second column is the data
  # test.stats options: norm, bern
  if(test.stats=='Normal meanvar' || test.stats=='Normal mean'){
    data2plot = dat[,1:2]
    ylim1 = min(data2plot[,2])
    ylim2 = max(data2plot[,2])
    ylabel = "Data"
  } else if(test.stats=='Bernoulli'){
    #need to plot the proportions instead of the actual data
    props = NULL
    for(i in 1:period.len){
      a = sum(dat[which(dat[,1]==i),2])/length(dat[which(dat[,1]==i),2])
      props = c(props, a)
      ylabel = "Proportions"
    }
    data2plot = cbind(1:period.len,props)
    ylim1 = 0
    ylim2 = 1
  }
  plot(data2plot, xlim=c(0,period.len+1), ylim=c(ylim1,ylim2),xlab='Fixed period',ylab=ylabel,main=title,pch=16,axes=defaultAxes)
  if(param.lines==TRUE){
    params = param.ests(data2plot, period.len, params2est=params2plot, cpts, circdata=FALSE)
    cpts = c(0,cpts)
    for(i in 1:(length(cpts)-1)){
      segments(x0=cpts[i],x1=cpts[i+1],y0=params[i],y1=params[i],col=param.est.col,lwd=2)
    }
    if(cpts[length(cpts)]!=period.len){segments(x0=cpts[length(cpts)],x1=period.len,y0=params[1],y1=params[1],col=param.est.col,lwd=2)}
  }
}



# calculates the estimated parameters in the estimated segments of the (circular time) data
param.ests <- function(data, period.len=96, params2est='Mean', cpts=c(period.len), circdata=TRUE){
  # circdata => does the data need to be put into circular time? If so, use data_circ
  # if not, assume data[,1] is the within period time index, data[,2] is the data
  if(circdata==TRUE){
    data = data_circ(data,period.len)
  }
  cpts = c(0,cpts)
  params.out = "No code to calculate this parameter estimate yet" # place here in case the params2est arg isn't an option below
  if(params2est=='Mean'){
    params.out = NULL
    if(cpts[length(cpts)]!=period.len){ # account for the wrapAround segment. Take this out and work this out accordingly.
      iterateVector = 2:(length(cpts)-1)
      wrapAround=TRUE # indicates if the segment does wrap around and if so, I'll need to work out the params for this separately
    } else{
      iterateVector = 1:(length(cpts)-1);
      wrapAround=FALSE
    }
    for(i in iterateVector){
      times = (cpts[i]+1):cpts[i+1] # all the time points in each segment
      subset.dat = NULL
      for(j in times){
        subset.dat = c(subset.dat,data[which(data[,1]==j),2]) # all the data within the segment
      }
      params.out[i] = mean(subset.dat)
    }
    if(isTRUE(wrapAround)){
      times = c(1:cpts[1],(cpts[length(cpts)]+1):period.len)
      subset.dat = NULL
      for(j in times){
        subset.dat = c(subset.dat,data[which(data[,1]==j),2]) # all the data within the segment
      }
      params.out[1] = mean(subset.dat)
    }
  }
  
  return(params.out)
}







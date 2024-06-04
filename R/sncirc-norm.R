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

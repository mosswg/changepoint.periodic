#' Run SNcirc on the provided data
#'
#' @param data CHANGEME
#' @param period.len CHANGEME
#' @param dist CHANGEME
#' @param max.cpts CHANGEME
#' @param minseglen CHANGEME
#' @param pen.val CHANGEME
#' @param cost CHANGEME
#' @param circData CHANGEME
#' @returns sncirc results
#'
#' @export
sncirc <- function(data, period.len=96, dist="Normal meanvar",max.cpts=5,minseglen=1,pen.val=3*log(length(data)),cost="Likelihood",circData=TRUE,useClass=TRUE){
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
    out = sncirc_c(data=data.input,period.len=period.len,max.cpts=max.cpts,minseglen=minseglen,pen.val=pen.val,dist_type=0)
  } else if(dist=="Normal meanvar"){
    out = sncirc_c(data=data.input,period.len=period.len,max.cpts=max.cpts,minseglen=minseglen,pen.val=pen.val,dist_type=1)
  } else if(dist=="Bernoulli"){
    out = sncirc_c(data=data.input,period.len=period.len,dist_type=2,max.cpts=max.cpts,minseglen=minseglen,pen.val=pen.val)
  }

  if (useClass) {
  }
  else {
    return(list(sncirc.results=out,period.len=period.len,dist=dist,max.cpts=max.cpts,minseglen=minseglen,pen.val=pen.val,cost=cost))
  }
}


# dist_type: 0 = normal mean, 1 = normal meanvar, 2 = bernoulli
#' @useDynLib changepoint.periodic, .registration = TRUE
sncirc_c <- function(data,period.len=96, dist_type=0, max.cpts=5, minseglen=1, pen.val=0) {
  # function that uses the PELT method to calculate changes in mean where the segments in the data are assumed to be Normal
  n = length(data)
  N = period.len
  M = max.cpts
  if(n<2){stop('Data must have atleast 2 observations to fit a changepoint model.')}

  error = 0
  op.cps = 0
  min_criterion = 0
  op.like = 0

  data.len = nrow(data)

  cptsout=rep(0,n) # sets up null vector for changepoint answer

  like.m <- array(0,c(M,N,N)) #in the format like.M[m,j,k]
  cp = array(-1,c(M,N,N)) # cp[m,j,k], last cpt location prior to j (mth changepoint), given m cpts and starting at k
  like.m.coll = matrix(ncol=N,nrow=M,0)
  cps.m=matrix(-1,ncol=M,nrow=M) #goes back and finds the optimal cpt locations for each m
  lv = array(-1, M)
  criterion = array(-1, M)
  f.cpts = array(-1, M)
  all.seg=matrix(-1,ncol=N,nrow=N) #ncol=N instead of N+(N-1)
  op.k = array(-1, M)

  answer=list()
  # answer[[6]]=1
  # on.exit(.C("FreePELT",answer[[6]]))

  storage.mode(data) = 'double'
  storage.mode(dist_type) = 'integer'
  storage.mode(data.len) = 'integer'

  storage.mode(period.len) = 'integer'
  storage.mode(max.cpts) = 'integer'
  storage.mode(minseglen) = 'integer'
  storage.mode(pen.val) = 'double'
  storage.mode(error) = 'integer'
  storage.mode(cp) = 'integer'
  storage.mode(op.cps) = 'integer'
  storage.mode(op.k) = 'integer'
  storage.mode(f.cpts) = 'integer'
  storage.mode(criterion) = 'double'
  storage.mode(lv) = 'integer'

  storage.mode(cptsout)='integer'
  storage.mode(cps.m) = 'integer'
  storage.mode(like.m) = 'double'
  storage.mode(like.m.coll) = 'double'
  storage.mode(op.like) = 'integer'

  answer = .C('sncirc', dist_type, data, data.len, period.len, max.cpts, minseglen, pen.val, cptsout, error, cp, cps.m, op.cps, op.k, f.cpts, criterion, like.m, like.m.coll, lv, all.seg, op.like)

  if(answer[[9]]>0){
    stop("C code error:",answer[[9]],call.=F)
  }

  answer[[10]][answer[[10]] == -1] = NA
  answer[[11]][answer[[11]]==-1] = NA

  answer[[11]] = cbind(answer[[14]], answer[[11]][,-ncol(answer[[11]])])

  if(answer[[12]]==(M)){warning('The number of segments identified is M, it is advised to increase M to make sure changepoints have not been missed.')}
  if(answer[[12]]==0){cpts=N}
  else{
    cpts = c(sort(answer[[11]][answer[[12]],][answer[[11]][answer[[12]],]>0])) #cpt locations for the optimal no. of cpts
  }

  return(list(cps=apply(answer[[11]],1,sort,na.last=TRUE),op.ncpts=answer[[12]],op.cpt.loc=cpts,op.like=answer[[20]], like.M=answer[[16]], period.len=period.len, pen.val=pen.val,
              cp=answer[[10]], like.M.coll=answer[[17]], op.k=answer[[13]], op.likes.m=answer[[18]], op.likes.m.pen=answer[[15]], max.cpts=max.cpts, all.seg=answer[[19]]))


  # return(list(cps=t(apply(answer[[8]],1,sort,na.last=TRUE)),cpts=sort(answer[[6]][answer[[6]]>0]),op.cpts=answer[[9]],pen=answer[[5]],like=answer[[10]],like.Q=-2*(answer[[11]])[,n]))
}



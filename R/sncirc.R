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


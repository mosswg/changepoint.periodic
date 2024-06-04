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

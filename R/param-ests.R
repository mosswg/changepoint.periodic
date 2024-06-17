# export

#' Calculate the estimated parameters in the estimated segments of the (circular time) data
#'
#' @param data CHANGEME
#' @param period.len CHANGEME
#' @param params2est CHANGEME
#' @param cpts CHANGEME
#' @param circdata CHANGEME
#' @returns estimated parameters
#'
#' @export
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


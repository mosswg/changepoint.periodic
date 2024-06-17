# rename and export

#' Plot period data
#'
#' @param data CHANGEME
#' @param period.ln CHANGEME
#' @param test.stats CHANGEME
#' @param title CHANGEME
#' @param param.lines CHANGEME
#' @param cpts CHANGEME
#' @param params2plot CHANGEME
#' @param circData CHANGEME
#' @param param.est.col CHANGEME
#' @param defaultAxes CHANGEME
#'
#' @export
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


#----------Functions----------#

getpeaks <- function(trace) {
  indexes <- which(diff(sign(diff(trace, na.pad = FALSE)), na.pad = FALSE) < 0) + 1              
  cbind(indexes, trace[indexes])
}

peakvalues <- function(x, pstart, pstop) {
  
  region <- x[x[,1] > pstart & x[,1] < pstop, ,drop=FALSE]
  
  if (length(region[,1]) == 0)
    return(c(0, NA))
  else
    return(c(max(region[,2], na.rm=TRUE), region[which.max(region[,2]),1]))
}

maxvalues <- function(x, pstart, pstop) {
  
  region <- x[(pstart+1) : (pstop-1)]
  
  if (length(region) == 0)
    return(NA)
  else
    return(max(region, na.rm=TRUE))
}

baseReturn <- function (index){
  if (index == 1) return("A")
  if (index == 2) return("C")
  if (index == 3) return("G")
  if (index == 4) return("T")
}

char2string<-function (chars){
  return(paste(chars, collapse = ""))
}

plotHist.fun<-function(data,main,y,line5,line3){
  x = seq(1:length(data))
  plot(data,
       type='h',
       #col='blue',
       col=ifelse((x < line5 | x > line3),"red","blue"),
       main=main,
       xlab="Bases",
       ylab=y)
  abline(h = NULL, v = line5, col = "red", lty=2)
  text(line5,max(data), as.character(line5), col = 2, adj = c(-.1, -.1))
  abline(h = NULL, v = line3, col = "red", lty=2)
  text(line3,max(data), as.character(line3), col = 2, adj = c(-.1, -.1))
}

dumbExtractor<- function(list){
  x<-NULL
  for (i in 1:length(list))
    x<-c(x,list[[i]])
  return(x)
}
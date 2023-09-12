mnntestimationresultantvector<-function(data,M=0,R=1){

  data <- as.matrix(data)
  nd <- nrow(data)
  if (R != length(M))
    return("Error: Length of M and number of dimensions are not equal")

  fork <- list(R)

  for (m in R:1) {
    fork[[m]] <- 0:M[m]
  }

  indk <- expand.grid(fork, KEEP.OUT.ATTRS = FALSE)
  indk <- as.matrix(indk)

  resultant<-rep(0+0i,nrow(indk))

  for (k in 1:nd){
    aux<-indk%*%(data[k,])
    resultant<-resultant+exp(-1i*aux)
  }

  resultant<-resultant/sqrt(sum(Mod(resultant)^2))

  resultant<-((1/(2*pi))^(R/2))*resultant
  res<-list()
  res[["cestimates"]]<-cbind.data.frame(indk,resultant)
  res
}

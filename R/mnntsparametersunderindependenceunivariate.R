mnntsparametersunderindependenceunivariate<-function(data,R,Mvector,cparlist){

  if (R != length(Mvector))
    return("Error: Length of M and number of dimensions are not equal")
  sec <- list(R)
  for (k in 1:R) {
    sec[[k]] <- 0:Mvector[k]
  }
  ind <- expand.grid(sec, KEEP.OUT.ATTRS = FALSE)
  ind <- as.matrix(ind)

  cparfinal<-rep(0+0i,prod(Mvector+1))
  cparaux<-cparlist[[R]]$cestimates[,2]

  for (k in (R-1):1){
    cparfinal<-kronecker(cparaux,cparlist[[k]]$cestimates[,2])
    cparaux<-cparfinal
  }
  cestimatesarray <- data.frame(cbind(ind,(cparfinal)))
  names(cestimatesarray)<-c(1:R,"cestimates")
  cestimatesarray[,1:R]<-as.integer(Re(as.matrix(cestimatesarray[,1:R])))

  loglik <- mnntsloglik(data, cparfinal, Mvector, R)
  AIC <- -2 * loglik + 2 * (2 * sum(Mvector))
  BIC <- -2 * loglik + (2 * sum(Mvector)) * log(nrow(data))

  res <- list(cestimates = cestimatesarray, loglik = loglik,
              AIC = AIC, BIC = BIC)
  return(res)
}

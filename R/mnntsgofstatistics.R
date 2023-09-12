mnntsgofstatistics<-function(data,charfunarray,R=1){

  if (ncol(data) != (ncol(charfunarray)-1))
    return("Error: ncol of data should be equal to ncol of charfunarray minus 1")

  if (ncol(charfunarray) != R+1)
    return("Error: ncol charfunarray should be equal to R+1")

  charfunarraynozero<-charfunarray[1:((nrow(charfunarray)-1)/2),]

  designmatrixaux<-as.matrix(data)%*%t(as.matrix(charfunarraynozero[,1:R]))

  designmatrixaux<-exp(1i*designmatrixaux)

  designmatrix<-cbind(Re(designmatrixaux),Im(designmatrixaux))

  empcharfunction<-apply(designmatrix,2,mean)
  empcharfunctionmatrix<-matrix(rep(empcharfunction,nrow(data)),nrow=nrow(data),byrow=TRUE)
  designmatrixf<-designmatrix-empcharfunctionmatrix

  omega<-matrix(0,nrow=2*nrow(charfunarraynozero),ncol=2*nrow(charfunarraynozero))

  for (k in 1:nrow(data)){
    omega<-omega + (designmatrixf[k,])%*%t(designmatrixf[k,])
  }

  omega<-(1/nrow(data))*omega

  #omegainv<-solve(omega)
  omegainv<-chol2inv(chol(omega))

  charfunmatrixaux<-cbind(matrix(rep(t(Re(charfunarraynozero[,R+1])),nrow(data)),nrow=nrow(data),byrow=TRUE),matrix(rep(t(Im(charfunarraynozero[,R+1])),nrow(data)),nrow=nrow(data),byrow=TRUE))

  designmatrixcharfun<-designmatrix-charfunmatrixaux

  designvector<-apply(designmatrixcharfun,2,mean)

  gofstat<-t(designvector)%*%omegainv%*%designvector
  gofstatnormal<-(nrow(data)*gofstat-2*nrow(charfunarraynozero))/(2*sqrt(nrow(charfunarraynozero)))

  res<-list()
  res$gofstat<-gofstat
  res$gofstatnormal<-gofstatnormal

  return(res)

}

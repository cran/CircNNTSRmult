mnntsgofdesignmatrix<-function(data,charfunarray,R=1){

  if (ncol(data) != (ncol(charfunarray)-1))
    return("Error: ncol of data should be equal to ncol of charfunarray minus 1")

  if (ncol(charfunarray) != R+1)
    return("Error: ncol charfunarray should be equal to R+1")

  charfunarraynozero<-charfunarray[1:((nrow(charfunarray)-1)/2),]

  designmatrixaux<-as.matrix(data)%*%t(as.matrix(charfunarraynozero[,1:R]))

  designmatrixaux<-exp(1i*designmatrixaux)
  designmatrix<-cbind(Re(designmatrixaux),Im(designmatrixaux))

  charfunmatrixaux<-cbind(matrix(rep(t(Re(charfunarraynozero[,R+1])),nrow(data)),nrow=nrow(data),byrow=TRUE),matrix(rep(t(Im(charfunarraynozero[,R+1])),nrow(data)),nrow=nrow(data),byrow=TRUE))

  designmatrix<-designmatrix-charfunmatrixaux

  designmatrix<-cbind(1,designmatrix)

  return(designmatrix)
}

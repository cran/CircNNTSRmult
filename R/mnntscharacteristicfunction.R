mnntscharacteristicfunction<-function(cestimatesarray=as.data.frame(matrix(c(0,1/(2*pi)),nrow=1,ncol=2)),M=0,R=1){

  if (R != length(M))
    return("Error: Length of M and number of dimensions are not equal")

  n <- nrow(cestimatesarray)

  cestimatesarrayres<-cestimatesarray

  cestimatesarrayres[,R+1]<-0+0i

  for (k1 in 1:n){
    for (k2 in 1:n){
      # minus the empirical characteristic function coefficient \underline{t}
      auxarray <- as.integer(cestimatesarray[k1,1:R] - cestimatesarray[k2,1:R])
      auxc<-((2*pi)^R)*(cestimatesarray[k1,R+1])*Conj(cestimatesarray[k2,R+1])
      aux<-c(auxarray,auxc)
      aux[1:R]<-as.integer(Re(aux[1:R]))
      cestimatesarrayres<-rbind(cestimatesarrayres,aux)
      cestimatesarrayres[,1:R]<-as.integer(Re(as.matrix(cestimatesarrayres[,1:R])))
    }
    cestimatesarrayres<-aggregate(cestimatesarrayres[,R+1],by=as.list(as.data.frame(cestimatesarrayres[,1:R])),FUN=sum,simplify=TRUE)
  }
  #cestimatesarrayres[1:R,]<--cestimatesarrayres[1:R,]
  cestimatesarrayres[,R+1]<-Conj(cestimatesarrayres[,R+1])
  names(cestimatesarrayres)<-c(1:R,"charfun")
  return(cestimatesarrayres)
}

mnntsmarginalgeneral<-function(cpars=as.data.frame(matrix(c(0,0,1/(2*pi)),nrow=1,ncol=3)),M=c(0,0),R=2,marginal=1){

  if (R != length(M))
    return("Error: Length of M and number of dimensions are not equal")

  if (nrow(cpars) != prod(M+1))
    return("Error: Number of rows of matrix cparameters not equal to prod(Mvector + 1)")

  cpars[,R+1] <- cpars[,R+1]/sqrt(sum(Mod(cpars[,R+1])^2))

  Rcomplete<-1:R
  Rmarginal<-marginal
  Rmarginalizing<-setdiff(Rcomplete,Rmarginal)

  sec <- list(length(Rmarginal))
  for (k in 1:length(Rmarginal)) {
    sec[[k]] <- 0:M[Rmarginal[k]]
  }
  ind.marginal <- expand.grid(sec, KEEP.OUT.ATTRS = FALSE)
  ind.marginal <- as.matrix(ind.marginal)

  sec <- list(length(Rmarginalizing))
  for (k in 1:length(Rmarginalizing)) {
    sec[[k]] <- 0:M[Rmarginalizing[k]]
  }
  ind.marginalizing <- expand.grid(sec, KEEP.OUT.ATTRS = FALSE)
  ind.marginalizing <- as.matrix(ind.marginalizing)

  Cmatrix<-matrix(0+0i,nrow=prod(M[Rmarginal]+1),ncol=prod(M[Rmarginal]+1))

  for (k in 1:nrow(ind.marginalizing)){
    aux <- cpars[apply(cpars[,Rmarginalizing]==matrix(ind.marginalizing[k,],nrow=nrow(cpars),ncol=length(Rmarginalizing),byrow=TRUE),1,sum)==length(Rmarginalizing),]
    Cmatrix <- Cmatrix + (aux[,R+1])%*%t(Conj(aux[,R+1]))
  }
  aux2<-eigen(Cmatrix,symmetric=TRUE)
  aux2vectors<-aux2$vectors
  aux2values<-aux2$values

  for (k in 1:prod(M[Rmarginal]+1)){
    if (Mod(aux2vectors[1,k]) < Mod(aux2vectors[prod(M[Rmarginal]+1),k])){
      aux2vectors[1:prod(M[Rmarginal]+1),k]<-aux2vectors[prod(M[Rmarginal]+1):1,k]
    }
    aux2vectors[,k]<-aux2vectors[,k]*rep(exp(-1i*Arg(aux2vectors[1,k])),prod(M[Rmarginal]+1))
    if (Re(aux2vectors[1,k])<0){
      aux2vectors[,k]<--aux2vectors[,k]
    }
  }
  aux2vectors<-((1/(2*pi))^(length(Rmarginal)/2))*aux2vectors

  res <- list()
  res[["index"]]<-ind.marginal
  res[["eigenvectors"]] <- aux2vectors
  res[["eigenvalues"]]<-aux2values
  res
}

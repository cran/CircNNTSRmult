mnntsmarginalgeneraldimension<-function(cpars=as.data.frame(matrix(c(0,0,1/(2*pi)),nrow=1,ncol=3)),M=c(0,0),R=2,marginal=1,theta=matrix(0,nrow=1,ncol=1)){

  if (ncol(theta) != length(marginal))
    return("Error: Dimension of the vector of theta angles not equal to dimension of the index to marginalize")

  if (R != length(M))
    return("Error: Length of M and number of dimensions are not equal")

  if (nrow(cpars) != prod(M+1))
    return("Error: Number of rows of matrix cparameters not equal to prod(Mvector + 1)")

  Rmarginal<-sort(marginal)

  sec <- list(length(Rmarginal))
  for (k in 1:length(Rmarginal)) {
    sec[[k]] <- 0:M[Rmarginal[k]]
  }
  ind.marginal <- expand.grid(sec, KEEP.OUT.ATTRS = FALSE)
  ind.marginal <- as.matrix(ind.marginal)

  auxsumvec <- rep(0+0i,prod(M[-Rmarginal]+1))

  for (k in 1:nrow(ind.marginal)){
    auxk <- cpars[apply(cpars[,Rmarginal]==matrix(ind.marginal[k,],nrow=nrow(cpars),ncol=length(Rmarginal),byrow=TRUE),1,sum)==length(Rmarginal),]
    auxsumvec <- auxsumvec + auxk[,R+1]*exp(1i*sum(ind.marginal[k,]*as.vector(theta)))
  }

  value <- sum(auxsumvec*Conj(auxsumvec))
  value <- ((2*pi)^(R-length(marginal)))*Re(value)
  return(value)
}

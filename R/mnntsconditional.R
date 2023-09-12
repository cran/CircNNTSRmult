mnntsconditional<-function(cpars=as.data.frame(matrix(c(0,0,1/(2*pi)),nrow=1,ncol=3)),
                           M=c(0,0),R=2,cond=1,cond.values=0){

  cond<-sort(cond)

  if (R != length(M))
    return("Error: Length of M and number of dimensions are not equal")

  if (sum(cond %in% 1:R) != length(cond)){
    return("Error: Conditional index must be a subset of 1:R")
  }

  if (length(cond.values) != length(cond)){
    return("Error: Number of conditional dimensions not equal to the number of conditional values")
  }

  uncond=setdiff(1:R,cond)

  cparsord <- cpars[,c(uncond,cond,R+1)]
  cparsord <- dfOrder(object=cparsord,columns=1:R)

  trigo.cond<-1
  for (k in 1:length(cond)){
    trigo.cond<-kronecker(trigo.cond,exp(1i*(0:M[cond[k]])*cond.values[k]))
  }

  cpar.cond <- kronecker(diag(prod(M[uncond]+1)),t((trigo.cond)))%*%cpars[,R+1]
  fmarg <- mnntsmarginalgeneraldimension(cpars,M,R,marginal=cond,theta=matrix(cond.values,nrow=1,ncol=length(cond.values)))
  cpar.cond<-cpar.cond/sqrt(fmarg)

  M.uncond<-M[uncond]

  sec <- list(length(uncond))
  for (k in 1:length(uncond)) {
    sec[[k]] <- 0:M.uncond[k]
  }

  ind.uncond <- expand.grid(sec, KEEP.OUT.ATTRS = FALSE)
  ind.uncond <- as.matrix(ind.uncond)

  param <- cbind.data.frame(ind.uncond,cpar.cond)

  param
}

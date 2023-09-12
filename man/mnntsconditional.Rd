\name{mnntsconditional}
\alias{mnntsconditional}
\title{Conditional MNNTS density}
\description{Computes the c parameters of a conditional MNNTS density at a particular value of the conditioning random vector
}
\usage{
mnntsconditional(cpars=as.data.frame(matrix(c(0,0,1/(2*pi)),nrow=1,ncol=3)),
M=c(0,0),R=2,cond=1,cond.values=0)
}
\arguments{
\item{cpars}{Matrix of parameters of an MNNTS density with the first R columns containing the index of the c parameter and the R+1 containing the complex parameter}

\item{M}{Vector of M parameters. A nonnegative integer number for each of the R components of the vector}

\item{R}{Number of dimensions}

\item{cond}{A subset of 1:R indicating the elements of the vector of variables to conditioning on}

\item{cond.values}{A vector of fixed values of the conditional elements of the random vector at which to conditioning on}
}

\value{
\item{param}{A matrix with the index and values of the c parameters for the MNNTS condtional density}
}
\references{
Fernandez-Duran and J. J. and Gregorio-Dominguez and M. M (2023). Multivariate Nonnegative Trigonometric Sums Distributions for High-Dimensional Multivariate Circular Data,
arXiv preprint arXiv:2301.03643v2
}
\author{
Juan Jose Fernandez-Duran and Maria Mercedes Gregorio-Dominguez
}

\examples{
# A univariate conditional from a bivariate joint

set.seed(200)
Mbiv<-c(2,3)
Rbiv<-length(Mbiv)
data(Nest)
data<-Nest*(pi/180)
est<-mnntsmanifoldnewtonestimation(data,Mbiv,Rbiv,100)
est
cpars2cond1<-mnntsconditional(cpars=est$cestimates,M=Mbiv,R=Rbiv,cond=1,cond.values=c(pi/2))
cpars2cond1
nntsplot(cpars2cond1$cpar.cond,M=Mbiv[2])

# A bivariate conditional from a trivariate joint

set.seed(200)
Mtriv<-c(2,3,3)
Rtriv<-length(Mtriv)
data(WindDirectionsTrivariate)
data<-WindDirectionsTrivariate
est<-mnntsmanifoldnewtonestimation(data,Mtriv,Rtriv,100)
est
cpars23cond1<-mnntsconditional(cpars=est$cestimates,M=Mtriv,R=Rtriv,cond=1,cond.values=pi/4)
cpars23cond1
mnntsplot(cpars23cond1,M=Mtriv[c(2,3)])
mnntsplotwithmarginals(cpars23cond1,M=Mtriv[c(2,3)])
}


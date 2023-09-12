\name{mnntsmarginalgeneraldimension}
\alias{mnntsmarginalgeneraldimension}
\title{Marginal Density Function at a Vector of Fixed Values}
\description{Computes the value of the marginal density function at a set of vector of angles}
\usage{
mnntsmarginalgeneraldimension(cpars=as.data.frame(matrix(c(0,0,1/(2*pi)),nrow=1,
ncol=3)),M=c(0,0),R=2,marginal=1,theta=matrix(0,nrow=1,ncol=1))
}

\arguments{
\item{cpars}{Matrix of parameters of an MNNTS density with the first R columns containing the index of the c parameter and the R+1 containing the complex parameter}

\item{M}{Vector of M parameters. A nonnegative integer number for each of the R components of the vector}

\item{R}{Number of dimensions}

\item{marginal}{A subset of 1:R indicating the elements of the vector of variables in the marginal}

\item{theta}{A vector of fixed values of the marginal elements of the random vector at which to obtain the value of the marginal density}
}

\value{A scalar with the value of the marginal density at the specified value of the marginal vector.
}

\references{
Fernandez-Duran and J. J. and Gregorio-Dominguez and M. M (2023). Multivariate Nonnegative Trigonometric Sums Distributions f
}
\author{
Juan Jose Fernandez-Duran and Maria Mercedes Gregorio-Dominguez
}

\examples{
# A univariate marginal from a bivariate joint

set.seed(200)
Mbiv<-c(2,3)
Rbiv<-length(Mbiv)
data(Nest)
data<-Nest
est<-mnntsmanifoldnewtonestimation(data,Mbiv,Rbiv,100)
est
marginal1value<-mnntsmarginalgeneraldimension(cpars=est$cestimates,
M=Mbiv,R=Rbiv,marginal=1,theta=matrix(c(pi/2),nrow=1,ncol=1))
marginal1value

# A bivariate marginal from a trivariate joint

set.seed(200)
Mtriv<-c(2,3,3)
Rtriv<-length(Mtriv)
data(WindDirectionsTrivariate)
data<-WindDirectionsTrivariate
est<-mnntsmanifoldnewtonestimation(data,Mtriv,Rtriv,100)
est
marginal12value<-mnntsmarginalgeneraldimension(cpars=est$cestimates,
M=Mtriv,R=Rtriv,marginal=c(1,2),theta=matrix(c(pi/4,pi/2),nrow=1,ncol=2))
marginal12value
}


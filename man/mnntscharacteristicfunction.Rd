\name{mnntscharacteristicfunction}
\alias{mnntscharacteristicfunction}
\title{Characteristic Function of an MNNTS Density}
\description{
Computes the characteristic function from the c parameters of an MNNTS density}
\usage{
mnntscharacteristicfunction(cestimatesarray=as.data.frame(matrix(c(0,1/(2*pi)),
nrow=1,ncol=2)),M=0,R=1)
}

\arguments{
\item{cestimatesarray}{output from mnntsmanifoldnewtonestimation function}

\item{M}{Vector of M parameters. A nonnegative integer number for each of the R components of the vector}

\item{R}{Number of dimensions}
}

\value{
A data frame (matrix) with the support and values of the characteristic function of the MNNTS density
}

\references{
Fernandez-Duran and J. J. and Gregorio-Dominguez and M. M (2023). Multivariate Nonnegative Trigonometric Sums Distributions for High-Dimensional Multivariate Circular Data,
arXiv preprint arXiv:2301.03643v2
}

\author{
Juan Jose Fernandez-Duran and Maria Mercedes Gregorio-Dominguez
}

\examples{
# A characteristic function from a bivariate MNNTS density

set.seed(200)
Mbiv<-c(2,3)
Rbiv<-length(Mbiv)
data(Nest)
data<-Nest*(pi/180)
est<-mnntsmanifoldnewtonestimation(data,Mbiv,Rbiv,50)
est
charfunbiv23<-mnntscharacteristicfunction(cestimatesarray=est$cestimates,M=Mbiv,R=Rbiv)
charfunbiv23

# A characteristic function from a trivariate MNNTS density

set.seed(200)
Mtriv<-c(2,3,3)
Rtriv<-length(Mtriv)
data(WindDirectionsTrivariate)
data<-WindDirectionsTrivariate
est<-mnntsmanifoldnewtonestimation(data,Mtriv,Rtriv,50)
est
charfuntriv233<-mnntscharacteristicfunction(cestimatesarray=est$cestimates,M=Mtriv,R=Rtriv)
charfuntriv233
}

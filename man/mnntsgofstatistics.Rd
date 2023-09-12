\name{mnntsgofstatistics}
\alias{mnntsgofstatistics}
\title{Statistics of the MNNTS Goodness of Fit Test}
\description{Computes the statistics of the goodness of fit test of an MNNTS density based on the estimated characteristic function
}
\usage{
mnntsgofstatistics(data,charfunarray,R=1)
}

\arguments{
\item{data}{Matrix of angles in radians (with R columns)}

\item{charfunarray}{A data frame (matrix) with the support and values of the characteristic function of the MNNTS density obtained by using the function mnntscharacteristicfunction with vector of parameters M of dimension R}

\item{R}{Number of dimensions}
}
\value{
\item{gofstat}{The value of the goodness of fit statistic}

\item{gofstatnormal}{The value of the normal approximation o fthe goodnes of fit statistic}
 }
\references{Fernandez-Duran and J. J. and Gregorio-Dominguez and M. M (2023). Multivariate Nonnegative Trigonometric Sums Distributions for High-Dimensional Multivariate Circular Data,arXiv preprint arXiv:2301.03643v2

Fan, Y. (1997). Goodness-of-fit tests for a multivariate distribution by the empirical characteristic function. Journal of Multivariate Analysis, 62, 36-63.

}
\author{
Juan Jose Fernandez-Duran and Maria Mercedes Gregorio-Dominguez
}

\examples{# A characteristic function from a bivariate MNNTS density

set.seed(200)
Mbiv<-c(2,3)
Rbiv<-length(Mbiv)
data(Nest)
data<-Nest*(pi/180)
est<-mnntsmanifoldnewtonestimation(data,Mbiv,Rbiv,70)
est
charfunbiv23<-mnntscharacteristicfunction(cestimatesarray=est$cestimates,M=Mbiv,R=Rbiv)
charfunbiv23
gofstats23<-mnntsgofstatistics(data,charfunbiv23,R=2)
gofstats23

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
gofstats233<-mnntsgofstatistics(data,charfuntriv233,R=3)
gofstats233
}


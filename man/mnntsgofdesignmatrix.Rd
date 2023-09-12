\name{mnntsgofdesignmatrix}
\alias{mnntsgofdesignmatrix}
\title{Design Matrix of the MNNTS Goodness of Fit Test}
\description{Computes the design matrix of the auxiliary regression for the goodness of fit test of an MNNTS density based on the estimated characteristic function}
\usage{
mnntsgofdesignmatrix(data,charfunarray,R=1)
}
\arguments{
\item{data}{Matrix of angles in radians (with R columns)}

\item{charfunarray}{A data frame (matrix) with the support and values of the characteristic function of the MNNTS density obtained by using the function mnntscharacteristic function with vector of parameters M of dimension R}

\item{R}{Number of dimensions}
}

\value{A matrix that is the design matrix to run the auxiliary regression for the goodness of fit test
}
\references{Fernandez-Duran and J. J. and Gregorio-Dominguez and M. M (2023). Multivariate Nonnegative Trigonometric Sums Distributions for High-Dimensional Multivariate Circular Data,
arXiv preprint arXiv:2301.03643v2

Fan, Y. (1997). Goodness-of-fit tests for a multivariate distribution by the empirical characteristic function. Journal of Multivariate Analysis, 62, 36-63.

}
\author{Juan Jose Fernandez-Duran and Maria Mercedes Gregorio-Dominguez
}

\examples{

# A characteristic function from a bivariate MNNTS density

set.seed(200)
Mbiv<-c(2,3)
Rbiv<-length(Mbiv)
data(Nest)
data<-Nest*(pi/180)
est<-mnntsmanifoldnewtonestimation(data,Mbiv,Rbiv,70)
est
charfunbiv23<-mnntscharacteristicfunction(cestimatesarray=est$cestimates,M=Mbiv,R=Rbiv)
charfunbiv23
designmatrix23<-mnntsgofdesignmatrix(data,charfunbiv23,R=2)
designmatrix23

# A characteristic function from a trivariate MNNTS density

set.seed(200)
Mtriv<-c(2,3,3)
Rtriv<-length(Mtriv)
data(WindDirectionsTrivariate)
data<-WindDirectionsTrivariate
est<-mnntsmanifoldnewtonestimation(data,Mtriv,Rtriv,40)
est
charfuntriv233<-mnntscharacteristicfunction(cestimatesarray=est$cestimates,M=Mtriv,R=Rtriv)
charfuntriv233
designmatrix233<-mnntsgofdesignmatrix(data,charfuntriv233,R=3)
designmatrix233
}


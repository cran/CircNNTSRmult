\name{mnntestimationresultantvector}
\alias{mnntestimationresultantvector}
\title{c Parameter Vector Estimate}
\description{Computes the c parameter vector estimate based on the mean resultant vector of the vectors of observed trigonometric moments
}
\usage{
mnntestimationresultantvector(data,M=0,R=1)
}

\arguments{

\item{data}{Data frame with the observed vectors of angles. The number of columns must be equal to R}

\item{M}{Vector of M parameters. A nonnegative integer number for each of the R components of the vector}

\item{R}{Number of dimensions}
}

\value{
\item{cestimates}{A matrix with the index and values of the c parameters estimates of the MNNTS density}
}

\references{
Fernandez-Duran and J. J. and Gregorio-Dominguez and M. M (2023). Multivariate Nonnegative Trigonometric Sums Distributions for High-Dimensional Multivariate Circular Data,
arXiv preprint arXiv:2301.03643v2
}
\author{
Juan Jose Fernandez-Duran and Maria Mercedes Gregorio-Dominguez
}

\examples{
# A bivariate dataset

Mbiv<-c(2,3)
Rbiv<-length(Mbiv)
data(Nest)
data<-Nest*(pi/180)
estmeanresultant<-mnntestimationresultantvector(data,M=Mbiv,R=Rbiv)
estmeanresultant

# A trivariate dataset

Mtriv<-c(2,3,3)
Rtriv<-length(Mtriv)
data(WindDirectionsTrivariate)
data<-WindDirectionsTrivariate
estmeanresultant<-mnntestimationresultantvector(data,M=Mtriv,R=Rtriv)
estmeanresultant

}


\name{mnntsparametersunderindependenceunivariate}
\alias{mnntsparametersunderindependenceunivariate}
\title{Marginal Density Function at a Vector of Fixed Values}
\description{
Computes the vector of c parameters of an MNNTS density from the vectors of c parameters of its independent marginals
}
\usage{
mnntsparametersunderindependenceunivariate(data,R,Mvector,cparlist)
}

\arguments{
\item{data}{Matrix of angles in radians (with R columns)}

\item{R}{Number of dimensions}

\item{Mvector}{Vector of M parameters. A nonnegative integer number for each of the R components of the vector}

\item{cparlist}{A list in which each element is a matrix containing the information of the vector of c parameters for each independent marginal component}
}

\value{
\item{cestimates}{Matrix of prod(M+1)*(R+1). The first R columns are the parameter number, and the last column is the c parameter's estimators}

\item{loglik}{Log-likelihood value}

\item{AIC}{Value of Akaike's Information Criterion}

\item{BIC}{Value of Bayesian Information Criterion}
}
\references{
Fernandez-Duran and J. J. and Gregorio-Dominguez and M. M (2023). Multivariate Nonnegative Trigonometric Sums Distributions for High-Dimensional Multivariate Circular Data, arXiv preprint arXiv:2301.03643v2
}
\author{
Juan Jose Fernandez-Duran and Maria Mercedes Gregorio-Dominguez
}

\examples{
# Bivariate MNNTS density from independent marginals

set.seed(200)
Mbiv<-c(2,3)
Rbiv<-length(Mbiv)
data(Nest)
data<-Nest*(pi/180)
est1<-nntsmanifoldnewtonestimation(data[,1],Mbiv[1])
est1
est2<-nntsmanifoldnewtonestimation(data[,2],Mbiv[2])
est2
est12independent<-mnntsparametersunderindependenceunivariate(data,R=Rbiv,
Mvector=Mbiv,cparlist=list(est1,est2))
est12independent

# Trivariate MNNTS density from independent marginals

set.seed(200)
Mtriv<-c(2,3,3)
Rtriv<-length(Mtriv)
data(WindDirectionsTrivariate)
data<-WindDirectionsTrivariate
est1<-nntsmanifoldnewtonestimation(data[,1],Mtriv[1],70)
est1
est2<-nntsmanifoldnewtonestimation(data[,2],Mtriv[2],70)
est2
est3<-nntsmanifoldnewtonestimation(data[,3],Mtriv[3],70)
est3
est123independent<-mnntsparametersunderindependenceunivariate(data,R=Rtriv,
Mvector=Mtriv,cparlist=list(est1,est2,est3))
est123independent
}


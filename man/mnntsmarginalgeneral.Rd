\name{mnntsmarginalgeneral}
\alias{mnntsmarginalgeneral}
\title{Mixing Probabilities of the Elements of the Mixture}
\description{Computes the mixing probabilities (eigenvalues) and parameter c vectors (eigenvectors) of the elements of the mixture defining a general MNNTS marginal of any dimension from an MNNTS density
}
\usage{
mnntsmarginalgeneral(cpars=as.data.frame(matrix(c(0,0,1/(2*pi)),nrow=1,ncol=3)),
M=c(0,0),R=2,marginal=1)
}
\arguments{

\item{cpars}{Matrix of parameters of an MNNTS density with the first R columns containing the index of the c parameter and the R+1 containing the complex parameter}

\item{M}{Vector of M parameters. A nonnegative integer number for each of the R components of the vector}

\item{R}{Number of dimensions}

\item{marginal}{A subset of 1:R indicating the elements of the random vector in the marginal}
}

\value{
\item{index}{Matrix of the index of the marginal MNNTS density}

\item{eigenvectors}{Matrix of the c parameter vectors of each element of the mixture. Each column is a parameter vector}

\item{eigenvalues}{The vector of mixing probabilities}
}

\references{
Fernandez-Duran and J. J. and Gregorio-Dominguez and M. M (2023). Multivariate Nonnegative Trigonometric Sums Distributions for High-Dimensional Multivariate Circular Data,
arXiv preprint arXiv:2301.03643v2
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
data<-Nest*(pi/180)
est<-mnntsmanifoldnewtonestimation(data,Mbiv,Rbiv,100)
est
cparsmarginal1<-mnntsmarginalgeneral(cpars=est$cestimates,M=Mbiv,R=Rbiv,marginal=1)
cparsmarginal1

# A bivariate marginal from a trivariate joint

set.seed(200)
Mtriv<-c(2,3,3)
Rtriv<-length(Mtriv)
data(WindDirectionsTrivariate)
data<-WindDirectionsTrivariate
est<-mnntsmanifoldnewtonestimation(data,Mtriv,Rtriv,100)
est
cparsmarginal12<-mnntsmarginalgeneral(cpars=est$cestimates,M=Mtriv,R=Rtriv,marginal=c(1,2))
cparsmarginal12
}

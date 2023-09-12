\name{CircNNTSRmult-package}
\Rdversion{0.1.0}
\alias{CircNNTSRmult-package}
\alias{CircNNTSRmult}
\docType{package}
\title{Multivariate Circular Data using MNNTS Models}
\description{A collection of utilities for the statistical analysis of multivariate circular data using distributions based on Multivariate Nonnegative Trigonometric Sums (MNNTS). The package includes functions for calculation of densities and distributions, for the estimation of parameters, and more.}

\details{\tabular{ll}{
Package: \tab CircNNTSRmult\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2023-09-09\cr
License: \tab GLP (>=2)\cr
Depends: \tab R (>= 3.5.0), stats, psychTools, CircNNTSR \cr
LazyLoad: \tab yes\cr
NeedsCompilation: \tab no \cr
}

The MNNTS (multivariate NNTS) density on a d-dimensional (\eqn{d>2}) hypertorus by Fernandez-Duran and Gregorio-Dominguez (2014) (see also Fernandez-Duran and Gregorio-Dominguez, 2016) for a vector of angles, \eqn{\underline{\Theta}=(\Theta_1,\Theta_2, \ldots, \Theta_d)^\top}, is defined as
\deqn{
f_{\underline{\Theta}}(\underline{\theta}) = \frac{1}{(2\pi)^d}\underline{c}_{12 \cdots d}^H\underline{e}\underline{e}^H\underline{c}_{12 \cdots d}}

\deqn{ = \frac{1}{(2\pi)^d}\sum_{k_1=0}^{M_1}\sum_{k_2=0}^{M_2} \cdots \sum_{k_d=0}^{M_d}\sum_{m_1=0}^{M_1}\sum_{m_2=0}^{M_2} \cdots \sum_{m_d=0}^{M_d} c_{k_1 k_2 \cdots k_d}\bar{c}_{k_1 k_2 \cdots k_d}e^{\sum_{r=1}^d i(k_r-m_r)\theta_r} \nonumber \\
}
where \eqn{\underline{c}_{12 \cdots d}} is a \eqn{d}-dimensional parameter vector of complex numbers of dimension \eqn{2\prod_{r=1}^{d}(M_r+1) - 1} with subindexes given for all the combinations (Kronecker products) of the \eqn{d} vectors \eqn{\underline{M}_r=(0,1, \ldots, M_r)^\top} for \eqn{r=1,2, \ldots, d} where \eqn{M_r} is the number of terms of the sum in the equation for the \eqn{r}-th component of the vector \eqn{\underline{\Theta}}. 
The vector \eqn{\underline{c}_{12 \cdots d}} must satisfy 
\eqn{\underline{c}_{12 \cdots d}^H\underline{c}_{12 \cdots d}=||\underline{c}_{12 \cdots d}||^2=\sum_{k_1=0}^{M_1}\sum_{k_2=0}^{M_2} \cdots \sum_{k_d=0}^{M_d} ||c_{k_1k_2 \cdots k_d}||^2=1}. 
For identifiabily, \eqn{c_{00 \cdots 0}} is a nonnegative real number. The vector \eqn{\underline{c}_{12 \cdots d}^H} is the Hermitian (conjugate and transpose) of vector \eqn{\underline{c}_{12 \cdots d}}. 
The MNNTS family has many desirable properties, the marginal and conditional densities of any order of an MNNTS density are also MNNTS densities and, independence among the elements of the vector \eqn{\underline{\Theta}} is 
translated into a Kronecker product decomposition in the parameter vector \eqn{\underline{c}_{12 \cdots d}}. For example, in the trivariate case 
\eqn{\underline{\Theta}=(\Theta_1, \Theta_2, \Theta_3)^\top}, if \eqn{\Theta_1}, \eqn{\Theta_2} and \eqn{\Theta_3} are joint independent then, \eqn{\underline{c}_{123}=\underline{c}_{1} \bigotimes \underline{c}_{2} \bigotimes \underline{c}_{3}} where \eqn{\underline{c}_1}, \eqn{\underline{c}_2} and \eqn{\underline{c}_3} are the parameter vectors of the NNTS marginal densities of \eqn{\Theta_1}, \eqn{\Theta_2} and \eqn{\Theta_3}, respectively. 
Similarly, if \eqn{\Theta_1} is groupwise independent of \eqn{(\Theta_2,\Theta_3)^\top} then, \eqn{\underline{c}_{123}=\underline{c}_{1} \bigotimes \underline{c}_{23}} where \eqn{\underline{c}_{23}} is the parameter vector 
of the bivariate MNNTS density of \eqn{(\Theta_2,\Theta_3)^\top}. These results apply to higher dimensions.
}

\author{
Juan Jose Fernandez-Duran and Maria Mercedes Gregorio-Dominguez

Maintainer: Maria Mercedes Gregorio Dominguez <mercedes@itam.mx>

}
\references{
Fernandez-Duran, J. J. and Gregorio-Dominguez M. M. (2014) Modeling angles in proteins and circular genomes using multivariate angular distributions based on nonnegative trigonometric sums. \emph{Statistical Applications in Genetics and Molecular Biology}, 13(1), 1-18.

Fernandez-Duran, J. J. and Gregorio-Dominguez, M. M. (2016). {CircNNTSR}: an {R} package for the statistical analysis of circular, multivariate circular, and spherical data using nonnegative trigonometric sums. \emph{Journal of Statistical Software}, 70, 1--19.

Fernandez-Duran, J. J. and Gregorio-Dominguez, M. M (2023). Multivariate Nonnegative Trigonometric Sums Distributions for High-Dimensional Multivariate Circular Data,
\emph{arXiv preprint arXiv:2301.03643v2}
}

\keyword{ package }


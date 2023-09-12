\name{EURUSDGBPBTCtimesminmax}
\alias{EURUSDGBPBTCtimesminmax}
\docType{data}
\title{Minimun and Maximun daily exchange rates}
\description{Minimun and maximun daily bid and ask exchange rates from Euro-US dollar, GB pound-US dollar, Bitcoin-US dollar from March 22, 2019 to March 22, 2023
}
\usage{data("EURUSDGBPBTCtimesminmax")}
\format{
  A data frame with 1048 observations on the following 14 variables.
  \describe{
    \item{\code{id}}{Observation number}
    \item{\code{day1}}{Date in format day/month/year}
    \item{\code{EURUSDAskMax}}{Daily maximum of ask Euro-US dollar exchange rate}
    \item{\code{EURUSDAskMin}}{Daily minimum of ask Euro-US dollar exchange rate}
    \item{\code{EURUSDBidMax}}{Daily maximum of bid Euro-US dollar exchange rate}
    \item{\code{EURUSDBidMin}}{Daily minimum of bid Euro-US dollar exchange rate}
    \item{\code{GBPUSDAskMax}}{Daily maximum of ask GB pound-US dollar exchange rate}
    \item{\code{GBPUSDAskMin}}{Daily minimum of ask GB pound-US dollar exchange rate}
    \item{\code{GBPUSDBidMax}}{Daily maximum of bid GB pound-US dollar exchange rate}
    \item{\code{GBPUSDBidMin}}{Daily minimum of bid GB pound-US dollar exchange rate}
    \item{\code{BTCUSDAskMax}}{Daily maximum of ask Bitcoin-US dollar exchange rate}
    \item{\code{BTCUSDAskMin}}{Daily minimum of ask Bitcoin-US dollar exchange rate}
    \item{\code{BTCUSDBidMax}}{Daily maximum of bid Bitcoin-US dollar exchange rate}
    \item{\code{BTCUSDBidMin}}{Daily minimum of bid Bitcoin-US dollar exchange rate}
  }
}

\source{Dukascopy publicly available tick-by-tick data
}

\keyword{datasets}

\name{rmev}
\alias{rmev}
\alias{pmev}
\alias{qmev}
\alias{dmev}

   

\title{The Metastatistical Extreme Value Distribution (MEVD)}

\description{Density, distribution function, quantile function and random generation for the MEV distribution with shape parameter 'w', scale parameter 'c' and mean number of wet days 'n'.}

\usage{
rmev(N, w, c, n)
pmev(q, w, c, n)
qmev(p, w, c, n)
dmev(x, w, c, n)
}

\arguments{
 \item{x,q}{vector or single values of quantiles for \code{dmev} and \code{pmev}.}
 \item{p}{vector or single value of probabilities for \code{qmev}.}
 \item{N}{Number of observations.}
 \item{w,c}{vector or single values of shape and scale parameter of the MEV distribution.}
 \item{n}{Mean number of wet days. Can be a vector for yearly MEV calculations (i.e. one value per year), or a single value for the simplified MEV (see details).}
}

\details{
\code{rmev,pmev,qmev,dmev}{ blabla}.} 

\value{
No default values are assumed for the shape and scale parameters.
}


\examples{
dmev(x=1000,w=0.7,c=70,n=100)
dmev(x=1000,w=c(0.7,0.8),c=c(70,80),n=c(100,120))

qmev(p=0.5,w=0.7,c=70,n=100)
qmev(p=0.5,w=c(0.7,0.8),c=c(70,80),n=c(100,120))

rmev(N=10,w=0.7,c=70,n=100)
rmev(N=10,w=c(0.7,0.8),c=c(70,80),n=c(100,120))
}



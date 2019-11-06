\name{plot.mevd}
\alias{plot.mevd}
\title{Plot return levels of MEVD fit}
\usage{
     \method{plot}{mevd}(x, q = c(2,5,10,20,50,100,200), 
     ci = FALSE, type=c("all","rl","qq"), ...)
}
\arguments{
     \item{x}{mevd object.}
     \item{q}{vector of return periods, \eqn{q > 1}.}
     \item{ci}{if \code{ci = TRUE}, confidence intervals will be computed.}
     \item{type}{if omitted a return level plot \code{type = 'rl'} and a qq-plot \code{type = 'qq'} are shown. }
     \item{...}{Further parameters may also be supplied as arguments. See e.g. \code{\link[graphics]{plot}}. }
}
\description{
     Plot return levels of an mevd object.
}
 
\examples{
w <- 0.7
c <- 48
data <- matrix(NA,365,50)
for(y in 1:50){
     data[,y] <- rweibull(365,shape=w,scale=c)
}

fit <- fmev(data)
plot(fit)

}
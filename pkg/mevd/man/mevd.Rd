\name{mevd}
\alias{mevd}
\alias{fmev}
\alias{fit_mev}
\alias{rlmev}
   

\title{Fitting the Metastatistical Extreme Value Distribution (MEV) to rainfall data}

\description{
  Fit a MEV distribution to rainfall observations using maximum likelihood or probability weighted moments.
}
\usage{
fmev(data, n, threshold = 0, type =c("simple","annual"), 
     method=c("pwm","mle"))

fit_mev(data,method)

rlmev(q, w, c, n)

}

\arguments{
     \item{data}{ data must be either numeric or matrix with years as columns and all values per year as rows (see details).}
     \item{n}{ mean number of wet days, i.e. precipitation values >= \code{threshold} per year (see details). Must be a vector for \code{type='annual'} or a single number for \code{type='simple'}. }
     \item{threshold}{ \eqn{data < threshold} is set to NA.}
     \item{type}{ \code{type='simple'} computes one single value for parameters C and w, whereas \code{type='annual'} calculates one C and w for each year in \code{data}. \code{type='annual'} can therefore only be used, when \code{data} is a matrix. In both cases, if not given explicitly, n will be computed as the mean number of wet days (i.e. days with precipitation >= \code{threshold}) from \code{data}. The \code{default} is \code{simple} (see details).}
     \item{method}{ MEV parameters C and w can be computed using probability weighted moments \code{method = 'pwm'} or maximum likelihood \code{method = 'mle'}. The \code{default} is \code{pwm}. (see details).}
     
     \item{q,w,c}{ \code{rlmev}: return periods \eqn{q > 1}, shape parameter \code{w} and scale parameter \code{c} for which return levels shall be computed.}
}

\details{
With the aim of weakening the requirement of an asymptotic assumption for the GEV distribution, a metastatistical approach was proposed Marani and Ignaccolo (2015). This is defined in terms of the distribution of the statistical parameters describing "ordinary" daily rainfall occurrence and intensity.The MEV accounts for the random process of event occurrence in each block and the possibly changing probability distribution of event magnitudes across different blocks, by recognizing the number of events in each block, n, and the values of the shape and scale parameters w and C of the parent Weibull distribution to be realisations of stochastic variables. The MEV can then be written as

\deqn{F = 1/T \sum(j=1,T) [1 - exp(-(x/C_j)^(w_j)]^(n_j)}

for \eqn{w > 0} and \eqn{C > 0}. With T fully recorded years, yearly C and w can be estimated by fitting a Weibull distribution to the values x of this year, and n is the number of ordinary events per year. In this case (\code{type = 'annual'}), single parameters w and C are estmiated for each year.

If the probability distribution of daily rainfall is assumed to be time-invariant, the MEV becomes simplified to 

\deqn{F = [1 - exp(-x/C)^w]^n}

with single values for the shape and scale parameters w and C. n is then the mean number of wet days at this location (\code{type = 'simple'}) (Marra et al., 2019; Schellander et al., 2019).

As is shown e.g. Schellander et al., 2019, probability weighted moments should be preferred over maximum likelihood for the estimation of the Weibull parameters w and C. Therefore \code{method = 'pwm'} is the default.


}

\value{
  An list of class \code{mevd} with components:
  \item{w}{ Shape parameter w of the MEVD, estimated with \code{method}. Single value if \code{type = 'simple'}, vector if \code{type = 'annual'}.}
  \item{C}{ Scale parameter C of the MEVD, estimated with \code{method}. Single value if \code{type = 'simple'}, vector if \code{type = 'annual'}.}
  \item{n}{ Number of wet days. Single value if \code{type = 'simple'}, vector if \code{type = 'annual'}.}
  \item{data}{ \eqn{data >= threshold} used to fit the MEVD.}
  \item{maxima}{ Yearly maxima estimated from data matrix.}
  \item{threshold}{ Threshold used to select data.}
  \item{method}{ Method used to fit the MEVD.}
  
  
  
  
}

\references{
 	Marani, M. and Ignaccolo, M. (2015) 'A metastatistical approach to rainfall extremes', Advances in Water Resources. Elsevier Ltd, 79(Supplement C), pp. 121-126. doi: 10.1016/j.advwatres.2015.03.001.
 	\cr\cr
 	Marra, F. et al. (2019) 'A simplified MEV formulation to model extremes emerging from multiple nonstationary underlying processes', Advances in Water Resources. Elsevier Ltd, 127(April), pp. 280-290. doi: 10.1016/j.advwatres.2019.04.002.
 	\cr\cr
 	Schellander, H., Lieb, A. and Hell, T. (2019) 'Error Structure of Metastatistical and Generalized Extreme Value Distributions for Modeling Extreme Rainfall in Austria', Earth and Space Science, 6, pp. 1616-1632. doi: 10.1029/2019ea000557.
 	\cr\cr
 	
}


\examples{
w <- 0.87
c <- 75
data <- matrix(NA,365,50)

# 50 years of raifall data
for(y in 1:50){
     data[,y] <- rweibull(365,shape=w,scale=c)
}

# mean number of wet days per year
threshold <- 20
n <- apply(data,2,function(x){
  length(which(x>threshold))
})

fit <- fmev(data,n=n,type="annual")
fit
plot(fit)

fit1 <- fmev(c(data),mean(n),type="simple")
plot(fit1, ci = TRUE, R = 150)
rlmev(2:50, fit1$w, fit1$c, fit1$n)
}
#------------------------------------------------------------------------------------------#
#  fit MEV based on Weibull variables 
#------------------------------------------------------------------------------------------#
fmev <- function(data, n = NULL, threshold = 0, type = c("simple","annual"), 
                 method=c("pwm","mle")){
  #Input:
  # data must be either a vector or a matrix of precipitation values
  if(!(is.matrix(data) | is.numeric(data)))
    stop("fmev: data must be either numeric or matrix with years as columns and all values per year as rows.")
  
  # if(!(type == "simple" | type == "annual"))
  #   stop("fmev: type must be one of 'simple' or 'annual'.")
  # 
  # if(!(method == "pwm" | method == "mle"))
  #   stop("fmev: method must be one of 'pwm' or 'mle'.")
  
  type <- match.arg(type)
  method <- match.arg(method)
  
  if(is.vector(data) && (type=="annual")) {
    stop("If 'data' is a vector type must be 'simple'")}
  
  # compute n from data
  if(is.null(n)){
    if(type == 'simple'){
      n <- mean(c(data[data>=threshold]),na.rm = TRUE)
    } else {
      n <- apply(data,2,mean,na.rm=TRUE)
    }
  } 
       
 
  # remove data smaller than threshold
  if(is.vector(data)){
    data <- data[data>=threshold]
  } else{
    data[data<threshold] <- NA
  }
  
  switch(type,
         # can only be used with data as matrix
         annual = {
           # data.frame for yearly parameters
           annual <- matrix(NA,nrow = ncol(data),ncol = 3)
           colnames(annual) <- c("w","c","n")
           annual <- as.data.frame(annual)
           
           for (i in 1:nrow(annual)){
             tryCatch({
               values <- as.numeric(na.omit(as.vector(data[,i])))
               fit <- fit_mev(values,method)
               annual$w[i] <- fit$w
               annual$c[i] <- fit$c
               #annual$n[i] <- n
             },error = function(e){cat("ERROR at year", paste(colnames(data)[i]),conditionMessage(e), "\n")})
           }
           w <- annual$w
           c <- annual$c
           #n <- annual$n
         },
         
         simple = {
           values <- as.numeric(na.omit(as.vector(data)))
           fit <- fit_mev(values,method)
           w <- fit$w
           c <- fit$c
         }
  )
  
  if(is.vector(data)){
      maxima <- c()
  } else{
    maxima <- apply(data,2,function(x) max(x, na.rm = TRUE)) 
  }
  
  res <- list(w=w,c=c,n=n, data=data, maxima=maxima, threshold=threshold, type=type, method=method)
  class(res) <- "mevd"
  return(res)
}




fit_mev <- function(data,method){
  
  switch(method,
  pwm = {
    data=sort(data)
    
    if(length(data)<2){
      stop("fit_mev: Too little data for fitting a Weibull-distribution (only one value greater or equal threshold)\n")
    }
    else{
      M0hat  = mean(data)
      M1hat  = 0
      N      = length(data) # sample size
      for (i in 1:N){
        M1hat   = M1hat + data[i]*(N - i)
      }
      
      M1hat = M1hat/(N*(N-1))
      c     = M0hat/gamma(log(M0hat/M1hat)/log(2)) # scale par
      w     = log(2)/log(M0hat/(2*M1hat)) # shape par
    }
    
    
  },
  
  mle = {
    
    if(length(data)<2){
      stop("fit_mev: Too little data for fitting a Weibull-distribution (only one value greater or equal threshold)\n")
    } else{
      estimation=eweibull(data, method = "mle")
      w=as.double(estimation$parameters[1])
      c=as.double(estimation$parameters[2])
    }
  }
         
         )
  
  if(length(data)<6){
    print("fit_mev: Maybe bad fitting results: less or equal 5 values used for fitting!")
  }
  
  return(list(w=w,c=c))
}





# x,w,C,n can be vectors or single values
#Distribution function
pmev <- function(q,w,c,n){
  nyears=length(n)
  ret=c()
  for(y in q){
    if(y>=0){
      val=sum((1-exp(-y^w/c^w))^n)/nyears
    }
    else{
      val=0
    }
    ret=c(ret,val)
  }
  return(ret)
}

#Density function
dmev <- function(x,w,c,n){
  nyears=length(n)
  ret=c()
  for(y in x){
    if(y>0){
      val=n*w*(y^(w-1)/(c^w))*exp(-y^w/c^w)*(1-exp(-y^w/c^w))^(n-1)
      val=sum(val)/nyears
    }
    else{
      val=0
    }
    ret=c(ret,val)
  }
  return(ret)
}

#Quantile function
qmev <- function(p,w,c,n){
  ret=rep(0,length(p))
  if(length(w)==1){
    for(i in 1:length(p)){
      if(p[i]==0){
        val=-Inf
      }
      else if(p[i]==1){
        val=Inf
      }
      else{
        val=c*(-log(1-p[i]^(1/n)))^(1/w)
      }
      ret[i]=val
    }
  }
  else if(length(w)>1){ #numerische Bestimmung
    for(i in 1:length(p)){
      if(p[i]==0){
        val=-Inf
      }
      else if(p[i]==1){
        val=Inf
      }
      else{
        min_fun=function(x){
          return(pmev(x,w,c,n)-p[i])
        }
        val=uniroot(min_fun,lower = 0, upper = 10^10)$root
      }
      ret[i]=val
    }
  }
  return(ret)
}

#Random generation
rmev <- function(N,w,c,n){
  x=runif(N)
  ret=qmev(x,w,c,n)
  return(ret)
}

#Return levels MEV
#q: return period(Vector)
rlmev <- function(q,w,c,n){
  if(all(q>1)){
    p=1/q
    ret=qmev(1-p,w,c,n)
    return(ret)
  }
  else{
    stop("return period 'q' has to be greater than 1")
  }
}

ci.mev <- function(x, alpha = 0.05, return.periods = c(2,10,20,30,50,75,100,150,200), R = 502){
  
  if(!inherits(x, "mevd"))
    stop("conf.int: x must be object of class 'mevd'")
  
  w <- x$w #shape
  c <- x$c #scale
  n <- x$n # number of wet days
  l <- length(x$data) 
  theta.hat <- c(shape=w,scale=c,n=n)
  
  
  # draw R samples of length l with w and C
  Z <- rweibull(n = l * R, shape = w, scale = c) 
  Z <- matrix(Z, l, R)
  
  
  # fit MEVD to simulated samples
  bfun <- function(z, n){
    #fit <- fitMEV_weibull(z, n, threshold = 0, type ="all_years", method="pwm") 
    fit <- fmev(z, n, threshold = x$threshold, type=x$type, method=x$method) 
    return(c(shape=fit$w,scale=fit$c,n=fit$n))
  }
  pars <- apply(Z, 2, bfun, n = n)
  shape <- pars["shape",]
  scale <- pars["scale",]
  n <- pars["n",]
  
  
  # compute return levels from R w and C parameters
  th <- rbind(shape, scale, n)
  th.est <- theta.hat
  rlfun <- function(theta, q) rlmev(q = q, w = theta[1], c = theta[2], n = theta[3])
  sam <- apply(th, 2, rlfun, q = return.periods)
  rownames(sam) <- paste0(return.periods, "-year")
  theta.hat <- rlmev(q = return.periods, w = th.est[1], c = th.est[2], n = th.est[3])
  
  
  
  # compute quantiles of simulated return levels
  out <- apply(sam, 1, quantile, probs = c(alpha/2,1 - alpha/2))
  out.names <- rownames(out)
  out <- rbind(out[1, ], theta.hat, out[2, ])
  rownames(out) <- c(out.names[1], "Estimate", out.names[2])
  colnames(out) <- rownames(sam)
  out <- t(out)
  
  return(out)
  
}





# Weibull Plotting-Position-Formula 
pp.weibull <- function(maxima){
  maxima=na.omit(maxima)
  xi=sort(maxima, decreasing = FALSE)
  Fi=(1:length(xi))/(length(xi)+1)
  TR=1/(1-Fi)
  return(TR)
}


# S3 method print
print.mevd <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("MEVD fitting\n\n")
  cat(paste0("type: ",x$type,"\n"))
  cat(paste0("Estimator: ",x$method,"\n"))

  cat("\nParameters:")
  cat("\nShape w:\n")
  shape <- x$w
  print.default(format(shape, digits = digits), print.gap = 2, quote = FALSE)
  cat("\nScale C:\n")
  scale <- x$c
  print.default(format(scale, digits = digits), print.gap = 2, quote = FALSE)
  cat("\nWet days n:\n")
  n <- x$n
  print.default(format(n, digits = digits), print.gap = 2, quote = FALSE)

  invisible(x)
}

 

# S3 method plot
plot.mevd <- function(x, q = c(2,5,10,20,50,100,200), ci = FALSE, type=c("all","rl","qq"), ...){
  
  type <- match.arg(type)
  
  if(type=="all" & !is.vector(x$data)){
    par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
  }
  
  
  if (is.element(type, c("all", "rl"))) {
    
    rls <- rlmev(q, x$w, x$c, x$n)
    plot(q, rls, type="l", xlab="Return Period [years]", ylab = "Return Level")
    
    # observations
    if(is.vector(x$data)){
      warning("plot.mevd: can't take yearly maxima of a vector")
    } else {
      obs.y <- x$maxima
      obs.x <- pp.weibull(obs.y)
      points(obs.x,sort(obs.y))
    }
    
    if(ci){
      bsd <- ci.mev(x, return.periods = q, ...)
      lines(q, bsd[,1],col="grey",lty=2)
      lines(q, bsd[,3],col="grey",lty=2)
    }
    
  }
  
  if (is.element(type, c("all", "qq"))) {
    if(!is.vector(x$data)){
      q.m <- rlmev(obs.x, x$w, x$c, x$n)
      plot(q.m,sort(obs.y), xlab="Model Quantiles", ylab="Empirical Quantiles")
      abline(c(1,1))
    }
  }
  
  title("MEVD", outer=TRUE)
  par(mfrow=c(1,1))
}

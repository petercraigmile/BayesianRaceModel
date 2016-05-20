
## ======================================================================
## Copyright 2005--2016, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================



## ======================================================================
## chain: a data structure used to store MCMC output.
##
## The collection of variables are stored in an R environment called
## 'the.chains' (to guard against global access).  For each variable
## the current state is stored in 'name' and the (partial) history of
## the variable is stored as a list of values in the 'name.chain'
## variable in the environment.
## ======================================================================



chain <- function (var.name, value=NULL, the.chains) {
    ## ======================================================================
    ## Create a new chain with a given name 'var.name' in the environment
    ## 'the.chains' and give it the initial value 'value'.
    ## ======================================================================

    assign(var.name, value, envir=the.chains)

    assign(paste(var.name, ".chain", sep=""), list(), envir=the.chains)

    invisible()
}


chain.update <- function (var.name, value, the.chains = .GlobalEnv, save = TRUE) {
    ## ======================================================================
    ## update the chain of name 'var.name' in the environment
    ## 'the.chains' with a new 'value'.  If 'save' is TRUE update the
    ## history to include the current value.
    ## ======================================================================
  
 assign(var.name, value, env=the.chains)

  ## if we want to append to the chain (when save is TRUE)
  if (save) {
    chain <- paste(var.name, ".chain", sep = "")
    m <- length(get(chain, env=the.chains)) + 1
    eval(substitute(a[[m]] <- av, list(a = as.name(chain), 
                                       av = value, m = m)), env=the.chains)
  }
}


chain.updates <- function (values, the.chains = .GlobalEnv, save = TRUE) {
    ## ======================================================================
    ## Expects 'value' to be a list of value, with names corresponding
    ## to chains stored in the environment 'the.chain'.  Update each
    ## of the variables in the list with the new values, and if 'save'
    ## is TRUE, update the history to include the current value.
    ## ======================================================================
    
    nms <- names(values)
    
    for (j in 1:length(values)) {

        var.j <- nms[[j]]
        val.j <- values[[j]]
    
        assign(var.j, val.j, env=the.chains)
        
        ## if we want to append to the chain (when save is TRUE)
        if (save) {
            chain <- paste(var.j, ".chain", sep = "")
            m <- length(get(chain, env=the.chains)) + 1
            eval(substitute(a[[m]] <- av, list(a = as.name(chain), 
                                               av = val.j, m = m)), env=the.chains)
        }
    }
}




get.chain <- function (var.name, index1, index2,
                       start=1, thin=1, transform,
                       the.chains) {
    ## ======================================================================
    ## From the environment 'the.chains', obtain the history of the
    ## variable of name 'var.name', subsetting the components of the
    ## variable by 'index1' and 'index2'.  You can subset which part
    ## of the history to return by changing the starting point
    ## 'start', and the 'thin' rate.  You can also transform the
    ## variable by including a 'transform' function.
    ## ======================================================================
     
    chain <- paste(var.name, ".chain", sep="")
    x <- get(chain, envir=the.chains)
    
    if (missing(index1)) {
        
        y <- unlist(x)
        label <- var.name    
    } else if (missing(index2)) {
        
        y <- sapply(x, function (z, index1) z[[index1]], index1=index1)
        label <- paste(var.name, "[", index1, "]", sep="")
    } else {
        
        y <- sapply(x, function (z, index1, index2) z[[index1]][[index2]],
                    index1=index1, index2=index2)
        label <- paste(var.name, index1, "[", index2, "]", sep="")   
    }
    
    sel <- seq(from=start, to=length(y), by=thin)
    
    if (missing(transform))
        x <- y[sel]
    else
        x <- transform(y[sel])
    
    attr(x, "label") <- label
    x
}




trace.plot <- function (var.name, index1, index2, start=1, thin=1, guide,
                        transform, the.chains) {
    ## ======================================================================
    ## Produce a trace plot -- the format of the function is for the
    ## mostpart the same as for get.chain.  The 'guide' argument when
    ## supplied adds a horizontal line at x='guide' (Can be useful for
    ## making comparisons.)
    ## ======================================================================

  y <- get.chain(var.name, index1, index2, start, thin, transform, the.chains)

  plot(ts(y), xlab="Iterations", ylab=attr(y, "label"))

  if (!missing(guide))
    abline(h=guide, col="blue", lwd=2)
  
  invisible()
}



my.save <- function (the.names, filename, the.path="output/") {
    ## ======================================================================
    ## Save the objects 'the.names' to the file 'filename', storing
    ## the files in the path 'the.path'.  The name of the machine is
    ## appended to the filename.
    ## ======================================================================

    hostname <- Sys.info()["nodename"]
    
    save(list=the.names, file=file.path(the.path,
                             paste(filename, "_", hostname, ".RData", sep="")))
}




## ======================================================================
## Copyright 2004-2015, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================


logit <- function (p) {
    ## ======================================================================
    ## The logit transformation of 'p'
    ## ======================================================================
  
  log(p/(1.0 - p))
}



inv.logit <- function (x) {
    ## ======================================================================
    ## The inverse logit transformation of 'x'
    ## ======================================================================
  
  exp(x) / (1.0 + exp(x))
}



## ======================================================================
## Copyright 2011-2015, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## Adapted from R code orginally written by Sungmin Kim.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================



sample.Gaussian.invGamma <- function (y, kappa) {
    ## ======================================================================
    ## Produce a posterior draw of the mean and variance for an IID
    ## N(mu, sigma^2) sample 'y', with a Gaussian-invGamma prior for
    ## mu and sigma^2, as specified in the kappa vector.
    ##
    ## kappa[1] = mu.0
    ## kappa[2] = nu
    ## kappa[3] = alpha
    ## kappa[4] = beta
    ##
    ## The Gaussian-invGamma prior is:
    ## mu | sigma^2 ~ N(mu.0, sigma^2 / nu)
    ## 1/sigma^2    ~ Gamma(alpha, beta)    
    ##
    ## Follows the parameterization in
    ## http://en.wikipedia.org/wiki/Normal-inverse_gamma_distribution
    ## ======================================================================

    n <- length(y)

    mu.0  <- kappa[1]
    nu    <- kappa[2]
    alpha <- kappa[3]
    beta  <- kappa[4]

    y.bar <- mean(y)
    sum.sq <- sum((y - y.bar)^2)
    
    sigma.sq <- 1.0 / rgamma(1, shape = alpha + 0.5 * n,
                             rate = beta + 0.5 * sum.sq + n*nu/(n + nu)* 0.5 * (y.bar - mu.0)^2)
    
    ## shape parameter of invGamma is rate parameter of Gamma.
    mu <- rnorm(1, (nu*mu.0 + n * y.bar)/(nu + n), sqrt(sigma.sq/(nu + n)))

    ## Return a vector -- first element is the mean parameter draw;
    ## the second element is the variance parameter draw.
    c(mu, sigma.sq)
}

 
r.Gaussian.invGamma <- function (n, kappa) {
    ## ======================================================================
    ## Produce 'n' draws from a Gaussian-invGamma distribution for
    ## mu and sigma^2, as specified in the 'kappa' vector.
    ##
    ## kappa[1] = mu.0
    ## kappa[2] = nu.0
    ## kappa[3] = alpha
    ## kappa[4] = beta
    ##
    ## The Gaussian-invGamma prior is:
    ## mu | sigma^2 ~ N(mu.0, sigma^2 / nu.0)
    ## 1/sigma^2    ~ Gamma(alpha, beta)    
    ##
    ## Follows the parameterization in
    ## http://en.wikipedia.org/wiki/Normal-inverse_gamma_distribution
    ## ======================================================================
    
    mu.0  <- kappa[1]
    nu.0  <- kappa[2]
    alpha <- kappa[3]
    beta  <- kappa[4]
    
    sigma.sq <- 1/rgamma(n, shape=alpha, scale=beta)
    
    mu <- rnorm(n, mu.0, sqrt(sigma.sq / nu.0))
    
    cbind(mu, sigma.sq)
}





## ======================================================================
## Copyright 2013-2014, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================



rmvnorm.cond.precision <- function (m, P, eigen.P) {
    ## ======================================================================
    ## Purpose: Sample from N( P^{-1} m, P^{-1} )
    ## The eigen decomposition of 'P' can be optionally provided in 'eigen.P'.
    ## Assumes: m and P are real-valued.
    ##          The length of m is equal to the dimension of the
    ##          square matrix P.
    ## ======================================================================

  if (missing(eigen.P)) {
    
    eigen.P <- eigen(P, symmetric=TRUE)
  }
  
  ev <- eigen.P$values

  if (!all(ev >= -sqrt(.Machine$double.eps) * abs(ev[1]))) {
    warning("P is numerically not positive definite")
  }
  
  A <- t(eigen.P$vectors) * sqrt(1/ev)
  
  drop(crossprod(A, A %*% m + rnorm(length(m))))
}



## ======================================================================
## Copyright 2014-16, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================


update.log.alpha.log.beta.logit.p <- function (log.alpha, log.beta, logit.p, T0, xi,
                                               log.alpha.pars,
                                               log.beta.mu, log.beta.sigma2,
                                               logit.p.mu, logit.p.sigma2,
                                               RTs, exposures, resp.new, block.sel,
                                               log.alpha.prop.sd,
                                               log.beta.prop.sd,
                                               logit.p.prop.sd, the.model) {
    ## ======================================================================
    ## Update log.alpha, log.beta and logit.p using random walk
    ## Metropolis-Hastings updates.
    ## The commented out code allows for proposal covariances to be
    ## provided for the three parameters by setting the global list objects
    ## La, Lb and Lp (each list item indexes the subjects).
    ## ======================================================================

    num.subjects <- length(RTs)
    num.blocks   <- ncol(log.alpha)
    
    the.log.alpha <- log.alpha
    the.log.beta  <- log.beta
    the.logit.p   <- logit.p

    ts <- 1:num.blocks

    jumps <- 0
    
    for (subj in 1:num.subjects) {

        RTs.subj <- RTs[[subj]]

        ## if (is.null(La)) {

        new.log.alpha <- rnorm(length(log.alpha[subj,]), log.alpha[subj,], log.alpha.prop.sd)
        ## } else {

        ##     new.log.alpha <- La[[subj]] %*% rnorm(length(log.alpha[subj,])) + log.alpha[subj,]
        ## }
        
        ## if (is.null(Lb)) {
        
        new.log.beta  <- rnorm(length(log.beta[subj,]),  log.beta[subj,],  log.beta.prop.sd)
        ## } else {

        ##     new.log.beta <- Lb[[subj]] %*% rnorm(length(log.beta[subj,])) + log.beta[subj,]
        ## }

        ## if (is.null(Lp)) {
            
        new.logit.p <- rnorm(length(logit.p[subj,]),   logit.p[subj,],   logit.p.prop.sd)
        ## } else {

        ##     new.logit.p <- Lp[[subj]] %*% rnorm(length(logit.p[subj,])) +  logit.p[subj,]
        ## }

        xi.sel <- which(xi[subj,]==1)
        
        log.alpha.mu <- log.alpha.pars[subj, 1] + ts * log.alpha.pars[subj, 2]

        log.numer <- sum(RT.log.likelihood(RTs.subj[xi.sel], T0[subj],
                                           resp.new[[subj]][xi.sel], exposures[[subj]][xi.sel],
                                           new.logit.p, new.log.alpha, new.log.beta,
                                           block.sel[xi.sel], the.model)) +
                         sum(dnorm(new.logit.p, logit.p.mu, sqrt(logit.p.sigma2), log=TRUE)) +
                             AR1LogLikelihood(new.log.alpha - log.alpha.mu,
                                              log.alpha.pars[subj,3],
                                              exp(log.alpha.pars[subj,4])) +
                                sum(dnorm(new.log.beta, log.beta.mu, sqrt(log.beta.sigma2), log=TRUE))

        log.denom <- sum(RT.log.likelihood(RTs.subj[xi.sel], T0[subj],
                                           resp.new[[subj]][xi.sel], exposures[[subj]][xi.sel],
                                           logit.p[subj,], log.alpha[subj,], log.beta[subj,],
                                           block.sel[xi.sel], the.model)) +
                           sum(dnorm(logit.p[subj,], logit.p.mu, sqrt(logit.p.sigma2), log=TRUE)) +
                               AR1LogLikelihood(log.alpha[subj,]-log.alpha.mu,
                                                log.alpha.pars[subj,3],
                                                exp(log.alpha.pars[subj,4])) +
                   sum(dnorm(log.beta[subj,], log.beta.mu, sqrt(log.beta.sigma2), log=TRUE))

        if (log(runif(1)) < (log.numer-log.denom)) {

            the.log.alpha[subj,] <- new.log.alpha
            the.log.beta[subj,]  <- new.log.beta
            the.logit.p[subj,]   <- new.logit.p
            jumps <- jumps + 1
        }
    }

    list(log.alpha = the.log.alpha,
         log.beta  = the.log.beta,
         logit.p   = the.logit.p,
         theta.jumps=jumps/num.subjects)
}


update.log.alpha.log.beta.W <- function (log.alpha, log.beta, T0, xi,
                                         log.alpha.pars,
                                         log.beta.mu, log.beta.sigma2,
                                         RTs, exposures, resp.new, block.sel,
                                         log.alpha.prop.sd,
                                         log.beta.prop.sd) {
    ## ======================================================================
    ## Updates for log.alpha and log.beta in the Weibull model.
    ## ======================================================================

    num.subjects <- length(RTs)
    num.blocks   <- ncol(log.alpha)
    
    the.log.alpha <- log.alpha
    the.log.beta  <- log.beta

    ts <- 1:num.blocks

    jumps <- 0
    
    for (subj in 1:num.subjects) {

        RTs.subj <- RTs[[subj]]

        new.log.alpha <- rnorm(length(log.alpha[subj,]), log.alpha[subj,], log.alpha.prop.sd)
        new.log.beta  <- rnorm(length(log.beta[subj,]),  log.beta[subj,],  log.beta.prop.sd)

        xi.sel <- which(xi[subj,]==1)
        
        log.alpha.mu <- log.alpha.pars[subj, 1] + ts * log.alpha.pars[subj, 2]

        log.numer <- sum(RT.log.likelihood.W(RTs.subj[xi.sel], T0[subj],
                                             resp.new[[subj]][xi.sel], exposures[[subj]][xi.sel],
                                             new.log.alpha, new.log.beta,
                                             block.sel[xi.sel])) +
                             AR1LogLikelihood(new.log.alpha-log.alpha.mu,
                                              log.alpha.pars[subj,3],
                                              exp(log.alpha.pars[subj,4])) +
                                sum(dnorm(new.log.beta, log.beta.mu, sqrt(log.beta.sigma2), log=TRUE))

        log.denom <- sum(RT.log.likelihood.W(RTs.subj[xi.sel], T0[subj],
                                             resp.new[[subj]][xi.sel], exposures[[subj]][xi.sel],
                                             log.alpha[subj,], log.beta[subj,],
                                             block.sel[xi.sel])) +
                               AR1LogLikelihood(log.alpha[subj,]-log.alpha.mu,
                                                log.alpha.pars[subj,3],
                                                exp(log.alpha.pars[subj,4])) +
                                   sum(dnorm(log.beta[subj,], log.beta.mu, sqrt(log.beta.sigma2), log=TRUE))

        if (log(runif(1)) < (log.numer-log.denom)) {

            the.log.alpha[subj,] <- new.log.alpha
            the.log.beta[subj,]  <- new.log.beta
            jumps <- jumps + 1
        }
    }

    list(log.alpha = the.log.alpha,
         log.beta  = the.log.beta,
         theta.jumps=jumps/num.subjects)
}




update.log.alpha.pars <- function (log.alpha.pars,
                                   log.alpha,
                                   phi.prop.sd, log.sigma.prop.sd,
                                   log.alpha.mu.prior,
                                   log.alpha.slope.prior,
                                   log.alpha.phi.prior,
                                   log.alpha.log.sigma.prior) {
    ## ======================================================================
    ## For each subject, update 'log.alpha.pars', the time series
    ## parameters for 'log.alpha' over blocks.
    ##
    ## phi.prop.sd       : proposal stdev for phi
    ## log.sigma.prop.sd : proposal stdev for log sigma
    ##
    ## log.alpha.mu.prior        : (mean, variance) for mu prior
    ## log.alpha.slope.prior     : (mean, variance) for slope prior
    ## log.alpha.phi.prior       : (mean, variance) for phi prior
    ## log.alpha.log.sigma.prior : (mean, variance) for log sigma prior
    ## ======================================================================

    num.subjects <- nrow(log.alpha)
    num.blocks   <- ncol(log.alpha)

    X  <- cbind(1, 1:num.blocks)    
    ts <- 2:num.blocks
 
    beta.mean <- c(log.alpha.mu.prior[1], log.alpha.slope.prior[1])
    beta.prec <- diag(1/c(log.alpha.mu.prior[2], log.alpha.slope.prior[2]))

    ## We store the new parameter values in this matrix.
    new.log.alpha.pars <- matrix(NA, num.subjects, 4)

    ## Keep track of how many jumps we make.
    jumps <- 0

    ## For each subject...
    for (subj in 1:num.subjects) {

        phi       <- log.alpha.pars[subj, 3]
        log.sigma <- log.alpha.pars[subj, 4]
        sigma2    <- exp(2*log.sigma)
        la        <- log.alpha[subj,]

        ## The update for the intercept and slope ('new.beta') is a
        ## multivariate normal draw.
        
        Z      <- rbind(X[1,] * sqrt(1-phi^2), X[ts,] - phi * X[ts-1,])
        innovs <-     c(la[1] * sqrt(1-phi^2), la[ts] - phi * la[ts-1])
        
        m <- crossprod(Z, innovs)/sigma2 + beta.prec %*% beta.mean
        P <- crossprod(Z)/sigma2 + beta.prec
        
        new.beta <- rmvnorm.cond.precision(m, P)
        
        z <- la - drop(X %*% new.beta)

        ## 'phi' and 'log.sigma' are sampled via Metropolis-Hastings.

        new.phi       <- rnorm.truncated(1, phi, phi.prop.sd, -1, 1) # a truncated normal random walk (not symmetric)
        new.log.sigma <- rnorm(1, log.sigma, log.sigma.prop.sd) # a random walk for log sigma

        log.numer <-  AR1LogLikelihood(z, new.phi, exp(new.log.sigma)) +
            dnorm.truncated(new.phi, log.alpha.phi.prior[1], sqrt(log.alpha.phi.prior[2]), -1, 1, log=TRUE) +
                dnorm(new.log.sigma, log.alpha.log.sigma.prior[1], sqrt(log.alpha.log.sigma.prior[2]), log=TRUE) -
                    dnorm.truncated(new.phi, phi, phi.prop.sd, -1, 1, log=TRUE)
        
        log.denom <-  AR1LogLikelihood(z, phi, sqrt(sigma2)) +
            dnorm.truncated(phi, log.alpha.phi.prior[1], sqrt(log.alpha.phi.prior[2]), -1, 1, log=TRUE) +
                dnorm(log.sigma, log.alpha.log.sigma.prior[1], sqrt(log.alpha.log.sigma.prior[2]), log=TRUE) -
                    dnorm.truncated(phi, new.phi, phi.prop.sd, -1, 1, log=TRUE)

        ## Do we accept these new values for phi and log sigma?
        if (log(runif(1)) < (log.numer-log.denom)) {

            new.log.alpha.pars[subj,] <- c(new.beta, new.phi, new.log.sigma)
            jumps <- jumps + 1

        } else {

            new.log.alpha.pars[subj,] <- c(new.beta, phi, log.sigma)
        }
    }

    list(log.alpha.pars=new.log.alpha.pars,
         log.alpha.pars.jumps=jumps/num.subjects)
}





update.logit.p.mu.sigma2 <- function (logit.p, logit.p.mu.sigma2.priors) {
    ## ======================================================================
    ## Produce a posterior draw for the logit.p hyperparameters under
    ## the assumption that column 'c' of logit.p is a Gaussian-Inverse
    ## Gamma distribution with parameters logit.p.mu.sigma2.prior[[c]].
    ## ======================================================================

    the.updates <- sapply(1:ncol(logit.p), function (j)
        sample.Gaussian.invGamma(logit.p[,j], logit.p.mu.sigma2.priors[[j]]))
    
    list(logit.p.mu=the.updates[1,],
         logit.p.sigma2=the.updates[2,])
}



update.log.beta.mu.sigma2 <- function (log.beta, log.beta.mu.sigma2.priors) {
    ## ======================================================================
    ## Produce a posterior draw for the log.beta hyperparameters under
    ## the assumption that column 'c' of log.beta is a Gaussian-Inverse
    ## Gamma distribution with parameters log.beta.mu.sigma2.prior[[c]].
    ## ======================================================================

    the.updates <- sapply(1:ncol(log.beta), function (j)
        sample.Gaussian.invGamma(log.beta[,j], log.beta.mu.sigma2.priors[[j]]))

    list(log.beta.mu=the.updates[1,],
         log.beta.sigma2=the.updates[2,])    
}




update.T0 <- function (T0, T0.mu, T0.sigma2,
                       log.alpha, log.beta, logit.p,
                       xi, muH, sigma2.H,
                       RTs, exposures, resp.new, block.sel,
                       T0.mu.sigma2.prior, the.model) {
    ## ======================================================================
    ## First update 'T0' for each subject assuming a prior of
    ## N(T0.mu, T0.sigma2) truncated to the domain [0, T0.max].
    ## Then update the hyperparameter T0.mu and T0.sigma2 assumption a
    ## Gaussian-Inverse Gamma prior with parameters 'T0.mu.sigma2.prior'.
    ##
    ## log.alpha, log.beta, logit.p: parameters of the decision process.
    ## xi : indicates for each subject and trial the status of
    ##      'abrupt' (0), 'decision process' (1) or 'too slow' (2).
    ## muH, sigma2.H : parameters of the too slow process.
    ## RTs           : the response time for each subject
    ##                 (one list item for each subject).
    ## exposure      : the number of exposures for each subject.
    ## resp.new      : indicates new response or not.
    ## block.sel     : indicates for each subject and trial the block.
    ##                 number that trial is in.
    ## ======================================================================

    num.subjects <- length(RTs)
    
    ## The new T0 values will be kept in this vector
    T0.vec <- rep(NA, num.subjects)

    ## How many jumps do we make?
    jumps <- 0

    ## For each subject...
    for (subj in 1:num.subjects) {

        RTs.subj <- RTs[[subj]]

        ## Which trials are decision ones or too slow ones?
        xi.sel1 <- which(xi[subj,]==1)
        xi.sel2 <- which(xi[subj,]==2)

        rt1 <- RTs.subj[xi.sel1]
        rt2 <- RTs.subj[xi.sel2]        

        ## T0 lives in the interval [0, T0.max] -- T0.max is the minimum RT of the decision and too slow RTs.
        T0.max <- min(c(rt1, rt2))

        ## Propose a new value for T0
        new.T0 <- T0.max * rbeta(1, 6, 1)

        ## Calculate the log prior density (ignoring the normalizing constant due to T0.max which cancels out).
        log.prior.old <-  dbeta(T0[subj]/T0.max, 6, 1, log=TRUE)
        
        if (log.prior.old==-Inf) {  ## if we are currently at somewhere impossible jump
            
            T0.vec[subj] <- new.T0
            jumps <- jumps + 1
            
        } else {

            ## the effect of the prior on T0 and proposal
            log.pi <- dnorm.truncated(new.T0, T0.mu, sqrt(T0.sigma2), 0, T0.max, log=TRUE) -
                dbeta(new.T0/T0.max, 6, 1, log=TRUE) -
                    dnorm.truncated(T0[subj], T0.mu, sqrt(T0.sigma2), 0, T0.max, log=TRUE) +
                        log.prior.old

            if (length(xi.sel1)>0) {

                ## if there is a contribution to the likelihood due to having decision trials add it here.
                if (the.model=="W") {

                    log.pi <- log.pi + sum(RT.log.likelihood.W(rt1, new.T0,
                                                               resp.new[[subj]][xi.sel1], exposures[[subj]][xi.sel1],
                                                               log.alpha[subj,], log.beta[subj,],
                                                               block.sel[xi.sel1])) -
                                                                   sum(RT.log.likelihood.W(rt1, T0[subj],
                                                                                           resp.new[[subj]][xi.sel1], exposures[[subj]][xi.sel1],
                                                                                           log.alpha[subj,], log.beta[subj,],
                                                                                           block.sel[xi.sel1]))
                    
                } else {

                    log.pi <- log.pi + sum(RT.log.likelihood(rt1, new.T0,
                                                             resp.new[[subj]][xi.sel1], exposures[[subj]][xi.sel1],
                                                             logit.p[subj,], log.alpha[subj,], log.beta[subj,],
                                                             block.sel[xi.sel1], the.model)) -
                                                                 sum(RT.log.likelihood(rt1, T0[subj],
                                                                                       resp.new[[subj]][xi.sel1], exposures[[subj]][xi.sel1],
                                                                                       logit.p[subj,], log.alpha[subj,], log.beta[subj,],
                                                                                       block.sel[xi.sel1], the.model))
                }
                
            }
            
            if (length(xi.sel2)>0) {

                ## if there is a contribution to the likelihood due to having too slow trials add it here.
                log.pi <- log.pi +
                    sum(dlnorm(rt2 - new.T0,  muH[subj],  sqrt(sigma2.H), log=TRUE)) -
                        sum(dlnorm(rt2 - T0[subj], muH[subj], sqrt(sigma2.H), log=TRUE))
            }


            ## Now see if the proposal is accepted.
            if (log(runif(1)) < log.pi) {
                
                T0.vec[subj] <- new.T0
                jumps <- jumps + 1
            } else {
                T0.vec[subj] <- T0[subj]
            }
        }
    }

    ## Now update the priors on T0.
    the.updates <- sample.Gaussian.invGamma(T0.vec, T0.mu.sigma2.prior)
    
    list(T0=T0.vec,
         T0.jumps=jumps/num.subjects,
         T0.mu=the.updates[1],
         T0.sigma2=the.updates[2])
}




update.xi.qvec <- function (qvec, T0,
                            logit.p, log.alpha, log.beta,
                            muG, sigma2.G, muH, sigma2.H,
                            RTs, exposures, resp.new, block.sel,
                            qvec.prior, the.model) {
    ## ======================================================================
    ## Update xi and the q.
    ## ======================================================================

    num.subjects <- length(RTs)
    num.trials <- length(RTs[[1]])
    num.blocks <- max(block.sel)

    new.xi <- t(sapply(1:num.subjects, function (subj) {

        RTs.subj <- RTs[[subj]]
        
        ## G
        w1 <- exp(log(qvec[subj, block.sel, 1]) +
                      dlnorm(RTs.subj, muG[subj], sqrt(sigma2.G), log=TRUE))
        
        ## C
        if (the.model=="W") {

            w2 <- exp(log(qvec[subj, block.sel, 2]) +
                          RT.log.likelihood.W(RTs.subj, T0[subj],
                                              resp.new[[subj]], exposures[[subj]],
                                              log.alpha[subj,], log.beta[subj,],
                                              block.sel))
            
        } else {
            
            w2 <- exp(log(qvec[subj, block.sel, 2]) +
                          RT.log.likelihood(RTs.subj, T0[subj],
                                            resp.new[[subj]], exposures[[subj]],
                                            logit.p[subj,], log.alpha[subj,], log.beta[subj,],
                                            block.sel, the.model))
        }
        
        ## H
        w3 <- exp(log(qvec[subj, block.sel, 3]) +
                      dlnorm(RTs.subj-T0[subj], muH[subj], sqrt(sigma2.H), log=TRUE))
        
        u <- runif(num.trials) * (w1+w2+w3)
        
        z <- rep(2, num.trials)
        z[u<(w1+w2)] <- 1
        z[u<w1] <- 0
        z
    }))

    ## update qvec
    new.qvec <- array(NA, c(num.subjects, num.blocks, 3))

    for (subj in 1:num.subjects) {

        z <- rep(0, num.blocks)

        new.q.subj <- .C("sample_q",
                         as.integer(new.xi[subj,]),
                         as.integer(num.trials),
                         as.integer(block.sel),
                         as.integer(num.blocks),
                         as.double(qvec.prior),
                         q0=as.double(z),
                         q1=as.double(z),
                         q2=as.double(z))

        new.qvec[subj,,1] <- new.q.subj$q0
        new.qvec[subj,,2] <- new.q.subj$q1
        new.qvec[subj,,3] <- new.q.subj$q2
    }

    list(xi=new.xi, qvec=new.qvec)
}



update.muG.muH <- function (xi, T0, RTs, 
                            muG.prior, sigma2.G,
                            muH.prior, sigma2.H) {
    ## ======================================================================
    ## Update 'muG' and 'muH' the parameters of the abrupt and too
    ## slow log normal distributions.
    ##
    ## RTs       : list of RTs for each subjects
    ## muG.prior : mean and variance for the normal prior for muG
    ## sigma2.G  : the variance of the abrupt time distribution.
    ## muH.prior : mean and variance for the normal prior for muH
    ## sigma2.H  : the variance of the too slow time distribution.
    ## ======================================================================

    num.subjects <- length(RTs)

    sum.G <- sum.H <- nG <- nH <- rep(NA, num.subjects)

    for (subj in 1:num.subjects) {

        xi.subj <- xi[subj,]
        RT.subj <- RTs[[subj]]

        xi.0 <- which(xi.subj==0)
        xi.2 <- which(xi.subj==2)

        nG[subj] <- length(xi.0)
        nH[subj] <- length(xi.2)

        if (nG[subj]>0) {
            sum.G[subj] <- sum(log(RT.subj[xi.0]))
        } else {
            sum.G[subj] <- 0
        }

        if (nH[subj] > 0) {
            sum.H[subj] <- sum(log(RT.subj[xi.2] - T0[subj]))
        } else {
            sum.H[subj] <- 0
        }
    }

    mG <- muG.prior[1]/muG.prior[2] + sum.G/sigma2.G
    mH <- muH.prior[1]/muH.prior[2] + sum.H/sigma2.H
    
    pG <- 1.0/muG.prior[2] + nG/sigma2.G
    pH <- 1.0/muH.prior[2] + nH/sigma2.H

    list(muG = rnorm(num.subjects, mG/pG, sqrt(1.0/pG)),
         muH = rnorm(num.subjects, mH/pH, sqrt(1.0/pH)))
}



## ======================================================================
## Copyright 2005--2016, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================

## Requires the truncatedNormals R library
## which is available from http://www.stat.osu.edu/~pfc/software/.
library(truncatedNormals)

## Utility code
source("functions/utility.R")

## The functions that calculate the RT likelihood
dyn.load("functions/RT_likelihood.so")
source("functions/RT_likelihood.R")

## Functions that calculate the AR(1) likelihood
library(Rcpp)
sourceCpp("functions/AR1.cpp")

## Read in the functions that carry out the updates for the MCMC algorithm
source("MCMC_updates.R")




hyperparameters <- function (the.model) {

    if (the.model=="MG3") {
        
        log.beta.mu.sigma2.priors <-
            list(c(log(400), 2, 20, 8),
                 c(log(600), 2, 20, 8),
                 c(log(600), 2, 20, 8),
                 c(log(400), 2, 20, 8))

        logit.p.mu.sigma2.priors <-
            list(c(logit(0.5), 1/4, 3, 5),
                 c(logit(0.5), 1/4, 3, 5),
                 c(logit(0.5), 1/4, 3, 5))
        
    } else if (the.model=="MG4") {
        
        log.beta.mu.sigma2.priors <-
            list(c(log(400), 2, 20, 8),
                 c(log(600), 2, 20, 8),
                 c(log(600), 2, 20, 8),
                 c(log(400), 2, 20, 8),
                 c(log(400), 2, 20, 8))

        ## Change 0.5 to 0.99 to move p_l closer to one.
        logit.p.mu.sigma2.priors <-
            list(c(logit(0.5), 1/4, 3, 5),
                 c(logit(0.5), 1/4, 3, 5),
                 c(logit(0.5), 1/4, 3, 5),
                 c(logit(0.5), 1/4, 3, 5))        

    } else if (the.model=="MG4B") {

        log.beta.mu.sigma2.priors <-
            list(c(log(400), 2, 20, 8),
                 c(log(600), 2, 20, 8),
                 c(log(600), 2, 20, 8),
                 c(log(400), 2, 20, 8),
                 c(log(400), 2, 20, 8))
        
        logit.p.mu.sigma2.priors <-
            list(c(logit(0.5), 1/4, 3, 5),
                 c(logit(0.5), 1/4, 3, 5),
                 c(logit(0.5), 1/4, 3, 5))
        
    } else if (the.model=="W") {

        log.beta.mu.sigma2.priors <-
            list(c(log(1),    2, 10, 8),
                 c(log(1500), 2, 20, 8),
                 c(log(2000), 2, 20, 8),
                 c(log(2000), 2, 20, 8),
                 c(log(1500), 2, 20, 8))

        logit.p.mu.sigma2.priors <- NULL
        
    }

    list(

        log.beta.mu.sigma2.priors = log.beta.mu.sigma2.priors,
        
        logit.p.mu.sigma2.priors = logit.p.mu.sigma2.priors,

        T0.mu.sigma2.prior = c(100, 1, 3, 20),
        
        muG.prior = c(log(50), 0.01^2),
        sigma2.G  = 1^2,
        
        muH.prior = c(log(7000), 0.01^2),
        sigma2.H  = 0.4^2,
        
        qvec.prior = c(0.1, 1, 0.1),
        
        log.alpha.mu.prior        = c(log(3), 0.5^2),
        log.alpha.slope.prior     = c(0, 10^2),
        log.alpha.phi.prior       = c(0, 0.15^2),
        log.alpha.log.sigma.prior = c(log(0.3), 0.1^2)
    )
}


init.chains <- function (the.data, hyper, the.model) {
    
    ## Posteriors are stored in 'chs'
    chs <- new.env()

    chs$the.data  <- the.data
    chs$hyper     <- hyper
    chs$the.model <- the.model
    
    chain("log.alpha", matrix(0, the.data$num.subjects, the.data$num.blocks), chs)

    if (the.model=="MG3") {
    
        ## ( \beta^{O_{stim=old}}, \beta^{N_{stim=old}}, \beta^{O_{stim=new}}, \beta^{N_{stim=new}} )^T
        chain("log.beta", t(sapply(the.data$RTs, function (x) rnorm(4, log(mean(x)), 0.1))), chs)
        
        chain("logit.p",t(matrix(rnorm(3*the.data$num.subjects,
                                       rep(sapply(hyper$logit.p.mu.sigma2.priors, function (x) x[1])), 0.01), 3)),
              chs)
    } else if (the.model=="MG4") {
        
        ## ( \beta^{O_{stim=old}}, \beta^{N_{stim=old}}, \beta^{O_{stim=new}}, \beta^{N_{stim=new}, \beta^{A_{stim=old}} )^T
        chain("log.beta", t(sapply(the.data$RTs, function (x) rnorm(5, log(mean(x)), 0.1))), chs)
        
        chain("logit.p",t(matrix(rnorm(4*the.data$num.subjects,
                                       rep(sapply(hyper$logit.p.mu.sigma2.priors, function (x) x[1])), 0.01), 4)),
              chs)
    
    } else if (the.model=="MG4B") {

        ## ( \beta^{O_{stim=old}}, \beta^{N_{stim=old}}, \beta^{O_{stim=new}}, \beta^{N_{stim=new}, \beta^{A_{stim=old}} )^T
        chain("log.beta", t(sapply(the.data$RTs, function (x) rnorm(5, log(mean(x)), 0.1))), chs)
        
        chain("logit.p",t(matrix(rnorm(3*the.data$num.subjects,
                                       rep(sapply(hyper$logit.p.mu.sigma2.priors, function (x) x[1])), 0.01), 3)),
              chs)
        
    } else if (the.model=="W") {

        chain("log.beta", t(sapply(the.data$RTs, function (x) c(-1, rnorm(4, log(mean(x)), 0.1)))), chs)
    }
    
    chain("theta.jumps", NULL, chs)
    
    chain("T0", as.numeric(drop(sapply(the.data$RTs, function (x) 0.4*min(x)))), chs)
    
    chain("muG", rep(log(200), the.data$num.subjects), chs)
    chain("muH", rep(log(1700), the.data$num.subjects), chs)
    
    ## for the log alpha processes ( beta0, beta1, phi, log.sd )^T
    chain("log.alpha.pars", matrix(0, the.data$num.subjects, 4), chs)
    chain("log.alpha.pars.jumps", NULL, chs)
    
    ## initialize the q vector
    qvec0 <- array(NA, c(the.data$num.subject, the.data$num.blocks, 3))
    for (subj in 1:the.data$num.subject)
        for (block in 1:the.data$num.blocks)
            qvec0[subj,block,] <- hyper$qvec.prior
    
    chain("qvec", qvec0, chs)
    rm(qvec0)
    
    chain("xi", matrix(1, the.data$num.subjects, the.data$num.trials), chs)
    
    log.beta.mu.sigma2.init <- update.log.beta.mu.sigma2(chs$log.beta, hyper$log.beta.mu.sigma2.priors)
    chain("log.beta.mu", log.beta.mu.sigma2.init$log.beta.mu, chs)
    chain("log.beta.sigma2", log.beta.mu.sigma2.init$log.beta.sigma2, chs)
    rm(log.beta.mu.sigma2.init)

    if (the.model!="W") {
        
        logit.p.mu.sigma2.init <- update.logit.p.mu.sigma2(chs$logit.p, hyper$logit.p.mu.sigma2.priors)
        chain("logit.p.mu", logit.p.mu.sigma2.init$logit.p.mu, chs)
        chain("logit.p.sigma2", logit.p.mu.sigma2.init$logit.p.sigma2, chs)
        rm(logit.p.mu.sigma2.init)
    }
    
    chain("T0.mu", hyper$T0.mu.sigma2.prior[1], chs)
    chain("T0.sigma2", hyper$T0.mu.sigma2.prior[4], chs)
    
    chain("T0.jumps", NULL, chs)

    chs
}


fit.model <- function (chs, N.ITER, thin=20, burn.in=FALSE, save.every=10000,
                       model.name="", path="output/", trace.plots=TRUE) {

    RTs <- chs$the.data$RTs
    block.sel <- chs$the.data$block.sel
    
    print(system.time(for (ITER in 1:N.ITER) {

        if (burn.in) {
            save <- FALSE
        } else {
            save <- ITER %% thin==0
        }

        chain.updates(update.T0(chs$T0, chs$T0.mu, chs$T0.sigma2,
                                chs$log.alpha, chs$log.beta, chs$logit.p, chs$xi,
                                chs$muH, chs$hyper$sigma2.H,
                                RTs, chs$the.data$exposures, chs$the.data$resp.new, block.sel,
                                chs$hyper$T0.mu.sigma2.prior, chs$the.model),
                      chs, save=save)
                
        chain.updates(update.log.beta.mu.sigma2(chs$log.beta, chs$hyper$log.beta.mu.sigma2.priors),
                       chs, save=save)

        chain.updates(update.muG.muH(chs$xi, chs$T0, RTs,
                                     chs$hyper$muG.prior, chs$hyper$sigma2.G,
                                     chs$hyper$muH.prior, chs$hyper$sigma2.H),
                      chs, save=save)
            
        chain.updates(update.log.alpha.pars(chs$log.alpha.pars, chs$log.alpha, 
                                            0.15, 0.20,
                                            chs$hyper$log.alpha.mu.prior,
                                            chs$hyper$log.alpha.slope.prior,
                                            chs$hyper$log.alpha.phi.prior,
                                            chs$hyper$log.alpha.log.sigma.prior),
                       chs, save=save)

        if (chs$the.model!="W") {
            
            chain.updates(update.logit.p.mu.sigma2(chs$logit.p, chs$hyper$logit.p.mu.sigma2.priors),
                          chs, save=save)
            
            chain.updates(update.xi.qvec(chs$qvec, chs$T0,
                                         chs$logit.p, chs$log.alpha, chs$log.beta,
                                         chs$muG, chs$hyper$sigma2.G, chs$muH, chs$hyper$sigma2.H,
                                         RTs, chs$the.data$exposures, chs$the.data$resp.new,
                                         block.sel, chs$hyper$qvec.prior, chs$the.model),
                          chs, save=save)

            
        chain.updates(update.log.alpha.log.beta.logit.p(chs$log.alpha, chs$log.beta, chs$logit.p,
                                                        chs$T0, chs$xi,
                                                        chs$log.alpha.pars,
                                                        chs$log.beta.mu, chs$log.beta.sigma2,
                                                        chs$logit.p.mu, chs$logit.p.sigma2,
                                                        RTs, chs$the.data$exposures,
                                                        chs$the.data$resp.new, block.sel,
                                                        0.03, 0.03, 0.9, chs$the.model),
                      chs, save=save)

            
        } else {

            chain.updates(update.xi.qvec(chs$qvec, chs$T0,
                                         NULL, chs$log.alpha, chs$log.beta,
                                         chs$muG, chs$hyper$sigma2.G, chs$muH, chs$hyper$sigma2.H,
                                         RTs, chs$the.data$exposures, chs$the.data$resp.new,
                                         block.sel, chs$hyper$qvec.prior, chs$the.model),
                          chs, save=save)

            
            chain.updates(update.log.alpha.log.beta.W(chs$log.alpha, chs$log.beta, 
                                                      chs$T0, chs$xi,
                                                      chs$log.alpha.pars,
                                                      chs$log.beta.mu, chs$log.beta.sigma2,
                                                      RTs, chs$the.data$exposures,
                                                      chs$the.data$resp.new, block.sel,
                                                      0.05, 0.05),
                          chs, save=save)

        }

        if (ITER%%save.every==0) {

            print(round(c(theta=mean(unlist(chs$theta.jumps.chain)),
                          T0=mean(unlist(chs$T0.jumps.chain)),
                          log.alpha.pars=mean(unlist(chs$log.alpha.pars.jumps.chain))),2))
            
            if (trace.plots) {

                trace.plot(chs, model.name)
            }

            my.save("chs", model.name, path)
        }
    }))
    
    chs    
}



trace.plot <- function (chs, stub) {

    if (chs$the.model!="W") {
        n <- ncol(chs$log.beta) + ncol(chs$logit.p) + 4 + 3 + 4
    } else {
        n <- ncol(chs$log.beta) + 4 + 3 + 4
    }
    
    for (k in 1:4) {
        
        png(file=paste("traces/", stub, "_traces", k, ".png", sep=""), width=2000, height=700)
        par(mfrow=c(5,n), cex=0.8, mar=c(1.5,1.5,0.7,0.5), mgp=c(2,0.5,0), bty="L")
        
        for (mm in 1:5) {
            
            subj <- (k-1)*4+mm
            
            for (l in 1:ncol(chs$log.beta))
                ts.plot( sapply(chs$log.beta.chain, function (x) x[subj,l]), xlab="", ylab="")
            ts.plot( sapply(chs$T0.chain, function (x) x[subj]), xlab="", ylab="")
            ts.plot(sapply(chs$muG.chain, function (x) x[subj]), xlab="", ylab="", col="blue")
            ts.plot(sapply(chs$muH.chain, function (x) x[subj]), xlab="", ylab="", col="blue")
            sapply(1:4, function (l)
                ts.plot(sapply(chs$log.alpha.pars.chain, function (x) x[subj,l]),
                        xlab="", ylab="", col="blue"))
            for (l in 1:4 )
                ts.plot(sapply(chs$log.alpha.chain, function (x) x[subj,l]), xlab="", ylab="")
            if (chs$the.model!="W") {
                for (l in 1:ncol(chs$logit.p))
                    ts.plot(sapply(chs$logit.p.chain, function (x) x[subj,l]), xlab="", ylab="")
            }
            
        }
        dev.off()
        
    }
}



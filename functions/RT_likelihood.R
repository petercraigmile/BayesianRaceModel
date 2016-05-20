
## ======================================================================
## Copyright 2013-2016, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================

## ======================================================================
## Remember
## $( \beta^{O_{stim=old}}, \beta^{N_{stim=old}},
##    \beta^{O_{stim=new}}, \beta^{N_{stim=new}} )^T$
## =====================================================================


RT.log.likelihood <- function (RT, T0,
                               resp.new, exposures,
                               logit.p, log.alpha, log.beta,
                               block.sel, the.model="MG3") {
    ## ======================================================================
    ## For the decision process, calculates the log probability
    ## density function value (equivalently the trial by trial log
    ## likelihood values) for each response time in the vector 'RT'
    ## with a minimum decision time 'T0' and
    ##
    ## resp.new  : vector of responses (TRUE: new; FALSE: old) for each trial.
    ## exposures : the number of expoures for each trial.
    ## logit.p   : the logit of p (vector of length 3).
    ## log.alpha : a log alpha value for each block (log shape parameters).
    ## log.beta  : the log beta value defined above (scale parameters).
    ## block.sel : indicate for each trial number the block that trial
    ##             is currently in (used for defining the alpha
    ##             values for each trial).
    ##
    ## The result is a vector of length equal to the RT vector.
    ##
    ## ASSUMES: The dynamical C code library 'RT_likelihood.so' or
    ## 'RT_likelihood.dll' is already loaded.
    ## ======================================================================
    
    n <- length(RT)
    
    the.alpha <- exp(log.alpha)[block.sel]

    if (length(the.alpha) != n) {
        stop("The alpha vector repeated over blocks must be of the same length as the RT vector.")
    }

    p <- exp(logit.p)/(1.0 + exp(logit.p))

    if (the.model=="MG3") {
        the.C.function <- "threep_RT_log_likelihood"
    } else if (the.model=="MG4") {
        the.C.function <- "fourp_RT_log_likelihood"
    } else if (the.model=="MG4B") {
        the.C.function <- "fourp_RT_log_likelihood"
        p <- c(p, p[3])
    } else if (the.model == "MN4") {
        the.C.function <- "fourp_RT_log_likelihood_L"
    }
    
    .C(the.C.function,
       ll=double(n), as.integer(n),
       as.double(RT), as.double(T0),
       as.integer(as.numeric(resp.new)),
       as.integer(exposures),
       as.double(p),
       as.double(the.alpha),
       as.double(exp(log.beta)))$ll
}


## ======================================================================
## For the Weibull model
## ======================================================================

RT.log.likelihood.W <- function (RT, T0,
                                 resp.new, exposures,
                                 log.alpha, log.beta,
                                 block.sel) {
    
    n <- length(RT)
    
    the.alpha <- exp(log.alpha)[block.sel]

    if (length(the.alpha) != n) {
        stop("The alpha vector repeated over blocks must be of the same length as the RT vector.")
    }

    .C("Weibull_RT_log_likelihood",
       ll=double(n), as.integer(n),
       as.double(RT), as.double(T0),
       as.integer(as.numeric(resp.new)),
       as.integer(exposures),
       as.double(the.alpha),
       as.double(exp(log.beta)))$ll
}




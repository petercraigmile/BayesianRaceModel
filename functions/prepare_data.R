
## ======================================================================
## Copyright 2005--2016, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================



calc.dprime <- function (are.old, are.new, respond.old) {

    ## Calculate the hit rate
    H  <- mean(respond.old[are.old])

    ## And the false discovery rate
    FA <- mean(respond.old[are.new])
    
    qnorm(H) - qnorm(FA)
}


prepare.data <- function (the.df, block.subset) {

    ## We exclude the first block from the analysis to ensure that the
    ## subjects had become familiar with the task. We also exclude
    ## responses to the first two and the last two items on the study list
    ## (which were exposed only once) to avoid primacy and recency effects.
    df <- the.df[the.df$block!=0 & the.df$exposures!=-1, ]

    if (!missing(block.subset)) {
        df <- df[as.numeric(t(sapply(block.subset, function (b) which(df$block==b)))), ]
    }

    ## Prepare the the.data (the data used in the MCMC)
    the.data <- list(    
        num.blocks   = max(df$block),
        num.subjects = max(df$subject),
        num.trials   = sum(df$subject==1),
        block.sel    = df$block[df$subject==1])
    
    the.data$block.indexes <- lapply(1:the.data$num.blocks, function (j) which(the.data$block.sel==j))
    
    the.data$RTs <- lapply(1:the.data$num.subjects, function (j) df$RT[df$subject==j])
    
    the.data$exposures <- lapply(1:the.data$num.subjects,
                             function (j) df$exposures[df$subject==j])
    
    the.data$resp.new <- lapply(1:the.data$num.subjects,
                            function (j) as.numeric(df$response[df$subject==j]=='m'))
    
    the.data$fastest.10pc <- lapply(1:the.data$num.subjects,
                                function (j) order(the.data$RTs[[j]])[1:floor(length(the.data$RTs[[j]])/10)])

    the.data$dprime <- sapply(1:the.data$num.subjects, function (subj) calc.dprime(are.old = the.data$exposures[[subj]]>0,
                                                                                   are.new = the.data$exposures[[subj]]==0,
                                                                                   respond.old = the.data$resp.new[[subj]]==0))
    
    the.data$dprime.exp <- sapply(1:4, function (e)
        sapply(1:the.data$num.subjects, function (subj) calc.dprime(are.old = the.data$exposures[[subj]]==e,
                                                                    are.new = the.data$exposures[[subj]]==0,
                                                                    respond.old = the.data$resp.new[[subj]]==0)))
    
    the.data
}

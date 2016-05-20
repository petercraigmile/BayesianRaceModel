


## ======================================================================
## Copyright 2005--2016, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================




## Read in the arguments from the command line.

args <- commandArgs(trailingOnly = TRUE)

the.model <- args[1]
the.num   <- args[2]

eighteen <- (if (length(args) > 2) TRUE else FALSE)

print(the.model)
print(the.num)
print(eighteen)


## Read in the data

source("functions/prepare_data.R")

main <- read.csv("../Data/main.csv", header=TRUE)

the.data   <- prepare.data(main)
the.data18 <- prepare.data(main, 1:18)



source("fit_models.R")

## Generate the hyperparameters
hyper <- hyperparameters(the.model)

## Initialize the model
if (eighteen) {
    
    model.name <- paste(the.model, "_18_chain", the.num, sep="")    
    chs <- init.chains(the.data18, hyper, the.model)
    
} else {

    model.name <- paste(the.model, "_chain", the.num, sep="")    
    chs <- init.chains(the.data, hyper, the.model)
}


## Burn in for 3000 iterations
chs <- fit.model(chs, 3000, burn.in=TRUE)


## Now run for 3000 * 100 iterations, thinning every 20 iterations.
chs <- fit.model(chs, 3000*100, save.every=120, thin=20,
                 model.name=model.name, path="output")

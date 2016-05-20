#!/bin/bash

R CMD BATCH --no-save --no-restore "--args $1 $2 $3" run_MCMC.R out$3_$1_$2.txt &

#!/bin/bash

## ======================================================================
## Copyright 2005--2016, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================


## Fit the minimum gamma model, generating four chains.

./run.sh MG4 1
./run.sh MG4 2
./run.sh MG4 3
./run.sh MG4 4

## Refit using only blocks 1..18

./run.sh MG4 1 18
./run.sh MG4 2 18
./run.sh MG4 3 18
./run.sh MG4 4 18



## Fit the minimum gamma model with p_0=1, generating four chains.

./run.sh MG3 1
./run.sh MG3 2
./run.sh MG3 3
./run.sh MG3 4

## Refit using only blocks 1..18

./run.sh MG3 1 18
./run.sh MG3 2 18
./run.sh MG3 3 18
./run.sh MG3 4 18



## Fit the minimum gamma model with p_2=p_3, generating four chains.

./run.sh MG4B 1
./run.sh MG4B 2
./run.sh MG4B 3
./run.sh MG4B 4

## Refit using only blocks 1..18

./run.sh MG4B 1 18
./run.sh MG4B 2 18
./run.sh MG4B 3 18
./run.sh MG4B 4 18



## Fit the Weibull model, generating four chains.

./run.sh W 1
./run.sh W 2
./run.sh W 3
./run.sh W 4

## Refit using only blocks 1..18

./run.sh W 1 18
./run.sh W 2 18
./run.sh W 3 18
./run.sh W 4 18

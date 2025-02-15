################################################################################
################################################################################
# MAIN FILE TO CONDUCT ANALYSIS FROM
# "Plant and animal functional diversity drive mutualistic network assembly
# across an elevational gradient"
# by Albrecht et al. Nature Communications 9:3177 (2018)
# https://doi.org/10.1038/s41467-018-05610-w
################################################################################
# clean work space
################################################################################
rm(list = ls(all = TRUE))
################################################################################
# load packages
################################################################################
library(vegan)
library(FD)
library(ade4)
library(bipartite)
library(rjags) # (JAGS needs to be installed on the system before installing rjags)
library(parallel)
library(coda)
library(qgraph)
################################################################################
# RUN SCRIPTS FROM SOURCE FILES
################################################################################
nRuns <- 8 # number of parallel chains
nAdpt <- 1e3 # adaptive iterations for burn-in
nThin <- 1e2 # thinning interval
nIter <- ceiling(nThin * 4e3 / nRuns) # number of iterations per chain
source("include/00_sourceFunctions.R")
source("include/01_prepareData.R")
source("include/02_setupModel.R")
source("include/03_fitModel.R")
source("include/04_inferenceModel.R")
################################################################################
# The script is compatible with linux, macOS and windows systems.
# If the code fails to run, check the setup of your machine using
# sessionInfo() and compare it to the session info provided in the README file
# that accompanies the data package
################################################################################
# END OF SCRIPT
################################################################################
################################################################################

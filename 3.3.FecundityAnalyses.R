### Elk Abundance - Fecundity Model
### Analysis script
### Last updated: Oct. 29, 2025
### Contact: xprockox@gmail.com

############################################################################################
### packages

library(nimble)
library(MCMCvis)
library(coda)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tibble)

dat <- read.csv('data/intermediate/productivity.csv')


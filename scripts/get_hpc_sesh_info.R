library(purrr)
library(dplyr)
library(tidyr)
library(tibble)
library(magrittr)

library(readr)
library(stringi)
library(stringr)
library(yaml)
library(git2r)
library(rstan)
library(loo)

library(splines)
library(matrixStats)
library(abind)

sesh <- sessionInfo()

saveRDS(sesh, "processed/hpc_sesh.Rds")
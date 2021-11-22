library(here)
library(tidyverse)
source("estimate-lrt.R")
source("tensor-processing.R")

phi <- readRDS(here("..", "data", "phi-tensor.rds"))
dict <- readRDS(here("..", "data", "dictionary.rds"))
theta <- readRDS(here("..", "data", "theta-coords.rds"))
y <- readRDS(here("..", "data", "y-vector.rds"))
meta <- readRDS(here("..", "data", "meta-info.rds"))

ijk <- as.matrix(phi[,1:3])
val <- pull(phi[,4])
rm(phi)
dims <- c(ncol(dict), meta$nvoxels, meta$nstreamlines)
Xunfold <- unfold(ijk, val, dims)
p <- meta$nstreamlines


fit <- big_lrt(Xunfold, y, theta, dict, p)

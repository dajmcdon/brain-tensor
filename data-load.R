# Packages and parameters -------------------------------------------------

library(tidyverse)
library(Matrix)
source("tensor-processing.R")
nangles <- 60


# Read in data ------------------------------------------------------------
# note, these files have been removed to save space
groups <- scan(
  "~/Downloads/sparsemat_streamlineGroups.txt",
  what = integer())

phi <- read_tsv("~/Downloads/sparsemat_phi_subs_vals.txt",
                col_names = c("row","column","slice","value"),
                col_types = "iiid") %>%
  arrange(slice, column, row)
# rows, columns, slice, value
# has 129241 rows, 148647 cols, and 25900 slices

Dmat <- scan("~/Downloads/sparsemat_DictSig.txt")
Dmat <- matrix(Dmat, nrow = nangles, byrow = TRUE) #  "!! rowise"
Y <- scan("~/Downloads/sparsemat_demeanedSignal.txt") # "rowwise"

nvoxels <- length(Y) / nangles


# Tests & Visual checks ----------------------------------------------------

testthat::expect_gte(ncol(Dmat),  max(phi$row))
testthat::expect_lte(max(phi$column), nvoxels)


# make Y contain: all angles then all voxels (angles vary fastest)
Y <- c(matrix(Y, nrow = nangles, byrow = TRUE)) # !! rowwise

stringr::str_glue("The tensor has {nrow} rows, ", "{ncol} columns ",
                  "and {nstr} streamlines.\n", "The image has {nangles} angles ",
                  "and {rem} voxels.",
                  nrow = ncol(Dmat),
                  ncol = length(Y) / nangles,
                  nstr = max(phi$slice),
                  rem = length(Y) / nangles)

theta <- read_tsv("~/Downloads/sparsemat_direction_vectors.txt", 
                  col_names = FALSE)

saveRDS(Dmat, here::here("..", "data", "dictionary.rds"))
saveRDS(Y, here::here("..", "data", "y-vector.rds"))
saveRDS(as.matrix(theta), here::here("..", "data", "theta-coords.rds"))
saveRDS(phi, here::here("..", "data", "phi-tensor.rds"))


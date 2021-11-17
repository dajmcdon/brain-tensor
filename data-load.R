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
Dmat <- matrix(Dmat, nrow = nangles, byrow = TRUE) #  "rowise"
Y <- scan("~/Downloads/sparsemat_demeanedSignal.txt") # (seemingly "rowwise")

nvoxels <- length(Y) / nangles


# Tests & Visual checks ----------------------------------------------------

testthat::expect_gte(ncol(Dmat),  max(phi$row))
testthat::expect_lte(max(phi$column), nvoxels)

matplot(matrix(Y[c(
  outer(0:(nangles - 1) * nvoxels + 1, 0:24, "+")
)], nangles), ty = "l", ylab = "", main =  "Am I somewhat periodic?")
matplot(Dmat[,1:100], ty = "l", xlab = "angles", main = "Strong pattern?")

# make Y contain: all angles then all voxels (angles vary fastest)
Y <- c(matrix(Y, nrow = nangles, byrow = TRUE))

stringr::str_glue("The tensor has {nrow} rows, ", "{ncol} columns ",
                  "and {nstr} streamlines.\n", "The image has {nangles} angles ",
                  "and {rem} voxels.",
                  nrow = ncol(Dmat),
                  ncol = length(Y) / nangles,
                  nstr = max(phi$slice),
                  rem = length(Y) / nangles)

# Unfold once
unfold <- function(ijk, val, dims, mode = 1) {
  assertthat::assert_that(mode %in% 1:3)
  assertthat::assert_that(dim(ijk)[2] == 3)
  assertthat::assert_that(length(dims) == 3)
  assertthat::assert_that(all(apply(ijk, 2, max) <= dims))
  
  i <- ijk[ ,mode, drop = TRUE]
  # useful if more than 3 dimensions, useless here, just need first dim
  J <- cumprod(dims[-mode])
  J <- c(1, J[-length(J)])
  ik <- ijk[ ,-mode, drop = FALSE]
  j <- drop(1 + (ik - 1) %*% J)
  uj <- unique.default(j) 
  # tells us the original column numbers with nonzero entries
  dense_jdx <- match(j, uj) # map to dense columns (drop rest)
  newdims <- c(dims[mode], max(dense_jdx))
  
  Xwide <- sparseMatrix(i = i, j = dense_jdx, x = val, dims = newdims)
  
  return(list(Xwide = Xwide, uj = uj - 1)) # 0 indexed original cols
}

refold <- function(mat, orig_j, nrows, ncols) {
  # input matrix is stacked with faces side by side
  # but only those columns that contain non-zeros
  # output stacks the faces vertically
  n <- nrow(mat)
  p <- max(orig_j)
  assertthat::assert_that(
    as.double(n) * as.double(p) <= as.double(nrows) * as.double(ncols),
    msg = "input matrix data is too big for output size.")
  assertthat::assert_that(
    nrows %% n == 0,
    msg = "nrows should be divisible by nrow(mat)")
  
  ij <- Matrix:::non0ind(mat) + 1 # row and column indices, 0 based so + 1
  
  ij[,1] <- ij[,1] + (orig_j[ij[,2]] %/% ncols) * n # new row index
  orig_j <- orig_j + 1
  ij[,2] <- orig_j[ij[,2]] %% ncols # new col index
  ij[ij[,2] == 0, 2] <- ncols # replace 0 with ncols
  
  sparseMatrix(i = ij[,1], j = ij[,2], x = mat@x, dims = c(nrows, ncols))
}

# Now we use this for all A_x1_Phi multiplications
# Xunfold is Phi stacked wide with some extra info
dense_sp_mult <- function(A, Xunfold, nrow_out, ncol_out) {
  # returns AX where X is sparse and A is dense, but avoids 0-dense columns
  # requires the output from unfold()
  assertthat::assert_that(ncol(A) == nrow(Xunfold$Xwide), 
                          msg = "Nonconformable matrices")
  # always dense, 
  # then drop 0 (potentially no effect),
  # then stack vertically rather than horizontally
  refold(drop0(A %*% Xunfold$Xwide), Xunfold$uj, nrow_out, ncol_out) 
}



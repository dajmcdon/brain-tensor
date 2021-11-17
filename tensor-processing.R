# Unfold once
unfold <- function(ijk, val, dims, mode = 1) {
  assertthat::assert_that(mode %in% 1:3)
  assertthat::assert_that(dim(ijk)[2] == 3)
  assertthat::assert_that(length(dims) == 3)
  assertthat::assert_that(all(apply(ijk, 2, max) <= dims))
  
  newdims <- c(dims[mode], prod(dims[-mode]))
  
  i <- ijk[ ,mode, drop = TRUE]
  # useful if more than 3 dimensions, useless here, just need first dim
  J <- cumprod(dims[-mode])
  J <- c(1, J[-length(J)])
  ik <- ijk[,-mode,drop = FALSE]
  j <- drop(1 + (ik - 1) %*% J)
  
  Xwide <- sparseMatrix(i = i, j = j, x = val, dims = newdims)
  dense_cols <- which(diff(Xwide@p) > 0) # find cols without 0
  ncols <- length(dense_cols)
  Xwide_max <- nnzero(Xwide)
  # this removes any columns that are all 0 (fewer columns)
  # Doesn't change storage, but tells me which columns in AXwide will be 0
  Xwide <- sparseMatrix(i = i, p = c(Xwide@p[dense_cols], Xwide_max),
                        x = Xwide@x, dims = c(nrow(Xwide), ncols),
                        index1 = FALSE)
  
  return(list(Xdense = Xwide, dense_cols = dense_cols, out_col = newdims[2]))
}

refold <- function(mat, ncols) {
  n <- nrow(mat)
  p <- ncol(mat)
  assertthat::assert_that(p %% ncols == 0, 
                          msg = "Nonconformable dimensions requested.")
  
  ij <- Matrix:::non0ind(mat) # row and column indices, 0 based
  ij[,1] <- ij[,1] + 1
  ij[,1] <- ij[,1] + (ij[,2] %/% ncols) * n # new row index
  ij[,2] <- (ij[,2] + 1) %% ncols # new col index
  ij[ij[,2] == 0, 2] <- ncols
  
  sparseMatrix(i = i, j = j, x = mat@x, dims = c(p / ncols * n, ncols))
}

# Now we use this for all A_x1_Phi multiplications
# Xunfold is Phi stacked wide with some extra info
dense_sp_mult <- function(A, Xunfold, ncol_out) {
  # returns AX where X is sparse and A is dense, but avoids 0-dense columns
  # requires the output from unfold()
  assertthat::assert_that(ncol(A) == nrow(Xunfold$Xdense), 
                          msg = "Nonconformable matrices")
  
  out_row <- nrow(A)
  out_col <- Xunfold$out_col
  ncols <- length(Xunfold$dense_cols)
  
  AX <- A %*% Xunfold$Xdense # always dense
  AX <- sparseMatrix(i = rep(seq(out_row), ncols),
                     j = rep(Xunfold$dense_cols, each = out_row),
                     x = AX, 
                     dims = c(out_row, out_col))
  AX <- drop0(AX)
  refold(AX, ncol_out) # stack tensor vertically rather than horizontally
}



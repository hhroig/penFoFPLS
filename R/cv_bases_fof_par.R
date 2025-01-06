
cv_bases_fof_par <- function(X,
                             Y,
                             center = TRUE,
                             argvals_X = NULL,
                             argvals_Y = NULL,
                             num_bases_X = 20,
                             num_bases_Y = 20,
                             fda_basis_func_X = fda::create.bspline.basis,
                             fda_basis_func_Y = fda::create.bspline.basis,
                             penalty_X = 0,
                             penalty_Y = 0,
                             ncomp = min(10, ncol(X)),
                             folds = 5,
                             verbose = TRUE,
                             stripped = TRUE,
                             ...) {

  tictoc::tic("Crossvalidation")


  # Initialize grid for the first component
  num_bases_grid <- expand.grid(numbases_X = num_bases_X,
                                numbases_Y = num_bases_Y)

  if (is.numeric(folds)) {

    num_folds <- folds
    folds <- caret::createFolds(1:nrow(Y), k = num_folds)

  }else if (is.list(folds)) {

    num_folds <- length(folds)

  }

  # Initialize CVEs:
  CVEs_ncomp <- array(data = NA, dim = ncomp) # averaged
  names(CVEs_ncomp) <- paste0("ncomp_", 1:ncomp)

  MSE_ncomp_fold <- matrix(data = NA,
                           nrow = ncomp,
                           ncol = num_folds) # MSE per component per fold
  colnames(MSE_ncomp_fold) <- paste0("fold_", 1:num_folds)
  rownames(MSE_ncomp_fold) <- paste0("ncomp_", 1:ncomp)

  # Best number of bases (unique) for each ncomp_i:
  best_num_bases <- matrix(data = NA,
                           nrow = ncomp,
                           ncol = 2)


  for (ncomp_i in 1:ncomp) {

    if (verbose) {
      cat("Component ", ncomp_i, "/", ncomp, "\n")
    }


    i <- row_num_bases <- NULL
    MSE_lambda_fold <- foreach::foreach (i = 1:num_folds,
                                         .packages = c("penFoFPLS"),
                                         .combine = "cbind") %:%
      foreach::foreach(row_num_bases = 1:nrow(num_bases_grid),
                       .packages = c("penFoFPLS"),
                       .combine = 'c' ) %dopar%
                       {

                         # MSE_lambda_fold <- matrix(NA, nrow = nrow(penalty_grid), ncol = num_folds)
                         # for (i in 1:num_folds) {
                         #   for (row_num_bases in 1:nrow(num_bases_grid)) {

                         # build train
                         Y_fold_train <- Y[-folds[[i]], , drop = F]
                         X_fold_train <- X[-folds[[i]], , drop = F]

                         # build test:
                         Y_fold_test <- Y[folds[[i]], , drop = F]
                         X_fold_test <- X[folds[[i]], , drop = F]

                         basisobj_X <- fda_basis_func_X(rangeval = range(argvals_X),
                                                        nbasis = num_bases_grid[row_num_bases, "numbases_X"])

                         basisobj_Y <- fda_basis_func_Y(rangeval = range(argvals_Y),
                                                        nbasis = num_bases_grid[row_num_bases, "numbases_Y"])

                         res_fpls <- ffpls_bs(X = X_fold_train,
                                              Y = Y_fold_train,
                                              argvals_X = argvals_X,
                                              argvals_Y = argvals_Y,
                                              ncomp = ncomp_i,
                                              center = center,
                                              basisobj_X = basisobj_X,
                                              basisobj_Y = basisobj_Y,
                                              penalty_X = penalty_X,
                                              penalty_Y = penalty_Y,
                                              verbose = FALSE,
                                              stripped = stripped
                                              ...      )


                         # MSE_lambda_fold[row_lambda , i] <-
                         num_int_1d(argvals = argvals_Y,
                                    f_obs = colMeans( (Y_fold_test -
                                                         stats::predict(object = res_fpls,
                                                                        newdata = X_fold_test)[, , ncomp_i]  )^2 ) )



                         #   } # loop row_lambda
                         # } # loop fold

                       } # nested loop parallel


    # Averaged MSE_fold:
    CVEs_ncomp_bases <- rowMeans(MSE_lambda_fold)

    # Best penalties per component:
    sel_num_bases <- which.min(CVEs_ncomp_bases)
    best_num_bases[ncomp_i, ] <- as.numeric(penalty_grid[sel_num_bases, ])

    # Save the folds-averaged CV error:
    CVEs_ncomp[ncomp_i] <- CVEs_ncomp_bases[sel_num_bases]

    # Save MSEs per fold, for the best lambda:
    MSE_ncomp_fold[ncomp_i, ] <- MSE_lambda_fold[sel_num_bases, ]


  } # loop in ncomp: number of components


  names(CVEs_ncomp) <- paste0("ncomp_", 1:ncomp)
  rownames(best_num_bases) <- paste0("ncomp_", 1:ncomp)
  colnames(MSE_ncomp_fold) <- paste0("fold_", 1:num_folds)
  colnames(best_num_bases) <- colnames(penalty_grid)
  rownames(MSE_ncomp_fold) <- paste0("ncomp_", 1:ncomp)

  if (stripped) {
    ret <- list(
      CVEs_ncomp = CVEs_ncomp,
      MSE_ncomp_fold = MSE_ncomp_fold,
      best_num_bases = best_num_bases,
      elapsed = tictoc::toc(quiet = !verbose)
    )
  }else {

    if (verbose) {
      cat("Fitting final model\n")
    }


    basisobj_X_best <- fda_basis_func_X(rangeval = range(argvals_X),
                                        nbasis = best_num_bases[ncomp, "numbases_X"])

    basisobj_Y_best <- fda_basis_func_Y(rangeval = range(argvals_Y),
                                        nbasis = best_num_bases[ncomp, "numbases_Y"])

    final_model <- ffpls_bs(X = X,
                            Y = Y,
                            argvals_X = argvals_X,
                            argvals_Y = argvals_Y,
                            ncomp = ncomp,
                            center = center,
                            basisobj_X = basisobj_X_best,
                            basisobj_Y = basisobj_Y_best,
                            penalty_X = penalty_X,
                            penalty_Y = penalty_Y,
                            verbose = FALSE,
                            stripped = stripped
                            ...)

    ret <- list(
      CVEs_ncomp = CVEs_ncomp,
      MSE_ncomp_fold = MSE_ncomp_fold,
      best_num_bases = best_num_bases,
      final_model = final_model,
      elapsed = tictoc::toc(quiet = !verbose)
    )

  }

  class(ret) <- "cv_fofr"

  return(ret)
}

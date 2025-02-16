#' Cross-Validation for Functional Partial Least Squares with Basis Selection
#'
#' Performs cross-validation to select the optimal number of basis functions for functional
#' predictors and responses in a functional partial least squares (FPLS) model.
#'
#' @param X Matrix of predictors. Rows represent observations, and columns represent variables.
#' @param Y Matrix of responses. Rows represent observations, and columns represent variables.
#' @param center Logical; if `TRUE`, the data is mean-centered before fitting the model. Default is `TRUE`.
#' @param argvals_X Vector of argument values corresponding to the functional predictors. Default is `NULL`.
#' @param argvals_Y Vector of argument values corresponding to the functional responses. Default is `NULL`.
#' @param num_bases_X Integer or vector specifying the number of basis functions for the predictors. Default is 20.
#' @param num_bases_Y Integer or vector specifying the number of basis functions for the responses. Default is 20.
#' @param fda_basis_func_X Function to create the basis functions for the predictors. Defaults to `fda::create.bspline.basis`.
#' @param fda_basis_func_Y Function to create the basis functions for the responses. Defaults to `fda::create.bspline.basis`.
#' @param penalty_X Numeric penalty for smoothing the predictor basis functions. Default is 0.
#' @param penalty_Y Numeric penalty for smoothing the response basis functions. Default is 0.
#' @param ncomp Integer specifying the maximum number of components to compute. Default is `min(10, ncol(X))`.
#' @param folds Integer or list specifying the cross-validation folds. If an integer, it indicates the number of folds. If a list, it provides predefined fold indices. Default is 5.
#' @param verbose Logical; if `TRUE`, progress messages are displayed during cross-validation. Default is `TRUE`.
#' @param stripped Logical; if `TRUE`, the final model is not included in the output, only cross-validation results. Default is `TRUE`.
#' @param RPhi Penalty matrix for the predictors. Default is `NULL`.
#' @param RPsi Penalty matrix for the responses. Default is `NULL`.
#' @param PX Penalty matrix for the predictors. Default is `NULL`.
#' @param PY Penalty matrix for the responses. Default is `NULL`.
#' @param harmaccelLfd_X Object of class `fd` representing the harmonic acceleration operator for the predictors. Default is `NULL`.
#' @param harmaccelLfd_Y Object of class `fd` representing the harmonic acceleration operator for the responses. Default is `NULL`.
#' @param ... Additional arguments passed to the `ffpls_bs` function.
#'
#' @return A list containing:
#' \describe{
#'   \item{CVEs_ncomp}{Vector of cross-validation errors for each number of components.}
#'   \item{MSE_ncomp_fold}{Matrix of mean squared errors (MSE) for each component and fold.}
#'   \item{best_num_bases}{Matrix indicating the best number of basis functions for predictors and responses for each component.}
#'   \item{final_model}{The final fitted FPLS model (only included if `stripped = FALSE`).}
#'   \item{elapsed}{Elapsed time for the cross-validation process.}
#' }
#'
#' @export
#'
#' @importFrom foreach %dopar%
#' @importFrom foreach %:%
#'
#' @examples
#' # 1D example:
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
                             RPhi = NULL,
                             RPsi = NULL,
                             PX = NULL,
                             PY = NULL,
                             harmaccelLfd_X = NULL,
                             harmaccelLfd_Y = NULL,
                             ...) {


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

  # Initialize computation times:
  cve_times_ncomp <- array(data = NA, dim = ncomp)
  names(cve_times_ncomp) <- paste0("ncomp_", 1:ncomp)

  # Initialize MSEs per fold:
  MSE_ncomp_fold <- matrix(data = NA,
                           nrow = ncomp,
                           ncol = num_folds) # MSE per component per fold
  colnames(MSE_ncomp_fold) <- paste0("fold_", 1:num_folds)
  rownames(MSE_ncomp_fold) <- paste0("ncomp_", 1:ncomp)

  # Best number of bases (unique) for each ncomp_i:
  best_num_bases <- matrix(data = NA,
                           nrow = ncomp,
                           ncol = 2)
  rownames(best_num_bases) <- paste0("ncomp_", 1:ncomp)
  colnames(best_num_bases) <- colnames(num_bases_grid)


  for (ncomp_i in 1:ncomp) {


    tictoc::tic(paste0("Crossvalidation component # ", ncomp_i))


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

        # MSE_lambda_fold <- matrix(NA, nrow = nrow(num_bases_grid), ncol = num_folds)
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

        if ((basisobj_X$type == "fourier" ) & (is.null(harmaccelLfd_X)) ) {

          return("Provide the harmonic accelaeration operator for X")

        }else if ((basisobj_X$type == "fourier" ) & (!is.null(harmaccelLfd_X)) ) {

          RPhi <- fda::fourierpen(basisobj = basisobj_X, Lfdobj = 0)

          #  compute the penalty matrix R
          PX = fda::eval.penalty(basisobj_X, harmaccelLfd_X)
        }

        if ((basisobj_Y$type == "fourier" ) & (is.null(harmaccelLfd_Y))) {

          return("Provide the harmonic accelaeration operator for Y")

        }else if ((basisobj_Y$type == "fourier" ) & (!is.null(harmaccelLfd_Y)) ) {

          RPsi <- fda::fourierpen(basisobj = basisobj_Y, Lfdobj = 0)

          #  compute the penalty matrix R
          PY = fda::eval.penalty(basisobj_Y, harmaccelLfd_Y)
        }

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
                             stripped = stripped,
                             RPhi = RPhi,
                             RPsi = RPsi,
                             PX = PX,
                             PY = PY,
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
    best_num_bases[ncomp_i, ] <- as.numeric(num_bases_grid[sel_num_bases, ])

    # Save the folds-averaged CV error:
    CVEs_ncomp[ncomp_i] <- CVEs_ncomp_bases[sel_num_bases]

    # Save MSEs per fold, for the best lambda:
    MSE_ncomp_fold[ncomp_i, ] <- MSE_lambda_fold[sel_num_bases, ]

    # Save times
    cve_times_ncomp[ncomp_i] <- tictoc::toc(quiet = !verbose)


  } # loop in ncomp: number of components


  # Transform to accumulated times:
  cve_times_ncomp <- cumsum(cve_times_ncomp)



  if (stripped) {
    ret <- list(
      CVEs_ncomp = CVEs_ncomp,
      MSE_ncomp_fold = MSE_ncomp_fold,
      best_num_bases = best_num_bases,
      elapsed = cve_times_ncomp
    )
  }else {

    if (verbose) {
      cat("Fitting final model\n")
    }


    basisobj_X_best <- fda_basis_func_X(rangeval = range(argvals_X),
                                        nbasis = best_num_bases[ncomp, "numbases_X"])

    basisobj_Y_best <- fda_basis_func_Y(rangeval = range(argvals_Y),
                                        nbasis = best_num_bases[ncomp, "numbases_Y"])


    if ((basisobj_X_best$type == "fourier" ) & (is.null(harmaccelLfd_X)) ) {

      return("Provide the harmonic accelaeration operator for X")

    }else if ((basisobj_X_best$type == "fourier" ) & (!is.null(harmaccelLfd_X)) ) {

      RPhi <- fda::fourierpen(basisobj = basisobj_X_best, Lfdobj = 0)

      #  compute the penalty matrix R
      PX = fda::eval.penalty(basisobj_X_best, harmaccelLfd_X)
    }

    if ((basisobj_Y_best$type == "fourier" ) & (is.null(harmaccelLfd_Y))) {

      return("Provide the harmonic accelaeration operator for Y")

    }else if ((basisobj_Y_best$type == "fourier" ) & (!is.null(harmaccelLfd_Y)) ) {

      RPsi <- fda::fourierpen(basisobj = basisobj_Y_best, Lfdobj = 0)

      #  compute the penalty matrix R
      PY = fda::eval.penalty(basisobj_Y_best, harmaccelLfd_Y)
    }


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
                            stripped = stripped,
                            RPhi = RPhi,
                            RPsi = RPsi,
                            PX = PX,
                            PY = PY,
                            ...)

    ret <- list(
      CVEs_ncomp = CVEs_ncomp,
      MSE_ncomp_fold = MSE_ncomp_fold,
      best_num_bases = best_num_bases,
      final_model = final_model,
      elapsed = cve_times_ncomp
    )

  }

  class(ret) <- "cv_fofr"

  return(ret)
}

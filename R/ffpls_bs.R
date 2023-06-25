#' Penalized functional PLS based on the basis representation of the data
#'
#' @param X a number-of-observations times  nodes-of-X matrix.
#' @param Y a number-of-observations times nodes-of-Y matrix.
#' @param center logical, indicating if data should be centered.
#' @param argvals_X a set of argument values for X.
#' @param argvals_Y a set of argument values for Y.
#' @param basisobj_X basis object for X.
#' @param basisobj_Y basis object for Y.
#' @param penalty_X penalty for the weights of X.
#' @param penalty_Y penalty for the weights of Y.
#' @param ncomp number of components, integer.
#' @param verbose logical, indicating if messages should be printed.
#' @param stripped logical.  If \code{TRUE} the calculations are stripped as
#' much as possible for speed. Particularly, if \code{FALSE} (default) it computes
#' the final models using the best combination of penalties.
#' Inspired by package \code{pls}.
#' @param RPhi matrix of inner products for the basis of X.
#' @param RPsi matrix of inner products for the basis of Y.
#' @param PX penalty matrix for the basis of X.
#' @param PY penalty matrix for the basis of Y.
#' @param ... extra arguments for \code{pls::plsr()} function.
#'
#' @return a fof_pls model.
#' @export
#'
#' @examples
#' # 1D example:
#' # library(fda)
ffpls_bs <- function(X,
                     Y,
                     center = TRUE,
                     argvals_X,
                     argvals_Y,
                     ncomp = 3,
                     basisobj_X,
                     basisobj_Y,
                     penalty_X = 0,
                     penalty_Y = 0,
                     verbose = TRUE,
                     stripped = FALSE,
                     RPhi = NULL,
                     RPsi = NULL,
                     PX = NULL,
                     PY = NULL,
                     ...
) {

  # tictoc::tic("fof_PLS")

  if (center) {

    # center:
    Xc <- scale(X, scale = F)
    Yc <- scale(Y, scale = F)

    # get mean of Y and X to add in the end:
    Y_mean <- attr(Yc, "scaled:center")
    X_mean <- attr(Xc, "scaled:center")


  }else {

    Xc <- X
    Yc <- Y
    X_mean <- rep(0, ncol(X))
    Y_mean <- rep(0, ncol(Y))

  } # center data

  # number of nodes and samples:
  n_samp <- nrow(Xc)
  n_nodes_X <- ncol(Xc)
  n_nodes_Y <- ncol(Yc)


  # Evaluate basis functions:
  # (rows corresponding to argument values and columns to basis functions)

  # Phi in the paper:
  Xbasis_eval <- fda::eval.basis(evalarg = argvals_X,
                                 basisobj = basisobj_X,
                                 Lfdobj=0, returnMatrix=TRUE)
  # tXbasis_eval <- Matrix::t(Xbasis_eval)

  # Psi in the paper:
  Ybasis_eval <- fda::eval.basis(evalarg = argvals_Y,
                                 basisobj = basisobj_Y,
                                 Lfdobj=0, returnMatrix=TRUE)
  # tYbasis_eval <- Matrix::t(Ybasis_eval)




  # Matrix of inner products (mass):
  if (is.null(RPhi)) {
    RPhi <- fda::bsplinepen(basisobj = basisobj_X, Lfdobj = 0)
  }
  inv_RPhi <- Matrix::solve(RPhi)


  if (is.null(RPsi)) {
    RPsi <- fda::bsplinepen(basisobj = basisobj_Y, Lfdobj = 0)
  }
  inv_RPsi <- Matrix::solve(RPsi)

  # Penalty matrix:
  dorder <- 2

  # Penalty X(t)
  if (is.null(PX)) {
    deltaX <- diff(diag(basisobj_X$nbasis),
                   differences = dorder) #d-order differences
    PX <- Matrix::t(deltaX) %*% deltaX
  }


  # Penalty Y(s)
  if (is.null(PY)) {
    deltaY <- diff(diag(basisobj_Y$nbasis),
                   differences = dorder) #d-order differences
    PY <- Matrix::t(deltaY) %*% deltaY
  }


  # Matrices LX and LY with penalty:
  # for X(t)
  LLXprim <- RPhi + penalty_X*PX
  LX <- expm::sqrtm(LLXprim)
  inv_LX <- Matrix::solve(LX)
  # for Y(s)
  LLYprim <- RPsi + penalty_Y*PY
  LY <- expm::sqrtm(LLYprim)
  inv_LY <- Matrix::solve(LY)

  # Represent X and Y using a B-spline basis:
  basis_repr_X <- fda::Data2fd(argvals = argvals_X,
                               y = Matrix::t(Xc),
                               basisobj = basisobj_X  )

  basis_repr_Y <- fda::Data2fd(argvals = argvals_Y,
                               y = Matrix::t(Yc),
                               basisobj = basisobj_Y  )

  # Coeff. matrix of size obs. times num. basis:
  A <- Matrix::t(basis_repr_X[["coefs"]]) # Alpha coeffs.
  G <- Matrix::t(basis_repr_Y[["coefs"]]) # Gamma coeffs.


  # left_PLS <- G %*% RPsi %*% Matrix::t(inv_LY)  # Pi matrix
  # right_PLS <- A %*% RPhi %*% Matrix::t(inv_LX) # Delta matrix

  left_PLS <- G %*% RPsi %*% inv_LY  # Pi matrix
  right_PLS <- A %*% RPhi %*% inv_LX # Delta matrix


  # PLS model:
  mvpls_model <- pls::plsr(left_PLS ~ right_PLS,
                           ncomp =  ncomp,
                           method = "oscorespls",
                           center = TRUE,
                           scale = FALSE,
                           ...)



  # Rertuns:
  if (stripped) { # fast return (for CV)

    ret <- list(argvals_X = argvals_X,
                # argvals_Y = argvals_Y,
                basisobj_X = basisobj_X,
                # basisobj_Y = basisobj_Y,
                mvpls_model = mvpls_model,
                RPhi = RPhi,
                LY = LY,
                Ybasis_eval = Ybasis_eval,
                inv_RPsi = inv_RPsi,
                inv_LX = inv_LX,
                G = G,
                ncomp = ncomp,
                X_mean = X_mean,
                Y_mean = Y_mean,
                elapsed = tictoc::toc(quiet = !verbose)
    )

    class(ret) <- "ffpls_bs"

  }else {         # full computations

    # Get MV-model components:
    # Notes:
    # Pi = left_PLS ~~ mvpls_model[["Yscores"]] %*% t(mvpls_model[["Yloadings"]])
    # Delta = right_PLS ~~ mvpls_model[["scores"]] %*% t(mvpls_model[["loadings"]])

    TT <- as.matrix(mvpls_model[["scores"]])

    V <-  as.matrix(mvpls_model[["Yscores"]]) # scores Y

    C <- as.matrix(mvpls_model[["loadings"]] )# regress. coeffs. Pi matrix

    D <-  as.matrix(mvpls_model[["Yloadings"]]) # regress. coeffs. Delta matrix

    W <- as.matrix( mvpls_model[["loading.weights"]] ) # weights Omega tilde


    # Coefficient function:

    Beta_hat <-  array(NA, dim = c(n_nodes_Y, n_nodes_X, ncomp))
    # Beta_hat <-  array(NA, dim = c(n_nodes_X, n_nodes_Y, ncomp))

    Yc_hat <- array(NA, dim = c(n_samp, n_nodes_Y, ncomp))

    for (h in 1:ncomp) {

      ## Coeffs. for the multivariate PLS:

      Theta <- as.matrix(mvpls_model[["coefficients"]][ , , h] )

      #############################################################################
      ## Coeffs. of the basis representation of Beta(p):
      B <- inv_LX %*% Theta %*% LY %*% inv_RPsi

      # Beta_hat observed in p_1, ..., p_mx and q_1, ..., q_my.
      Beta_hat[ , , h] <- Ybasis_eval %*% Matrix::t(B) %*% Matrix::t(Xbasis_eval )
      # Beta_hat[ , , h] <- Xbasis_eval %*% B %*% Matrix::t( Ybasis_eval )
      #############################################################################



      # Estimated Y(q_1, ..., q_my)
      hat_left_PLS <- mvpls_model$fitted.values[, , h]
      Yc_hat[, , h] <- hat_left_PLS %*% LY %*% inv_RPsi %*% Matrix::t(Ybasis_eval)

    }




    ret <- list(argvals_X = argvals_X,
                argvals_Y = argvals_Y,
                basisobj_X = basisobj_X,
                basisobj_Y = basisobj_Y,
                RPhi = RPhi,
                RPsi = RPsi,
                inv_RPsi = inv_RPsi,
                inv_RPhi = inv_RPhi,
                inv_LY = inv_LY,
                inv_LX = inv_LX,
                LX = LX,
                LY = LY,
                G = G,
                Ybasis_eval = Ybasis_eval,
                Xbasis_eval = Xbasis_eval,
                mvpls_model = mvpls_model,
                ncomp = ncomp,
                V = V,
                TT = TT,
                C = C,
                W = W,
                D = D,
                X_mean = X_mean,
                Y_mean = Y_mean,
                fitted.values = Yc_hat + Y_mean,
                coefficient_function = Beta_hat,
                elapsed = tictoc::toc(quiet = !verbose)
    )

    class(ret) <- "ffpls_bs"

  }

  return(ret)


} # end: FPLS function

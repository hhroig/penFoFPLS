#' Predict new data for FFoFPLS-FDA model (penalized FPLS algorithm using an fda basis,
#' as in Aguilera et al. 2023)
#'
#' @param object an \code{ffpls_bs} model.
#' @param newdata new data matrix.
#' @param mode if \code{response} (default) it returns the predicted Y curves. If
#'             \code{coefficients} it returns the basis coefficients for the predicted
#'             curves.
#' @param ... further arguments.  Currently not used
#'
#' @return predicted values
#' @export
#'
#' @examples
#' # 1D example:
predict.ffpls_bs <-  function(object, newdata, mode = "response", ...){

  newXc <- scale(newdata, center = object$X_mean, scale = FALSE)

  # Represent X_test and Y_test using a B-spline basis:
  bsplineXtest <- fda::Data2fd(argvals = object$argvals_X,
                               y = Matrix::t(newXc),
                               basisobj = object$basisobj_X  )

  # Coeff. matrix of size obs. times num. basis:
  A_test <- Matrix::t(bsplineXtest[["coefs"]]) # Alpha coeffs.

  # Input for mv-pls:
  new_right_PLS <- A_test %*% object$RPhi %*% object$inv_LX # Delta (test) matrix

  if (mode == "response") { # Estimated Y

    Yc_test_hat <- array(NA, dim = c(nrow(newdata), nrow(object$Ybasis_eval), object$ncomp))

    for (h in 1:object$ncomp) {
      # Estimated Y(q_1, ..., q_my)
      hat_left_PLS <- stats::predict(object$mvpls_model, new_right_PLS)[, , h]
      Yc_test_hat[, , h] <- hat_left_PLS %*% object$LY %*% object$inv_RPsi %*% Matrix::t(object$Ybasis_eval)
    }
    pred <- Yc_test_hat + object$Y_mean


  }else if (mode == "coefficients") { # Estimated basis coefficients for Y

    G_hat <- array(NA, dim = c(nrow(newdata), ncol(object$inv_RPsi), object$ncomp))

    for (h in 1:object$ncomp) {
      # Estimated Y(q_1, ..., q_my)
      hat_left_PLS <- stats::predict(object$mvpls_model, new_right_PLS)[, , h]
      G_hat[, , h] <- hat_left_PLS %*% object$LY %*% object$inv_RPsi
    }
    pred <- G_hat

  }



  return(pred)

}

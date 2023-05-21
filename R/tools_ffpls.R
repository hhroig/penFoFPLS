# 1D numerical integration
#
num_int_1d <- function(argvals, f_obs) {

  f <- stats::splinefun(x = argvals, y = f_obs)
  # f <- stats::approxfun(x = argvals, y = f_obs)

  subds <- 100

  res <- "no_convergence"

  count_iter <- 0

  while (!inherits(res, "numeric") && count_iter < 20) {

    res <- try( stats::integrate(f,
                                 min(argvals),
                                 max(argvals),
                                 subdivisions = 10*subds,
                                 stop.on.error = FALSE)$value, silent = TRUE)

    count_iter <- count_iter + 1

  }



  if (!inherits(res, "numeric")) {

    message("Integral did not converge with 10^", count_iter+2, "subdivisions.\n")
    message(res)

  }

  # stopifnot( inherits(res, "numeric") )

  return(res)


}





# Integrated mean square error for the evaluation of the estimated beta(q,p) FoF regression.
imse_beta_ffpls_bs <- function(beta_true,
                               beta_hat,
                               argvals_X,
                               argvals_Y){

  D2_area <- (max(argvals_Y) - min(argvals_Y))*(max(argvals_X) - min(argvals_X))

  beta_diff_sq <- (beta_true - beta_hat)^2

  int_beta_in_q <- array(NA, dim = length(argvals_Y))

  for (q_i in 1:length(argvals_Y)) {

    int_beta_in_q[q_i] <- num_int_1d(argvals = argvals_X,
                                     f_obs = beta_diff_sq[q_i, ]  )
  }



  imse <- num_int_1d(argvals = argvals_Y,
                     f_obs = int_beta_in_q  )

  imse <- imse/D2_area

  return(imse)

}

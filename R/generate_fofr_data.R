#' Generate artificial data for function-on-function regression
#'
#' @param nbasisX number of basis Bspline basis to use for X
#' @param nbasisY number of basis Bspline basis to use for Y
#' @param nbeta number identifying which beta(q, p) to generate
#' @param nnodesX number of nodes to use in the generation of X
#' @param nnodesY number of nodes to use in the generation of Y
#' @param Rsq R^2 between the generated Y(t) and perturbed Y(t) + error.
#' @param n sample size.
#' @param n_validation additional data for validation purposes.
#'
#' @return a list with the generate X, Y, and true beta
#' @export
#'
#' @examples
#' generate_fofr_data(nbasisX = 7, nbasisY = 5, nbeta = 1,
#'                    nnodesX = 100, nnodesY = 100)
generate_fofr_data <- function(nbasisX = 7, nbasisY = 5, nbeta = 1,
                               nnodesX = 99, nnodesY = 98, Rsq = 0.9,
                               n = 100, n_validation = 50) {

  # number of basis to use:
  K <- nbasisX
  L <- nbasisY

  # intervals are defined in [0,1]:
  TX <- 1
  TY <- 1


  # Cubic B spline basis :
  phiX <- fda::create.bspline.basis(c(0, TX), nbasis=K)
  phiY <- fda::create.bspline.basis(c(0, TY), nbasis=L)


  # Generate X:
  n <- n + n_validation

  # random coefficients for X:
  alpha <- matrix(stats::runif(n*K, -1, 1), nrow=n)

  X <- fda::fd(t(alpha), phiX)


  # define Beta:
  p <- seq(0, TX, length.out = nnodesX) #intervalle [0,TX]
  q <- seq(0, TY, length.out = nnodesY) #intervalle [0,TY]

  if (nbeta == 1) {

    # As in Preda, C., & Schiltz, J. (2011) http://orbilu.uni.lu/bitstream/10993/5912/2/2457_001.pdf
    f <- function(q, p) {
      return((q-p)^2)
    }

  }else if (nbeta == 2) {

    # single exponential top right corner:
    f <- function(q, p){
      return(5*exp(-((p - 0.75)^2 + (q - 0.75)^2)/( 2*0.2^2 ))    )
    }

  }else if (nbeta == 3) {

    # monkey saddle:
    f <- function(q, p){
      return(((p*4 - 2)^3 - 3*(p*4 - 2)*((q*4 - 2)^2))  )
    }

  }else if (nbeta == 4) {

    # double exponential top right and bottom left:
    f <- function(q, p){
      return( 5*exp(-((p - 0.75)^2 + (q - 0.75)^2)/( 2*0.25^2 )) +
                5*exp(-((p - 0.1)^2 + (q - 0.1)^2)/( 2*0.25^2 )) )
    }


  }


  beta <- outer(q, p, f)

  beta_fd <- fda::Data2fd(t(beta), argvals = p, phiX )

  # Center data X(p):
  Xc_fd <- fda::center.fd(X)
  Xmean_fd <- fda::mean.fd(X)


  # Generate Y clean (no noise):
  Y_clean <- fda::inprod(Xc_fd, beta_fd)

  # Noisy Y , approx. R^2(t) ~~ Rsq:
  Y <- Y_clean

  for (q_ind in 1:length(q)) {

    # variance of residuals at point "q[q_ind]"
    var_e <- (1/Rsq - 1)*stats::var(Y_clean[, q_ind])

    # Noisy Y:
    Y[ , q_ind] <- Y_clean[ , q_ind] +
      as.matrix(stats::rnorm(n = nrow(Y_clean), mean = 0, sd = sqrt(var_e)))

  }

  # Savings:

  X <- t( fda::eval.fd(evalarg = p, fdobj = X) )
  beta_true <- beta
  Xc <- t( fda::eval.fd(evalarg = p, fdobj = Xc_fd) )

  n_traintest = n - n_validation

  out <- list(Xc = Xc[1:n_traintest, ],
              X = X[1:n_traintest, ],
              Y = Y[1:n_traintest, ],
              beta_true = beta_true,
              Xc_val = Xc[(n_traintest+1):nrow(Xc), ],
              X_val = X[(n_traintest+1):nrow(X), ],
              Y_val = Y[(n_traintest+1):nrow(X), ]
              )

  return(out)

}


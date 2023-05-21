#' Generate artificial data for function-on-function regression
#'
#' @param nbasisX number of basis Bspline basis to use for X
#' @param nbasisY number of basis Bspline basis to use for Y
#' @param nbeta number identifying which beta(q, p) to generate
#' @param nnodesX number of nodes to use in the generation of X
#' @param nnodesY number of nodes to use in the generation of Y
#'
#' @return a list with the generate X, Y, and true beta
#' @export
#'
#' @examples
#' # To reproduce the paper by Preda & Schiltz (2011):
#' generate_fofr_data(nbasisX = 7, nbasisY = 5, nbeta = 1,
#'                    nnodesX = 100, nnodesY = 100)
generate_fofr_data <- function(nbasisX = 7, nbasisY = 5, nbeta = 1,
                               nnodesX = 99, nnodesY = 98) {

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
  n <- 100

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


  # Generate Y:
  Y <- fda::inprod(Xc_fd, beta_fd)

  varthY <- rep(0, length(q))
  somme <- 0

  for(k in 1:length(q)) {

    somme <- 0
    coef_phi <- rep(0,K)

    for(i in 1:K) {
      coef_phi[i] <- 1
      somme <- somme + fda::inprod(  fda::fd(coef_phi, phiX), beta_fd[k]  )^{2}
      coef_phi[i] <- 0
    }

    varthY[k] <- (1/3)*somme
  }


  # add error to Y(q)
  tauxbruit <- 0.1
  bruit <- rep(0,length(q))

  Yb = Y

  for (i in 1:n) {

    for (k in 1:length(q)) {
      bruit[k] <- stats::rnorm(1 , 0, sqrt(tauxbruit*varthY[k]))
    }

    Yb[i, ] <- Y[i, ] + bruit
  }

  # Savings:

  X <- t( fda::eval.fd(evalarg = p, fdobj = X) )
  beta_true <- beta
  Xc <- t( fda::eval.fd(evalarg = p, fdobj = Xc_fd) )

  out <- list(Xc = Xc, X = X, Y = Y, beta_true = beta_true)

  return(out)

}


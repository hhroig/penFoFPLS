% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_ffpls.R
\name{predict.ffpls_bs}
\alias{predict.ffpls_bs}
\title{Predict new data for FFoFPLS-FDA model (penalized FPLS algorithm using an fda basis,
as in Aguilera et al. 2023)}
\usage{
\method{predict}{ffpls_bs}(object, newdata, mode = "response", ...)
}
\arguments{
\item{object}{an \code{ffpls_bs} model.}

\item{newdata}{new data matrix.}

\item{mode}{if \code{response} (default) it returns the predicted Y curves. If
\code{coefficients} it returns the basis coefficients for the predicted
curves.}

\item{...}{further arguments.  Currently not used}
}
\value{
predicted values
}
\description{
Predict new data for FFoFPLS-FDA model (penalized FPLS algorithm using an fda basis,
as in Aguilera et al. 2023)
}
\examples{
# 1D example:
}

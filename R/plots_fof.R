#' Plot Methods for cv_sofr objects
#'
#' @details The function is simply a wrapper for the underlying plot functions used to
#' make the selected plots.  See \code{box_plot} and \code{means_plot}
#' for details.
#'
#' @param x an object of class \code{cv_sofr} resulting from a cross-validation procedure.
#' @param plot_type character.  What kind of plot to plot.
#' @param ... additional parameters that depends on the \code{plot_type}.
#'
#' @return returns whatever the underlying plot function.
#'
#' @export
#'
#' @examples
#' # 1D example:
#' \dontrun{
#' }
plot.cv_fofr <- function(x, plot_type = c("means", "boxplot"), ...)
{
  plot_type <- match.arg(plot_type)
  plotFunc <- switch(plot_type,
                     means = means_plot,
                     boxplot = box_plot )
  plotFunc(x, ...)
}


#' Validation plot (MSEP/RMSEP/AUC)
#'
#' @param object results from a cross-validation procedure.
#' @param val_type character.  Any of "RMSE", "MSE", or "AUC".
#'
#' @export
#'
#'@importFrom ggplot2 ggplot aes geom_boxplot xlab ylab theme_bw
#'
#' @examples
##' # 1D example:

#' \dontrun{
#'
#' }
box_plot <- function(object, val_type = c("MSE", "RMSE") )
{

  MSE <- ncomp <- NULL

  val_type <- match.arg(val_type)

  if (val_type == "MSE") {

    df <- t(object[["MSE_ncomp_fold"]]) %>%
      dplyr::as_tibble() %>%
      tidyr::pivot_longer(cols = rownames(object[["MSE_ncomp_fold"]]),
                          names_to = "ncomp",
                          values_to = "MSE") %>%
      dplyr::mutate(as.factor(ncomp))

    p <- ggplot(data = df, aes(x = ncomp, y = MSE)) +
      geom_boxplot() +
      xlab("number of components") + ylab("MSE") +
      theme_bw()

  }else if (val_type == "RMSE") {

    df <- t(object[["MSE_ncomp_fold"]]) %>%
      dplyr::as_tibble() %>%
      tidyr::pivot_longer(cols = rownames(object[["MSE_ncomp_fold"]]),
                          names_to = "ncomp",
                          values_to = "MSE") %>%
      dplyr::mutate(as.factor(ncomp))

    p <- ggplot(data = df, aes(x = ncomp, y = sqrt(MSE))) +
      geom_boxplot() +
      xlab("number of components") + ylab("RMSE") +
      theme_bw()

  }

  return(p)


}


#' Mean validation plot (mean of MSEP/RMSEP/AUC)
#'
#' @param object results from a cross-validation procedure.
#' @param val_type character.  Any of "RMSE", "MSE", or "AUC".
#'
#' @export
#'
#'@importFrom ggplot2 ggplot aes geom_line geom_point xlab ylab theme_bw
#'
#' @examples
##' # 1D example:
#'
#' \dontrun{
#' }
means_plot <- function(object, val_type = c("MSE", "RMSE") )
{

  MSE <- ncomp <- NULL

  val_type <- match.arg(val_type)

  if (val_type == "MSE") {

    df <- t(object[["MSE_ncomp_fold"]]) %>%
      dplyr::as_tibble() %>%
      tidyr::pivot_longer(cols = rownames(object[["MSE_ncomp_fold"]]),
                          names_to = "ncomp",
                          values_to = "MSE") %>%
      dplyr::mutate(ncomp = as.numeric(as.factor(ncomp))) %>%
      dplyr::group_by(ncomp) %>%
      dplyr::summarise(MSE = mean(MSE))

    p <- ggplot(data = df, aes(x = ncomp, y = MSE)) +
      geom_point(size = 2) +
      geom_line() +
      xlab("number of components") + ylab("MSE") +
      theme_bw()

  }else if (val_type == "RMSE") {

    df <- t(object[["MSE_ncomp_fold"]]) %>%
      dplyr::as_tibble() %>%
      tidyr::pivot_longer(cols = rownames(object[["MSE_ncomp_fold"]]),
                          names_to = "ncomp",
                          values_to = "MSE") %>%
      dplyr::mutate(ncomp = as.numeric(as.factor(ncomp))) %>%
      dplyr::group_by(ncomp) %>%
      dplyr::summarise(MSE = mean(MSE))

    p <- ggplot(data = df, aes(x = ncomp, y = sqrt(MSE))) +
      geom_point(size = 2) +
      geom_line() +
      xlab("number of components") + ylab("RMSE") +
      theme_bw()

  }

  return(p)


}

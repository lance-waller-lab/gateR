#' Calculate p-value corrections
#'
#' Internal function to calculate various p-value corrections including a correlated and uncorrelated Bonferroni correction for use within the \code{\link{rrs}} and \code{\link{lotrrs}} function.
#'
#' @param input An object of class 'rrs' from the \code{\link{rrs}} or \code{\link{lotrrs}} function.
#' @param alpha Numeric. The two-tailed alpha level for significance threshold (default in \code{\link{rrs}} and \code{\link{lotrrs}} functions is 0.05).
#' @param nbc Integer. The number of bins. Similar to \code{nbclass} argument in \code{\link[pgirmess]{correlog}} function. The default is the average number of gridded knots in one-dimension (i.e., x-axis). 
#' 
#' @details This function provides functionality for multiple testing correction in two ways:
#' 
#' \enumerate{
#' \item Computes a conventional Bonferroni correction ("uncorrelated") by dividing the \code{alpha} level by the number of gridded knots across the estimated surface. The default in the \code{\link[sparr]{risk}} function is a resolution of 128 x 128 or n = 16,384 knots. 
#' \item Computes a correlated Bonferroni correction ("correlated") by taking in account the spatial correlation of the relative risk surface values (if using the \code{rrs} function for a single condition gate) or the ratio of relative risk surfaces values (if using the \code{lotrrs} function for a two condition gate). The \code{alpha} level is divided by the minimum number of knots that are not spatially correlated. The minimum number of knots that are not spatially correlated is computed by counting the knots that are a distance apart that exceeds the minimum distance of non-significant spatial correlation based on a correlogram using the \code{\link[pgirmess]{correlog}} function. 
#' }
#'
#' @return An object of class 'list'. This is a named list with the following components:
#'
#' \describe{
#' \item{\code{uncorrected}}{Numeric. Returns the uncorrected p-value.}
#' \item{\code{correlated}}{Numeric. Returns the correlated Bonferroni corrected p-value.}
#' \item{\code{uncorrelated}}{Numeric. Returns the uncorrelated Bonferroni corrected p-value.}
#' }
#'
#' @importFrom pgirmess correlog
#' @export
#'
#' @keywords internal
#' 
pval_correct <- function(input,
                         alpha = 0.05,
                         nbc = NULL) {

# Inputs
  if (class(input) != "im") {
    stop("The 'input' argument must be of class 'im' from an 'rrs' object.")
  }

# Calculate correlated Bonferroni correction
## Prepare data
### Coordinates of grid points within input 'im'
  rx <- rep(input$xcol, length(input$yrow))
  for(i in 1:length(input$yrow)) {
    if (i == 1) { ry <- rep(input$yrow[i], length(input$xcol)) }
    if (i != 1) { ry <- c(ry, rep(input$yrow[i], length(input$xcol))) }
  }
### Create data frame of values to estimate spatial correlation
  out <- data.frame("x" = rx,
                    "y" = ry,
                    "v" = as.vector(t(input$v)))
  out$v <- ifelse(is.infinite(out$v), NA, out$v) # omit inifite values
  out <- out[!is.na(out$v),] # omit NA values

  if(is.null(nbc)) { nbc <- round(mean(length(input$xcol), length(input$yrow)))}

## Estimate spatial correlation
  out_cor <- pgirmess::correlog(coords = out[ , 1:2],
                                z = out[ , 3],
                                method = "Moran",
                                randomisation = FALSE,
                                nbclass = nbc)

## Find shortest distance without significant correlation
  corr_eff <- out_cor[ , 1][min(which(out_cor[ , 3] >= alpha))]

## Bonferroni
  x_kern <- length(input$xcol) # number of kernels in x-axis
  y_kern <- length(input$yrow) # number of kernels in y-axis
  x_dist_kern <- input$ystep # distance between kernels in x-axis
  y_dist_kern <- input$ystep # distance between kernels in y-axis
  n_x_kern <- corr_eff / x_dist_kern # number of x-axis kernels within correlation distance
  n_y_kern <- corr_eff / y_dist_kern # number of y- axis kernels within correlation distance
  x_uncorr <- x_kern / n_x_kern # number of uncorrelated x-axis kernels
  y_uncorr <- y_kern / n_y_kern # number of uncorrelated y-axis kernels
  n_uncorr <- x_uncorr * y_uncorr # total number of uncorrelated kernels

  alpha <- alpha # uncorrected
  alpha_correlated <- alpha / n_uncorr # correlated Bonferroni
  alpha_bonferroni <- alpha / (x_kern * y_kern) # uncorrelated Bonferroni (including tests outside window)

  out_alpha <- list("uncorrected" = alpha,
                    "correlated" = alpha_correlated,
                    "uncorrelated" = alpha_bonferroni)
  return(out_alpha)
}

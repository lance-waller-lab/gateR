#' Calculate p-value corrections
#'
#' Internal function to calculate various p-value corrections for use within the \code{\link{rrs}} and \code{\link{lotrrs}} functions.
#'
#' @param input An object of class 'rrs' from the \code{\link{rrs}} or \code{\link{lotrrs}} function.
#' @param type Character string specifying which correction for multiple comparisons. Options include a False Discovery Rate \code{p_correct = "FDR"}, a spatiall dependent Sidak correction \code{p_correct = "correlated Sidak"}, a spatially dependent Bonferroni correction \code{p_correct = "correlated Bonferroni"}, an independent Sidak correction \code{p_correct = "uncorrelated Sidak"}, and an independent Bonferroni correction \code{p_correct = "uncorrelated Bonferroni"}.
#' @param alpha Numeric. The two-tailed alpha level for significance threshold (default in \code{\link{rrs}} and \code{\link{lotrrs}} functions is 0.05).
#' @param nbc Integer. The number of bins. Similar to \code{nbclass} argument in \code{\link[pgirmess]{correlog}} function. The default is the average number of gridded knots in one-dimension (i.e., x-axis). 
#' 
#' @details This function provides functionality for multiple testing correction in five ways:
#' 
#' \enumerate{
#' \item Computes a False Discovery Rate by Benjamini and Hochberg \doi{10.1111/j.2517-6161.1995.tb02031.x} (\code{p_correct = "FDR"}) by: 1) sorting the p-values (p_i) of each knot in ascending order (p_1 <= p_2 <= ... <= p_m), 2) starting from p_m find the first p_i for which p_i <= (i/m) * alpha.
#' \item Computes an independent Sidak correction \doi{10.2307/2283989} (\code{p_correct = "uncorrelated Sidak"}) by 1 - (1 - \code{alpha}) ^ (1 / total number of gridded knots across the estimated surface). The default in the \code{\link[sparr]{risk}} function is a resolution of 128 x 128 or n = 16,384 knots and a custom resolution can be specified using the \code{resolution} argument within the \code{\link[sparr]{risk}} function.
#' \item Computes an independent Bonferroni correction (\code{p_correct = "uncorrelated Bonferroni"}) by \code{alpha} / total number of gridded knots across the estimated surface. The default in the \code{\link[sparr]{risk}} function is a resolution of 128 x 128 or n = 16,384 knots and a custom resolution can be specified using the \code{resolution} argument within the \code{\link[sparr]{risk}} function.
#' \item Computes a spatially dependent Sidak correction (\code{p_correct = "correlated Sidak"}) by taking into account the spatial correlation of the relative risk surface values (if using the \code{rrs} function for a single condition gate) or the ratio of relative risk surfaces values (if using the \code{lotrrs} function for a two condition gate). The correction use the minimum number of knots that are not spatially correlated instead of the total number of knots.The minimum number of knots that are not spatially correlated is computed by counting the knots that are a distance apart that exceeds the minimum distance of non-significant spatial correlation based on a correlogram using the \code{\link[pgirmess]{correlog}} function. 
#' \item Computes a spatially dependent Bonferroni correction (\code{p_correct = "correlated Bonferroni"}) by taking into account the spatial correlation of the relative risk surface values (if using the \code{rrs} function for a single condition gate) or the ratio of relative risk surfaces values (if using the \code{lotrrs} function for a two condition gate). The correction use the minimum number of knots that are not spatially correlated instead of the total number of knots.The minimum number of knots that are not spatially correlated is computed by counting the knots that are a distance apart that exceeds the minimum distance of non-significant spatial correlation based on a correlogram using the \code{\link[pgirmess]{correlog}} function. 
#' }
#' 
#' @return An object of class 'numeric' with the corrected alpha level.
#'
#' @importFrom pgirmess correlog
#' @export
#'
#' @keywords internal
#' 
pval_correct <- function(input,
                         type = c("FDR", "correlated Sidak", "correlated Bonferroni", "uncorrelated Sidak", "uncorrelated Bonferroni"),
                         alpha = 0.05,
                         nbc = NULL) {

# False Discovery Rate
  if (type == "FDR") {
  sort_pvals <- sort(as.vector(input$P$v), decreasing = TRUE)
  
  fdr <- function(pvals, alpha) {
    m <- length(pvals)
    for (i in 1:length(pvals)) {
      if (pvals[i] <= (i/m) * alpha) { return(pvals[i]) } 
    }
  }
  
  out_alpha <- fdr(sort_pvals, alpha)
  return(out_alpha)
  }

# Correlated Bonferroni correction
  if (type == "correlated Bonferroni" | 
      type == "correlated Sidak") {
## Prepare data
### Coordinates of grid points within input 'im'
  rx <- rep(input$lrr$xcol, length(input$lrr$yrow))
  for(i in 1:length(input$lrr$yrow)) {
    if (i == 1) { ry <- rep(input$lrr$yrow[i], length(input$lrr$xcol)) }
    if (i != 1) { ry <- c(ry, rep(input$lrr$yrow[i], length(input$lrr$xcol))) }
  }
### Create data frame of values to estimate spatial correlation
  out <- data.frame("x" = rx,
                    "y" = ry,
                    "v" = as.vector(t(input$lrr$v)))
  out$v <- ifelse(is.infinite(out$v), NA, out$v) # omit inifite values
  out <- out[!is.na(out$v),] # omit NA values

  if(is.null(nbc)) { nbc <- round(mean(length(input$lrr$xcol), length(input$lrr$yrow)))}

## Estimate spatial correlation
  out_cor <- pgirmess::correlog(coords = out[ , 1:2],
                                z = out[ , 3],
                                method = "Moran",
                                randomisation = FALSE,
                                nbclass = nbc)

## Find shortest distance without significant correlation
  corr_eff <- out_cor[ , 1][min(which(out_cor[ , 3] >= alpha))]

## Bonferroni
  x_kern <- length(input$lrr$xcol) # number of kernels in x-axis
  y_kern <- length(input$lrr$yrow) # number of kernels in y-axis
  x_dist_kern <- input$lrr$xstep # distance between kernels in x-axis
  y_dist_kern <- input$lrr$ystep # distance between kernels in y-axis
  n_x_kern <- corr_eff / x_dist_kern # number of x-axis kernels within correlation distance
  n_y_kern <- corr_eff / y_dist_kern # number of y- axis kernels within correlation distance
  x_uncorr <- x_kern / n_x_kern # number of uncorrelated x-axis kernels
  y_uncorr <- y_kern / n_y_kern # number of uncorrelated y-axis kernels
  n_uncorr <- x_uncorr * y_uncorr # total number of uncorrelated kernels

  # Correlated Sidak correction
  if (type == "correlated Sidak") {
    out_alpha <- 1 - (1 - alpha) ^ (1 / n_uncorr ) 
  }
  
  # Correlated Bonferroni correction
  if (type == "correlated Bonferroni") {
    out_alpha <- alpha / n_uncorr
  }
  
  return(out_alpha)
  }
  
  # Uncorrelated Bonferroni correction
  if (type == "uncorrelated Sidak") {
    out_alpha <- 1 - (1 - alpha) ^ (1 / prod(input$lrr$dim) )
    return(out_alpha)
  }
  
  # Uncorrelated Bonferroni correction
  if (type == "uncorrelated Bonferroni") {
    out_alpha <- alpha / prod(input$lrr$dim)
    return(out_alpha)
  }
  
  
}

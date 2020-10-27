#' Prepare significant p-values for plotting
#'
#' Internal function to convert an object of class 'im' to values readable by \code{\link[fields]{image.plot}} function within the \code{\link{rrs}}, \code{\link{lotrrs}}, and \code{\link{gating}} functions.
#'
#' @param input An object of class 'rrs' from the \code{\link{rrs}} or \code{\link{lotrrs}} function.
#' @param alpha Numeric. The two-tailed alpha level for significance threshold (default in \code{\link{rrs}}, \code{\link{lotrrs}}, and \code{\link{gating}} functions is 0.05).
#'
#' @return An object of class 'raster' with categorical values:
#' 
#' \itemize{
#' \item A value of 1: Significant numerator.
#' \item A value of 2: Insignificant.
#' \item A value of 3: Significant denominator.
#' }
#'
#' @importFrom raster cut raster
#' @importFrom spatstat as.im
#' @export
#'
#' @keywords internal
#' 
pval_plot <- function(input,
                      alpha) {

  # Inputs
  if (class(input) != "im") {
    stop("The 'input' argument must be an object of class 'im'")
  }

  out <- raster::raster(spatstat::as.im(input))  # create raster
  out <- raster::cut(out,
                     breaks = c(-Inf, alpha / 2, 1 - alpha / 2, Inf),
                     right = FALSE)
  return(out)
}

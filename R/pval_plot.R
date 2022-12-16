#' Prepare significant p-values for plotting
#'
#' Internal function to convert an object of class 'im' to values readable by \code{\link[fields]{image.plot}} function within the \code{\link{rrs}}, \code{\link{lotrrs}}, and \code{\link{gating}} functions.
#'
#' @param input An object of class 'rrs' from the \code{\link{rrs}} or \code{\link{lotrrs}} function.
#' @param alpha Numeric. The two-tailed alpha level for significance threshold (default in \code{\link{rrs}}, \code{\link{lotrrs}}, and \code{\link{gating}} functions is 0.05).
#'
#' @return An object of class 'SpatRaster' with categorical values:
#' 
#' \itemize{
#' \item A value of 1: Significant numerator.
#' \item A value of 2: Insignificant.
#' \item A value of 3: Significant denominator.
#' }
#'
#' @importFrom spatstat.geom as.im
#' @importFrom terra rast values
#' @export
#'
#' @keywords internal
#' 
pval_plot <- function(input,
                      alpha) {

  # Inputs
  if (!inherits(input, "im")) {
    stop("The 'input' argument must be an object of class 'im'")
  }

  out <- terra::rast(spatstat.geom::as.im(input))  # create SpatRaster
  terra::values(out) <- cut(terra::values(out),
                            breaks = c(-Inf, alpha / 2, 1 - alpha / 2, Inf))
  return(out)
}

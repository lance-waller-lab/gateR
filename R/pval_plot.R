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
#' @importFrom sp coordinates gridded
#' @importFrom stats na.omit
#' @export
#'
#' @keywords internal
#' 
pval_plot <- function(input,
                      alpha) {

  # Inputs
  if (class(input) != "im") {
    stop("The 'input' argument must be of class 'im' from an 'rrs' object.")
  }

  # Coordinates of grid points within input 'im'
  rx <- rep(input$xcol, length(input$yrow))
  for(i in 1:length(input$yrow)) {
    if (i == 1) { ry <- rep(input$yrow[i], length(input$xcol)) }
    if (i != 1) { ry <- c(ry, rep(input$yrow[i], length(input$xcol))) }
  }

  out <- data.frame("x" = rx,
                    "y" = ry,
                    "v" = as.vector(t(input$v)))
  out$v <- ifelse(is.infinite(out$v), NA, out$v)
  out <- stats::na.omit(out) # remove NAs
  sp::coordinates(out) <- ~ x + y # convert to spatialpixelsdataframe
  suppressMessages(suppressWarnings(sp::gridded(out) <- TRUE)) # gridded
  out <- raster::raster(out)  # create raster
  out <- raster::cut(out,
                     breaks = c(-Inf, alpha / 2, 1 - alpha / 2, Inf),
                     right = FALSE)
  return(out)
}

#' A single gate for a single condition
#'
#' Estimates a relative risk surface and computes the asymptotic p-value surface for a single gate with a single condition. Includes features for basic visualization. This function is used internally within the \code{\link{gating}} function to extract the points within the significant areas.  This function can also be used as a standalone function.
#'
#' @param dat Input data frame flow cytometry data with four (4) features (columns): 1) ID, 2) Condition A ID, 3) Marker A as x-coordinate, 4) Marker B as y-coordinate.
#' @param alpha Numeric. The two-tailed alpha level for significance threshold (default is 0.05).
#' @param p_correct Character string specifying whether to apply a correction for multiple comparisons including a Bonferroni correction \code{p_correct = "uncorrelated"} or a correlated Bonferroni correction \code{p_correct = "correlated"}. If \code{p_correct = "none"} then no correction is applied. 
#' @param nbc Optional. An integer for the number of bins when \code{p_correct = "correlated"}. Similar to \code{nbclass} argument in \code{\link[pgirmess]{nbclass}}. The default is the average number of gridded knots in one-dimension (i.e., x-axis). 
#' @param doplot Logical. If \code{TRUE}, the output includes basic data visualizations.
#' @param rcols Character string of length three (3) specifying the colors for: 1) Group A, 2) Neither, and 3) Group B designations. The defaults are \code{c("#FF0000", "#cccccc", "#0000FF")} or \code{c("red", "grey80", "blue")}.
#' @param win Optional. Object of class \code{owin} for a custom two-dimensional window within which to estimate the surfaces. The default is NULL and calculates a convex hull around the data. 
#' @param verbose Logical. If \code{TRUE} will print function progress during execution. If \code{FALSE} (the default), will not print.
#' @param ... Arguments passed to \code{\link[sparr]{risk}} to select bandwidth, edge correction, and resolution.
#'
#' @details This function estimates a relative risk surface and computes the asymptotic p-value surface for a single gate and single condition using the \code{\link[sparr]{risk}} function. Bandwidth is fixed across both layers (numerator and demoninator spatial densities). Basic visualization is available if \code{doplot = TRUE}. 
#' 
#' Provides functionality for a correction for multiple testing. If \code{p_correct = "uncorrelated"}, then a conventional Bonferroni correction is calculated by dividing the \code{alpha} level by the number of gridded knots across the estimated surface. The default in the \code{\link[sparr]{risk}} function is a resolution of 128 x 128 or n = 16,384 knots and a custom resolution can be specified using the \code{resolution} argument within the \code{\link[sparr]{risk}} function. If \code{p_correct = "correlated"} (NOTE: May take a considerable amount of computation resources and time), then a Bonferroni correction that takes into account the spatial correlation of the surface is calculated within the internal \code{pval_correct} function. The \code{alpha} level is divided by the minimum number of knots that are not spatially correlated. The minimum number of knots that are not spatially correlated is computed by counting the knots that are a distance apart that exceeds the minimum distance of non-significant spatial correlation based on a correlogram using the \code{\link[pgirmess]{correlog}} function. If \code{p_correct = "none"}, then the function does not account for multiple testing and uses the uncorrected \code{alpha} level. See the internal \code{pval_correct} function documentation for more details.
#'
#' @return An object of class 'list' where each element is a object of class 'rrs' created by the \code{\link[sparr]{risk}} function with two additional components:
#' 
#' \describe{
#' \item{\code{rr}}{An object of class 'im' with the relative risk surface.}
#' \item{\code{f}}{An object of class 'im' with the spatial density of the numerator.}
#' \item{\code{g}}{An object of class 'im' with the spatial density of the denominator.}
#' \item{\code{P}}{An object of class 'im' with the asymptotic p-value surface.}
#' \item{\code{lrr}}{An object of class 'im' with the log relative risk surface.}
#' \item{\code{alpha}}{A numeric value for the alpha level used within the gate.}
#' }
#'
#' @importFrom fields image.plot
#' @importFrom graphics close.screen par screen split.screen 
#' @importFrom grDevices chull
#' @importFrom raster extent values
#' @importFrom spatstat.core owin ppp
#' @importFrom sparr OS risk
#' @export 
#'
#' @examples
#'   library(flowWorkspaceData)
#'   library(ncdfFlow)
#' 
#' # Use 'extdata' from the {flowWorkspaceData} package
#'   flowDataPath <- system.file("extdata", package = "flowWorkspaceData")
#'   fcsFiles <- list.files(pattern = "CytoTrol", flowDataPath, full = TRUE)
#'   ncfs  <- ncdfFlow::read.ncdfFlowSet(fcsFiles)
#'   fr1 <- ncfs[[1]]
#'   fr2 <- ncfs[[2]]
#' 
#' ## Comparison of two samples (single condition) "g1"
#' ## One gate (two markers) "CD4", "CD38"
#' ## Log10 Transformation for both markers
#' ## Remove cells with NA and Inf values
#' 
#' # First sample
#'   obs_dat1 <- data.frame("id" = seq(1, nrow(fr1@exprs), 1),
#'                          "g1" = rep(1, nrow(fr1@exprs)),
#'                          "log10_CD4" = log(fr1@exprs[ , 5], 10),
#'                          "log10_CD38" = log(fr1@exprs[ , 6], 10))
#' # Second sample
#'   obs_dat2 <- data.frame("id" = seq(1, nrow(fr2@exprs), 1),
#'                          "g1" = rep(2, nrow(fr2@exprs)),
#'                          "log10_CD4" = log(fr2@exprs[ , 5], 10),
#'                          "log10_CD38" = log(fr2@exprs[ , 6], 10))
#' # Full set
#'   obs_dat <- rbind(obs_dat1, obs_dat2)
#'   obs_dat <- obs_dat[complete.cases(obs_dat), ] # remove NAs
#'   obs_dat <- obs_dat[is.finite(rowSums(obs_dat)), ] # remove Infs
#'   obs_dat$g1 <- as.factor(obs_dat$g1) # set "g1" as binary factor
#' 
#' # Run rrs() function
#'   test_rrs <- rrs(dat = obs_dat, p_correct = "none")
#'   
rrs <- function(dat, 
                alpha = 0.05, 
                p_correct = "none",
                nbc = NULL,
                doplot = FALSE, 
                rcols = c("#FF0000", "#cccccc", "#0000FF"),
                win = NULL,
                verbose = FALSE,
                ...) {
  
  `%!in%` <- function(x, y)!(`%in%`(x, y)) # helpful custom function
  
  # Checks
  ## dat
  if ("data.frame" %!in% class(dat)) { stop("'dat' must be class 'data.frame'") }
  
  ## group
  if (nlevels(dat[ , 2]) != 2) { stop("The second feature of 'dat' must be a binary factor.") }
  
  ## p_correct
  match.arg(p_correct, choices = c("none", "correlated", "uncorrelated"))
  
  ## alpha
  if (alpha <= 0 | alpha >= 1 ) {
    stop("The argument 'alpha' must be a numeric value between zero (0) and one (1).")
  }
  
  ## rcols
  if (length(rcols) != 3) {
    stop("The argument 'rcols' must be a vector of length three (3).")
  }
  
  ## win
  if (!is.null(win) & class(win) != "owin") { stop("'win' must be class 'owin'") }
  if (is.null(win)) {
    dat <- as.data.frame(dat)
    dat <- dat[!is.na(dat[ , 4]) & !is.na(dat[ , 5]) , ]
    chul <- grDevices::chull(dat[ , 4:5])
    chul_coords <- dat[ , 4:5][c(chul, chul[1]), ]
    win <- spatstat.core::owin(poly = list(x = rev(chul_coords[ , 1]),
                                           y = rev(chul_coords[ , 2])))
  }
  
  Vnames <- names(dat) # axis names
  names(dat) <- c("id", "G1", "G2", "V1", "V2")
  levels(dat$G1) <- c("case", "control")
  
  # Create PPP
  suppressMessages(suppressWarnings(c1_ppp <- spatstat.core::ppp(x = dat$V1,
                                                                 y = dat$V2,
                                                                 marks = dat$G1,
                                                                 window = win)))
  
  # Estimate SRR and p-values
  c1_h0 <- sparr::OS(c1_ppp, "geometric")
  suppressMessages(suppressWarnings(out <- sparr::risk(f = c1_ppp,
                                                       h0 = c1_h0,
                                                       tolerate = TRUE,
                                                       edge = "diggle",
                                                       verbose = verbose,
                                                       log = FALSE,
                                                       ...)))
  
  if (all(is.na(out$rr$v))) { 
    cat("relative risk unable to be estimated")
    return(out)
  }
  
  ## Calculate log RR
  suppressMessages(suppressWarnings(out$lrr <- log(out$rr)))
  
  # Alpha level
  if (p_correct == "none") { out$alpha <- alpha }
  if (p_correct == "correlated") {
    cat("\nPlease be patient... Calculating correlated Bonferroni correction\n")
    alpha_correct <- pval_correct(input = out$lrr, alpha = alpha, nbc = nbc)
    out$alpha <- alpha_correct$correlated 
  } 
  if (p_correct == "uncorrelated") { out$alpha <- alpha / prod(out$rr$dim) }
  
  if (doplot == TRUE) {
    # Graphics
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op))
    # Plotting inputs
    ## Colorkeys
    c1_plot <- lrr_plot(input = out$lrr, cols = rcols, midpoint = 0)
    ## Extent
    blim <- as.vector(raster::extent(c1_plot$v))
    xlims <- blim[1:2]
    ylims <- blim[3:4]
    
    # Plot of ratio and p-values
    #Vnames <- names(dat) # axis names
    Ps <- pval_plot(out$P, alpha = out$alpha)
    if (all(raster::values(Ps)[!is.na(raster::values(Ps))] == 2) | all(is.na(raster::values(Ps)))) { 
      pcols <- rcols[2]
      brp <- c(1, 3)
      atp <- 2
      labp <- "No"
    } else {
      pcols <- rcols
      brp <- c(1, 1.67, 2.33, 3)
      atp <- c(1.33, 2, 2.67)
      labp <- c("Numerator", "Insignificant", "Denominator")
    }
    
    graphics::par(pty = "s", bg = "white")
    invisible(graphics::split.screen(matrix(c(0, 0.45, 0.55, 1, 0.14, 0.14, 0.86, 0.86),
                                            ncol = 4)))
    graphics::screen(1)
    fields::image.plot(c1_plot$v, 
                       breaks = c1_plot$breaks,
                       col = c1_plot$cols,
                       xlim = xlims,
                       ylim = ylims,
                       axes = TRUE,
                       cex.lab = 1,
                       xlab = paste(Vnames[4], "\n", sep = ""),
                       ylab = Vnames[5],
                       cex = 1,
                       horizontal = TRUE,
                       axis.args = list(at = c1_plot$at,
                                        labels = c1_plot$labels,
                                        cex.axis = 0.67),
                       main = paste("log relative risk\n(bandwidth:",
                                    round(c1_h0, digits = 3),
                                    "units)",
                                    sep = " "))
    par(bg = "transparent")
    graphics::screen(2)
    fields::image.plot(Ps, 
                       breaks = brp,
                       col = pcols,
                       xlim = xlims,
                       ylim = ylims, 
                       xlab = paste(Vnames[4], "\n", sep = ""),
                       ylab = "",
                       cex = 1,
                       axes = TRUE,
                       horizontal = TRUE,
                       axis.args = list(at = atp,
                                        labels = labp,
                                        cex.axis = 0.67),
                       main = paste("significant p-values\n(alpha = ",
                                    formatC(out$alpha, format = "e", digits = 2),
                                    ")",
                                    sep = ""))
    graphics::close.screen(all = TRUE) # exit split-screen mode
  }
  
  return(out)
}

#' A single gate for two conditions
#'
#' Estimates a ratio of relative risk surfaces and computes the asymptotic p-value surface for a single gate with two conditions. Includes features for basic visualization. This function is used internally within the \code{\link{gating}} function to extract the points within the significant areas. This function can also be used as a standalone function.
#'
#' @param dat Input data frame flow cytometry data with five (5) features (columns): 1) ID, 2) Condition A ID, 3) Condition B ID, 4) Marker A as x-coordinate, 5) Marker B as y-coordinate.
#' @param alpha Numeric. The two-tailed alpha level for significance threshold (default is 0.05).
#' @param p_correct Character string specifying whether to apply a correction for multiple comparisons including a Bonferroni correction \code{p_correct = "uncorrelated"} or a correlated Bonferroni correction \code{p_correct = "correlated"}. If \code{p_correct = "none"} then no correction is applied. 
#' @param nbc Optional. An integer for the number of bins when \code{p_correct = "correlated"}. Similar to \code{nbclass} argument in \code{\link[pgirmess]{correlog}}. The default is the average number of gridded knots in one-dimension (i.e., x-axis). 
#' @param doplot Logical. If \code{TRUE}, the output includes basic data visualizations.
#' @param rcols Character string of length three (3) specifying the colors for: 1) group A (numerator), 2) neither, and 3) group B (denominator) designations. The defaults are \code{c("#FF0000", "#cccccc", "#0000FF")} or \code{c("red", "grey80", "blue")}.
#' @param lower_lrr Optional, numeric. Lower cut-off value for the log relative risk value in the color key (typically a negative value). The default is no limit and the color key will include the minimum value of the log relative risk surface. 
#' @param upper_lrr Optional, numeric. Upper cut-off value for the log relative risk value in the color key (typically a positive value). The default is no limit and the color key will include the maximum value of the log relative risk surface.
#' @param win Optional. Object of class \code{owin} for a custom two-dimensional window within which to estimate the surfaces. The default is NULL and calculates a convex hull around the data. 
#' @param verbose Logical. If \code{TRUE} will print function progress during execution. If \code{FALSE} (the default), will not print.
#' @param ... Arguments passed to \code{\link[sparr]{risk}} to select bandwidth, edge correction, and resolution.
#'
#' @details This function estimates a ratio of relative risk surfaces and computes the asymptotic p-value surface for a single gate with two conditions using three successive \code{\link[sparr]{risk}} functions. A relative risk surface is estimated for Condition A at each level of Condition B and then a ratio of the two relative risk surfaces is computed. 
#' 
#' \deqn{RR_{Condition B1} = \frac{Condition A2 of B1}{Condition A1 of B1}}
#' \deqn{RR_{Condition B2} = \frac{Condition A2 of B2}{Condition A1 of B2}}
#' \deqn{ln(rRR) = ln\left (\frac{RR_{Condition B2}}{CRR_{Condition B2}}\right )}
#' 
#' The p-value surface of the ratio of relative risk surfaces is estimated assuming asymptotic normality of the ratio value at each gridded knot. The bandwidth is fixed across all layers. Basic visualization is available if \code{doplot = TRUE}. 
#' 
#' Provides functionality for a correction for multiple testing.  If \code{p_correct = "uncorrelated"}, then a conventional Bonferroni correction is calculated by dividing the \code{alpha} level by the number of gridded knots across the estimated surface. The default in the \code{\link[sparr]{risk}} function is a resolution of 128 x 128 or n = 16,384 knots and a custom resolution can be specified using the \code{resolution} argument within the \code{\link[sparr]{risk}} function. If \code{p_correct = "correlated"} (NOTE: May take a considerable amount of computation resources and time), then a Bonferroni correction that takes into account the spatial correlation of the surface is calculated within the internal \code{pval_correct} function. The \code{alpha} level is divided by the minimum number of knots that are not spatially correlated. The minimum number of knots that are not spatially correlated is computed by counting the knots that are a distance apart that exceeds the minimum distance of non-significant spatial correlation based on a correlogram using the \code{\link[pgirmess]{correlog}} function. If \code{p_correct = "none"}, then the function does not account for multiple testing and uses the uncorrected \code{alpha} level. See the internal \code{pval_correct} function documentation for more details.
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
#' @importFrom spatstat owin ppp
#' @importFrom sparr OS risk
#' @export 
#'
#' @examples
#' if (interactive()) {
#' # Use 'extdata' from the {flowWorkspaceData} package
#'   flowDataPath <- system.file("extdata", package = "flowWorkspaceData")
#'   fcsFiles <- list.files(pattern = "CytoTrol", flowDataPath, full = TRUE)
#'   ncfs  <- ncdfFlow::read.ncdfFlowSet(fcsFiles)
#'   fr1 <- ncfs[[1]]
#'   fr2 <- ncfs[[2]]
#' 
#' # Comparison of two samples at two time points (two conditions) "g1" and "g2"
#' ## (Create a random binary variable for "g2")
#' ## One gate (two markers) "CD4", "CD38"
#' ## Log10 Transformation for both markers
#' ## Remove cells with NA and Inf values
#' 
#' # First sample
#'   obs_dat1 <- data.frame("id" = seq(1, nrow(fr1@exprs), 1),
#'                          "g1" = rep(1, nrow(fr1@exprs)),
#'                          "g2" = stats::rbinom(nrow(fr1@exprs), 1, 0.5),
#'                          "log10_CD4" = log(fr1@exprs[ , 5], 10),
#'                          "log10_CD38" = log(fr1@exprs[ , 6], 10))
#' # Second sample
#'   obs_dat2 <- data.frame("id" = seq(1, nrow(fr2@exprs), 1),
#'                          "g1" = rep(2, nrow(fr2@exprs)),
#'                          "g2" = stats::rbinom(nrow(fr2@exprs), 1, 0.5),
#'                          "log10_CD4" = log(fr2@exprs[ , 5], 10),
#'                          "log10_CD38" = log(fr2@exprs[ , 6], 10))
#' # Full set
#'   obs_dat <- rbind(obs_dat1, obs_dat2)
#'   obs_dat <- obs_dat[complete.cases(obs_dat), ] # remove NAs
#'   obs_dat <- obs_dat[is.finite(rowSums(obs_dat)), ] # remove Infs
#'   obs_dat$g1 <- as.factor(obs_dat$g1) # set "g1" as binary factor
#'   obs_dat$g2 <- as.factor(obs_dat$g2) # set "g2" as binary factor
#' 
#' # Run lotrrs() function
#'   test_lotrrs <- lotrrs(dat = obs_dat, p_correct = "none")
#' }   
#' 
lotrrs <- function(dat, 
                   alpha = 0.05, 
                   p_correct = "none",
                   nbc = NULL,
                   doplot = FALSE, 
                   rcols = c("#FF0000", "#cccccc", "#0000FF"),
                   lower_lrr = NULL,
                   upper_lrr = NULL,
                   win = NULL, 
                   verbose = FALSE, 
                   ...) {
  
  # Checks
  ## dat
  if ("data.frame" %!in% class(dat)) { stop("'dat' must be class 'data.frame'") }
  
  ## group
  if (nlevels(dat[ , 2]) != 2) { stop("The second feature of 'dat' must be a binary factor.") }
  
  if (nlevels(dat[ , 3]) != 2) { stop("The third feature of 'dat' must be a binary factor.") }
  
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
    win <- spatstat::owin(poly = list(x = rev(chul_coords[ , 1]),
                                           y = rev(chul_coords[ , 2])))
  }
  
  Vnames <- names(dat) # axis names
  names(dat) <- c("id", "G1", "G2", "V1", "V2")
  levels(dat$G1) <- c("case", "control")
  levels(dat$G2) <- c("0", "1")
  
  # Create two PPPs
  t0_df <- dat[dat$G2 == "0", ]
  t1_df <- dat[dat$G2 == "1", ]
  suppressMessages(suppressWarnings(t_ppp <- spatstat::ppp(x = dat$V1,
                                                                y = dat$V2,
                                                                marks = dat$G1,
                                                                window = win)))
  suppressMessages(suppressWarnings(t0_ppp <- spatstat::ppp(x = t0_df$V1,
                                                                 y = t0_df$V2,
                                                                 marks = t0_df$G1,
                                                                 window = win)))
  suppressMessages(suppressWarnings(t1_ppp <- spatstat::ppp(x = t1_df$V1,
                                                                 y = t1_df$V2,
                                                                 marks = t1_df$G1,
                                                                 window = win)))
  
  # Estimate two SRRs
  t_h0 <- sparr::OS(t_ppp, "geometric")
  rm(t_ppp) # conserve memory
  t0_rr <- sparr::risk(t0_ppp,
                       h0 = t_h0,
                       log = FALSE,
                       edge = "diggle",
                       verbose = verbose,
                       ...)
  t1_rr <- sparr::risk(t1_ppp,
                       h0 = t_h0,
                       log = FALSE,
                       edge = "diggle",
                       verbose = verbose,
                       ...)
  
  # Create objects of class 'bivden'
  t0_bd <- sparr::bivariate.density(t0_ppp,
                                    h0 = t_h0,
                                    edge = "diggle",
                                    verbose = FALSE,
                                    ...)
  t1_bd <- sparr::bivariate.density(t1_ppp,
                                    h0 = t_h0,
                                    edge = "diggle",
                                    verbose = FALSE,
                                    ...)
  t0_rr_bd <- list("z" = t0_rr$rr,
                   "h0" = t0_rr$f$h0,
                   "hp" = t0_rr$f$hp,
                   "h" = t0_rr$f$h,
                   "him" = t0_rr$f$him,
                   "q" = t0_bd$q,
                   "gamma" = t0_rr$f$gamma,
                   "geometric" = t0_rr$f$geometric,
                   "pp" = t0_ppp)
  t1_rr_bd <- list("z" = t1_rr$rr,
                   "h0" = t1_rr$f$h0,
                   "hp" = t1_rr$f$hp,
                   "h" = t1_rr$f$h,
                   "him" = t1_rr$f$him,
                   "q" = t1_bd$q,
                   "gamma" = t1_rr$f$gamma,
                   "geometric" = t1_rr$f$geometric,
                   "pp" = t1_ppp)
  class(t0_rr_bd) <- "bivden"
  class(t1_rr_bd) <- "bivden"
  rm(t0_ppp, t1_ppp, t0_bd, t1_bd) # conserve memory
  
  # Estimate the spatial relative risk surface of the two spatial relative risks
  ## Also estimate the asymptotic p-value surface
  suppressMessages(suppressWarnings(out <- sparr::risk(f = t1_rr_bd,
                                                         g = t0_rr_bd,
                                                         h0 = t_h0,
                                                         tolerate = TRUE,
                                                         edge = "diggle",
                                                         verbose = verbose,
                                                         log = FALSE,
                                                         ...)))
  
  suppressMessages(suppressWarnings(out$rr <- t1_rr_bd$z/t0_rr_bd$z))
  suppressMessages(suppressWarnings(out$lrr <- log(t1_rr_bd$z) - log(t0_rr_bd$z)))
  rm(t0_rr_bd, t1_rr_bd) # conserve memory
  
  if (all(is.na(out$rr$v))) { 
    cat("relative risk unable to be estimated")
    return(out)
  }
  
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
    tr_plot <- lrr_plot(input = out$lrr,
                        cols = rcols,
                        midpoint = 0,
                        upper_lrr = upper_lrr,
                        lower_lrr = lower_lrr)
    ## Extent
    blim <- as.vector(raster::extent(tr_plot$v))
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
      labp <- c("numerator", "insignificant", "denominator")
    }

    graphics::par(pty = "s", bg = "white")
    invisible(graphics::split.screen(matrix(c(0, 0.45, 0.55, 1, 0.14, 0.14, 0.86, 0.86),
                                            ncol = 4)))
    graphics::screen(1)
    fields::image.plot(tr_plot$v,
                       breaks = tr_plot$breaks,
                       col = tr_plot$cols,
                       xlim = xlims,
                       ylim = ylims,
                       axes = TRUE,
                       cex.lab = 1,
                       xlab = paste(Vnames[4], "\n", sep = ""),
                       ylab = Vnames[5],
                       cex = 1,
                       horizontal = TRUE,
                       axis.args = list(at = tr_plot$at,
                                        labels = tr_plot$labels,
                                        cex.axis = 0.67),
                       main = paste("log of the\nrelative risk surfaces\n(bandwidth:",
                                    round(t_h0, digits = 3),
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

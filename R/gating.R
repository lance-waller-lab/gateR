#' Gating strategy for mass cytometry data using spatial relative risk functions
#' 
#' Extracts cells within statistically significant combinations of fluorescent markers, successively, for a set of markers. Statistically significant combinations are identified using two-tailed p-values of a relative risk surface assuming asymptotic normality. This function is currently available for two-level comparisons of a single group (e.g., case/control) or two groups (e.g., case/control at time 1 and time 2). Provides functionality for basic visualization and multiple testing correction.
#'
#' @param dat Input data frame flow cytometry data with the following features (columns): 1) ID, 2) Condition A ID, 3) Condition B ID (optional), and a set of markers.
#' @param vars A vector of characters with the name of features (columns) within \code{dat} to use as markers for each gate. See details below.
#' @param n_condition A numeric value of either 1 or 2 designating if the gating is performed with one condition or two conditions. 
#' @param numerator Logical. If \code{TRUE} (the default), cells will be extracted within all statistically significant numerator (i.e., case) clusters. If \code{FALSE}, cells will be extracted within all statistically significant denominator (i.e., control) clusters. 
#' @param alpha Numeric. The two-tailed alpha level for significance threshold (default is 0.05).
#' @param p_correct Character string specifying whether to apply a correction for multiple comparisons including a Bonferroni correction \code{p_correct = "uncorrelated"} or a correlated Bonferroni correction \code{p_correct = "correlated"}. If \code{p_correct = "none"} (the default), then no correction is applied. 
#' @param nbc Optional. An integer for the number of bins when \code{p_correct = "correlated"}. Similar to \code{nbclass} argument in \code{\link[pgirmess]{correlog}}. The default is the average number of gridded knots in one-dimension (i.e., x-axis).
#' @param doplot Logical. If \code{TRUE}, the output includes basic data visualizations.
#' @param rcols Character string of length three (3) specifying the colors for: 1) Group A, 2) Neither, and 3) Group B designations. The defaults are \code{c("#FF0000", "#cccccc", "#0000FF")} or \code{c("red", "grey80", "blue")}.
#' @param lower_lrr Optional, numeric. Lower cut-off value for the log relative risk value in the color key (typically a negative value). The default is no limit and the color key will include the minimum value of the log relative risk surface. 
#' @param upper_lrr Optional, numeric. Upper cut-off value for the log relative risk value in the color key (typically a positive value). The default is no limit The default is no limit and the color key will include the maximum value of the log relative risk surface.
#' @param win Optional. Object of class \code{owin} for a custom two-dimensional window within which to estimate the surfaces. The default is NULL and calculates a convex hull around the data. 
#' @param verbose Logical. If \code{TRUE} will print function progress during execution. If \code{FALSE} (the default), will not print.
#' @param ... Arguments passed to \code{\link[sparr]{risk}} to select bandwidth, edge correction, and resolution.
#'
#' @details This function performs a sequential gating strategy for mass cytometry data comparing two levels with one or two conditions. Gates are typically two-dimensional space comprised of two fluorescent markers. The two-level comparison allows for the estimation of a spatial relative risk function and the computation of p-value based on an assumption of asymptotic normality. Cells within statistically significant areas are extracted and used in the next gate. This function relies heavily upon the \code{\link[sparr]{risk}} function. Basic visualization is available if \code{doplot = TRUE}. 
#' 
#' The \code{vars} argument must be a vector with an even-numbered length where the odd-numbered elements are the markers used on the x-axis of a gate and the even-numbered elements are the markers used on the y-axis of a gate. For example, if \code{vars = c("V1", "V2", "V3", and "V4")} then the first gate is "V1" on the x-axis and "V2" on the y-axis and then the second gate is V3" on the x-axis and "V4" on the y-axis. Makers can be repeated in successive gates. 
#' 
#' The \code{n_condition} argument specifies if the gating strategy is performed for one condition or two conditions. If \code{n_condition = 1}, then the function performs a one condition gating strategy using the internal \code{rrs} function, which computes the statistically significant areas (clusters) of a relative risk surface at each gate and selects the cells within the clusters specified by the \code{numerator} argument. If \code{n_condition = 2}, then the function performs a two conditions gating strategy using the internal \code{lotrrs} function, which computes the statistically significant areas (clusters) of a ratio of relative risk surfaces at each gate and selects the cells within the clusters specified by the \code{numerator} argument. See the documentation for the internal \code{rrs} and \code{lotrrs} functions for more details.
#' 
#' The p-value surface of the ratio of relative risk surfaces is estimated assuming asymptotic normality of the ratio value at each gridded knot. The bandwidth is fixed across all layers.
#' 
#' Provides functionality for a correction for multiple testing.  If \code{p_correct = "uncorrelated"}, then a conventional Bonferroni correction is calculated by dividing the \code{alpha} level by the number of gridded knots across the estimated surface. The default in the \code{\link[sparr]{risk}} function is a resolution of 128 x 128 or n = 16,384 knots and a custom resolution can be specified using the \code{resolution} argument within the \code{\link[sparr]{risk}} function. If \code{p_correct = "correlated"} (NOTE: May take a considerable amount of computation resources and time), then a Bonferroni correction that takes into account the spatial correlation of the surface is calculated within the internal \code{pval_correct} function. The \code{alpha} level is divided by the minimum number of knots that are not spatially correlated. The minimum number of knots that are not spatially correlated is computed by counting the knots that are a distance apart that exceeds the minimum distance of non-significant spatial correlation based on a correlogram using the \code{\link[pgirmess]{correlog}} function. If \code{p_correct = "none"} (the default), then the function does not account for multiple testing and uses the uncorrected \code{alpha} level. See the internal \code{pval_correct} function documentation for more details.
#' 
#'
#' @return An object of class \code{list}. This is a named list with the following components:
#' 
#' #' \describe{
#' \item{\code{obs}}{An object of class 'tibble' of the same features as \code{dat} that includes the information for the cells extracted with significant clusters in the final gate.}
#' \item{\code{gate}}{An object of class 'list' of 'rrs' objects from each gate.}
#' }
#' 
#' The objects of class 'rrs' is similar to the output of the \code{\link[sparr]{risk}} function with two additional components:
#' \describe{
#' \item{\code{rr}}{An object of class 'im' with the relative risk surface.}
#' \item{\code{f}}{An object of class 'im' with the spatial density of the numerator.}
#' \item{\code{g}}{An object of class 'im' with the spatial density of the denominator.}
#' \item{\code{P}}{An object of class 'im' with the asymptotic p-value surface.}
#' \item{\code{lrr}}{An object of class 'im' with the log relative risk surface.}
#' \item{\code{alpha}}{A numeric value for the alpha level used within the gate.}
#' }
#'
#' @importFrom grDevices chull
#' @importFrom maptools unionSpatialPolygons
#' @importFrom raster rasterToPolygons values
#' @importFrom sp coordinates over
#' @importFrom spatstat.core owin
#' @importFrom tibble add_column
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
#' ## Two gates (four markers) "CD4", "CD38", "CD8", and "CD3"
#' ## Log10 Transformation for both markers
#' ## Remove cells with NA and Inf values
#' 
#' # First sample
#'   obs_dat1 <- data.frame("id" = seq(1, nrow(fr1@exprs), 1),
#'                          "g1" = rep(1, nrow(fr1@exprs)),
#'                          "log10_CD4" = log(fr1@exprs[ , 5], 10),
#'                          "log10_CD38" = log(fr1@exprs[ , 6], 10),
#'                          "log10_CD8" = log(fr1@exprs[ , 7], 10),
#'                          "log10_CD3" = log(fr1@exprs[ , 8], 10))
#' # Second sample
#'   obs_dat2 <- data.frame("id" = seq(1, nrow(fr2@exprs), 1),
#'                          "g1" = rep(2, nrow(fr2@exprs)),
#'                          "log10_CD4" = log(fr2@exprs[ , 5], 10),
#'                          "log10_CD38" = log(fr2@exprs[ , 6], 10),
#'                          "log10_CD8" = log(fr2@exprs[ , 7], 10),
#'                          "log10_CD3" = log(fr2@exprs[ , 8], 10))
#' # Full set
#'   obs_dat <- rbind(obs_dat1, obs_dat2)
#'   obs_dat <- obs_dat[complete.cases(obs_dat), ] # remove NAs
#'   obs_dat <- obs_dat[is.finite(rowSums(obs_dat)), ] # remove Infs
#'   obs_dat$g1 <- as.factor(obs_dat$g1) # set "g1" as binary factor
#' 
#' # Run gating() function
#' ## Single condition, no multiple testing correction
#'   test_gate <- gateR::gating(dat = obs_dat,
#'                              vars = c("log10_CD4", "log10_CD38",
#'                                       "log10_CD8", "log10_CD3"),
#'                              n_condition = 1,
#'                              p_correct = "none")
#' 
gating <- function(dat,
                   vars,
                   n_condition = c(1, 2),
                   numerator = TRUE,
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
  
  ## vars
  if (!all(vars %in% names(dat))) {
    stop("All elements in the argument 'vars' must match named features of 'dat'." )
  }
  
  if ((length(vars) %% 2) != 0 ) {
    stop("The argument 'vars' must be a character vector with an even-numbered length.")
  }
  
  ## n_condition
  match.arg(as.character(n_condition), choices = 1:2)
  
  ## p_correct
  match.arg(p_correct, choices = c("none", "correlated", "uncorrelated"))
  
  ## numerator
  if (numerator == TRUE) { type_cluster <- "numerator" 
  } else { 
    type_cluster <- "denominator" 
  }
  
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
  
  # Format data input
  dat <- dat[!is.na(dat[ , which(colnames(dat) %in% vars[1])]) &
               !is.na(dat[ , which(colnames(dat) %in% vars[2])]), ]

  if (n_condition == 1 & nlevels(dat[ , 2]) != 2) { stop("The second feature of 'dat' must be a binary factor.") }
  
  if (n_condition == 2 & nlevels(dat[ , 2]) != 2) { stop("The second feature of 'dat' must be a binary factor.") }
  
  if (n_condition == 2 & nlevels(dat[ , 3]) != 2) { stop("The third feature of 'dat' must be a binary factor.") }
  
  if (n_condition == 2) { dat_gate <- dat[ , c(1:3, which(colnames(dat) %in% vars))]
  } else {
    dum_condition <- rep(1, nrow(dat))
    dat_gate <- dat[ , c(1:2, which(colnames(dat) %in% vars))]
    dat_gate <- tibble::add_column(dat_gate, condition = dum_condition, .after = 2)
  }
  dat <- as.data.frame(dat)
  dat_gate <- as.data.frame(dat_gate)
  n_gate <- length(vars) / 2 # calculate the number of gates

  # Create empty list to save output of each gate
  list_gate <- vector('list', length(n_gate))

  # Gating
  for (k in 1:n_gate) {

  if (k == 1) { j <- 1 } else { j <- k * 2 - 1 } # vars indexing per gate

    # format data by selecting the desired vars
    dat_gate <- dat_gate[!is.na(dat_gate[ , which(colnames(dat_gate) %in% vars[j])]) &
                           !is.na(dat_gate[ , which(colnames(dat_gate) %in% vars[j + 1])]), ]
    xvar <-  which(colnames(dat_gate) %in% vars[j])
    yvar <-  which(colnames(dat_gate) %in% vars[j + 1])
    df <- dat_gate[ , c(1:3, xvar, yvar)]
    df <- df[!is.na(df[ , 4]) & !is.na(df[ , 5]), ] # remove NAs

    # Create custom window with a convex hull
    chul <- grDevices::chull(df[ , 4:5])
    chul_coords <- df[ , 4:5][c(chul, chul[1]), ]
    win_gate <- spatstat.core::owin(poly = list(x = rev(chul_coords[ , 1]),
                                                y = rev(chul_coords[ , 2])))

    # Estimate significant relative risk areas
    ## Bonferroni correction only necessary in first gate
    if (k == 1) { p_correct <- p_correct } else { p_correct <- "none"}

    if (n_condition == 2) {
    out <- lotrrs(dat = df,
                  win = win_gate,
                  doplot = doplot,
                  alpha = alpha,
                  p_correct = p_correct,
                  nbc = nbc,
                  rcols = rcols,
                  lower_lrr = lower_lrr,
                  upper_lrr = upper_lrr,
                  verbose = verbose,
                  ...)
    } else {
    out <- rrs(dat = df,
               win = win_gate,
               doplot = doplot,
               alpha = alpha,
               p_correct = p_correct,
               nbc = nbc,
               rcols = rcols,
               lower_lrr = lower_lrr,
               upper_lrr = upper_lrr,
               verbose = verbose,
               ...)
    }

    # Convert p-value surface into a categorized raster
    ## v == 2 : significant T1 numerator;  v == 1: not
    Ps <- pval_plot(input = out$P, alpha = out$alpha)
    list_gate[[k]] <- out # save for output
    rm(out, df, win_gate) # conserve memory

    # Go back one gate if current gate has no significant area and produce output of previous gate
    if (all(raster::values(Ps)[!is.na(raster::values(Ps))] == 2) | all(is.na(raster::values(Ps)))) {
      cat(paste("\nGate", k, "yeilded no significant", type_cluster, "cluster(s)...",
                "Returning results from previous gate\n",
                sep = " "))
      output <- dat[which(dat[ , 1] %in% dat_gate[ , 1]), ]
      out_list <- list("obs" = output, "gate" = list_gate)
      return(out_list)
    }

    # convert categorized raster to gridded polygons
    out_pol <- raster::rasterToPolygons(Ps)
    rm(Ps) # conserve memory

    if (numerator == TRUE) { v <- 1 } else { v <- 3 }
    # combine gridded polygons of similar categories
    pols <- try(maptools::unionSpatialPolygons(out_pol[out_pol$layer == v, ],
                                               IDs = rep(1, length(out_pol[out_pol$layer == v, ]))), silent = TRUE)
    if("try-error" %in% class(pols)) {
      cat(paste("Gate", k, "yeilded no significant", type_cluster, "cluster(s)...",
                "Returning results from previous gate\n",
                sep = " "))
      output <- dat[which(dat[ , 1] %in% dat_gate[ , 1]), ]
      out_list <- list("obs" = output, "gate" = list_gate)
      return(out_list)
    }
    rm(out_pol) # conserve memory

    # Overlay points
    sp::coordinates(dat_gate) <- ~ cbind(dat_gate[,which(colnames(dat_gate) %in% vars[j])],
                                         dat_gate[,which(colnames(dat_gate) %in% vars[j + 1])])

    # extract points within significant cluster
    dat_gate <- as.data.frame(dat_gate[!is.na(sp::over(dat_gate, pols)), ])
    rm(pols) # conserve memory

    # Output for the final gate
    if (k == n_gate) {
      cat(paste("\nObservations within significant", type_cluster, "cluster(s) of Gate", k, "\n",sep = " "))
      output <- dat[which(dat[ , 1] %in% dat_gate[ , 1]), ]
      out_list <- list("obs" = output, "gate" = list_gate)
      return(out_list)
    }
  }
}

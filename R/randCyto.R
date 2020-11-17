#' Subset of the 'extdata' data in the 'flowWorkspaceData' package
#'
#' A sample dataset containing information about flow cytometry data with two binary conditions and four markers. The data are a random subset of the 'extdata' data in the 'flowWorkspaceData' package found on Bioconductor \url{http://bioconductor.org/packages/release/data/experiment/html/flowWorkspaceData.html} and formated for 'gateR' input. The selected markers are arcsinh transformed.
#'
#' @format A data frame with 11763 rows and 7 variables:
#' \describe{
#'   \item{id}{cell ID number}
#'   \item{g1}{binary condition #1}
#'   \item{g2}{binary condition #2}
#'   \item{arcsinh_CD4}{arcsinh-transformed CD4}
#'   \item{arcsinh_CD38}{arcsinh-transformed CD38}
#'   \item{arcsinh_CD8}{arcsinh-transformed CD8}
#'   \item{arcsinh_CD3}{arcsinh-transformed CD3}
#' }
#' @examples
#' head(randCyto)
#'
#' @source \url{https://github.com/Waller-SUSAN/gateR/blob/master/README.md}
"randCyto"

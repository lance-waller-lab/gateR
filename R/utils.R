# Not within. Opposite of %in%
`%!in%` <- function(x, y)!(`%in%`(x, y))

# False Discovery Rate (Benjamini & Hochberg)
fdr <- function(pvals, alpha) {
  pcrit <- NULL
  m <- length(pvals)
  for (i in 1:m) {
    if (pvals[i] <= (i/m) * alpha) { 
      pcrit <- pvals[i]
    }
  }
  return(max(pcrit, pvals[1]))
}

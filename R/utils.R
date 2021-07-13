# Not within. Opposite of %in%
`%!in%` <- function(x, y)!(`%in%`(x, y))

# 
# https://stackoverflow.com/a/28017199/6784787
move_to_start <- function(x, to_move) { x[ , c(to_move, setdiff(colnames(x), to_move))] } 

x <- matrix(data = rpois(3500 * 100000, 1000), nrow = 3500, ncol = 100000)
v <- apply(X = x, MARGIN = 2, FUN = function(col) which.max(col))

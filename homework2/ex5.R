# Approach (1): No integration

b_grid <- seq(0.01, 1, by = 0.00001)
med <- qgamma(0.5, 100*b_grid^2, b_grid)

b <- b_grid[which.min(abs(med - 75))]
a <- 100 * b^2

# Approach (2): Using integration

b_grid <- seq(0.01, 1, by = 0.0001)

int <- sapply(b_grid, function(b) {
  1 - integrate(dgamma, shape = 100*b^2, rate = b, lower = 0, upper = 75)$value
})

b <- b_grid[which.min(abs(int - 0.5))]
a <- 100 * b^2

# Validate using Monte Carlo sampling

s <- 1e+7
lambda <- rgamma(n = s, shape = a, rate = b)

hist(lambda, main = NA, xlab = "λ", ylab = "Frequency")

mean(lambda > 75)

var(lambda)
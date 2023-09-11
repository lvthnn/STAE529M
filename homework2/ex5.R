# Exercise 2.4.5

# Method (*)

b_grid <- seq(0.01, 1, by = 0.00001)
med <- qgamma(0.5, 100*b_grid^2, b_grid)

b <- b_grid[which.min(abs(med - 75))]
a <- 100 * b^2

print(c(a, b))

# Method (**)

b_grid <- seq(0.01, 1, by = 0.0001)

int <- sapply(b_grid, function(b) {
  integrate(dgamma, shape = 100*b^2, rate = b, lower = 0, upper = 75)$value
})

b <- b_grid[which.min(abs(int - 0.5))]
a <- 100 * b^2

print(c(a, b))

# Validate using Monte Carlo sampling

s <- 1e+7
lambda <- rgamma(n = s, shape = a, rate = b)

hist(lambda, main = NA, xlab = "λ", ylab = "Frequency")

mean(lambda > 75)

var(lambda)

# Visualise estimate distribution

t <- seq(0, 150, by = 0.01)
ft <- dgamma(x = t, shape = a, rate = b)

plot(t, ft, 'l', xlab = latex2exp::TeX("\\lambda"), ylab = "Density")

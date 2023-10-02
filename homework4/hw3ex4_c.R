library(invgamma)
library(latex2exp)

set.seed(42)

S <- 1e+4 + 2000
Y <- 1:10
n <- 10
a <- 1

theta           <- matrix(NA, nrow = S, ncol = 11)
theta[1,]       <- c(rep(var(Y), 10), 1)
colnames(theta) <- c(paste0("s_", 1:n), "b")


for (s in 2:S) {
  for (i in 1:10) {
    A <- a + 1/2
    B <- 1/2 * Y[i]^2 + theta[s - 1,"b"]

    si <- rinvgamma(1, shape = A, rate = B)
    theta[s,i] <- si
  }

  A <- n*a + 1
  B <- sum(1/theta[s,1:10]) + 1

  b <- rgamma(1, shape = A, rate = B)
  theta[s,11] <- b
}

theta <- theta[2001:S,]

apply(theta, 2, mean)

# VISUALISATIONS

labs <- gsub("s_([0-9]*)", "$\\\\\\sigma_{\\1}^2$", colnames(theta))

mas <- matrix(NA, nrow = 1e+4, ncol = 11)

for (i in 1:11) {
  mas[,i] <- cumsum(theta[,i]) / seq_along(theta[,i])
}

pdf("trace_plots_vars_skewed.pdf", width = 9)
par(mfrow = c(3, 4))
for (i in 1:11) {
  plot(theta[,i], ylab = TeX(labs[i]))
  lines(mas[,i], col = "red")
}
dev.off()

pdf("densities_vars_skewed.pdf", width = 9)
par(mfrow = c(3, 4))
for (i in 1:11) {
  plot(density(theta[,i]), xlab = TeX(labs[i]), main = NA)
}
dev.off()

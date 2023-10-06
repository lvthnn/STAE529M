library(cubature)

Y <- c(12, 10, 22)
S <- c(3, 3, 10)

m  <- 10000
mu <- seq(-28, 50, length.out = m)

D <- 1 # prior µ = 1
for (i in 1:length(Y)) {
  D <- D * dnorm(Y[i], mu, S[i]) # likelihood
}
W <- D/sum(D) # standardise to form weights on partition

plot(mu, D, type='l')

# find standardising constant for integration
post <- function(mu, Y) { prod(dnorm(Y, mu, S)) }

g0 <- function(mu, Y) { post(mu,Y) }
g1 <- function(mu, Y) { mu * post(mu, Y) }

m0 <- adaptIntegrate(g0, -28, 50, Y = Y)$int
m1 <- adaptIntegrate(g1, -28, 50, Y = Y)$int

post_mean <- m1/m0
map_est <- sum(Y/(S^2))/sum(1/(S^2))

# alternative method is
sum(mu * W)

pdf("posterior_density_map_mean.pdf")
plot(mu, D/m0, type = 'l', xlab = expression(mu), ylab = "Posterior Density", xlim = c(-2, 25))
abline(v = post_mean, lty = 2, col = "blue", lwd = 1.5)
abline(v = map_est, lty = 3, col = "orange", lwd = 1.5)
legend("topright", inset = c(0.05, 0.05), legend = c("Posterior", "Posterior Mean", "MAP"), col = c("black", "blue", "orange"), lty = c(1, 2, 3))
dev.off()

# data
Y1 <- c( 2.0, -3.1, -1.0,  0.2,  0.3, 0.4)
Y2 <- c(-3.5, -1.6, -4.6, -0.9, -5.1, 0.1)
n1 <- n2 <- 6

Y1_bar <- mean(Y1)
Y2_bar <- mean(Y2)

var1 <- mean((Y1 - Y1_bar)^2)
var2 <- mean((Y2 - Y2_bar)^2)
sp <- sqrt(var1/2 + var2/2) # pooled sd

# posterior parameters
delta   <- Y2_bar - Y1_bar
sd_post <- sp * sqrt(1/n1 + 1/n2)

ci <- delta + sd_post * qt(c(0.025, 0.975), df = n1 + n2)

t  <- seq(delta - 5, delta + 5, length.out = 1000)
ft <- dnorm(t, delta, sd_post)

pdf("ex1-posterior.pdf")
plot(t, ft, type = "l", xlab = expression(delta), ylab = "Posterior Density", lwd = 1.5)
segments(x0 = ci[1], x1 = ci[1], y0 = 0, y1 = dnorm(ci[1], delta, sd_post))
segments(x0 = ci[2], x1 = ci[2], y0 = 0, y1 = dnorm(ci[2], delta, sd_post))
segments(x0 = delta, x1 = delta, y0 = 0, y1 = dnorm(delta, delta, sd_post), col = 2)
legend("topright", inset = c(0.05, 0.05), legend = c("Posterior mean", "95% credible set"), col = c(2,1), lty = c(1, 1))
dev.off()


# jags to analyse sensititity

library(rjags)

model_init <- textConnection("model{
  for (i in 1:n) {
    Y1[i] ~ dnorm(mu, tau)
    Y2[i] ~ dnorm(mu + delta, tau)
  } 

  mu    ~ dnorm(0, 0.0001)
  delta ~ dnorm(0, 0.0001)
  tau   ~ dgamma(0.1, 0.1)
  sigma <- 1/sqrt(tau)
}")

data  <- list(Y1 = Y1, Y2 = Y2, n = n1)

model <- jags.model(model_init, data = data, n.chains = 2, quiet = TRUE)

update(model, 1e+4, progress.bar = "none")

params  <- c("delta")
samples <- coda.samples(model, variable.names = params, n.iter = 2e+4, progress.bar = "none")

pdf("ex1-jags.pdf", width = 12)
plot(samples)
dev.off()

# delta is negative in plug-in calculations with 95% credible set disjoint from zero,
# so this suggests effectiveness 

# however jags returns result slightly different result with 95% cs including zero, so
# the analysis seems to be sensitive to prior

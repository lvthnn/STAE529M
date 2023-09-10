# Draw the posterior distribution
m <- 5
Y <- 10
a <- b <- 1

t <- seq(0, 1, by = 0.001)
ft <- dbeta(x = t, shape1 = a + m, shape2 = b + Y)

pdf("posterior_negbinom_beta_density.pdf")
plot(t, ft, type="l", xlab = latex2exp::TeX("\\theta | Y"), ylab = "Posterior Density")
abline(v = qbeta(p = 0.025, shape1 = a + m, shape2 = b + Y), col = "red", lty = 2, lwd = 1.5)
abline(v = qbeta(p = 0.975, shape1 = a + m, shape2 = b + Y), col = "red", lty = 2, lwd = 1.5)
legend("topright", inset = c(0.05, 0.05), cex = 0.8, legend = c("95% credible interval"), col = "red", lty = 2, lwd = 1.5)
dev.off()

# Monte Carlo validation of credible interval
s <- 1e+7
l <- qbeta(p = 0.025, shape1 = a + m, shape2 = b + Y)
u <- qbeta(p = 0.975, shape1 = a + m, shape2 = b + Y)

samp <- rbeta(n = s, shape1 = a + m, shape2 = b + Y)
mean(samp > l & samp < u)

# Plot the priot

library(latex2exp)

t <- seq(0.001, 10, by = 0.001)
ft <- 1/t

pdf("jeffreys_prior_lambda.pdf")
plot(t, ft, type = 'l', ylim = c(0, 10), xlab = TeX("\\lambda"), ylab = TeX("\\pi(\\lambda)"))
dev.off()

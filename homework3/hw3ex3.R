library(lattice)

cs <- 1:10
ys <- c(12, 90, 80, 5, 63, 15, 67, 22, 56, 33)
ns <- c(50, 150, 63, 10, 63, 8, 56, 19, 63, 19)

dat <- data.frame(cs, ys, ns)
colnames(dat) <- c("county", "approve", "disapprove")

dat$prop_approve <- dat$approve / (dat$approve + dat$disapprove)
dat$n <- dat$approve + dat$disapprove

# (a) analyse using jeffrey's prior

cis <- lapply(seq_len(nrow(dat)), function(i) {
  yi <- dat[i,][[2]]
  ni <- dat[i,][[5]]

  l <- qbeta(p = 0.025, shape1 = yi + 1/2, shape2 = ni - yi + 1/2)
  u <- qbeta(p = 0.975, shape1 = yi + 1/2, shape2 = ni - yi + 1/2)
  return(c(l, u))
})

cis <- do.call(rbind, cis)

dat$ci_l <- cis[,1]
dat$ci_u <- cis[,2]

# Monte Carlo simulation to verify validity of intervals
sapply(seq_len(nrow(dat)), function(i) {
  s <- 1e+6
  yi <- dat[i,][[2]]
  ni <- dat[i,][[5]]
  ci_l <- dat[i,][[6]]
  ci_u <- dat[i,][[7]]

  r <- rbeta(n = s, shape1 = yi + 1/2, shape2 = ni - yi + 1/2)
  return(mean(r > ci_l & r < ci_u))
})

# (b) select a and b so that the E() and Var() of the Beta(a, b) dist.
# match mean and variance of sample proportions p1,...,p10

prop_mean <- mean(dat$prop_approve)
prop_var <- var(dat$prop_approve)

cost <- function(pars) {
  a <- pars[1]
  b <- pars[2]

  eq1 <- a/(a + b) - prop_mean 
  eq2 <- (a*b)/((a + b)^2 * (a + b + 1)) - prop_var

  return(eq1^2 + eq2^2)
}

res <- optim(par = c(1, 1), cost)
a <- res$par[1]
b <- res$par[2]

t <- seq(0, 1, by = 1e-4)
ft <- dbeta(t, a, b)
plot(t, ft, type = "l")

# Monte Carlo for verification
s <- 1e+6
r <- rbeta(s, a, b)
hist(r)

mean(r)
var(r)

# Empirical Bayes using a and b from last term

cis_emp_bayes <- lapply(seq_len(nrow(dat)), function(i) {
  row <- dat[i, ]
  yi <- row[[2]]
  ni <- row[[5]]

  l <- qbeta(p = 0.025, shape = yi + a, shape2 = ni - yi + b)
  u <- qbeta(p = 0.975, shape = yi + a, shape2 = ni - yi + b)

  return(c(l, u))
})

cis_emp_bayes <- do.call(rbind, cis_emp_bayes)

dat$ci_emp_l <- cis_emp_bayes[,1]
dat$ci_emp_u <- cis_emp_bayes[,2]

pdf("comparison_pci_identity.pdf", width = 12)
par(mfrow = c(1, 2))
with(dat, plot(ci_l, ci_emp_l, xlab = "Lower JP 95% PCI", ylab = "Lower EB 95% PCI"))
abline(a = 0, b = 1, lty = 2)
with(dat, plot(ci_u, ci_emp_u, xlab = "Upper JP 95% PCI", ylab = "Upper EB 95% PCI"))
abline(a = 0, b = 1, lty = 2)
dev.off()

pdf("comparison_pci_interval.pdf", width = 12)
par(mfrow = c(1, 2))
plot(dat$county, dat$prop_approve, xlab = "County", ylab = "Sample approval proportion", main = "JP 95% PCIs", ylim = range(dat$ci_l, dat$ci_u))
arrows(x0 = dat$county, y0 = dat$ci_l, x1 = dat$county, y1 = dat$ci_u, code = 3, angle = 90)

plot(dat$county, dat$prop_approve, xlab = "County", ylab = "Sample approval proportion", main = "EB 95% PCIs", ylim = range(dat$ci_l, dat$ci_u))
arrows(x0 = dat$county, y0 = dat$ci_emp_l, x1 = dat$county, y1 = dat$ci_emp_u, code = 3, angle = 90)
dev.off()

# Empiric Bayes credible intervals don't seem to be symmetric -- to be expected?
# Some cases where variance is significantly underestimated on one end by Bayes,
# but similar in general.

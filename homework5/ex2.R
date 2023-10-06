library(MASS)
library(rjags)

data(Boston)

# (a) fit bayesian lm with uninformative normal priors for reg. coef.
# verify convergence of MCMC and summarise posterior dists. of all b_i's

n <- nrow(Boston)
p <- ncol(Boston)

X <- cbind(1, scale(Boston[,1:p-1])) # design matrix
Y <- Boston[,p]

data  <- list(X = X, Y = Y, n = n, p = p)

model_init <- textConnection("model{
  for (i in 1:n) {
    Y[i] ~ dnorm(inprod(X[i,], beta), taue)
  } 

  for (i in 1:p) {
    beta[i] ~ dnorm(0, 0.001)
  }

  taue ~ dgamma(0.1, 0.1)
}")

model <- jags.model(model_init, data = data, n.chains = 2, quiet = TRUE)
update(model, 1e+4)

samples <- coda.samples(model, variable.names = c("beta"), n.iter = 2e+4)

sum <- summary(samples)
names <- c("(Intercept)", colnames(Boston)[1:p - 1])
rownames(sum$statistics) <- names
rownames(sum$quantiles)  <- names

# drop intercept term as R linear model already has it
lin <- lm(Y ~ X[,-1])
names(lin$coef) <- names

# lin and our estimates coincide very well! let's have a look visually
# why are sd and se so similar between lin and blr? ask birgir
dat <- as.data.frame(cbind(0:13, sum$statistics[,c(1,2)], summary(lin)$coef[,c(1,2)]))
rownames(dat) <- NULL
colnames(dat) <- c("beta", "blr_beta", "blr_sd", "lm_beta", "lm_se")

sep <- 0.15
plot(dat$beta - sep, dat$blr_beta, pch = 21, xlab = expression(beta), ylab = "Estimate")
points(dat$beta + sep, dat$lm_beta, pch = 24)
arrows(x0 = dat$beta - sep, x1 = dat$beta - sep, y0 = dat$blr_beta - dat$blr_sd, y1 = dat$blr_beta + dat$blr_sd, code = 3, angle = 90, length = 0.05)
arrows(x0 = dat$beta + sep, x1 = dat$beta + sep, y0 = dat$lm_beta - dat$lm_se, y1 = dat$lm_beta + dat$lm_se, code = 3, angle = 90, length =0.05)

legend("topright", inset = c(0.05, 0.05), legend = c("BLR", "OLS"), pch = c(21, 24))


library(MASS)
library(rjags)

data(Boston)


# (A)
n <- nrow(Boston)
p <- ncol(Boston)

X <- cbind(1, scale(Boston[,1:p-1])) # design matrix
Y <- Boston[,p]

data  <- list(X = X, Y = Y, n = n, p = p)

model_init_uninf <- textConnection("model{
  for (i in 1:n) {
    Y[i] ~ dnorm(inprod(X[i,], beta), tau)
  } 

  for (i in 1:p) {
    beta[i] ~ dnorm(0, 0.0001)
  }

  tau ~ dgamma(0.1, 0.1)
}")

model_uninf <- jags.model(model_init_uninf, data = data, n.chains = 2, quiet = TRUE)
update(model_uninf, 1e+4)

samples_uninf <- coda.samples(model_uninf, variable.names = c("beta"), n.iter = 2e+4)

sum_uninf <- summary(samples_uninf)
names <- c("(Intercept)", colnames(Boston)[1:p - 1])
rownames(sum_uninf$statistics) <- names
rownames(sum_uninf$quantiles)  <- names




# (B)
lin <- lm(Y ~ X[,-1])
names(lin$coef) <- names

# lin and our estimates coincide very well! let's have a look visually
# why are sd and se so similar between lin and blr? ask birgir
dat <- as.data.frame(cbind(0:13, sum_uninf$statistics[,c(1,2)], summary(lin)$coef[,c(1,2)]))
rownames(dat) <- NULL
colnames(dat) <- c("beta", "blr_beta", "blr_sd", "lm_beta", "lm_se")

sep <- 0.15

pdf("ex2-estimates-blr-lm.pdf", width = 12)
plot(dat$beta - sep, dat$blr_beta, pch = 21, xlab = "Covariate", ylab = "Estimate", xaxt = "none")
axis(1, at=0:13, labels=names)
abline(h = 0, lty = 2)
points(dat$beta + sep, dat$lm_beta, pch = 24)
arrows(x0 = dat$beta - sep, x1 = dat$beta - sep, y0 = dat$blr_beta - dat$blr_sd, y1 = dat$blr_beta + dat$blr_sd, code = 3, angle = 90, length = 0.05)
arrows(x0 = dat$beta + sep, x1 = dat$beta + sep, y0 = dat$lm_beta - dat$lm_se, y1 = dat$lm_beta + dat$lm_se, code = 3, angle = 90, length = 0.05)
legend("topright", inset = c(0.05, 0.05), legend = c("BLR", "OLS"), pch = c(21, 24))
dev.off()




# (C)
model_init_lasso <- textConnection("model{
  for (i in 1:n) { 
    Y[i] ~ dnorm(inprod(X[i,], beta), taue) 
  } 

  beta[1] ~ dnorm(0, 0.0001)
  for (i in 2:p) {
    beta[i] ~ ddexp(0, taub * taue)
  }

  taue ~ dgamma(0.1, 0.1)
  taub ~ dgamma(0.1, 0.1)
}")

model_lasso <- jags.model(model_init_lasso, data = data, n.chains = 2, quiet = TRUE)
update(model_lasso, 1e+4)

samples_lasso <- coda.samples(model_lasso, variable.names = c("beta"), n.iter = 2e+4)

sum_lasso <- summary(samples_lasso)
rownames(sum_lasso$statistics) <- names
rownames(sum_lasso$quantiles)  <- names

dat <- as.data.frame(cbind(0:13, sum_uninf$statistics[,c(1,2)], sum_lasso$statistics[,c(1,2)]))
rownames(dat) <- NULL
colnames(dat) <- c("beta", "blr_beta", "blr_sd", "blasso_beta", "blasso_sd")

pdf("ex2-estimates-lasso-blr.pdf", width = 12)
plot(dat$beta - sep, dat$blr_beta, pch = 21, xlab = "Covariate", xaxt = "none", ylab = "Estimate")
axis(1, at=0:13, labels=names)
abline(h = 0, lty = 2)
points(dat$beta + sep, dat$blasso_beta, pch = 24)
arrows(x0 = dat$beta - sep, x1 = dat$beta - sep, y0 = dat$blr_beta - dat$blr_sd, y1 = dat$blr_beta + dat$blr_sd, code = 3, angle = 90, length = 0.05)
arrows(x0 = dat$beta + sep, x1 = dat$beta + sep, y0 = dat$blasso_beta - dat$blasso_sd, y1 = dat$blasso_beta + dat$blasso_sd, code = 3, angle = 90, length = 0.05)
legend("topright", inset = c(0.05, 0.05), legend = c("BLR", "BLASSO"), pch = c(21, 24))
dev.off()




# (D)

n_tr <- 500

X_tr <- X[1:n_tr,]
Y_tr <- Y[1:n_tr]

X_pr  <- X[(n_tr + 1):n,]
Y_pr  <- Y[(n_tr + 1):n]

n_pr <- nrow(X_pred)

data <- list(n_tr = n_tr, n_pr = n_pr, p = p, X_tr = X_tr, X_pr = X_pr, Y_tr = Y_tr)

model_init_pred <- textConnection("model{
  # Likelihood 
  for (i in 1:n_tr) {
    Y_tr[i]  ~ dnorm(inprod(X_tr[i,], beta), tau)
  } 
  
  # Priors
  for (i in 1:p) {
    beta[i] ~ dnorm(0, 0.0001)
  }

  tau ~ dgamma(0.1, 0.1)

  # Predictions
  for (i in 1:n_pr) {
    Y_pr[i] ~ dnorm(inprod(X_pr[i,], beta), tau)
  }
}")

model_pred <- jags.model(model_init_pred, data = data, n.chains = 2)

update(model_pred, 1e+4)

params <- c("beta", "Y_pr")

samples_pred <- coda.samples(model_pred, variable.names = params, n.iter = 2e+4)

summary(samples_pred)

par(mfrow = c(2, 3))

for (i in 1:6) {
  PPD <- samples_pred[[1]][,i]
  plot(density(PPD), main = NA)
  abline(v = Y_pr[i], lwd = 1.5, col = 2)
  abline(v = mean(PPD), lwd = 1.5, lty = 2)
}

#  _____________
# < Bayesically >
#  -------------
#         \   ^__^
#          \  (oo)\_______
#             (__)\       )\/\
#                 ||----w |
#                 ||     ||

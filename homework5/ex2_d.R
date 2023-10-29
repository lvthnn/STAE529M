library(MASS)
library(rjags)

data(Boston)

n <- nrow(Boston)
p <- ncol(Boston)

X <- cbind(1, scale(Boston[,1:p-1])) # design matrix
Y <- Boston[,p]

n_tr <- 500

X_tr <- X[1:n_tr,]
Y_tr <- Y[1:n_tr]

X_pr  <- X[(n_tr + 1):n,]
Y_pr  <- Y[(n_tr + 1):n]

n_pr <- nrow(X_pr)

data <- list(n_tr = n_tr, n_pr = n_pr, p = p, X_tr = X_tr, X_pr = X_pr, Y_tr = Y_tr)

model_init_pred_uninf <- textConnection("model{
  # Likelihood 
  for (i in 1:n_tr) {
    Y_tr[i]  ~ dnorm(inprod(X_tr[i,], beta[]), tau)
  } 

  # Priors
  for (i in 1:p) {
    beta[i] ~ dnorm(0, 0.0001)
  }

  tau ~ dgamma(0.1, 0.1)

  # Predictions
  for (i in 1:n_pr) {
    Y_pr[i] ~ dnorm(inprod(X_pr[i,], beta[]), tau)
  }

}")

model_pred_uninf <- jags.model(model_init_pred_uninf, data = data, n.chains = 2, quiet = TRUE)

update(model_pred_uninf, 1e+4, progress.bar = "none")

params <- c("beta", "Y_pr")

samples_pred_uninf <- coda.samples(model_pred_uninf, variable.names = params, n.iter = 2e+4, progress.bar = "none")

sum_pred_uninf <- summary(samples_pred_uninf)

# LASSO

model_init_pred_lasso <- textConnection("model{
  for (i in 1:n_tr) { 
    Y_tr[i] ~ dnorm(inprod(X_tr[i,], beta), taue) 
  } 


  # Predictions
  for (i in 1:n_pr) {
    Y_pr[i] ~ dnorm(inprod(X_pr[i,], beta), taue)
  }

  beta[1] ~ dnorm(0, 0.0001)
  for (i in 2:p) {
    beta[i] ~ ddexp(0, taub * taue)
  }

  taue ~ dgamma(0.1, 0.1)
  taub ~ dgamma(0.1, 0.1)

}")

model_pred_lasso <- jags.model(model_init_pred_lasso, data = data, n.chains = 2, quiet = TRUE)

update(model_pred_lasso, 1e+4, progress.bar = "none")

samples_pred_lasso <- coda.samples(model_pred_lasso, variable.names = params, n.iter = 2e+4,
                                   progress.bar = "none")

sum_pred_lasso <- summary(samples_pred_lasso)


pdf("ex2-ppd-uninf.pdf")
par(mfrow = c(2, 3))
for (i in 1:6) {
  PPD_uninf    <- samples_pred_uninf[[1]][,i] 
  PPD_CI_uninf <- sum_pred_uninf$quantiles[i, c(1, 5)]

  plot(density(PPD_uninf), main = NA, xlab = "medv")
  abline(v = Y_pr[i], col = 2, lwd = 1.5)
  abline(v = PPD_CI_uninf[1], lty = 2)
  abline(v = PPD_CI_uninf[2], lty = 2)
}
dev.off()

dat_tr <- cbind(X_tr, Y_tr)
dat_tr <- data.frame(dat_tr[,-1])

lin <- lm(Y_tr ~ ., data = dat_tr)

predict(lin, data.frame(X_pr[,-1]))

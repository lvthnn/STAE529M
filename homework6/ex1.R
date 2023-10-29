library(rjags)

# load data
tmp <- tempfile()
download.file(
  url = "https://www4.stat.ncsu.edu/~bjreich/BSMdata/guns.RData",
  destfile = tmp
)

load(tmp)
rm(tmp)

data_guns <- list(X = X, N = N, Y = Y, Z = Z)
n_reg     <- rowSums(data_guns$X)
X         <- cbind(1, data_guns$Z, n_reg)

n <- nrow(X)
p <- ncol(X)

data <- list(X = X, Y = Y, N = N, n = n, p = p)


# (a) fit poisson model


model_init <- textConnection("model{
  for (i in 1:n) {
    Y[i] ~ dpois(N[i] * lambda[i])
    lambda[i] <- exp(inprod(X[i,], beta))
  }

  for (j in 1:p) { beta[j] ~ dnorm(0, 0.0001) }
}")

model_poisson <- jags.model(model_init, data = data, n.chains = 2, quiet = TRUE)

update(model_poisson, 1e+4, progress.bar = "none")

samples_poisson <- coda.samples(model_poisson, variable.names = c("beta"), 
                                n.iter = 2e+4, progress.bar = "none")

sum_poisson <- summary(samples_poisson)

# (b) fit negative binomial model

model_init <- textConnection("model{
  for (i in 1:n) {
    Y[i]      ~ dnegbin(q[i], m)
    q[i]      <- m / (m + N[i] * lambda[i])
    lambda[i] <- exp(inprod(X[i,], beta))
  }

  for (j in 1:p) { beta[j] ~ dnorm(0, 0.0001) }
  m ~ dgamma(0.1, 0.1)
}")

model_negbin <- jags.model(model_init, data = data, n.chains = 2, quiet = TRUE)

update(model_negbin, 1e+4, progress.bar = "none")

samples_negbin <- coda.samples(model_negbin, variable.names = c("beta"),
                               n.iter = 2e+4, progress.bar = "none")

sum_negbin <- summary(samples_negbin)

# higher variance in negative binomial model. in, fact sqrt(sd) in the poisson
# model is approximately the sd in the negative binomial model. this is also true
# in the data:

mean(Y)
var(Y) # ~ E(Y)^2
sd(Y)

# some plots :D
posterior_means_poisson <- sum_poisson$statistics[,1]
posterior_means_negbin  <- sum_negbin$statistics[,1]

posterior_ci_poisson <- sum_poisson$quantiles[,c(1, 5)]
posterior_ci_negbin  <- sum_negbin$quantiles[,c(1, 5)]

pdf("ex1-comparison-poisson-negbin.pdf", width = 10)
sep <- 0.15
plot(1:7 - sep, posterior_means_poisson, pch = 1, xlim = c(1, 7), ylim = range(posterior_ci_poisson, posterior_ci_negbin), xaxt = "none",  xlab = "Covariate",
     ylab = "Posterior Mean")
axis(1, at = 1:7, colnames(X))
points(1:7 + sep, posterior_means_negbin, pch = 2)
arrows(x0 = 1:7 - sep, x1 = 1:7 - sep, y0 = posterior_ci_poisson[,1], y1 = posterior_ci_poisson[,2], code = 3, angle = 90, length = 0.15)
arrows(x0 = 1:7 + sep, x1 = 1:7 + sep, y0 = posterior_ci_negbin[,1],  y1 = posterior_ci_negbin[,2],  code = 3, angle = 90, length = 0.15)
abline(h = 0, lty = 2)
legend("bottomright", inset = c(0.05, 0.05), legend = c("Poisson", "Negative Binomial"), pch = 1:2)
dev.off()


# (c) 

# create prediction sets
Xn     <- X
Xn[,7] <- 0

Xf     <- X
Xf[,7] <- 25

data_ppd <- list(X = X, Xn = Xn, Xf = Xf, Y = Y, N = N, n = n, p = p)

model_init <- textConnection("model{

  for (i in 1:n) {
    Y[i] ~ dpois(N[i] * lambda[i])
    lambda[i] <- exp(inprod(X[i,], beta))
  }

  for (j in 1:p) { beta[j] ~ dnorm(0, 0.0001) }

  for (i in 1:n) { 
    Yn[i] ~ dpois(N[i] * lambda_n[i])
    lambda_n[i] <- exp(inprod(Xn[i,], beta))

    Yf[i] ~ dpois(N[i] * lambda_f[i])
    lambda_f[i] <- exp(inprod(Xf[i,], beta))
  }

}")

model_poisson <- jags.model(model_init, data = data_ppd, n.chains = 2, quiet = TRUE)

update(model_poisson, 1e+4, progress.bar = "none")

samples_poisson <- coda.samples(model_poisson, variable.names = c("beta", "Yf", "Yn", "Y"), 
                                n.iter = 2e+4, progress.bar = "none")

sum_homicide <- summary(samples_poisson)

pdf("lol.pdf", width = 20, height = 10)
par(mfrow = c(5, 10))
for (i in 1:50) {
  Yi  <- samples_poisson[[1]][,i]
  Yif <- samples_poisson[[1]][,50+i]
  Yin <- samples_poisson[[1]][,100+i]

  xl <- range(Yin, Yif)
 
  plot(density(Yin), xlim = xl, col = "red", lwd = 1.5, main = NA, xlab = "No. homicides")
  lines(density(Yif), lty = 1, col = "blue", lwd = 1.5)
  lines(density(Yi), lty = 2)

  abline(v = Y[i])
}
dev.off()

# Problematic fit in the model -- no gun law regulations seems to
# decrease gun homicides in some cases?

# Hmm. Look at code more.

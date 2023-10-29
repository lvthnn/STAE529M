library(rjags)
library(scales)

set.seed(42)

# load in data
load("gambia.rda")
plot(gambia.borders, type = "l", axes = F, xlab = "", ylab = "")


# set up data set
X <- cbind(
  1,
  scale(gambia$x),
  scale(gambia$y),
  scale(gambia$age),
  scale(gambia$green),
  gambia[,c("netuse", "treated", "phc")]
)

colnames(X) <- c("intercept", "x", "y", "age", "green", "netuse", "treated", "phc")
Y <- gambia[,"pos"]

loc <- as.numeric(factor(paste0(X$x, "|", X$y), labels = 1:65))
n <- nrow(X)
p <- ncol(X)
m <- length(unique(loc))
data <- list(X = X, Y = Y, loc = loc, n = n, p = p, m = m)

# logistic regression model
model_init <- textConnection("model{
  for (i in 1:n) {
    Y[i] ~ dbern(pr[i])
    logit(pr[i]) <- inprod(X[i,], beta)
  }

  for (j in 1:p) { beta[j] ~ dnorm(0, 0.01) }
}")

model <- jags.model(model_init, data = data, n.chains = 2, quiet = TRUE)

# train model
update(model, n.iter = 1e+4)
samples     <- coda.samples(model, variable.names = c("beta"), n.chains = 2, 
                            n.iter = 2e+4)
sum_samples <- summary(samples)

# verify convergence
gelman.diag(samples)

# draw some pretty significance pictures
pdf("ex3-significance-covariaties-non-random.pdf", width = 12)
par(mfrow = c(1, 2))
posterior_means <- sum_samples$statistics[,1]
posterior_cis   <- sum_samples$quantiles[,c(1,5)]

plot(posterior_means, xaxt = "none", ylim = range(posterior_cis),
     xlab = "Covariate", ylab = "Posterior Mean")
axis(1, at = 1:8, labels = colnames(X))
arrows(x0 = 1:8, x1 = 1:8, y0 = posterior_cis[,1], y1 = posterior_cis[,2],
       code = 3, angle = 90)
abline(h = 0, lty = 2)

vals <- samples[[1]][,1:8]
plot(density(vals[,1]), xlim = range(vals), ylim = c(0, 8), col = 1, lty = 1, main = NA,
     xlab = "Estimate", "Density")
for (i in 2:8) { lines(density(vals[,i]), lty = i, col = i) }
abline(v = 0, lty = 2)
legend("topleft", inset = c(0.05, 0.05), legend = colnames(X), lty = 1:8, col = 1:8)
dev.off()


# summarise odds ratios

# scale back
qnt      <- sum_samples$quantiles
sd_x     <- sd(gambia$x)
sd_y     <- sd(gambia$y)
sd_age   <- sd(gambia$age)
sd_green <- sd(gambia$green)

rbind(
  exp(qnt[2:5,] / c(sd_x, sd_y, sd_age/10, sd_green)),
  exp(qnt[6:8,])
)

pdf("ex3-jags-norandomeffects.pdf")
plot(samples)
dev.off()

# fit random effect model

model_re <- textConnection("model{
  for (i in 1:n) {
    Y[i] ~ dbern(q[i])
    logit(q[i]) <- inprod(X[i,], beta) + RE[loc[i]]
  }
 
  for (j in 1:m) { RE[j]   ~ dnorm(0, tau) }
  for (j in 1:p) { beta[j] ~ dnorm(0, 0.0001) }

  tau ~ dgamma(0.1, 0.1)
}")

model_REs <- jags.model(model_re, data = data, n.chains = 2, quiet = TRUE)

update(model_REs, n.iter = 1e+4)

params <- c("beta", "RE", "tau")

samples_re     <- coda.samples(model_REs, variable.names = params, n.iter = 2e+4)
sum_samples_re <- summary(samples_re)

# verify convergence
gelman.diag(samples_re)

sum_samples_re <- summary(samples_re)
rownames(sum_samples_re$statistics)[66:73] <- colnames(X)
rownames(sum_samples_re$quantiles)[66:73] <- colnames(X)

# compare covariate effects
effects_no_re <- sum_samples$statistics[,1]
effects_re    <- sum_samples_re$statistics[66:73,1]

cis_no_re <- sum_samples$quantiles[,c(1,5)]
cis_re <- sum_samples_re$quantiles[66:73,c(1,5)]

sep <- 0.15

pdf("ex3-comparison-re-no-re.pdf", width = 10)
plot(1:8 - sep, effects_no_re, xaxt = "none", pch = 1, ylim = range(cis_no_re, cis_re), xlim = c(1, 8),
     xlab = "Covariate", ylab = "Estimate")
points(1:8 + sep, effects_re, pch = 2)
arrows(x0 = 1:8 - sep, x1 = 1:8 - sep, y0 = cis_no_re[,1], y1 = cis_no_re[,2], code = 3, angle = 90, length = 0.15)
arrows(x0 = 1:8 + sep, x1 = 1:8 + sep, y0 = cis_re[,1], y1 = cis_re[,2], code = 3, angle = 90, length = 0.15)
abline(h = 0, lty = 2)
legend("topright", inset = c(0.05, 0.05), legend = c("No random effects", "Random effects"), pch = c(1, 2))
axis(1, at = 1:8, colnames(X))
dev.off()

RE_post_means <- sum_samples_re$statistics[1:65,1]

spatial <- data.frame(1:65, unique(gambia[,c("x", "y")]), RE_post_means)
colnames(spatial) <- c("loc", "x", "y", "re")

pdf("ex3-gambia-regional-plot.pdf", width = 12)
par(mfrow = c(1, 2))
plot(gambia.borders, type = 'l', axes = F, xlab = NA, ylab = NA)
points(spatial$x, spatial$y, cex = 2.5 * exp(spatial$re / max(spatial$re)), pch = 21, col = "red",
       bg = alpha("red", 0.4))

plot(spatial$x, spatial$y,cex = 2.5 * exp(spatial$re / max(spatial$re)), pch = 21, col = "red", 
     bg = alpha("red", 0.4), xlab = "x-coordinate", ylab = "y-coordinate", xlim = range(gambia.borders$x, na.rm = TRUE),
     ylim = range(gambia.borders$y, na.rm = TRUE))
dev.off()

# load in samples to avoid regeneration

save(samples, samples_re, file = "jags-samples.rda")


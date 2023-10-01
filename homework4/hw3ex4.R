S <- 1e+6
Y <- 1:10
n <- a <- 10

theta           <- matrix(NA, nrow = S, ncol = 11)
theta[1,]       <- c(rep(var(Y), 10), 1)
colnames(theta) <- c(paste0("s_", 1:n), "b")


for (s in 2:S) {
  vals <- c()

  for (i in 1:10) {
    A <- a + 1/2
    B <- Y[i]^2 + theta[s-1,"b"]

    si <- rinvgamma(1, A, B)
    theta[s,i] <- si
  }

  A <- n*a + 1
  B <- sum(1/theta[s,1:10]) + 1

  b <- rgamma(1, A, B)
  theta[s,11] <- b
}

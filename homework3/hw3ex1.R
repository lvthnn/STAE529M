# 2.4.6

library(invgamma)

n <- 20
ssy <- 15
sy <- -2
k <- 15

as <- c(0.1, 1.0)
bs <- c(0.1, 1.0)
cs <- rep(1:k, each = 2)

params <- cbind(data.frame(as, bs), cs)
colnames(params) <- c("a", "b", "c")

probs <- sapply(1:nrow(params), function(i) {
  a <- params[i,"a"]
  b <- params[i,"b"]
  c <- params[i,"c"]

  pinvgamma(q = c^2, shape = 0.5 * n + a, rate = 0.5 * 15 + b, lower.tail = FALSE)
})

sens <- sapply(seq(1, nrow(params), by = 2), function(i) probs[i]/probs[i + 1])

# probs[1]/probs[2] = 0.9507509
# probs[3]/probs[4] = 1.7706740

# Results with higher values of c are more sensitive to the prior 
# see by increasing value of k in code


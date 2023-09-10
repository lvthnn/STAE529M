# Exercise 2.4.2

n1 <- 2820
n2 <- 27

y1 <- 563
y2 <- 10

a <- b <- 0.1

s <- 1e+5
lambda1 <- rgamma(s, shape = y1 + a, rate = n1 + b)
lambda2 <- rgamma(s, shape = y2 + a, rate = n2 + b)

plot(density(lambda1))
lines(density(lambda2))

print(mean(lambda2 > lambda1))

# Visualizations for report

xs <- seq(0, 1, by = 0.0001)
p1 <- dgamma(xs, shape = y1 + a, rate = n1 + b)
p2 <- dgamma(xs, shape = y2 + a, rate = n2 + b)

plot(xs, p1, type = 'l')

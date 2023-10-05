library(MASS)
library(rjags)

data(Boston)

head(Boston, 10)

# (a) fit bayesian lm with uninformative normal priors for reg. coef.
# verify convergence of MCMC and summarise posterior dists. of all b_i's



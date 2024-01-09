# Final exam fall 2023 - STÆ529M Bayesian Data Analysis
#
# Problem 5
# Load the data
data5 <- read.table("problem5.txt", header=T)
id <- data5[,1] 
time <- data5[,2]
height <- data5[,3]

T <- 9
n <- 26

X <- matrix(nrow = n, ncol = T)
Y <- matrix(nrow = n, ncol = T)

for(i in 1:n){
  for(j in 1:T){
    k <- j + T*(i - 1)
    X[i,j] <- time[k]  
    Y[i,j] <- height[k]
  }
} 


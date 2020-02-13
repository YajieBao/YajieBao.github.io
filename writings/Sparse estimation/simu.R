########################################
##### Sparse mean vector estimation ####
########################################

### normal case ###

# covariance matrix
p = 2500
A <- matrix(runif(p^2)*2-1, ncol = p) 
Sigma <- t(A) %*% A

# sample generation
library(MASS)
n= 500
s = 30 #sparsity
beta = c(runif(s, min = 1, max = 50), rep(0, p - s)) #true mean vector
beta = sample(beta)
X = mvrnorm(n, mu = beta, Sigma = Sigma)

# mean estimation
delta = 0.01
thres = max(sqrt(diag(Sigma))) * sqrt(2 * log(2 * p / delta) / n)
mean_1 = apply(X, 2, mean)
mean_t = sapply(mean_1, function(x) ifelse(abs(x) >= thres, x, 0))

beta[setdiff(which(beta != 0), which(mean_t != 0))]
all(which(mean_t != 0) %in% which(beta != 0))
min(beta[which(beta!=0)]) #less than thresholding

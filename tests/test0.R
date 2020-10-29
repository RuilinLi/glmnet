library(glmnet)
library(pgenlibr)
library(profvis)
rm(list=ls())
set.seed(1)
n = 1200
p = 1200
sample_subset = 1:n
vsubset = 1:p
sample_subset = sample_subset[-15]
vsubset = vsubset[-9]
n = n - 1
p = p - 1

beta = rbinom(p,1,0.3) * rnorm(p)

pgen <- pgenlibr::NewPgen("data/toy_data.pgen", pvar = NULL, 
                          sample_subset =sample_subset)
X <- pgenlibr::ReadList(pgen, vsubset, meanimpute=T)
y = X %*% beta 
y = y/sd(y)

# Reference
profvis({
fit_ref = glmnet(X, y, family = 'gaussian',  standardize = T, intercept = T)
})

# Dense
profvis({
  fit_dense = glmnet(X, y, family = gaussian(), standardize = T, intercept = T, thresh = 1e-7)
  
})


y2 = rep(1,n)
y2[y < median(y)] = 0
fit_dense = glmnet(X, y2, family = "binomial", standardize = T, intercept = T, thresh=1e-7)
fit_dense2 = glmnet(X, y2, family = binomial(), standardize = T, intercept = T, thresh=1e-7)

# Sparse
Xs = as(X, 'sparseMatrix')
fit_sparse = glmnet(Xs, y, family = gaussian(), standardize = T, intercept = T)


# Plink Format
Xp = glmnet::PlinkMatrix("data/toy_data.pgen", sample_subset, vsubset)


fit_plink = glmnet(Xp, y, family = gaussian(),  standardize = T, intercept = T)

# Test matrix-vector multiplication
# Notice that, after calling glmnet on Xp with standardize and intercept
# the value has been modified
y2 = rnorm(p)
z = rnorm(n)

X = center(X, rep(1.0, n))
X = X$x
X = standardize(X,rep(1.0, n))
X = X$x
max(abs(drop(X%*% y2) - Xp %*% y2))
max(abs(drop(z %*% X) - z %*% Xp))



me_dense = max(abs(fit_dense$beta -  fit_ref$beta))
me_sparse = max(abs(fit_sparse$beta -  fit_ref$beta))
me_plink = max(abs(fit_plink$beta - fit_ref$beta))




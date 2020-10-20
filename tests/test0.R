library(glmnet)
library(pgenlibr)
rm(list=ls())
set.seed(1)
n = 100
p = 100
sample_subset = 1:n
vsubset = 1:p
sample_subset = sample_subset[-15]
vsubset = vsubset[-9]
n = n - 1
p = p - 1

beta = rbinom(p,1,0.3) * rnorm(p)

pgen <- pgenlibr::NewPgen("/Users/ruilinli/plink-ng/toy_data.pgen", pvar = NULL, 
                          sample_subset =sample_subset)
X <- pgenlibr::ReadList(pgen, vsubset, meanimpute=F)
y = X %*% beta 

lam = 0.12190
# Reference
fit_ref = glmnet(X, y, family = 'gaussian', lambda=lam, standardize = T, intercept = T)

# Dense
#fit_dense = glmnet(X, y, family = gaussian(), lambda=lam, standardize = T, intercept = T)
fit_dense = glmnet(X, y, family = gaussian(), standardize = T, intercept = T)

# Sparse
Xs = as(X, 'sparseMatrix')
fit_sparse = glmnet(Xs, y, family = gaussian(), lambda=lam, standardize = T, intercept = T)


# Plink Format
a = glmnet::PlinkMatrix("/Users/ruilinli/plink-ng/toy_data.pgen", sample_subset, vsubset)


fit_plink = glmnet(a, y, family = gaussian(),  standardize = T, intercept = T)

y2 = rnorm(p)
z = rnorm(n)

X = center(X, rep(1.0, n))
X = X$x
X = standardize(X,rep(1.0, n))
X = X$x
max(abs(drop(X%*% y2) - a %*% y2))
max(abs(drop(z %*% X) - z %*% a))



me_dense = max(abs( fit_dense$beta -  fit_ref$beta))
me_sparse = max(abs(fit_sparse$beta -  fit_ref$beta))
me_plink = max(abs(fit_plink$beta - fit_ref$beta))




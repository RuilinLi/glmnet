library(glmnet)
set.seed(1)
n = 800
p = 400
X = matrix(rnorm(n*p),nrow = n, ncol = p)
beta = rbinom(p,1,0.3) * rnorm(p)
y = X %*% beta 

fit = glmnet(X, y, family = gaussian())


# fit2 = glmnet(X, y, family = 'gaussian', lambda=0.02)
# 
# 
# diff = abs(fit$beta - fit2$beta)
# denom = abs(beta)
# diff[which(denom != 0)]/ denom[which(denom != 0)]

.Call('wls_dense2', alm0,almc,alpha,m,nobs,nvars,
      x,r,v,intr,ju,vp,cl,nx,thr,
      maxit,a,aint,g,ia,iy,iz,mm,
      nino,rsqc,nlp,jerr)

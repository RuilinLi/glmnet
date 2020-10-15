library(glmnet)
library(pgenlibr)
rm(list=ls())
n = 100
p = 200
set.seed(1)
sample_subset = 1:n
vsubset = 1:p
beta = rbinom(p,1,0.3) * rnorm(p)

typeof(vsubset)

pgen <- pgenlibr::NewPgen("/Users/ruilinli/plink-ng/toy_data.pgen", pvar = NULL, 
                          sample_subset =sample_subset)

X <- pgenlibr::ReadList(pgen, vsubset, meanimpute=F)
y = X %*% beta 
# Dense example
fit = glmnet(X, y, family = 'gaussian', lambda=0.12190, standardize = F, intercept = F)

# plink example

x2 = list()
x2[[1]] = "/Users/ruilinli/plink-ng/toy_data.pgen"
x2[[2]] = sample_subset
x2[[3]] = vsubset
class(x2) <- c(class(x2), 'plink2')
a=wls_plink_cpp(alm0=alm0,almc=almc,alpha=alpha,m=m,no=nobs,ni=nvars,
               x=x2,r=r,v=v,intr=intr,ju=ju,vp=vp,cl=cl,nx=nx,thr=thr,
               maxit=maxit,a=a,aint=aint,g=g,ia=ia,iy=iy,iz=iz,mm=mm,
               nino=nino,rsqc=rsqc,nlp=nlp,jerr=jerr)


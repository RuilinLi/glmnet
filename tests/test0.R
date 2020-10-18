library(glmnet)
library(pgenlibr)
rm(list=ls())
n = 300
p = 200
set.seed(1)
sample_subset = 1:n
#sample_subset = sample_subset[-3]
vsubset = 1:p
#vsubset = vsubset[-12]
n = n-1
p=p-1
beta = rbinom(p,1,0.3) * rnorm(p)
res = glmnet::xptrtest("/Users/ruilinli/plink-ng/toy_data.pgen", sample_subset, vsubset)

typeof(vsubset)

pgen <- pgenlibr::NewPgen("/Users/ruilinli/plink-ng/toy_data.pgen", pvar = NULL, 
                          sample_subset =sample_subset)

X <- pgenlibr::ReadList(pgen, vsubset, meanimpute=F)
y = X %*% beta 
# Dense example\
#fit = glmnet(X, y, family = 'gaussian', lambda=0.12190, standardize = F, intercept = F)

fit = glmnet(X, y, family = gaussian(), lambda=0.12190, standardize = F, intercept = F)

# plink example

x2 = list()
x2[[1]] = "/Users/ruilinli/plink-ng/toy_data.pgen"
x2[[2]] = sample_subset
x2[[3]] = vsubset
a=wls_plink_cpp(alm0=alm0,almc=almc,alpha=alpha,m=m,no=nobs,ni=nvars,
               x=x2,r=r,v=v,intr=intr,ju=ju,vp=vp,cl=cl,nx=nx,thr=thr,
               maxit=maxit,a=a,aint=aint,g=g,ia=ia,iy=iy,iz=iz,mm=mm,
               nino=nino,rsqc=rsqc,nlp=nlp,jerr=jerr)

glmnet::xptrtest("/Users/ruilinli/plink-ng/toy_data.pgen", sample_subset, vsubset)


setClass("PlinkMatrix", representation(samples = "integer", variants="integer", fname="character"),
         contains = "Matrix")
PlinkMatrix <- function(fname, samples, variants)
{
  samples =as.integer(sort(unique(samples)))
  new("PlinkMatrix", samples=samples, variants = variants, fname=fname, Dim=c(length(samples),length(variants)))
}
a = glmnet::PlinkMatrix("/Users/ruilinli/plink-ng/toy_data.pgen", 1:100, 1:200)



setGeneric("rowmax", function(object) {
  standardGeneric("rowmax")
})

setMethod("rowmax", signature(object = "matrix"), function(object) {
  apply(object, 2, max)
})
x = matrix(rnorm(20),4,5)

rowmax(x)

a = glmnet::PlinkMatrix("/Users/ruilinli/plink-ng/toy_data.pgen", 1:300, 1:200)



# w = rep(1/300, 300)
# 
# apply(X, 2, function(r) sqrt(weighted.mean(r^2, w) -
#                                weighted.mean(r, w)^2))

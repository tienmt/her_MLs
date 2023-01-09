

library(glmnet)
library(HiLMM)
require(doMC)
registerDoMC(cores=10)
library(RSpectra)
noise = 1

load('/data2/thetm/MA/new.masnps.rda')
masnps.filter = masnps.filter/2

p = ncol(x.penici)
n = nrow(x.penici)
source('/data2/thetm/review_HERI/mainfunctions.R')

load('/data2/thetm/review_HERI/penicilin_genes.rda')
x.penici = masnps.filter[,is.element(colnames(masnps.filter),penicilin.genes)]
load('/data2/thetm/MA/ypenincillin.rda')


# Enet method
glmnetcv = cv.glmnet(x.penici, ypenicilin, alpha = 0.01,parallel = T)
tam = predict(glmnetcv, type = 'nonzero', s = 'lambda.min')[[1]] 
enet = coef( glmnetcv,  s='lambda.min')[-1]
enet[tam] %*% cov2(x.penici[,tam]) %*% enet[tam]/ var(ypenicilin)


# convex optimization method
EigenPrism(ypenicilin, x.penici, target = 'heritability')

# MLE method
M = tcrossprod(x.penici)
O = eigen(M)$vectors
lambda = eigen(M)$values
yTilde = t(O) %*% ypenicilin
mle.r =ridge.mle(yTilde,lambda,1.e-6)
eta_chap = 1 - mle.r$sig2hat / var(ypenicilin)
s_eta = mle.r$sig2hat/sqrt(2)
eta_chap + qnorm(1-0.05/2)* s_eta*c(-1,1) 


# moment method
mm2dicker(ypenicilin, x.penici,CI = T )

# scaled lasso method
scaled_lasso <- scaled.lasso(x.penici , ypenicilin)
Slasso = 1-scaled_lasso$sigma2^2/var(ypenicilin)
k = scaled_lasso$non.zeros
Slasso- c(-1,1)*log(.5)*(k*log(p)/n + 1/sqrt(n))


# boostHER
xdat2 = x.penici 
for (itrrs in 1:(log2(p/n)+1) ) {
  model <- cv.glmnet(xdat2 , ypenicilin , alpha=0, parallel = T)
  coef0 <- abs(predict(model,type="coefficients",s='lambda.min')[-1] )
  xdat2 <- xdat2[, coef0 > quantile(coef0,.5)]
  if(ncol(xdat2) < n ) break
}
outlist <- foreach(jj = seq(100) ) %dopar%  {
  foreach(i = seq(2), .packages=c("glmnet")) %dopar%  {
    which <- sample(rep(seq(2), length = n)) == i
    xols1 = xdat2[, predict(cv.glmnet(xdat2[which,], ypenicilin[which], parallel = T,alpha = 1),
                            type = 'nonzero', s = 'lambda.min')[[1]] ]
    a = NA
    tryCatch({a <-(summary(lm(ypenicilin[!which] ~ xols1[!which,] ))$sigma)^2 
    },  error=function(e){ } )
  }}
1 - mean( unlist(outlist), na.rm=T)/var(ypenicilin)
quantile(1 - unlist(outlist)/var(ypenicilin), probs = c(0.025,0.975))
sd(1 -  unlist(outlist), na.rm=T)/var(ypenicilin)






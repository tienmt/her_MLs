

library(glmnet)
library(HiLMM)
require(doMC)
registerDoMC(cores=10)
library(RSpectra)
noise = 1

load('/data2/thetm/meala/newsnps.rda')
load('/data2/thetm/heritability-LM/zone_size_realdata/zone_size_dat.rda')
load('/data2/thetm/review_HERI/trimethoprim_causal_gene.rda')
source('/data2/thetm/review_HERI/mainfunctions.R')
Trimethoprim_causal_gene = paste('coor',sep='_',Trimethoprim_causal_gene)
xpeni = mesnp[new.sam,is.element(colnames(mesnp),tetracylin_causal_gene)]

p = ncol(xpeni)
n = nrow(xpeni)

# Enet method
glmnetcv = cv.glmnet(xpeni, zone_size_trimetho, alpha = 0.01,parallel = T)
tam = predict(glmnetcv, type = 'nonzero', s = 'lambda.min')[[1]] 
enet = coef( glmnetcv,  s='lambda.min')[-1]
enet[tam] %*% cov2(xpeni[,tam]) %*% enet[tam]/ var(zone_size_trimetho)


# convex optimization method
EigenPrism(zone_size_trimetho, xpeni, target = 'heritability')

# MLE method
M = tcrossprod(xpeni)
O = eigen(M)$vectors
lambda = eigen(M)$values
yTilde = t(O) %*% zone_size_trimetho
mle.r =ridge.mle(yTilde,lambda,1.e-6)
1 - mle.r$sig2hat / var(zone_size_trimetho)
eta_chap = 1 -mle.r$sig2hat / var(zone_size_trimetho)
eta_chap + qnorm(1-0.05/2)*c(-1,1) *mle.r$sig2hat/sqrt(n)


# moment method
mm2dicker(zone_size_trimetho, xpeni,CI = T )

# scaled lasso method
scaled_lasso <- scaled.lasso(xpeni , zone_size_trimetho)
Slasso = 1-scaled_lasso$sigma2^2/var(zone_size_trimetho)
k = scaled_lasso$non.zeros
Slasso- c(-1,1)*log(.5)*(k*log(p)/n + 1/sqrt(n))
Slasso

# boostHER
xdat2 = xpeni 
for (itrrs in 1:(log2(p/n)+1) ) {
  if(ncol(xdat2) < n ) break
  model <- cv.glmnet(xdat2 , zone_size_trimetho , alpha=0, parallel = T)
  coef0 <- abs(predict(model,type="coefficients",s='lambda.min')[-1] )
  xdat2 <- xdat2[, coef0 > quantile(coef0,.5)]
}
#xdat2 = xpeni 
outlist <- foreach(jj = seq(100) ) %dopar%  {
  foreach(i = seq(2), .packages=c("glmnet")) %dopar%  {
    which <- sample(rep(seq(2), length = n)) == i
    xols1 = xdat2[, predict(cv.glmnet(xdat2[which,], zone_size_trimetho[which], parallel = T,alpha = 1),
                            type = 'nonzero', s = 'lambda.min')[[1]] ]
    a = NA
    tryCatch({a <-(summary(lm(zone_size_trimetho[!which] ~ xols1[!which,] ))$sigma)^2 
    },  error=function(e){ } )
  }}
1 - mean( unlist(outlist), na.rm=T)/var(zone_size_trimetho)
quantile(1 - unlist(outlist)/var(zone_size_trimetho), probs = c(0.025,0.975))
sd(1 -  unlist(outlist), na.rm=T)/var(zone_size_trimetho)


outlist <- foreach(jj = seq(100) ) %dopar%  {
  foreach(i = seq(2) ) %dopar%  {
    which <- sample(rep(seq(2), length = n)) == i
   a = NA
    tryCatch({a <-(summary(lm(zone_size_trimetho[!which] ~ xdat2[!which,] ))$sigma)^2 
    },  error=function(e){ } )
  }}
1 - mean( unlist(outlist), na.rm=T)/var(zone_size_trimetho)
quantile(1 - unlist(outlist)/var(zone_size_trimetho), probs = c(0.025,0.975))







library(glmnet)
require(doMC)
registerDoMC(cores=10)

load('/data2/thetm/meala/newsnps.rda')
load('/data2/thetm/review_HERI/penicilin_genes.rda')
source('/data2/thetm/review_HERI/mainfunctions.R')

tam = colnames(mesnp)[!is.element(colnames(mesnp),penicilin.genes)]
tam1 = sample(tam,5000 - 827)
x = mesnp[,is.element(colnames(mesnp),penicilin.genes)]
x = cbind(x,mesnp[,tam1])

peni.gene = colnames(mesnp)[is.element(colnames(mesnp),penicilin.genes)]

p = ncol(x)
n = nrow(x)
noise = 1
h2 = 0.8
beta0 = rep(0,p)
names(beta0) = colnames(x)
beta0[peni.gene] = rnorm(length(peni.gene))
beta0 <- beta0 * sqrt((noise^2)*h2/(1-h2) / as.numeric(beta0[peni.gene]%*%cov2(x[,peni.gene])%*%beta0[peni.gene]) )
y = c(x%*%beta0 + rnorm(n)*noise)
1 - noise^2/var(y)

oracle = elasnet = eprism = mle = mm = slasso = bher = c()
for(ss in 1:50){
  y = c(x%*%beta0 + rnorm(n)*noise)
  #oracle
  oracle[ss] = 1 - noise^2/var(y)
  
  # Enet method
  glmnetcv = cv.glmnet(x, y, alpha = 0.01,parallel = T)
  tam = predict(glmnetcv, type = 'nonzero', s = 'lambda.min')[[1]] 
  enet = coef( glmnetcv,  s='lambda.min')[-1]
  elasnet[ss] = enet[tam] %*% cov2(x[,tam]) %*% enet[tam]/ var(y)
  
  # convex optimization method
  eprism[ss] = EigenPrism(y, x, target = 'heritability')$estimate
  
  # MLE method
  M = tcrossprod(x)
  O = eigen(M)$vectors
  lambda = eigen(M)$values
  yTilde = t(O) %*% y
  mle.r = ridge.mle(yTilde,lambda,1.e-6)
  mle[ss] = 1 - mle.r$sig2hat / var(y)
  
  # moment method
  mm[ss] = mm2dicker(y, x,CI = T )$est
  
  # scaled lasso method
  scaled_lasso <- scaled.lasso(x , y)
  slasso[ss] = 1-scaled_lasso$sigma2^2/var(y)

  # boostHER
  xdat2 = x 
  for (itrrs in 1:(log2(p/n)+1) ) {
    model <- cv.glmnet(xdat2 , y , alpha=0, parallel = T)
    model <- abs(predict(model,type="coefficients",s='lambda.min')[-1] )
    xdat2 <- xdat2[, model > quantile(model,.5)]
    if(ncol(xdat2) < n ) break
  }
  outlist <- foreach(jj = seq(50) ) %dopar%  {
    foreach(i = seq(2), .packages=c("glmnet")) %dopar%  {
      which <- sample(rep(seq(2), length = n)) == i
      xols1 = xdat2[, predict(cv.glmnet(xdat2[which,], y[which], parallel = T,alpha = .01),
                              type = 'nonzero', s = 'lambda.min')[[1]] ]
      a = NA
      tryCatch({a <-(summary(lm(y[!which] ~ xols1[!which,] ))$sigma)^2 
      },  error=function(e){ } )
    }}
  bher[ss] = 1 - mean( unlist(outlist), na.rm=T)/var(y)
  print(ss)
}








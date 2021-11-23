
library(glmnet)
require(doMC)
registerDoMC(cores=10)

load('/data2/thetm/meala/newsnps.rda')
load('/data2/thetm/review_HERI/penicilin_genes.rda')
source('/data2/thetm/review_HERI/mainfunctions.R')

tam = colnames(mesnp)[!is.element(colnames(mesnp),penicilin.genes)]
tam1 = sample(tam,5000 - 827)
x = mesnp[,is.element(colnames(mesnp),penicilin.genes)]
x1 = cbind(x,mesnp[,tam1])

peni.gene = colnames(mesnp)[is.element(colnames(mesnp),penicilin.genes)]

p = ncol(x1)
n = nrow(x1)
noise = 1
h2 = 0.8
beta0 = rep(0,p)
names(beta0) = colnames(x)
beta0[peni.gene] = rnorm(length(peni.gene))
beta0 <- beta0 * sqrt((noise^2)*h2/(1-h2) / as.numeric(beta0[peni.gene]%*%cov2(x[,peni.gene])%*%beta0[peni.gene]) )
y1 = c(x1%*%beta0 + rnorm(nrow(x1))*noise)
rd_samples = sample(n,500)
1 - noise^2/var(y1[rd_samples])

oracle = elasnet = eprism = mle = mm = slasso = bher = c()
for(ss in 1:50){
  rd_samples = sample(nrow(x1),2000)
  y = y1[rd_samples]
  x = x1[rd_samples,]
  p = ncol(x)
  n = nrow(x)
  #oracle
  oracle[ss] = 1 - noise^2/var(y)
  
  # Enet method
  glmnetcv = cv.glmnet(x, y, alpha = 0.01,parallel = T)
  tam = predict(glmnetcv, type = 'nonzero', s = 'lambda.min')[[1]] 
  enet = coef( glmnetcv,  s='lambda.min')[-1]
  elasnet[ss] = enet[tam] %*% cov2(x[,tam]) %*% enet[tam]/ var(y)
  
  # convex optimization method
  tryCatch({
    eprism[ss] = EigenPrism(y, x, target = 'heritability')$estimate
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  
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
    coef0 <- abs(predict(model,type="coefficients",s='lambda.min')[-1] )
    xdat2 <- xdat2[, coef0 > quantile(coef0,.5)]
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








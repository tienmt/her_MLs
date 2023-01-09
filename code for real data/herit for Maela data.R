
library(glmnet)
library(HiLMM)
require(doMC)
registerDoMC(cores=10)
library(RSpectra)
noise = 1

load('/data2/thetm/meala/newsnps.rda')
load('/data2/thetm/heritability-LM/zone_size_realdata/zone_size_dat.rda')

p = ncol(mesnp)
n = nrow(mesnp)
source('/data2/thetm/review_HERI/mainfunctions.R')
zsmaela = mesnp[new.sam,]

# convex optimization method
start_time <- Sys.time()
eprsm <- EigenPrism(zone_size_trimetho, zsmaela, target = 'heritability')
eprsm
end_time <- Sys.time()
end_time - start_time

# MLE method
start_time <- Sys.time()
M = tcrossprod(zsmaela)
O = eigen(M)$vectors
lambda = eigen(M)$values
yTilde = t(O) %*% zone_size_trimetho
mle.r =ridge.mle(yTilde,lambda,1.e-6)
eta_chap = 1 - mle.r$sig2hat / var(zone_size_trimetho)
eta_chap
eta_chap + qnorm(1-0.05/2)* c(-1,1)/sqrt(2*n)
end_time <- Sys.time()
end_time - start_time

# moment method
start_time <- Sys.time()
mmdker = mm2dicker(zone_size_trimetho, zsmaela,CI = T )
mmdker
end_time <- Sys.time()
end_time - start_time

# scaled lasso method
start_time <- Sys.time()
scaled_lasso <- scaled.lasso(zsmaela , zone_size_trimetho)
Slasso = 1-scaled_lasso$sigma2^2/var(zone_size_trimetho)
k = scaled_lasso$non.zeros
Slasso
Slasso- c(-1,1)*log(.5)*(k*log(p)/n + 1/sqrt(n))
end_time <- Sys.time()
end_time - start_time

# Enet method
start_time <- Sys.time()
glmnetcv = cv.glmnet(zsmaela, zone_size_trimetho, alpha = 0.01,parallel = T)
tam = predict(glmnetcv, type = 'nonzero', s = 'lambda.min')[[1]] 
enet = coef( glmnetcv,  s='lambda.min')[-1]
Enet = enet[tam] %*% cov2(zsmaela[,tam]) %*% enet[tam]/ var(zone_size_trimetho)
Enet
end_time <- Sys.time()
end_time - start_time


n = nrow(zsmaela)
splitdata <- sample(rep(seq(2), length = n)) == 1
glmnetcv = cv.glmnet(zsmaela[!splitdata,], zone_size_trimetho[!splitdata], alpha = 0.01,parallel = T)
tam = predict(glmnetcv, type = 'nonzero', s = 'lambda.min')[[1]] 


mmdker = mm2dicker(zone_size_trimetho[splitdata], zsmaela[splitdata,tam],CI = T )
mmdker


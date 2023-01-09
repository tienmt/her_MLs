
library(glmnet)
library(HiLMM)
require(doMC)
registerDoMC(cores=10)
noise = 1
source('/data2/thetm/review_HERI/mainfunctions.R')
load('/data2/thetm/co_Heritablity/John sharing data/Nei.G.data.rda')

p = ncol(snp.NeisGo.filtered)
n = nrow(snp.NeisGo.filtered)


ng.pheno.tet = ng.pheno.tet[!is.na(ng.pheno.tet)]
matSNPS = snp.NeisGo.filtered[names(ng.pheno.tet),] 
# Enet method
start_time <- Sys.time()
glmnetcv = cv.glmnet(matSNPS, ng.pheno.tet, alpha = 0.01,parallel = T)
tam = predict(glmnetcv, type = 'nonzero', s = 'lambda.min')[[1]] 
enet = coef( glmnetcv,  s='lambda.min')[-1]
Enet = enet[tam] %*% cov2(matSNPS[,tam]) %*% enet[tam]/ var(ng.pheno.tet)
Enet
end_time <- Sys.time()
end_time - start_time


# convex optimization method
start_time <- Sys.time()
eprsm <- EigenPrism(ng.pheno.tet, matSNPS, target = 'heritability')
eprsm
end_time <- Sys.time()
end_time - start_time

# MLE method
start_time <- Sys.time()
M = tcrossprod(matSNPS)
O = eigen(M)$vectors
lambda = eigen(M)$values
yTilde = t(O) %*% ng.pheno.tet
mle.r =ridge.mle(yTilde,lambda,1.e-6)
eta_chap = her.trunc(1 - mle.r$sig2hat / var(ng.pheno.tet) )
eta_chap
eta_chap + qnorm(1-0.05/2)* c(-1,1)/sqrt(2*n)
end_time <- Sys.time()
end_time - start_time

# moment method
start_time <- Sys.time()
mmdker = mm2dicker(ng.pheno.tet, matSNPS,CI = T )
mmdker
end_time <- Sys.time()
end_time - start_time

# scaled lasso method
start_time <- Sys.time()
scaled_lasso <- scaled.lasso(matSNPS , ng.pheno.tet)
Slasso = 1-scaled_lasso$sigma2^2/var(ng.pheno.tet)
k = scaled_lasso$non.zeros
Slasso
Slasso- c(-1,1)*log(.5)*(k*log(p)/n + 1/sqrt(n))
end_time <- Sys.time()
end_time - start_time





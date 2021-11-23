
library(glmnet)
library(HiLMM)
require(doMC)
registerDoMC(cores=10)
library(RSpectra)
noise = 1
source('/data2/thetm/review_HERI/mainfunctions.R')
load('/data2/thetm/Ecoli_data/newsnps.rda')
load('/data2/thetm/Ecoli_data/phenos.map.sampl.name.rda')
colnames(snp.ecoli) <- 1:121779

177858	180383	424550	425707 
674206	675417	890044	891267	
1965977	1968025 2637747	2638571 
2822543	2824819 3812429	3815005

tam = as.numeric(test)
aa = tam [tam >= 177858 & tam <= 180383]
aa = c(aa,tam [tam >= 3812429 & tam <= 3815005])


y.gentamicin = y.gentamicin[!is.na(y.gentamicin)]
matSNP = snp.ecoli[names(y.gentamicin),] 
p = ncol(matSNP)
n = nrow(matSNP)
# convex optimization method
start_time <- Sys.time()
eprsm <- EigenPrism(y.gentamicin, matSNP, target = 'heritability')
eprsm
end_time <- Sys.time()
end_time - start_time

# MLE method
start_time <- Sys.time()
M = tcrossprod(matSNP)
O = eigen(M)$vectors
lambda = eigen(M)$values
yTilde = t(O) %*% y.gentamicin
mle.r =ridge.mle(yTilde,lambda,1.e-6)
eta_chap = her.trunc(1 - mle.r$sig2hat / var(y.gentamicin) )
eta_chap
eta_chap + qnorm(1-0.05/2)* c(-1,1)/sqrt(2*n)
end_time <- Sys.time()
end_time - start_time

# moment method
start_time <- Sys.time()
mmdker = mm2dicker(y.gentamicin, matSNP,CI = T )
mmdker
end_time <- Sys.time()
end_time - start_time

# scaled lasso method
start_time <- Sys.time()
scaled_lasso <- scaled.lasso(matSNP , y.gentamicin)
Slasso = 1-scaled_lasso$sigma2^2/var(y.gentamicin)
k = scaled_lasso$non.zeros
Slasso
Slasso- c(-1,1)*log(.5)*(k*log(p)/n + 1/sqrt(n))
end_time <- Sys.time()
end_time - start_time

# Enet method
start_time <- Sys.time()
glmnetcv = cv.glmnet(matSNP, y.gentamicin, alpha = 0.01,parallel = T)
tam = predict(glmnetcv, type = 'nonzero', s = 'lambda.min')[[1]] 
enet = coef( glmnetcv,  s='lambda.min')[-1]
Enet = enet[tam] %*% cov2(matSNP[,tam]) %*% enet[tam]/ var(y.gentamicin)
Enet
end_time <- Sys.time()
end_time - start_time




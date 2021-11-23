

cov2 <- function ( x ) {
  crossprod ( scale (x , TRUE , FALSE ) )/( NROW ( x ) -1)
}
cov2 <- compiler::cmpfun(cov2)
cor.scale = function(snp,y){
  a = snp -mean(snp)
  b = y - mean(y)
  abs(a%*%b / sqrt(sum(a^2)*sum(b^2)) )
}
cor.scale <- compiler::cmpfun(cor.scale)


library(glmnet)
require(doMC)
registerDoMC(cores=10)

load('/data2/thetm/MA/new.masnps.rda')
masnps.filter = masnps.filter/2
#

zsmaela = mesnp[new.sam,]
p = ncol(zsmaela)
n = nrow(zsmaela)
 #SCREENING

  sam.cor <- apply(zsmaela, 2,function(snp) cor.scale(snp,zone_size_tetra) )
  xsele = zsmaela[,sam.cor> quantile(sam.cor, 0.6) ]
  p = ncol(xsele)
  xdat2 = xsele 
  for (itrrs in 1:(log2(p/n)+1) ) {
    model <- cv.glmnet(xdat2 , zone_size_tetra , alpha=0, parallel = T)
    coef0 <- abs(predict(model,type="coefficients",s='lambda.min')[-1] )
    xdat2 <- xdat2[, coef0 > quantile(coef0,.5)]
    if(ncol(xdat2) < n ) break
  }
  
  #boosting HERRA
  outlist <- foreach(jj = seq(100) ) %dopar%  {
    foreach(i = seq(2), .packages=c("glmnet")) %dopar%  {
      which <- sample(rep(seq(2), length = n)) == i
      xols1 = xdat2[, predict(cv.glmnet(xdat2[which,], zone_size_tetra[which], parallel = T,alpha = 1),
                              type = 'nonzero', s = 'lambda.min')[[1]] ]
      a = NA
      tryCatch({a <-(summary(lm(zone_size_tetra[!which] ~ xols1[!which,] ))$sigma)^2 
      },  error=function(e){ } )
    }}
 1 - mean( unlist(outlist), na.rm=T)/var(zone_size_tetra)
  
 sd(1 -  unlist(outlist), na.rm=T)/var(zone_size_tetra)

 
 
 
 
 z.s.list = list(ng.pheno.azi = ng.pheno.azi,
                 ng.pheno.cfx = ng.pheno.cfx,
                 ng.pheno.cip = ng.pheno.cip,
                 ng.pheno.pen = ng.pheno.pen,
                 ng.pheno.tet = ng.pheno.tet)
 
 

 herra.outlist <- function(x,y){
   y = y[!is.na(y)]
   matSNP = x[names(y),] 
   p = ncol(matSNP)
   n = nrow(matSNP)
   xdat2 = matSNP 
   for (itrrs in 1:(log2(p/n)+1) ) {
     model <- cv.glmnet(xdat2 , y , alpha=0, parallel = T)
     coef0 <- abs(predict(model,type="coefficients",s='lambda.min')[-1] )
     xdat2 <- xdat2[, coef0 > quantile(coef0,.5)]
     if(ncol(xdat2) < n ) break
   }
   foreach(jj = seq(100) ) %dopar%  {
     foreach(i = seq(2), .packages=c("glmnet")) %dopar%  {
       which =  sample(rep(seq(2), length = n)) == i
       xols1 = xdat2[, predict(cv.glmnet(xdat2[which,], y[which], parallel = T,alpha = 1),
                               type = 'nonzero', s = 'lambda.min')[[1]] ]
       a = NA
       tryCatch({a <-(summary(lm(y[!which] ~ xols1[!which,] ))$sigma)^2 
       },  error=function(e){ } )
     }}
 }
 herra.outlist <- compiler::cmpfun(herra.outlist)
 
 start_time <- Sys.time()
 herra.boost <- mclapply(z.s.list ,function(pheno) unlist(herra.outlist(snp.NeisGo.filtered,pheno) ) ,mc.cores = 5)
 end_time <- Sys.time()
 end_time - start_time
 
 1 - mean(unlist(herra.boost$trime.sulfa), na.rm=T)/var((trime.sulfa == 'R')*1)
 
 1 - mean(unlist(herra.boost$trime.sulfa), na.rm=T)/var((trime.sulfa == 'R')*1) + 
   c(-1,1)*sd(unlist(herra.boost$trime.sulfa), na.rm=T)/var((trime.sulfa == 'R')*1)
 
 
 1 - mean(unlist(herra.boost$ng.pheno.tet), na.rm=T)/var(ng.pheno.tet, na.rm=T)
 1 - mean(unlist(herra.boost$ng.pheno.tet), na.rm=T)/var(ng.pheno.tet, na.rm=T) + 
   c(-1,1)*sd(unlist(herra.boost$ng.pheno.tet), na.rm=T)/var(ng.pheno.tet, na.rm=T)
 
 
 
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

library(corrplot)
par(mfrow= c(1,2),mar = c(.1,.1,.1,.1))
#samples correlation
corrplot(cor(t(x[sample(n,100),])),method = "shade",tl.pos='n',order = "hclust",col = gray.colors(100))
#snp correlation
corrplot(cor(x[,sample(p,100)]),method = "shade",tl.pos='n',order = "hclust",col = gray.colors(100))



x = 1:7
avg = c(0.71846,0.7763, 0.8166, 0.9999, 0.7656, 0.6461, 0.8048)
lower = c(0.71846,0.5136, 0.7572, 0.6547, 0.5679, 0.5679, 0.7757)
upper = c(0.71846,1.0000, 0.8761, 1.0000, 0.9642, 0.7244, 0.8338)

causal = c(.7650, .2422, .8252, .9999, .7757, .8561, .8030)
cau.low = c(.7650, .0000, .7685, .9644, .6426, .8230, .7391)
cau.up = c(.7650, .5387, .8818, 1, .9088, .8892, .8610)

par(mar=c(2,4,2,.5))
plot(x, avg, lwd= 2, ylim= c(0, 1),
     pch= 0, col= "#56B4E9", xlab="",  ylab="heritability", main="", xaxt = "n" )
arrows(x, lower, x, upper, length=0.09, angle=90, code=3,
       lwd=3, col= "#56B4E9")
grid(lwd=2)
points(x+0.1,causal,col= "#E69F00",pch=2, lwd =2)
arrows(x+0.1, cau.low, x+0.1, cau.up, length=0.09, angle=90, code=3,
       lwd=3, col= "#E69F00")
axis(1, at=1:7, labels=c('Enet','Eprism','MLE','Moment','SLasso','GCTA','BoostHER'))
add_legend("top",inset=.03,
           c("whole genomes      ",'causal genes'),
           lty=c(1,1),lwd=c(2,2), col = c("#56B4E9", "#E69F00"),
           pch = c(0,2), horiz=TRUE, bty='n', cex=1)


setwd("~/Dropbox/ongoing_works/Heritability/review Heritablity/R codes")

wholegenes = cbind(oracle , elasnet, eprism , mle , mm , slasso , bher)
wholegenes = cbind(wholegenes,'wholegenes')
wholegenes = as.data.frame(wholegenes)

causalgenes = cbind(oracle , elasnet, eprism , mle , mm , slasso , bher)
causalgenes = cbind(causalgenes,"V8"='causalgenes')
causalgenes = as.data.frame(causalgenes)

subsample1500 = cbind(oracle , elasnet, eprism , mle , mm , slasso , bher)
subsample1500 = cbind(subsample1500,'subsample1500')
subsample1500 = as.data.frame(subsample1500)

subsample500 = cbind(oracle , elasnet, eprism , mle , mm , slasso , bher)
subsample500 = cbind(subsample500,'subsample500')
subsample500 = as.data.frame(subsample500)

mydat = gtools::smartbind(wholegenes,subsample1500,subsample500)

mydat[,1:7] = sapply(mydat[,1:7], as.numeric)
mydat$V8 <- factor(mydat$V8, levels = c('wholegenes','subsample1500','subsample500'),ordered = TRUE)
mydat = reshape2::melt(mydat, id = "V8") 
library(ggplot2)
bp<-ggplot(data =mydat , aes(x = variable, y = value,color = V8)) + 
 geom_boxplot() + 
  xlab('') + ylab('heritability')+
  labs(color = 'settings\n')
bp +scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))







add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


avg = c(0.7974,0.5846, 0.8638, 0.9999, 0.8206, 0.7123, 0.8462)
lower = c(0.7974,0.4438, 0.8387, 0.9109, 0.4784, 0.6727, 0.8330)
upper = c(0.7974,0.7254, 0.8889, 1.0000, 1.0000, 0.7518, 0.8594)

causal = c(0.8012, NA, 0.8325, 0.8595, 0.7912, 0.8394, 0.8215)
cau.low = c(0.8012, NA, 0.5881, 0.8530, 0.7264, 0.8120, 0.7983)
cau.up = c(0.8012, NA, 1.0000, 0.8661, 0.8560, 0.8667, 0.8423)
x = 1:7
par(mar=c(2,4,2,.5))
plot(x, avg, lwd= 2, ylim= c(0.44, 1),
     pch= 0, col= "#56B4E9", xlab="",  ylab="heritability", main="", xaxt = "n" )
arrows(x, lower, x, upper, length=0.09, angle=90, code=3,
       lwd=3, col= "#56B4E9")
grid(lwd=2)
points(x+0.1,causal,col= "#E69F00",pch=2, lwd =2)
arrows(x+0.1, cau.low, x+0.1, cau.up, length=0.09, angle=90, code=3,
       lwd=3, col= "#E69F00")
axis(1, at=1:7, labels=c('Enet','Eprism','MLE','Moment','SLasso','GCTA','BoostHER'))
add_legend("top",inset=.03,
           c("whole genomes      ",'causal genes'),
           lty=c(1,1),lwd=c(2,2), col = c("#56B4E9", "#E69F00"),
           pch = c(0,2), horiz=TRUE, bty='n', cex=1)










### Tetracycline
avg = c(0.7342,0.9295, 0.8302, 0.9999, 0.7593, 0.7514, 0.8435)
lower = c(0.7342, 0.7887, 0.8051, 0.9109, 0.3851, 0.7165, 0.8300)
upper = c(0.7342, 1.0000 , 0.8554, 1.0000, 1.0000, 0.7863, 0.8570)

causal = c(0.7362, NA, NA, 0.3612, NA, 0.7869, 0.7458)
cau.low = c(0.7362, NA, NA, 0.3587, NA, 0.7280, 0.7201)
cau.up = c(0.7362, NA, NA, 0.3638, NA, 0.8459, 0.7746)
x = 1:7
par(mar=c(2,4,2,.5))
plot(x, avg, lwd= 2, ylim= c(0, 1),
     pch= 0, col= "#56B4E9", xlab="",  ylab="heritability", main="", xaxt = "n" )
arrows(x, lower, x, upper, length=0.09, angle=90, code=3,
       lwd=3, col= "#56B4E9")
grid(lwd=2)
points(x+0.1,causal,col= "#E69F00",pch=2, lwd =2)
arrows(x+0.1, cau.low, x+0.1, cau.up, length=0.09, angle=90, code=3,
       lwd=3, col= "#E69F00")
axis(1, at=1:7, labels=c('Enet','Eprism','MLE','Moment','SLasso','GCTA','BoostHER'))


### Co-trimoxazole
avg = c(0.7362,0.9999, 0.8188, 0.9999, 0.7406, 0.7826, 0.7571)
lower = c(0.7362, 0.8582, 0.7936, 0.9109, 0.4887, 0.7518, 0.7361)
upper = c(0.7362, 1.0000 , 0.8439, 1.0000, 0.9924, 0.8134, 0.7781)

causal = c(0.6499, NA, 0.5873, 0.4084, 0.5082, 0.7228, 0.6652)
cau.low = c(0.6499, NA, 0.0480, 0.4062, 0.4842, 0.6428, 0.6414)
cau.up = c(0.6499, NA, 1, 0.4105, 0.5321, 0.8029, 0.6927)
x = 1:7
par(mar=c(2,4,2,.5))
plot(x, avg, lwd= 2, ylim= c(0, 1),
     pch= 0, col= "#56B4E9", xlab="",  ylab="heritability", main="", xaxt = "n" )
arrows(x, lower, x, upper, length=0.09, angle=90, code=3,
       lwd=3, col= "#56B4E9")
grid(lwd=2)
points(x+0.1,causal,col= "#E69F00",pch=2, lwd =2)
arrows(x+0.1, cau.low, x+0.1, cau.up, length=0.09, angle=90, code=3,
       lwd=3, col= "#E69F00")
axis(1, at=1:7, labels=c('Enet','Eprism','MLE','Moment','SLasso','GCTA','BoostHER'))



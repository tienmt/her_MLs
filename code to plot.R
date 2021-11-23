
Eprism[Eprism>1] = 1
epr.boost[epr.boost>1] = 1
Eprism[Eprism < 0] = 0
epr.boost[epr.boost < 0] = 0

boxplot(cbind(h2aprx,Enet,
              Eprism,B_Eprism=epr.boost,
              MLE=Re(MLE),B_MLE=mle.boost,
              Moment,B_mm = mm.boost,
              HERRA, B_herra = herra.boost),
main=bquote('screening 90%, 20 random SNPs,'~ sigma[epsilon]^2==1~","~ h^2==0.) )
abline(h=mean(h2aprx))
~","~ h^2==0.8
1015x600



boxplot(cbind(Chloramphenical= 1-unlist(herra.boost$y.chl)/var(phenotypes$y.chl),
              Erythromycin= 1-unlist(herra.boost$y.eryth)/var(phenotypes$y.eryth),
              Tetracycline= 1-unlist(herra.boost$y.tetracy)/var(phenotypes$y.tetracy),
              penicillin= 1-unlist(herra.boost$ylactam)/var(phenotypes$ylactam),
              'Co-trimoxazole'= 1-unlist(herra.boost$ytrime)/var(phenotypes$ytrime))
        ,staplelwd = 2 , boxfill = NA, border = NA)
boxplot(cbind(Chloramphenical= 1-unlist(herra.boost$y.chl)/var(phenotypes$y.chl),
              Erythromycin= 1-unlist(herra.boost$y.eryth)/var(phenotypes$y.eryth),
              Tetracycline= 1-unlist(herra.boost$y.tetracy)/var(phenotypes$y.tetracy),
              penicillin= 1-unlist(herra.boost$ylactam)/var(phenotypes$ylactam),
              'Co-trimoxazole'= 1-unlist(herra.boost$ytrime)/var(phenotypes$ytrime))
        , xaxt = "n", add = TRUE,  boxwex=0.25, 
        boxfill="black", at = 1:5 - 0.15)
boxplot(cbind(Chloramphenical= zs.chl.b.herra,
              Erythromycin= zs.eryth.b.herra,
              Tetracycline= zs.tetra.b.herra,
              penicillin= zs.oxacilliin.b.herra,
              'Co-trimoxazole'= zs.trimtho.b.herra)
        , xaxt = "n", add = TRUE,  boxwex=0.25,
        boxfill="grey",  at = 1:5+ 0.15)
legend("topright",inset=.01, cex = 1,
       c("binary phenotypes","continuous zone sizes"),bg="grey96",
       horiz=F, lty=c(5,5),lwd=c(30,30),
       col = c("white","grey"))



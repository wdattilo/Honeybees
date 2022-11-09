SCRIPT
###GAMS
library(ggplot2)
library(mgcv)
DATOS <- read.csv("METRICAS_APIS_2_1-112.csv", header = TRUE)
DATOS <- DATOS[-(108:110), ]
names(DATOS)
dim(DATOS)
attach(DATOS)


######## Diferencias PC1 exotica / nativa
library(lmerTest)
library(car)
lmm <- glmer(PC1_postive~status+(1|id), family=Gamma)
summary(lmm)
Anova(lmm,type="III")
boxplot(PC1_postive~status, outline=F,col=c('red', 'green'))




library(gridExtra)

a = ggplot(DATOS, aes(x=degree)) + geom_density(color="darkblue", fill="lightblue") +  xlab("Species degree")

b = ggplot(DATOS, aes(x=betweenness)) + geom_density(color="darkblue", fill="lightblue") +  xlab("Betweenness centrality")

c = ggplot(DATOS, aes(x=strength)) + geom_density(color="darkblue", fill="lightblue") +  xlab("Species strength")

d = ggplot(DATOS, aes(x=katzz)) + geom_density(color="darkblue", fill="lightblue") +  xlab("Katz centrality")


grid.arrange(a,b,c,d, ncol = 2)


par(mfrow=c(3,2)) 



for(i in seq(13, 41, 4)[-c(2,3)]){
  EFFCT <-  ((DATOS[,i+1] - DATOS[,i+2])/DATOS[,i+3]) #)EFCTOS DE APIS -  EFCTO PROMDIO ) / STDV
  M <- gamm(EFFCT ~  s(PC1_postive, k=3, bs="cr"), random=list(id=~1),
, data=DATOS)
  #plot(fitted(M$gam), residuals(M$gam))
  if( summary(M$gam)$s.table[1,4] <= 0.05){
    print(colnames(DATOS)[i])
    print(paste("R2 = ", round(summary(M$gam)$r.sq, 2)))
    print(summary(M$gam))
    
    SEQ <- seq(0, 13, length=30)
    PRED <- predict(M$gam, data.frame(PC1_postive=SEQ))
    plot(DATOS$PC1_postive, EFFCT, pch = 21, col="gray40", bg="gray90", ylab=c(paste("Standarized effect in", 
    colnames(DATOS)[i]), "after excluding Apis mellifera"), 
    xlab="Positive PC1 scores", )
    lines(SEQ, PRED, col="blue", lwd=1.5)
  }
}










library(mgcv)
DATOS <- read.csv("METRICAS_APIS_2_1-112.csv", header = TRUE)
DATOS <- DATOS[-(108:110), ]
names(DATOS)
dim(DATOS)

par(mfrow=c(3,1))

  i<-13 #"H2"
  colnames(DATOS)[i]
  EFFCT <-  ((DATOS[,i+1] - DATOS[,i+2])/DATOS[,i+3])
  M <- gamm(EFFCT ~  s(PC1_postive, k=3, bs="cr"), random=list(id=~1), data=DATOS)
  summary(M$gam)
  SEQ <- seq(0, 13, length=30)
  PRED <- predict(M$gam, data.frame(PC1_postive=SEQ))
  plot(DATOS$PC1_postive, EFFCT, pch = 21, col="gray40", bg="gray90", cex=1, ylab=c(paste("Standarized effect in", colnames(DATOS)[i]), "after excluding Apis mellifera"), xlab="Interactive role of Apis mellifera (PC1)", )
  lines(SEQ, PRED, col="blue", lwd=1.5)
    
  
  
  i<-37 #"niche.plants"
  colnames(DATOS)[i]
  EFFCT <-  ((DATOS[,i+1] - DATOS[,i+2])/DATOS[,i+3])
  M <- gamm(EFFCT ~  s(PC1_postive, k=3, bs="cr"), random=list(id=~1), correlation=corAR1(), data=DATOS)
  #plot(fitted(M$gam), residuals(M$gam))
  summary(M$gam)
  SEQ <- seq(0, 13, length=30)
  PRED <- predict(M$gam, data.frame(PC1_postive=SEQ))
  plot(DATOS$PC1_postive, EFFCT,pch = 21, col="gray40", bg="gray90", cex=1,  ylab=c(paste("Standarized effect in", colnames(DATOS)[i]), "after excluding Apis mellifera"), xlab="Interactive role of Apis mellifera (PC1", )
  lines(SEQ, PRED, col="blue", lwd=1.5)
  
  

  i<-41 #"robustness"
  colnames(DATOS)[i]
  EFFCT <-  ((DATOS[,i+1] - DATOS[,i+2])/DATOS[,i+3])
  M <- gamm(EFFCT ~  s(PC1_postive, k=3, bs="cr"), random=list(id=~1), correlation=corAR1(), data=DATOS)
  #plot(fitted(M$gam), residuals(M$gam))
  summary(M$gam)
  SEQ <- seq(0, 13, length=30)
  PRED <- predict(M$gam, data.frame(PC1_postive=SEQ))
  plot(DATOS$PC1_postive, EFFCT, pch = 21, col="gray40", bg="gray90", cex=1, ylab=c(paste("Standarized effect in", colnames(DATOS)[i]), "after excluding Apis mellifera"), xlab="Interactive role of Apis mellifera (PC1", )
  lines(SEQ, PRED, col="blue", lwd=1.5)
  
  
  
  
######## Autocorrelacion espacial
library(letsR)
library(ape)

dists <- as.matrix(dist(cbind(DATOS$long, DATOS$lat)))
dists.inv <- 1/dists
diag(dists.inv) <- 0
dists.inv
R5 = residuals(M$gam, type="pearson")

Moran.I(R5, dists.inv)


### Correlation plot
library(corrplot)
library(RColorBrewer)
attach(DATOS)
DATOS <- DATOS[,c('degree','betweenness','strength','katzz')]


M <-cor(DATOS)
corrplot(M, type="upper", order="hclust",col=brewer.pal(n=8, name="RdYlBu"))




#FIGURAS FINAIS APENAS NO WEIGHTED..

par(mfrow=c(3,1))

#"niche.plants"
  i <- which(colnames(DATOS)=="niche.plants")
  colnames(DATOS)[i]
  #EFFECT <-  ((DATOS[,i+1] - DATOS[,i+2])/DATOS[,i+3])
  EFFECT <-  ((DATOS[,i+1] - DATOS[,i+2]) / (DATOS[,i+1] + DATOS[,i+2]))
  
  M <- gamm(EFFECT ~  s(PC1_postive, k=3, bs="tp"), random=list(id=~1), correlation=corAR1(), data=DATOS)
  
  #plot(fitted(M$gam), residuals(M$gam))
  #qqnorm(residuals(M$gam))
  #qqline(residuals(M$gam))
  #gam.check(M$gam, plat=FALSE)
  summary(M$gam)
  
  plot(PC1_postive, EFFECT, col="gray40", bg="gray90", cex=1.2,  ylab="", xlab="", pch=21, axes=F)
  axis(2, labels= round(seq(-1.2, 0.2, 0.2),2), at=seq(-1.2, 0.2, 0.2))
  axis(1, labels=rep("",9), at=seq(-2,14,2))
  axis(3, at=c(-100,100), label=c(-100,100))
  axis(4, at=c(-100,100), label=c(-100,100))
  mtext("Weighted impact on", side=2, line=3, at=-0.5, cex=0.8)
  mtext(expression(paste("niche overlap (", italic("NO"),")")), side=2, line=1.8, at=-0.5, cex=0.8)
  mtext("(A)", side=3, at=12.5, line=-1.2, cex=0.8)
  
  SEQ <- seq(min(PC1_postive), max(PC1_postive), length=30)
  PRED <- predict(M$gam, data.frame(PC1_postive=SEQ), se=TRUE)
  lines(SEQ, PRED$fit, col="blue", lwd=1.5)
  lines(SEQ, PRED$fit+PRED$se.fit*2, col="cornflowerblue", lwd=1.5, lty=2)
  lines(SEQ, PRED$fit-PRED$se.fit*2, col="cornflowerblue", lwd=1.5, lty=2)
  
  
  #"H2"
  i <- which(colnames(DATOS)=="H2")
  colnames(DATOS)[i]
  #EFFECT <-  ((DATOS[,i+1] - DATOS[,i+2])/DATOS[,i+3])
  EFFECT <-  ((DATOS[,i+1] - DATOS[,i+2]) / (DATOS[,i+1] + DATOS[,i+2]))
  M <- gamm(EFFECT ~  s(PC1_postive, k=3, bs="tp"), random=list(id=~1), correlation=corAR1(), data=DATOS)
  #plot(fitted(M$gam), residuals(M$gam))
  #qqnorm(residuals(M$gam))
  #qqline(residuals(M$gam))
  #gam.check(M$gam)
  
  summary(M$gam)
  plot(PC1_postive, EFFECT, col="gray40", bg="gray90", cex=1,  ylab= "", xlab="", pch=21, axes=F)
  axis(2, labels= round(seq(-0.2, 1.2, 0.1),2), at= round(seq(-0.2, 1.2, 0.1),2))
  axis(1, labels=rep("",9), at=seq(-2,14,2))
  axis(3, at=c(-100,100), label=c(-100,100))
  axis(4, at=c(-100,100), label=c(-100,100))
  mtext("Weighted impact on", side=2, line=3, at=0.18, cex=0.8)
  mtext(expression(paste("network specialization (", italic("H"[2]),")")), side=2, line=1.8, at=0.18, cex=0.8)
  mtext("(B)", side=3, at=12.5, line=-1.2, cex=0.8)
  legend(9, -5, c("Native", "Exotic"), pch=c(19, 15), col=paste(c("#7F7F7F", "#8B0000"),95, sep=""), bg=c("#7F7F7F", "#8B0000"), bty="n")
  
  SEQ <- seq(min(PC1_postive), max(PC1_postive), length=30)
  PRED <- predict(M$gam, data.frame(PC1_postive=SEQ), se=TRUE)
  lines(SEQ, PRED$fit, col="blue", lwd=1.5)
  lines(SEQ, PRED$fit+PRED$se.fit*2, col="cornflowerblue", lwd=1.5, lty=2)
  lines(SEQ, PRED$fit-PRED$se.fit*2, col="cornflowerblue", lwd=1.5, lty=2)
  
  
  #"robustness"
  i<- which(colnames(DATOS)=="robustness")
  colnames(DATOS)[i]
  #shapiro.test(DATOS[,i]^(0.99))
  #hist(DATOS[,i]^0.99)
  
  #EFFECT <-  ((DATOS[,i+1] - DATOS[,i+2])/DATOS[,i+3])
  EFFECT <-  ((DATOS[,i+1] - DATOS[,i+2]) / (DATOS[,i+1] + DATOS[,i+2]))
  M <- gamm(EFFECT ~  s(PC1_postive, k=3, bs="tp"), random=list(id=~1), correlation=corAR1(), data=DATOS)
  #plot(fitted(M$gam), residuals(M$gam))
  summary(M$gam)
  
  plot(PC1_postive, EFFECT, col="gray40", bg="gray90",  cex=1, ylab="", xlab="", pch=21, axes=T)
  mtext("Interactive role of", side=1, at=6, line=2, cex=0.8)
  mtext(expression(paste(italic("A. mellifera"), "(NMDS 1)")), side=1, at=6, line=3.5, cex=0.8)
  #axis(2, labels=seq(-5, 20, 3), at=seq(-5, 20, 3))
  #axis(1, labels=seq(-2, 14, 2), at=seq(-2, 14, 2))
  #axis(3, at=c(-100,100), label=c(-100,100))
  #axis(4, at=c(-100,100), label=c(-100,100))
  mtext("Weighted impact on", side=2, line=3, at=0.13, cex=0.8)
  mtext(expression(paste("network robustness (", italic("R"),")")), side=2, line=1.8, at=0.13, cex=0.8)
  mtext("(C)", side=3, at=12.5, line=-1.2, cex=0.8)
  
  SEQ <- seq(min(PC1_postive), max(PC1_postive), length=30)
  PRED <- predict(M$gam, data.frame(PC1_postive=SEQ), se=TRUE)
  lines(SEQ, PRED$fit, col="blue", lwd=1.5)
  lines(SEQ, PRED$fit+PRED$se.fit*2, col="cornflowerblue", lwd=1.5, lty=2)
  lines(SEQ, PRED$fit-PRED$se.fit*2, col="cornflowerblue", lwd=1.5, lty=2)


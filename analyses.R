library(mgcv)
library(geodist)
library(car)
library(fields)
library(lmerTest)

DATOS <- read.csv("variables.csv")
DATOS$site <- as.factor(DATOS$site)
DATOS <- DATOS[-c(82, 93), ]
DATOS$site <- as.factor(DATOS$site)
DATOS$influence <- DATOS$influence*100
DATOS$PC1positivo <- DATOS$PC1+abs(min(DATOS$PC1))+0.0001
#DATOS$latjitter <- jitter(DATOS$lat, amount = .005)
#DATOS$longjitter <- jitter(DATOS$long, amount = .005)
#LATjitter <- DATOS$latjitter
#LONGjitter <- DATOS$longjitter
#write.csv(cbind(LONGjitter, LATjitter), "~/Downloads/COORDSjittered.csv", row.names = FALSE)
attach(DATOS)

##DISTANCIA DE DESPLAZAMIENTO ENTRE COORDENADAS ORIGINALES Y CORDENADAS CON RUIDO
DESPLZAMIENTO <- numeric()
for(i in 1:nrow(DATOS))DESPLZAMIENTO <- c(DESPLZAMIENTO, geodist(cbind(lon=c(long[i],LONGjitter[i]), lat=c(lat[i],LATjitter[i]), measure="haversine"))[1,2])
mean(DESPLZAMIENTO)
sd(DESPLZAMIENTO)
DIST <- geodist(cbind(lon=LONGjitter, lat=LATjitter), measure="haversine")
library(MASS)
NDMS <- isoMDS(dist(DIST))
DATOS$PCA1_DIST <- NDMS$points[,1]
DATOS$absLAT <- abs(DATOS$LATjitter)
attach(DATOS)


###PC1
###############################################
##ESTE SALE BIEN, Y ES EL QUE HAY QUE REPORTAR, Tiene 26% de pseudo R2 y buen ajuste
Mnull <- gamm(PC1positivo ~ 1, random=list(site = ~ PCA1_DIST|site), correlation=corGaus(1,form=~PCA1_DIST), family=Gamma("log"), data=DATOS)
M <- gamm(PC1positivo ~ s(precipitation, k=4, bs="tp")  + s(absLAT, k=4, bs="tp") + te(precipitation, absLAT, k=8) +  s(influence, k=4, bs="tp")  + s(influence, absLAT, k=5) + s(temperature, k=4, bs="tp"), random=list(site = ~ PCA1_DIST), correlation=corGaus(form=~PCA1_DIST|site), family=Gamma("log"), data=DATOS)
gam.check(M$gam)#INDICA QUE LOS AJUSTES SON ADECUADO (NO ES LA SIGNIFICACINA DE LAS VARIABLES EXPLICATIVAS)
summary(M$gam)
plot(M$gam$fitted.values, residuals(M$gam, type="pearson")) #Medio chafa
abline(h=0)
qqnorm(M$gam$residuals)
qqline(M$gam$residuals)
CONCURVITY<- concurvity(M$gam,full=FALSE)
VIF <- as.dist(1/(1-CONCURVITY$estimate^2))
VIF#Eliminamos s(precipitation, absLAT). SON COLINEALES

M <- gamm(PC1positivo ~ s(precipitation, k=4, bs="tp")  + s(absLAT, k=4, bs="tp")  +  s(influence, k=4, bs="tp")  + s(influence, absLAT, k=5) + s(temperature, k=4, bs="tp"), random=list(site = ~ PCA1_DIST), correlation=corGaus(form=~PCA1_DIST|site), family=Gamma("log"), data=DATOS)
gam.check(M$gam)#INDICA QUE LOS AJUSTES SON ADECUADO (NO ES LA SIGNIFICACINA DE LAS VARIABLES EXPLICATIVAS)

plot(M$gam$fitted.values, residuals(M$gam, type="pearson"))#Medio chafa
abline(h=0)
qqnorm(M$gam$residuals)
qqline(M$gam$residuals)

CONCURVITY<- concurvity(M$gam,full=FALSE)
VIF <- as.dist(1/(1-CONCURVITY$estimate^2))
VIF#TODO CHIDO!

#pseudoR2
DN <- sum(residuals(Mnull$gam, type = "deviance")^2) # Deviance from NULL GAMM
DR <- sum(residuals(M$gam, type = "deviance")^2) # Deviance from our GAMM
DEV_EXP <- (DN-DR)/DN # Deviance Explained from our GAMM model
DEV_EXP

summary(M$gam)
round(gam.vcomp(M$gam,rescale=F)^2,3)

#COMPONENTES DE DEVIANZA EXPPLICAOS POR CADA VARIABLE: TOTAL ES 0.2598981
#+ s(influence, absLAT, k=5)
Mparcial <- gamm(PC1positivo ~ s(precipitation, k=4, bs="tp")  + s(absLAT, k=4, bs="tp")  +  s(influence, k=4, bs="tp")  + s(temperature, k=4, bs="tp"), random=list(site = ~ PCA1_DIST), correlation=corGaus(form=~PCA1_DIST|site), family=Gamma("log"), data=DATOS)
DN <- sum(residuals(Mnull$gam, type = "deviance")^2) # Deviance from NULL GAMM
DRparcial <- sum(residuals(Mparcial$gam, type = "deviance")^2) # Deviance from our GAMM
WEIGHT_LAT_INFLUE <- (DN-DRparcial)/DN
WEIGHT_LAT_INFLUE


#+ s(temperature, k=4, bs="tp") + s(influence, absLAT, k=5)
Mparcial <- gamm(PC1positivo ~ s(precipitation, k=4, bs="tp")  + s(absLAT, k=4, bs="tp")  +  s(influence, k=4, bs="tp")  , random=list(site = ~ PCA1_DIST), correlation=corGaus(form=~PCA1_DIST|site), family=Gamma("log"), data=DATOS)
DN <- sum(residuals(Mnull$gam, type = "deviance")^2) # Deviance from NULL GAMM
DRparcial <- sum(residuals(Mparcial$gam, type = "deviance")^2) # Deviance from our GAMM
WEIGHT_TEMP <- (DN-DRparcial)/DN
WEIGHT_TEMP


# s(precipitation, k=4, bs="tp")+ s(influence, absLAT, k=5)
Mparcial <- gamm(PC1positivo ~  s(absLAT, k=4, bs="tp")  +  s(influence, k=4, bs="tp")  + s(temperature, k=4, bs="tp") , random=list(site = ~ PCA1_DIST), correlation=corGaus(form=~PCA1_DIST|site), family=Gamma("log"), data=DATOS)
DN <- sum(residuals(Mnull$gam, type = "deviance")^2) # Deviance from NULL GAMM
DRparcial <- sum(residuals(Mparcial$gam, type = "deviance")^2) # Deviance from our GAMM
WEIGHT_PRECIP <- (DN-DRparcial)/DN
WEIGHT_PRECIP


# + s(absLAT, k=4, bs="tp")   + s(influence, absLAT, k=5)
Mparcial <- gamm(PC1positivo ~  s(precipitation, k=4, bs="tp") +  s(influence, k=4, bs="tp")+ s(temperature, k=4, bs="tp") + s(influence, absLAT, k=5), random=list(site = ~ PCA1_DIST), correlation=corGaus(form=~PCA1_DIST|site), family=Gamma("log"), data=DATOS)
DN <- sum(residuals(Mnull$gam, type = "deviance")^2) # Deviance from NULL GAMM
DRparcial <- sum(residuals(Mparcial$gam, type = "deviance")^2) # Deviance from our GAMM
WEIGHT_LAT <-  (DN-DRparcial)/DN
WEIGHT_LAT

# + s(influence, k=4, bs="tp") + s(influence, absLAT, k=5)
Mparcial <- gamm(PC1positivo ~  s(precipitation, k=4, bs="tp") +  s(absLAT, k=4, bs="tp") + s(temperature, k=4, bs="tp") + s(influence, absLAT, k=5), random=list(site = ~ PCA1_DIST), correlation=corGaus(form=~PCA1_DIST|site), family=Gamma("log"), data=DATOS)
DN <- sum(residuals(Mnull$gam, type = "deviance")^2) # Deviance from NULL GAMM
DRparcial <- sum(residuals(Mparcial$gam, type = "deviance")^2) # Deviance from our GAMM
WEIGHT_INFL <-  (DN-DRparcial)/DN
WEIGHT_INFL

WEIGHTS <- 1-( c(PRECIP=WEIGHT_PRECIP, LAT=WEIGHT_LAT, INFL=WEIGHT_INFL,LAT_INLF=WEIGHT_LAT_INFLUE, TEM=WEIGHT_TEMP)/((DN-DR)/DN))

cbind(summary(M$gam)$s.table, WEIGHTS)

#PRECIPITACIÓN
SEQprec <- seq(min(precipitation), max(precipitation), length=30)
PRED <- predict(M$gam, data.frame(precipitation=SEQprec, temperature=mean(temperature), absLAT=mean(absLAT), influence=mean(influence)))
plot((PC1positivo) ~ precipitation, las=1, ylab=expression(paste("PC1 (importance of " ,italic("Apis mellifera"), ")")), xlab="Precipitation", pch=19, col="#AAAAAA99", cex=0.75, ylim=c(0,12))
lines(SEQprec, exp(PRED))


#LATITUDE
SEQlat <- seq(min(absLAT), max(absLAT), length=30)
PRED <- predict(M$gam, data.frame(precipitation=mean(precipitation), temperature=mean(temperature), absLAT=SEQlat, influence=mean(influence)))
plot((PC1positivo) ~ absLAT, las=1, ylab=expression(paste("PC1 (importance of " ,italic("Apis mellifera"), ")")), xlab="Latitude", pch=19, col="#AAAAAA99", cex=0.75, ylim=c(0,12))
lines(SEQlat, exp(PRED))

#TEMPERATURA
SEQtemp <- seq(min(temperature), max(temperature), length=30)
PRED<- predict(M$gam, data.frame(precipitation=mean(precipitation), temperature=SEQtemp, absLAT=mean(absLAT), influence=mean(influence)))
plot((PC1positivo)~temperature, las=1, ylab=expression(paste("PC1 (importance of " ,italic("Apis mellifera"), ")")), xlab="Temperature", pch=19, col="#AAAAAA99", cex=0.75, ylim=c(0,12))
lines(SEQtemp, exp(PRED))


###INFLUENCIA:LATITUD
EXPDATA <- expand.grid(Influence=seq(min(influence), max(influence), length.out = 20), ABS_latitude=seq(min(absLAT), max(absLAT), length.out = 20))
PRED_PC1<- exp(predict(M$gam, data.frame(precipitation=mean(precipitation), temperature=mean(temperature), absLAT=EXPDATA[,2], influence=EXPDATA[,1])))

COL <- paste(colorRampPalette(c( "BLUE", "#388EE9", "#28ACEA", "#19CEEB","#09EBEE","#00FEEF"))(20), "55", sep="")
plot(c(0,100), c(0, 4), type="n", las=1, ylab=expression(paste("PC1 (importance of " ,italic("Apis mellifera"), ")")), xlab="Global human influence index (%)")
k<-0
kk<-0
for(i in 1:20){
  kk<- kk+1
polygon(c(EXPDATA[(381:400)-k, 1], EXPDATA[(400:381)-k, 1]), c(PRED_PC1[(400:381)-k],rep(0, 20)), las=1, col=COL[kk], border=COL[kk])
  k<-k+20
}
#points(influence, (PC1positivo), pch=19, col="black", cex=0.5, line.col="white")
COL <- paste(colorRampPalette(c("BLUE", "#388EE9", "#28ACEA", "#19CEEB","#09EBEE","#00FEEF"))(20), "99", sep="")
colorbar.plot(100, 3.7, strip=seq(0,53, length=20), adj.x = 1, adj.y=1, strip.width=0.025, strip.length=0.4, col=COL[20:1], horizontal=TRUE)
text(c(65, 102), c(3.9, 3.9), c(0, 55), cex=1)
mtext("Latitude", side=3, line=-3, at=84, cex=1)

######## Autocorrelacion espacial
library(letsR)
library(ape)

dists <- as.matrix(dist(cbind(DATOS$LONGjitter, DATOS$LATjitter)))
dists.inv <- 1/dists
diag(dists.inv) <- 0
dists.inv
R5 = residuals(M$gam, type="pearson")

Moran.I(R5, dists.inv)



######## Diferencias PC1 exotica / nativa
library(lmerTest)
lmm <- glmer(PC1positivo~exotica+(1|site), family=Gamma)
summary(lmm)
Anova(lmm,type="III")
plot(lmm)
boxplot(PC1positivo~exotica, outline=F,col=c('darkgoldenrod1', 'brown2'))



######## MAPA

library(maps)
library(scales)
space <- read.csv("variables.csv", header=TRUE) 
attach(space)
map("world", fill = TRUE, col = "antiquewhite")
map.axes()
map.scale(-40, -50, ratio = F, cex = 0.4) 



rbPal <- colorRampPalette(c("darkgoldenrod1", "brown2"))

space$Col <- rbPal(10)[as.numeric(cut(exotica2,breaks = 10))]


points(x=space$long, y=space$lat, pch = c(21, 24)[as.numeric(exotica)], cex=1.5, col="white", bg=space$Col) 

compassRose<-function(x,y,rot=0,cex=1) {
  oldcex<-par(cex=cex)
  mheight<-strheight("M")
  xylim<-par("usr")
  plotdim<-par("pin")
  xmult<-(xylim[2]-xylim[1])/(xylim[4]-xylim[3])*plotdim[2]/plotdim[1]
  point.angles<-seq(0,2*pi,by=pi/4)+pi*rot/180
  crspans<-rep(c(mheight*3,mheight/2),length.out=9)
  xpoints<-cos(point.angles)*crspans*xmult+x
  ypoints<-sin(point.angles)*crspans+y
  for(point in 1:8) {
    pcol<-ifelse(point%%2,"black","white")
    polygon(c(xpoints[c(point,point+1)],x),c(ypoints[c(point,point+1)],y),col=pcol)
  }
  txtxpoints<-cos(point.angles[c(1,3,5,7)])*1.2*crspans[1]*xmult+x
  txtypoints<-sin(point.angles[c(1,3,5,7)])*1.2*crspans[1]+y
  text(txtxpoints,txtypoints,c("E","N","W","S"))
  par(oldcex)
}
compassRose(-10,-30, rot=0, cex=0.5)




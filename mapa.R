library(maps)
library(scales)
space <- read.csv("METRICAS_APIS_2_1-112.csv", header=TRUE) 
space <- space[-(108:110), ]
attach(space)
map("world", fill = TRUE, interior = FALSE, col = "snow2", bg="white")
map.axes()
map.scale(-40, -50, ratio = F, cex = 0.4) 



rbPal <- colorRampPalette(c("darkgoldenrod1", "brown2"))

space$Col <- rbPal(10)[as.numeric(cut(size,breaks = 10))]


points(x=space$long, y=space$lat, pch = c(21), cex=log(space$size)/2+1, col="white", bg="brown1") 

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

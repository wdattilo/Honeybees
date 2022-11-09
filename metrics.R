SCRIPT PARA CALCULAR TODAS LAS METRICAS

library(bipartite)
library(stringr)
library(readr)
source("~/Desktop/2021/INECOL/funciones.r")
##

FILES <- list.files("~/Downloads/Networks/", ".txt")
NETWORKS <- list()
id.redes <- character()
for(i in 1:length(FILES)){
  ENCODING <- as.character(guess_encoding(paste("~/Downloads/Networks/", FILES[i], sep=""), n_max = 1000)[1,1])
  NETWORKS[[i]] <- read.delim(paste("~/Downloads/Networks/", FILES[i], sep=""), row.names = 1, fileEncoding = ENCODING)
  id.redes <- c(id.redes, paste("R", unlist(strsplit(unlist(lapply(strsplit(FILES[i], "red"), "[",2)), ".txt")),sep="_"))
}

names(NETWORKS) <- id.redes

##METRICAS DE CENTRALIDAD
resultados<-data.frame(rede=numeric(),degree=numeric(), katzz=numeric(), katzscore=numeric(), betweenness=numeric(), closeness=numeric(), strength=numeric())

for(i in 1:112){
  net <- as.matrix(NETWORKS[[i]])
  n_p = nrow(net)  #informe o número de linhas 
  n_a = ncol(net)  #informe o numero de colunas
  
  # coloca ela em formato em quadrado
  A <- rbind(cbind(matrix(0, n_p, n_p), net), cbind(t(net), matrix(0, n_a, n_a)))
  
  net2 = katz(A)
  
  a = nrow(net) + 1
  b = length(net2)
  c = net2[a:b]
  d = mean(c)
  e = sd(c)
  katzz = net2[["Apis_mellifera"]]
  katzscore = (katzz-d)/e
  bet = specieslevel(net, index="betweenness")$`higher level`
  betweenness = bet["Apis_mellifera",2]
  clos = specieslevel(net, index="closeness")$`higher level`
  closeness = clos["Apis_mellifera",2]
  stre = specieslevel(net, index="species strength")$`higher level`
  strength= stre["Apis_mellifera",]
  deg = specieslevel(net, index="degree")$`higher level`
  degree= deg["Apis_mellifera",]
  resultados[i,1]<-id.redes[[i]]
  resultados[i,2]<-degree
  resultados[i,3]<-katzz
  resultados[i,4]<-katzscore
  resultados[i,5]<-betweenness
  resultados[i,6]<-closeness
  resultados[i,7]<-strength
}
resultados

contribution2 = data.frame(resultados$degree, resultados$strength,resultados$katzz, resultados$betweenness)
contribution2[is.na(contribution2)] <- 0 #Limpar los NAs

#### Calcular el PCA
model<- prcomp(contribution2, center=TRUE, scale=TRUE)
summary(model)  #Proportion of Variance= % de explicacion de cada eje
plot(model, col="blue") #Importancia de cada PC... PC1, PC2, y etc.
biplot(model) #Importancia de los atributos para cada especie
PC1 <- predict(model)[,1] #Resumen de las especies en cada componente.
PC1_postive =  PC1+abs(min(PC1))+0.0001
write.csv(cbind(resultados, PC1, PC1_postive), "~/Downloads/CENTRALITY.csv", row.names=FALSE)



###OTRAS MÉTRICAS
COL_Apis_mellifera <- unlist(lapply(NETWORKS, function(x){match("Apis_mellifera", colnames(x))}))#Column named Apis_mellifera in each network
NCOL_NETWORK <- lapply(NETWORKS, ncol)#No. of columns in each network


#METRICAS <- matrix(0, length(NETWORKS), 32)
METRICAS <- matrix(0, length(NETWORKS), 4)
rownames(METRICAS) <- names(NETWORKS)
colnames(METRICAS) <- c("robust", "robust_no_Apis.mel", "robust_mean", "robust_std")
#colnames(METRICAS) <- c("H2", "H2_no_Apis.mel", "H2_mean","H2_std", "numberH", "numberH_no_Apis.mel", "numberH_mean", "numberH_std", "numberL", "numberL_no_Apis.mel", "numberL_mean", "numberL_std", "vulne", "vulne_no_Apis.mel", "vulne_mean", "vulne_std", "diversity", "diversity_no_Apis.mel", "diversity_mean", "diversity_std", "niche1", "niche1_no_Apis.mel", "niche1_mean", "niche1_std", "niche2", "niche2_no_Apis.mel", "niche2_mean", "niche2_std", "robust", "robust_no_Apis.mel", "robust_mean", "robust_std")

for(j in 1:112){
  NET <- list()
for(i in (1:NCOL_NETWORK[[j]])){
  NET[[i]] <- NETWORKS[[j]][,-i]
} 
  NET <- lapply(NET, function(x){no <- which(rowSums(x)==0);if(length(no)>0){x[-no,]} else x})
  #metrics, H2
# EXT <- unlist(lapply(NET, function(x){networklevel(x, index="H2")}))
# H2 <- networklevel(NETWORKS[[j]], index="H2")
# H2_no_Apis.mel <- EXT[COL_Apis_mellifera[j]]
# H2_mean <- mean(EXT[-COL_Apis_mellifera[j]])
# H2_std <- sd(EXT[-COL_Apis_mellifera[j]])
# 
# EXT <- unlist(lapply(NET, function(x) grouplevel(x, index="number of species", level="higher")))
# numberH <- grouplevel(NETWORKS[[j]], index="number of species", level="higher")
# numberH_no_Apis.mel <- EXT[COL_Apis_mellifera[j]]
# numberH_mean <- mean(EXT[-COL_Apis_mellifera[j]])
# numberH_std <- sd(EXT[-COL_Apis_mellifera[j]])
# 
# EXT <- unlist(lapply(NET, function(x) grouplevel(x, index="number of species", level="lower")))
# numberL <- grouplevel(NETWORKS[[j]], index="number of species", level="lower")
# numberL_no_Apis.mel <- EXT[COL_Apis_mellifera[j]]
# numberL_mean <- mean(EXT[-COL_Apis_mellifera[j]])
# numberL_std <- sd(EXT[-COL_Apis_mellifera[j]])
# 
# 
# EXT <- unlist(lapply(NET, function(x) grouplevel(x, index="vulnerability", level="lower")))
# vulne <- grouplevel(NETWORKS[[j]], index="vulnerability", level="lower")
# vulne_no_Apis.mel <- EXT[COL_Apis_mellifera[j]]
# vulne_mean <- mean(EXT[-COL_Apis_mellifera[j]])
# vulne_std <- sd(EXT[-COL_Apis_mellifera[j]])
# 
# EXT <- unlist(lapply(NET, function(x) networklevel(x, index="Shannon diversity")))
# diversity <- networklevel(NETWORKS[[j]], index="Shannon diversity")
# diversity_no_Apis.mel <- EXT[COL_Apis_mellifera[j]]
# diversity_mean <- mean(EXT[-COL_Apis_mellifera[j]])
# diversity_std <- sd(EXT[-COL_Apis_mellifera[j]])
# 
# EXT <- unlist(lapply(NET, function(x) grouplevel(x, index="niche overlap", level="higher")))
# niche1 <- grouplevel(NETWORKS[[j]], index="niche overlap",level="higher")
# niche1_no_Apis.mel <- EXT[COL_Apis_mellifera[j]]
# niche1_mean <- mean(EXT[-COL_Apis_mellifera[j]])
# niche1_std <- sd(EXT[-COL_Apis_mellifera[j]])
# 
# EXT <- unlist(lapply(NET, function(x) grouplevel(x, index="niche overlap", level="lower")))
# niche2 <- grouplevel(NETWORKS[[j]], index="niche overlap",level="lower")
# niche2_no_Apis.mel <- EXT[COL_Apis_mellifera[j]]
# niche2_mean <- mean(EXT[-COL_Apis_mellifera[j]])
# niche2_std <- sd(EXT[-COL_Apis_mellifera[j]])
# 
EXT <- unlist(lapply(NET, function(x) robustness(second.extinct(x, participant="higher", method="degree", nrep=1, details=FALSE))))
robust <- robustness(second.extinct(NETWORKS[[j]], participant="higher", method="degree", nrep=1, details=FALSE))
robust_no_Apis.mel <- EXT[COL_Apis_mellifera[j]]
robust_mean <- mean(EXT[-COL_Apis_mellifera[j]])
robust_std <- sd(EXT[-COL_Apis_mellifera[j]])

METRICAS[j, ] <- c(robust, robust_no_Apis.mel, robust_mean, robust_std)
# METRICAS[j, ] <- c(H2, H2_no_Apis.mel, H2_mean,H2_std, numberH, numberH_no_Apis.mel, numberH_mean, numberH_std, numberL, numberL_no_Apis.mel, numberL_mean, numberL_std, vulne, vulne_no_Apis.mel, vulne_mean, vulne_std, diversity, diversity_no_Apis.mel, diversity_mean, diversity_std, niche1, niche1_no_Apis.mel, niche1_mean, niche1_std, niche2, niche2_no_Apis.mel, niche2_mean, niche2_std, robust, robust_no_Apis.mel, robust_mean, robust_std)
print(j)
}
write.csv(METRICAS, "~/Downloads/METRICAS_APIS_2_1-112.csv", row.names = TRUE)
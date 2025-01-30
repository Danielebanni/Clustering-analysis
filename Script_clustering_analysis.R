library(ggplot2)
library(GGally)
library(mclust)
library(cluster)
library(factoextra)
library(ape)
library(datasetsICR)
library(NbClust)
library(stats)
library(StatMatch)
library(clusterCrit) 
library(teigen)
library(mhca)
library(FSAdata)
library(devtools)
library(gridExtra)
library(grid)
library(mixtools)
library(MASS)
library(ContaminatedMixt)
library(flexmix)
library(biotools)
library(mnormt)
library(fpc)

ggpairs(wreath)
wreath2 <- as.data.frame(wreath)
wreath1 <- as.data.frame(wreath) 
plot(wreath)

res.wreath <- NbClust(data = wreath2, method ="ward.D" ) #serve solo per capire il numero ottimale di cluster
summary(res.wreath)
res.wreath$Best.partition
table(res.wreath$Best.partition)
res.wreath$Best.nc #comparare con pam_clustering
res.wreath$All.index 

par(mfrow=c(2,2))
ch_values <- res.wreath$All.index[,"CH"]  # Colonna dell'indice CH
num_clusters <- 2:15  # Range dei numeri di cluster valutati

# Grafico dell'indice CH
plot(num_clusters, ch_values, type = "b", pch = 19, col = "blue",
     xlab = "Number of clusters", ylab = "CH index",
     main = "Calinski-Harabasz index per number of clusters")

# Evidenzia il numero ottimale di cluster
optimal_clusters <- which.max(ch_values) + 1  # Aggiungi 1 perché l'indice parte da min.nc
points(optimal_clusters, max(ch_values), col = "red", pch = 19, cex = 2)
text(optimal_clusters, max(ch_values), labels = paste("Opt:", optimal_clusters), pos = 3)

# Estrai i valori dell'indice di Calinski-Harabasz (CH)
sl_values <- res.wreath$All.index[,"Silhouette"]  # Colonna dell'indice CH
num_clusters <- 2:15  # Range dei numeri di cluster valutati

# Grafico dell'indice CH
plot(num_clusters, sl_values, type = "b", pch = 19, col = "blue",
     xlab = "Number of clusters", ylab = "Silhouette index",
     main = "Silhouette index per number of clusters")

# Evidenzia il numero ottimale di cluster
optimal_clusters <- which.max(sl_values) + 1  # Aggiungi 1 perché l'indice parte da min.nc
points(optimal_clusters, max(sl_values), col = "red", pch = 19, cex = 2)
text(optimal_clusters, max(sl_values), labels = paste("Opt:", optimal_clusters), pos = 3)

# Estrai i valori dell'indice di Calinski-Harabasz (CH)
KL_values <- res.wreath$All.index[,"KL"]  # Colonna dell'indice CH
num_clusters <- 2:15  # Range dei numeri di cluster valutati

# Grafico dell'indice CH
plot(num_clusters, KL_values, type = "b", pch = 19, col = "blue",
     xlab = "Number of clusters", ylab = "KL index",
     main = "Krzanowski-Lai index per number of clusters")

# Evidenzia il numero ottimale di cluster
optimal_clusters <- which.max(KL_values) + 1  # Aggiungi 1 perché l'indice parte da min.nc
points(optimal_clusters, max(KL_values), col = "red", pch = 19, cex = 2)
text(optimal_clusters, max(KL_values), labels = paste("Opt:", optimal_clusters), pos = 3)

# Estrai i valori dell'indice di Calinski-Harabasz (CH)
dunn_values <- res.wreath$All.index[,"Dunn"]  # Colonna dell'indice CH
num_clusters <- 2:15  # Range dei numeri di cluster valutati

# Grafico dell'indice CH
plot(num_clusters, dunn_values, type = "b", pch = 19, col = "blue",
     xlab = "Number of clusters", ylab = "Dunn index",
     main = "Dunn index per number of clusters")

# Evidenzia il numero ottimale di cluster
optimal_clusters <- which.max(dunn_values) + 1  # Aggiungi 1 perché l'indice parte da min.nc
points(optimal_clusters, max(dunn_values), col = "red", pch = 19, cex = 2)
text(optimal_clusters, max(dunn_values), labels = paste("Opt:", optimal_clusters), pos = 3)

ggh1 <- ggplot(wreath2, aes(x = x1)) +
  geom_histogram(aes(y = ..density..), bins = 100, fill = "skyblue", color = "black", alpha = 0.7) +
  scale_x_continuous(limits = c(-20, 20)) +  # Impostare limiti per l'asse x
  scale_y_continuous(limits = c(0, 0.15)) +  # Impostare limiti per l'asse y
  labs(title = "Distribution of x1", x = "x1", y = "Density") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Istogramma per x2
ggh2 <- ggplot(wreath2, aes(x = x2)) +
  geom_histogram(aes(y = ..density..), bins = 100, fill = "lightcoral", color = "black", alpha = 0.7) +
  scale_x_continuous(limits = c(-20, 20)) +  # Impostare limiti per l'asse x
  scale_y_continuous(limits = c(0, 0.15)) +  # Impostare limiti per l'asse y
  labs(title = "Distribution of x2", x = "x2", y = "Density") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
grid.arrange(ggh1, ggh2, nrow = 1) 

par(mfrow=c(1,1))
BIC_clust1 <- mclustBIC(wreath, G = 1:20)
plot(BIC_clust1)
summary(BIC_clust1)

ICL_clust1 <- mclustICL(wreath, G = 1:20) #usato per cluster ben separati e facili da interpretare
summary(ICL_clust1)
plot(ICL_clust1)

LRT_clust1 <- mclustBootstrapLRT(wreath, modelName = "EEV", maxG = 20)
LRT_clust1
summary(LRT_clust1)
plot(LRT_clust1)
par(mfrow=c(2,2))
plot(LRT_clust1,G=7, hist.col = "grey", hist.border = "lightgrey", breaks = "Scott", 
     col = "forestgreen", lwd = 2, lty = 3, main = NULL)
plot(LRT_clust1,G=8, hist.col = "grey", hist.border = "lightgrey", breaks = "Scott", 
     col = "forestgreen", lwd = 2, lty = 3, main = NULL)
plot(LRT_clust1,G=13, hist.col = "grey", hist.border = "lightgrey", breaks = "Scott", 
     col = "forestgreen", lwd = 2, lty = 3, main = NULL)
plot(LRT_clust1,G=14, hist.col = "grey", hist.border = "lightgrey", breaks = "Scott", 
     col = "forestgreen", lwd = 2, lty = 3, main = NULL)

par(mfrow=c(1,1))
mod2 <- Mclust(wreath, G=14, modelNames = "EEV")
summary(mod2, parameters = TRUE)
plot(mod2)
mod2$loglik
dist_matrix <- dist(wreath2)
cluster.stats(d= dist_matrix, clustering= mod2$classification)

#mahalanobis distances matrix
D <- sqrt( D2.dist(wreath2, cov(wreath2)) )

h1 <- hclust(D, method = 'average')
h2 <- hclust(D, method = 'ward.D2')
par(mfrow = c(1, 2), cex = 0.7)
plot(h1, hang = -1, main = 'UPGMA')
abline(h=0.75, col=1)
plot(h2, hang = -1, main = 'Ward')
abline(h=5, col=1)

h1clust = cutree(h1, k=14)
h2clust = cutree(h2, k=14)

plot(wreath2[,1:2], col= h1clust)
plot(wreath2[,1:2], col= h2clust)

#PAM
par(mfrow = c(1,1))
pam_clustering <- pam(D, k=14, diss=TRUE) 
clusplot(pam_clustering, main = "Cluster PAM", color=TRUE, shade = TRUE, labels = 2, lines = 0)

data_clustered <- data.frame(
  wreath,
  Cluster = as.factor(pam_clustering$clustering)  # Aggiunge i cluster
)
plot(data_clustered[,1:2], col=data_clustered$Cluster)
# Visualizzazione con ggplot2
ggplot(data_clustered, aes(x = x1, y = x2, color = Cluster)) +
  geom_point(size = 3) +  # Punti colorati per cluster
  labs(
    title = "Clustering with PAM ",
    x = "x1",
    y = "x2"
  ) +
  theme_minimal() #comparazione con mod2 ha solo una osservazione classificata diversamente
table(data_clustered$Cluster)

#mhca package
n<- nrow(wreath2)
mh<-mhclust(wreath2 )

cmh<-cutree(mh,k=14)
cluster.stats(d= dist_matrix, clustering= cmh)

# Genera una palette di k colori distinti

col_palette <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", 
                 "#800000", "#808000", "#008000", "#800080", "#008080", "#C0C0C0", 
                 "#FF6347", "#D2691E") 

# Assegna un colore unico a ciascun cluster
col_cluster <- col_palette[cmh]

#plot(wreath2[,1],wreath2[,2],col=col_cluster,main='Mahalanobis HCA',frame=FALSE)
plot(wreath2[,1],wreath2[,2],col=cmh,main='Mahalanobis HCA',frame=FALSE)
plot(mh,labels=FALSE,main='Dendrogram of MHCA')
y<-min(mh$height)-diff(range(mh$height))/10 
text(1:n,y,(1:n)[mh$order],col=cmh[mh$order],srt=90)

plot(wreath2[,1],wreath2[,2],col=cmh,main='Mahalanobis HCA',frame=FALSE)


#hierarchical clustering
hc.single1 = hclust(x.dist1, method="single")
hc.complete1 = hclust(x.dist1, method="complete")
hc.ward1 =hclust(x.dist1, method= "ward.D")
hc.average1 =hclust(x.dist1, method= "average")

par(mfrow= c(2,2))
plot(hc.single1)
plot(hc.complete1)
plot(hc.ward1)
plot(hc.average1)

hc.single2 = hclust(x.dist2, method="single")
hc.complete2 = hclust(x.dist2, method="complete")
hc.ward2 =hclust(x.dist2, method= "ward.D")
hc.average2 =hclust(x.dist2, method= "average")

par(mfrow= c(2,2))
plot(hc.single2)
plot(hc.complete2)
plot(hc.ward2)
plot(hc.average2)

hc.single3 = hclust(x.dist3, method="single")
hc.complete3 = hclust(x.dist3, method="complete")
hc.ward3 =hclust(x.dist3, method= "ward.D")
hc.average3 =hclust(x.dist3, method= "average")

par(mfrow=c(2,2))
plot(hc.single3)
plot(hc.complete3)
plot(hc.ward3)
plot(hc.average3)


hiclust11 = cutree(hc.single1, k=14)
hiclust12 = cutree(hc.complete1, k=14)
hiclust13 = cutree(hc.ward1, k=14)
hiclust14 = cutree(hc.average1, k=14)

par(mfrow=c(2,2))
plot(wreath2, col= hiclust11, main= "Single")
plot(wreath2, col =hiclust12, main="Complete")
plot(wreath2 , col= hiclust13, main="Ward.D")
plot(wreath2 , col= hiclust14, main= "Average")
table(hiclust11)
table(hiclust12)
table(hiclust13)
table(hiclust14)

cluster.stats(d= dist_matrix, clustering= hiclust11)
cluster.stats(d= dist_matrix, clustering= hiclust12)
cluster.stats(d= dist_matrix, clustering= hiclust13)
cluster.stats(d= dist_matrix, clustering= hiclust14)


hiclust21 = cutree(hc.single2, k=14)
hiclust22 = cutree(hc.complete2, k=14)
hiclust23 = cutree(hc.ward2, k=14)
hiclust24 = cutree(hc.average2, k=14)

par(mfrow=c(2,2))
plot(wreath2, col= hiclust21, main= "Single")
plot(wreath2, col =hiclust22, main= "Complete")
plot(wreath2 , col= hiclust23, main= "Ward.D")
plot(wreath2 , col= hiclust24, main= "Average")
table(hiclust21)
table(hiclust22)
table(hiclust23)
table(hiclust24)
cluster.stats(d= dist_matrix, clustering= hiclust21)
cluster.stats(d= dist_matrix, clustering= hiclust22)
cluster.stats(d= dist_matrix, clustering= hiclust23)
cluster.stats(d= dist_matrix, clustering= hiclust24)


hiclust31 = cutree(hc.single3, k=14)
hiclust32 = cutree(hc.complete3, k=14)
hiclust33 = cutree(hc.ward3, k=14)
hiclust34 = cutree(hc.average3, k=14)

par(mfrow=c(2,2))
plot(wreath2, col= hiclust31, main= "Single")
plot(wreath2, col =hiclust32, main= "Complete")
plot(wreath2 , col= hiclust33, main= "Ward.D")
plot(wreath2 , col= hiclust34, main= "Average")
table(hiclust31)
table(hiclust32)
table(hiclust33)
table(hiclust34)
cluster.stats(d= dist_matrix, clustering= hiclust31)
cluster.stats(d= dist_matrix, clustering= hiclust32)
cluster.stats(d= dist_matrix, clustering= hiclust33)
cluster.stats(d= dist_matrix, clustering= hiclust34)

#teigen
res2 <- teigen(wreath2, Gs=c(8,14), models = "all", verbose= FALSE)

res2$allbic
summary(res2)
res2$bestmodel
res2$classification
table(res2$classification)
res2$modelname


res2$logl

CUCC_14 <- teigen(wreath2, 14, models = "CUCC", verbose = FALSE)
plot(CUCC_14)
cluster.stats(d= D, clustering= CUCC_14$classification)
cluster.stats(d= dist_matrix, clustering= CUCC_14$classification)


#Matrice di covarianza e centroide
cov_matrix <- cov(wreath)          
center <- colMeans(wreath)         

# hierarchical clustering using mahalanobis distance
mahal_dist_matrix <- as.matrix(dist(wreath, method = "mahalanobis", diag = TRUE)) 
mahal_dist <- as.dist(mahal_dist_matrix)

ha1 = hclust(mahal_dist, method="single")
ha2 = hclust(mahal_dist, method="complete")
ha3 =hclust(mahal_dist, method= "ward.D")
ha4 = hclust(mahal_dist, method="average")

par(mfrow =c(2,2))
plot(ha1)
plot(ha2)
plot(ha3)
plot(ha4)

haclust1 = cutree(ha1, k=14)
haclust2 = cutree(ha2, k=14)
haclust3 = cutree(ha3, k=14)
haclust4 = cutree(ha4, k=14)
par(mfrow =c(2,2))
plot(wreath2, col= haclust1, main= "Single")
plot(wreath2, col =haclust2, main= "Complete")
plot(wreath2 , col= haclust3, main= "Ward.D")
plot(wreath2 , col= haclust4, main= "Average")
table(haclust1)
table(haclust2)
table(haclust3)
table(haclust4)

cluster.stats(d= dist_matrix, clustering= haclust1)
cluster.stats(d= dist_matrix, clustering= haclust2)
cluster.stats(d= dist_matrix, clustering= haclust3)
cluster.stats(d= dist_matrix, clustering= haclust4)

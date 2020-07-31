#Librerias
library(cluster)
library(dplyr)
library(ggplot2)
library(readr)
library(Rtsne)
library(gridExtra)
library(factoextra)
library(NbClust)


############################## Parte 1: Lectura de dataset ############################## 
#Lectura del data
ruta <- file.choose()
hepatitis <- read.table(ruta,sep=",", na.strings = c("?"))

#Asignación de los nombres de cada columa, respecto a las variables indicadas.
clase.variables <- c("CLASS","AGE","SEX","STEROID","ANTIVIRALS","FATIGUE","MALAISE",
                     "ANOREXIA","LIVER_BIG","LIVER_FIRM","SPLEEN_PALPABLE","SPIDERS",
                     "ASCITES","VARICES","BILIRUBIN","ALK_PHOSPHATE","SGOT","ALBUMIN",
                     "PROTIME","HISTOLOGY")

colnames(hepatitis) <- clase.variables
cantidad.nulos <- sapply(hepatitis, function(x) sum(is.na(x))) #Contando cantidad de NA por variable

############################## Parte 2: Filtro y eliminación de variables ############################## 

#Parte 2.1: Eliminando variables que no afectan el estudio
hepatitis$VARICES <- NULL #Eliminando columna VARICES
hepatitis$CLASS <- NULL #Eliminando columna CLASS

hepatitis.datos.original.kmeans <- hepatitis                    #Guardando datos originales para uso de kmeans
hepatitis.datos.original.pam <- hepatitis                       #Guardando datos originales para uso de pam

#Parte 2.2: Modificando datos de variables con valor NA (para el uso de kmeans):

#Varibles discreta (3 = no hay información)
hepatitis$STEROID[is.na(hepatitis$STEROID)] <- 3
hepatitis$FATIGUE[is.na(hepatitis$FATIGUE)] <- 3
hepatitis$MALAISE[is.na(hepatitis$MALAISE)] <- 3
hepatitis$ANOREXIA[is.na(hepatitis$ANOREXIA)] <- 3
hepatitis$LIVER_BIG[is.na(hepatitis$LIVER_BIG)] <- 3
hepatitis$LIVER_FIRM[is.na(hepatitis$LIVER_FIRM)] <- 3
hepatitis$SPLEEN_PALPABLE[is.na(hepatitis$SPLEEN_PALPABLE)] <- 3
hepatitis$SPIDERS[is.na(hepatitis$SPIDERS)] <- 3
hepatitis$ASCITES[is.na(hepatitis$ASCITES)] <- 3

#Variables continuas con NA tienen el minimo valor posible
hepatitis$BILIRUBIN[is.na(hepatitis$BILIRUBIN)] <- 0
hepatitis$ALK_PHOSPHATE[is.na(hepatitis$ALK_PHOSPHATE)] <- 0
hepatitis$SGOT[is.na(hepatitis$SGOT)] <- 0
hepatitis$ALBUMIN[is.na(hepatitis$ALBUMIN)] <- 0
hepatitis$PROTIME[is.na(hepatitis$PROTIME)] <- 0

############################## Parte 3: Determinando n° Cluster optimo ############################## 

#Parte 3.1: cluster optimo para kmeans

optimo <- kmeans(hepatitis, center = 1)$betweenss
for(i in 1:10) optimo[i] <- kmeans(hepatitis, centers = i)$betweenss
plot(1:10, optimo, type = "b", xlab = "número de cluster", ylab = "suma de cuadrados inter grupos")
fviz_nbclust(hepatitis,kmeans,method = "wss")
fviz_nbclust(hepatitis,kmeans,method = "silhouette")
fviz_nbclust(hepatitis,kmeans,method = "gap_stat")
resnumclust <- NbClust(hepatitis,distance="euclidean",min.nc=2,max.nc=10,method="kmeans")
fviz_nbclust(resnumclust)

#Parte 3.2: cluster optimo para PAM
sil_width <- c(NA)
for(i in 1:8){  
  pam_fit <- pam(gower_dist, diss = TRUE, k = i)  
  sil_width[i] <- pam_fit$silinfo$avg.width  
}
plot(1:8, sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(1:8, sil_width)


############################## Parte 4: Aplicando K-Means ##################################### 

#Parte 4.1: utilización de kmeans con columna PROTIME
hepatitis.datos.original.kmeans <- hepatitis
hepatitis.datos.original.kmeans <- scale(hepatitis.datos.original.kmeans)           #Escalando datos de la tabla
kmedia.contenido <- kmeans(hepatitis.datos.original.kmeans,2)
#Graficando resultados
clusplot(hepatitis.datos.original.kmeans, kmedia.contenido$cluster)
grafica1 <- fviz_cluster(kmedia.contenido,data=hepatitis.datos.original.kmeans,ellipse.type = "euclid",repel=TRUE,star.plot = TRUE)
grafica2 <- fviz_cluster(kmedia.contenido,data=hepatitis.datos.original.kmeans,ellipse.type = "norm")

#Parte 4.2: utilización de kmeans sin columna PROTIME
hepatitis.datos.original.kmeans <- hepatitis
hepatitis.datos.original.kmeans$PROTIME <- NULL
hepatitis.datos.original.kmeans <- scale(hepatitis.datos.original.kmeans)           #Escalando datos de la tabla
kmedia.contenido <- kmeans(hepatitis.datos.original.kmeans,2)
#Graficando resultados
clusplot(hepatitis.datos.original.kmeans, kmedia.contenido$cluster)
grafica3 <- fviz_cluster(kmedia.contenido,data=hepatitis.datos.original.kmeans,ellipse.type = "euclid",repel=TRUE,star.plot = TRUE)
grafica4 <- fviz_cluster(kmedia.contenido,data=hepatitis.datos.original.kmeans,ellipse.type = "norm")
#### Preguntar si colores de grafico influyen al saber si un paciente vive o muere

#Parte 4.3: utilización de pam

#Obteniendo las distancias gower
gower_dist <- daisy(hepatitis.datos.original.pam, metric = "gower")
gower_mat <- as.matrix(gower_dist)

#Calculando kmeans
k <- 2
pam_fit <- pam(gower_dist, diss = TRUE, k)
pam_results <- hepatitis.datos.original.pam %>%
  mutate(cluster = pam_fit$clustering) %>%
  group_by(cluster) %>%
  do(the_summary = summary(.))
pam_results$the_summary

tsne_obj <- Rtsne(gower_dist, is_distance = TRUE)
tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering))
ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = cluster))





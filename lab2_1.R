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

hepatitis <- subset( hepatitis, select = -c(VARICES,LIVER_FIRM,CLASS ) ) #quitando variables VARICES, LIVER_FIRM y CLASS


hepatitis.datos.original.kmeans <- hepatitis                    #Guardando datos originales para uso de kmeans
hepatitis.datos.original.pam <- hepatitis                       #Guardando datos originales para uso de pam

#Parte 2.2: Modificando datos de variables con valor NA (para el uso de kmeans):

#Varibles discreta (3 = no hay información)
hepatitis$STEROID[is.na(hepatitis$STEROID)] <- 3
hepatitis$FATIGUE[is.na(hepatitis$FATIGUE)] <- 3
hepatitis$MALAISE[is.na(hepatitis$MALAISE)] <- 3
hepatitis$ANOREXIA[is.na(hepatitis$ANOREXIA)] <- 3
hepatitis$LIVER_BIG[is.na(hepatitis$LIVER_BIG)] <- 3

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

optimoHepatitisKmeans <- function(hepatitis) {
  optimo <- kmeans(hepatitis, center = 1)$betweenss
  for(i in 1:10) optimo[i] <- kmeans(hepatitis, centers = i)$betweenss
  plot(1:10, optimo, type = "b", xlab = "número de cluster", ylab = "suma de cuadrados inter grupos")
  fviz_nbclust(hepatitis,kmeans,method = "wss")
  fviz_nbclust(hepatitis,kmeans,method = "silhouette")
  fviz_nbclust(hepatitis,kmeans,method = "gap_stat")
  resnumclust <- NbClust(hepatitis,distance="euclidean",min.nc=2,max.nc=10,method="kmeans")
  
}


#Parte 3.2: cluster optimo para PAM
optimoHepatitisPam <- function(hepatitis){
  fviz_nbclust(x = hepatitis.datos.original.pam, FUNcluster = pam, method = "wss", k.max = 10,
               diss = dist(hepatitis.datos.original.pam, method = "manhattan"))
}

############################## Parte 4: Aplicando K-Means ##################################### 


#Parte 4.1: utilización de kmeans con columna todas las columnas
hepatitis.datos.original.kmeans <- hepatitis
#hepatitis.datos.original.kmeans <- scale(hepatitis.datos.original.kmeans)           #Escalando datos de la tabla
cat("Obteniendo Clúster óptimo: todas las variables \n")
optimoHepatitisKmeans(hepatitis.datos.original.kmeans)                               #Obteniendo clúster óptimo
k<-2
kmedia.contenido <- kmeans(hepatitis.datos.original.kmeans,k)
#Graficando resultados
grafica1 <- fviz_cluster(kmedia.contenido,data=hepatitis.datos.original.kmeans,ellipse.type = "euclid",repel=TRUE,star.plot = TRUE)
grafica2 <- fviz_cluster(kmedia.contenido,data=hepatitis.datos.original.kmeans,ellipse.type = "norm")
cat("#############################################################################################\n")

#Parte 4.2: utilización de kmeans sin columna PROTIME
hepatitis.datos.original.kmeans.sin.protime <- hepatitis
hepatitis.datos.original.kmeans.sin.protime$PROTIME <- NULL
hepatitis.datos.original.kmeans.sin.protime <- scale(hepatitis.datos.original.kmeans.sin.protime)           #Escalando datos de la tabla
cat("Obteniendo Clúster óptimo: sin variable PROTIME \n")
optimoHepatitisKmeans(hepatitis.datos.original.kmeans.sin.protime)                                          #Obteniendo clúster óptimo
k<-3
kmedia.contenido <- kmeans(hepatitis.datos.original.kmeans.sin.protime,k)
#Graficando resultados
grafica3 <- fviz_cluster(kmedia.contenido,data=hepatitis.datos.original.kmeans.sin.protime,ellipse.type = "euclid",repel=TRUE,star.plot = TRUE)
grafica4 <- fviz_cluster(kmedia.contenido,data=hepatitis.datos.original.kmeans.sin.protime,ellipse.type = "norm")
cat("#############################################################################################\n")

#Parte 4.4: utilización de kmeans sin columna Alk phosphate
hepatitis.datos.original.kmeans.sin.alk <- hepatitis
hepatitis.datos.original.kmeans.sin.alk$ALK_PHOSPHATE <- NULL
hepatitis.datos.original.kmeans.sin.alk <- scale(hepatitis.datos.original.kmeans.sin.alk)           #Escalando datos de la tabla
cat("Obteniendo Clúster óptimo: sin variable ALK PHOSPHATE \n")
optimoHepatitisKmeans(hepatitis.datos.original.kmeans.sin.protime)                                  #Obteniendo clúster óptimo
k<-3
kmedia.contenido <- kmeans(hepatitis.datos.original.kmeans.sin.alk,k)
#Graficando resultados
grafica5 <- fviz_cluster(kmedia.contenido,data=hepatitis.datos.original.kmeans.sin.alk,ellipse.type = "euclid",repel=TRUE,star.plot = TRUE)
grafica6 <- fviz_cluster(kmedia.contenido,data=hepatitis.datos.original.kmeans.sin.alk,ellipse.type = "norm")
cat("#############################################################################################\n")

#Parte 4.5: utilización de kmeans sin columna SEX
hepatitis.datos.original.kmeans.sin.sex <- hepatitis
hepatitis.datos.original.kmeans.sin.sex$SEX <- NULL
hepatitis.datos.original.kmeans.sin.sex <- scale(hepatitis.datos.original.kmeans.sin.sex)           #Escalando datos de la tabla
print("Obteniendo Clúster óptimo: sin variable SEX \n")
optimoHepatitisKmeans(hepatitis.datos.original.kmeans.sin.protime)                                  #Obteniendo clúster óptimo
k<-3
kmedia.contenido <- kmeans(hepatitis.datos.original.kmeans.sin.sex,k)
#Graficando resultados
grafica7 <- fviz_cluster(kmedia.contenido,data=hepatitis.datos.original.kmeans.sin.sex,ellipse.type = "euclid",repel=TRUE,star.plot = TRUE)
grafica8 <- fviz_cluster(kmedia.contenido,data=hepatitis.datos.original.kmeans.sin.sex,ellipse.type = "norm")
cat("#############################################################################################\n")

#Parte 4.6: utilización de kmeans sin columnas PROTIME, ALK PHOSPHATE y sex
hepatitis.datos.original.kmeans.sin.col <- hepatitis
hepatitis.datos.original.kmeans.sin.col$PROTIME <- NULL
hepatitis.datos.original.kmeans.sin.col$ALK_PHOSPHATE <- NULL
hepatitis.datos.original.kmeans.sin.col$SEX <- NULL
#hepatitis <- subset( hepatitis.datos.original.kmeans.sin.col, select = -c(PROTIME,ALK_PHOSPHATE,SEX ) ) #quitando variables PROTIME, ALK_PHOSPHATE y SEX
hepatitis.datos.original.kmeans.sin.col <- scale(hepatitis.datos.original.kmeans.sin.col)                #Escalando datos de la tabla
cat("Obteniendo Clúster óptimo: sin variable PROTIME,ALK PHOSPHATE y SEX \n")
optimoHepatitisKmeans(hepatitis.datos.original.kmeans.sin.protime)                                        #Obteniendo clúster óptimo
k<-3
kmedia.contenido <- kmeans(hepatitis.datos.original.kmeans.sin.col,k)
#Graficando resultados
grafica9 <- fviz_cluster(kmedia.contenido,data=hepatitis.datos.original.kmeans.sin.col,ellipse.type = "euclid",repel=TRUE,star.plot = TRUE)
grafica10 <- fviz_cluster(kmedia.contenido,data=hepatitis.datos.original.kmeans.sin.col,ellipse.type = "norm")
cat("#############################################################################################\n")


############################## Parte 5: Aplicando pam ##################################### 

k <- 2


#Parte 5.1: utilización de pam con todas las columnas

#Obteniendo las distancias gower
optimoHepatitisPam(hepatitis.datos.original.pam)
set.seed(123)
gower_dist <- daisy(hepatitis.datos.original.pam, metric = "gower")
gower_mat <- as.matrix(gower_dist)
#Calculando kmeans
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


#Parte 5.1: utilización de pam sin PROTIME

hepatitis.datos.original.pam.sin.protime <- hepatitis
hepatitis.datos.original.pam.sin.protime$PROTIME <- NULL
optimoHepatitisPam(hepatitis.datos.original.pam.sin.protime)    #Obteniendo nuevo cluste optimo

#Obteniendo las distancias gower
set.seed(123)
gower_dist <- daisy(hepatitis.datos.original.pam.sin.protime, metric = "gower")
gower_mat <- as.matrix(gower_dist)
#Calculando kmeans

pam_fit <- pam(gower_dist, diss = TRUE, k)
pam_results <- hepatitis.datos.original.pam.sin.protime %>%
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

#Parte 5.1: utilización de pam sin PROTIME

hepatitis.datos.original.pam.sin.protime <- hepatitis
hepatitis.datos.original.pam.sin.protime$PROTIME <- NULL
optimoHepatitisPam(hepatitis.datos.original.pam.sin.protime)    #Obteniendo nuevo cluste optimo

#Obteniendo las distancias gower
set.seed(123)
gower_dist <- daisy(hepatitis.datos.original.pam.sin.protime, metric = "gower")
gower_mat <- as.matrix(gower_dist)
#Calculando kmeans

pam_fit <- pam(gower_dist, diss = TRUE, k)
pam_results <- hepatitis.datos.original.pam.sin.protime %>%
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


#Parte 5.2: utilización de pam sin ALK PHOSPHATE

hepatitis.datos.original.pam.sin.alk <- hepatitis
hepatitis.datos.original.pam.sin.alk$ALK_PHOSPHATE <- NULL
optimoHepatitisPam(hepatitis.datos.original.pam.sin.alk)    #Obteniendo nuevo cluste optimo

#Obteniendo las distancias gower
set.seed(123)
gower_dist <- daisy(hepatitis.datos.original.pam.sin.alk, metric = "gower")
gower_mat <- as.matrix(gower_dist)
#Calculando kmeans

pam_fit <- pam(gower_dist, diss = TRUE, k)
pam_results <- hepatitis.datos.original.pam.sin.alk %>%
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

#Parte 5.3: utilización de pam sin SEX

hepatitis.datos.original.pam.sin.sex <- hepatitis
hepatitis.datos.original.pam.sin.sex$SEX <- NULL
optimoHepatitisPam(hepatitis.datos.original.pam.sin.sex)    #Obteniendo nuevo cluste optimo

#Obteniendo las distancias gower
set.seed(123)
gower_dist <- daisy(hepatitis.datos.original.pam.sin.sex, metric = "gower")
gower_mat <- as.matrix(gower_dist)
#Calculando kmeans

pam_fit <- pam(gower_dist, diss = TRUE, k)
pam_results <- hepatitis.datos.original.pam.sin.sex %>%
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


#Parte 5.4: utilización de pam sin PROTIME, ALK PHOSPHATE y SEX

hepatitis.datos.original.pam.sin.col <- hepatitis
hepatitis.datos.original.pam.sin.col$PROTIME <- NULL
hepatitis.datos.original.pam.sin.col$ALK_PHOSPHATE <- NULL
hepatitis.datos.original.pam.sin.col$SEX <- NULL
optimoHepatitisPam(hepatitis.datos.original.pam.sin.col)    #Obteniendo nuevo cluste optimo

#Obteniendo las distancias gower
set.seed(123)
gower_dist <- daisy(hepatitis.datos.original.pam.sin.col, metric = "gower")
gower_mat <- as.matrix(gower_dist)
#Calculando kmeans

pam_fit <- pam(gower_dist, diss = TRUE, k)
pam_results <- hepatitis.datos.original.pam.sin.col %>%
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








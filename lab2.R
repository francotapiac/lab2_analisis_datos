library(ggplot2)
library(gridExtra)
library(cluster)
library(factoextra)
library(NbClust)
library(dplyr)
############################## Parte 1: Lectura de dataset ############################## 
#Lectura del data
ruta <- file.choose()
hepatitis <- read.table(ruta,sep=",", na.strings = c("?"))

#Asignación de los nombres de cada columa, respecto a las variables indicadas.
clase_variables <- c("CLASS","AGE","SEX","STEROID","ANTIVIRALS","FATIGUE","MALAISE",
                     "ANOREXIA","LIVER_BIG","LIVER_FIRM","SPLEEN_PALPABLE","SPIDERS",
                     "ASCITES","VARICES","BILIRUBIN","ALK_PHOSPHATE","SGOT","ALBUMIN",
                     "PROTIME","HISTOLOGY")

colnames(hepatitis) <- clase_variables
sapply(hepatitis, function(x) sum(is.na(x))) #Contando cantidad de NA por variable
hepatitis.datos.original <- hepatitis        #Guardando datos originales 

############################## Parte 2: Filtro y eliminación de variables ############################## 

#Eliminando variables que no afectan el estudio
hepatitis$VARICES <- NULL #Eliminando columna VARICES
hepatitis$CLASS <- NULL #Eliminando columna CLASS

#Modificando datos de variables importantes con NA:

#Varibles discreta (3 = no se información)
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

hepatitis <- scale(hepatitis)

############################## Parte 3:Analizando los datos ##############################


#Determinando n° cluster óptimo
set.seed(80)
optimo <- kmeans(hepatitis, center = 1)$betweenss
for(i in 1:10) optimo[i] <- kmeans(hepatitis, centers = i)$betweenss
plot(1:10, optimo, type = "b", xlab = "número de cluster", ylab = "suma de cuadrados inter grupos")

fviz_nbclust(hepatitis,kmeans,method = "wss")
fviz_nbclust(hepatitis,kmeans,method = "silhouette")
fviz_nbclust(hepatitis,kmeans,method = "gap_stat")

resnumclust <- NbClust(hepatitis,distance="euclidean",min.nc=2,max.nc=10,method="kmeans")
fviz_nbclust(resnumclust)

g.dist = daisy(hepatitis, metric="gower", type=list(logratio=2))
############################## Parte 4: Obteniendo k-media ##############################
#Escalado de las variables
#hepatitis <- as.data.frame(scale(hepatitis[,1:18]))

#hepatitis.sin.acceso.nulo <- na.omit(hepatitis) #Eliminación de cada fila con uno o más atributos NA
kmedia.contenido <- kmeans(hepatitis,2)

############################## Parte 5: Gráficando resultWados ##############################
clusplot(hepatitis, kmedia.contenido$cluster)
grafica1 <- fviz_cluster(kmedia.contenido,data=hepatitis,ellipse.type = "euclid",repel=TRUE,star.plot = TRUE)
grafica2 <- fviz_cluster(kmedia.contenido,data=hepatitis,ellipse.type = "norm")

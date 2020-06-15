library(ggplot2)
library(gridExtra)
library(cluster)

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

############################## Parte 1: Filtro y eliminación de variables ############################## 

#Eliminando variables que no afectan el estudio
hepatitis$PROTIME <- NULL #Eliminando columna PROTIME

#Modificando datos de variables importantes con NA:

#Varibles discreta (3 = no se tienen datos)
hepatitis$STEROID[is.na(hepatitis$STEROID)] <- 3
hepatitis$FATIGUE[is.na(hepatitis$FATIGUE)] <- 3
hepatitis$MALAISE[is.na(hepatitis$MALAISE)] <- 3
hepatitis$ANOREXIA[is.na(hepatitis$ANOREXIA)] <- 3
hepatitis$LIVER_BIG[is.na(hepatitis$LIVER_BIG)] <- 3

#Variables continuas con NA tienen el maximo valor posible

#Escalado de las variables
hepatitis.escalado <- as.data.frame(scale(hepatitis[,1:19]))
hepatitis.sin.acceso.nulo <- na.omit(hepatitis) #Eliminación de cada fila con uno o más atributos NA


############################## Parte 2: ##############################


#Determinando n° cluster óptimo
set.seed(80)
optimo <- kmeans(hepatitis.sin.acceso.nulo, center = 1)$betweenss
for(i in 2:10) optimo[i] <- kmeans(hepatitis.sin.acceso.nulo, centers = i)$betweenss
plot(1:10, optimo, type = "b", xlab = "número de cluster", ylab = "suma de cuadrados inter grupos")




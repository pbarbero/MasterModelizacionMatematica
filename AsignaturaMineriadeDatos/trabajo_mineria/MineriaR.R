library(TH.data)
library(ggplot2)
library(gridExtra)
data("wpbc")
View(wpbc)

summary(wpbc)

p1 <- qplot(as.factor(status),data=wpbc,geom="bar",
            fill=I("black"), alpha=I(0.5),color=I("black"),
            main="Estado",xlab="",ylab="")
p2 <- qplot(as.factor(time),data=wpbc,geom="bar",
            fill=I("blue"), alpha=I(0.5),color=I("black"),
            main="Tiempo",xlab="",
            ylab="")

grid.arrange(p1,p2,nrow=1)

qplot(wpbc$time, wpbc$status, geom="line")


# Creamos la nueva variable. 
wpbc$clase = factor(wpbc$clase, labels = c("Positivo", "Negativo"))
attach(wpbc)
wpbc$clase[status == "R" & time <= 24] <- "Positivo"
wpbc$clase[status == "R" & time > 24] <- "Negativo"
wpbc$clase[status == "N"] <- "Negativo"
detach(wpbc)
View(wpbc)

qplot(wpbc$clase)

summary(wpbc)


# Empezamos a clasificar. Hay 4 NA's en pnodes, si más adelante dan problemas los 
# obviaremos. 1 en Positivo y 3 en Negativo.

# Si tenemos que obviarlos usaremos:
# use = "pairwise.complete.obs"

# Análisis descriptivo e importancia de cada variable.

# Normalidad de las variables
shapiro.test(wpbc$time)
shapiro.test(wpbc$mean_radius)
shapiro.test(wpbc$mean_texture)
shapiro.test(wpbc$mean_perimeter)
shapiro.test(wpbc$mean_area)
shapiro.test(wpbc$mean_smoothness)
shapiro.test(wpbc$mean_compactness)
shapiro.test(wpbc$mean_concavity)
shapiro.test(wpbc$mean_concavepoints)
shapiro.test(wpbc$mean_symmetry)
shapiro.test(wpbc$mean_fractaldim)

shapiro.test(wpbc$SE_radius)
shapiro.test(wpbc$SE_texture)
shapiro.test(wpbc$SE_perimeter)
shapiro.test(wpbc$SE_area)
shapiro.test(wpbc$SE_smoothness)
shapiro.test(wpbc$SE_compactness)
shapiro.test(wpbc$SE_concavity)
shapiro.test(wpbc$SE_concavepoints)
shapiro.test(wpbc$SE_symmetry)
shapiro.test(wpbc$SE_fractaldim)

shapiro.test(wpbc$worst_radius)
shapiro.test(wpbc$worst_texture)
shapiro.test(wpbc$worst_perimeter)
shapiro.test(wpbc$worst_area)
shapiro.test(wpbc$worst_smoothness)
shapiro.test(wpbc$worst_compactness)
shapiro.test(wpbc$worst_concavity)
shapiro.test(wpbc$worst_concavepoints)
shapiro.test(wpbc$worst_symmetry)
shapiro.test(wpbc$worst_fractaldim)

# Shapiro Wilk no nos asegura la normalidad de la variable worst_texture,
# sin embargo al realizar el histograma y la qqplot observamos que sí que lo es.
library(nortest)
hist(wpbc$worst_texture)
qqnorm(wpbc$worst_texture)

# Todas las variables son normales, podemos obtener las correlaciones que hay entre
# ellas.

library("psych")

wpbc1<-wpbc[,-c(1,2,35)]
x <- subset(wpbc1) # Omit missing values
r <- cor(x,  use = "pairwise.complete.obs") # Correlation matrix
r2 <- r^2 # Squared correlation coefficients
i <- solve(r) # Inverse matrix of correlation matrix
d <- diag(i) # Diagonal elements of inverse matrix
p2 <- (-i/sqrt(outer(d, d)))^2 # Squared partial correlation coefficients
diag(r2) <- diag(p2) <- 0 # Delete diagonal elements
KMO <- sum(r2)/(sum(r2)+sum(p2))
MSA <- colSums(r2)/(colSums(r2)+colSums(p2))

MSA
# Son candidatos a salir del análisis las variables SE_texture y pnodes
# ya que son poco relevantes.
KMO
# el KMO es bien

wpbc2<-wpbc1[,-c(2,12,22,32)]
x <- subset(wpbc2) # Omit missing values
r <- cor(x,  use = "pairwise.complete.obs") # Correlation matrix
r2 <- r^2 # Squared correlation coefficients
i <- solve(r) # Inverse matrix of correlation matrix
d <- diag(i) # Diagonal elements of inverse matrix
p2 <- (-i/sqrt(outer(d, d)))^2 # Squared partial correlation coefficients
diag(r2) <- diag(p2) <- 0 # Delete diagonal elements
KMO <- sum(r2)/(sum(r2)+sum(p2))
MSA <- colSums(r2)/(colSums(r2)+colSums(p2))

# Quitamos todas las texture y pnodes y el KMO mejora un poquito.
MSA
KMO


# Matriz de correlaciones
matriz.correlaciones<-cor(wpbc2, use = "pairwise.complete.obs")
matriz.correlaciones

library(reshape2)
qplot(x=Var1, y=Var2, data=melt(cor(wpbc2, use = "pairwise.complete.obs")), fill=abs(value), geom="tile", 
      main="Mapa de calor de la matriz de correlación")

det<-det(matriz.correlaciones)
det
# Están poco correladas entre ellas

# Prueba de esfericidad de Bartlett
n<-198 #numero de casos 
p<-28 #numero de variables
est<- (n-1-((2*p+5)/6))*log(det)
est
gl<-(p^2-p)/2
gl
pchisq(c(est),df=gl, lower.tail=FALSE)

# No podemos aplicar AF porque la matriz de correlaciones es la identidad
# Tenemos que seguir quitando variables
# Quitamos tsize, no nos aporta nada

wpbc3<-wpbc2[,-c(28)]
# Matriz de correlaciones
matriz.correlaciones<-cor(wpbc3, use = "pairwise.complete.obs")
matriz.correlaciones

library(reshape2)
qplot(x=Var1, y=Var2, data=melt(cor(wpbc3, use = "pairwise.complete.obs")), fill=abs(value), geom="tile", 
      main="Mapa de calor de la matriz de correlación")

det<-det(matriz.correlaciones)
det
n<-198 #numero de casos 
p<-27 #numero de variables
est<- (n-1-((2*p+5)/6))*log(det)
est
gl<-(p^2-p)/2
gl
pchisq(c(est),df=gl, lower.tail=FALSE)

# Observamos que Area, Perimetro y Radio están totalmente correlados
# En el resto de variables no se observa esa correlación
# No podemos aplicar Análisis Factorial para preveer la estructura de asociación,
# están las variables tan poco correladas que la mayoría de factores que obtengamos
# estaran formados por una variable sólo.




# El conjunto de datos no aconseja realizar un conjunto de entrenamiento y otro de
# validación, por tanto veremos qué es más aconsejable en nuestro problema si
# realizar CV ó Bootstrap. 

# Creamos un nuevo conjunto de datos con las variables continuas y la nueva
# variable clase, que son las que nos interesan.

datos <- wpbc[,-c(1,2)]
summary(datos)
# Haremos las técnicas de clasificación con ambos procedimientos y observaremos
# cual de ellos es mejor para luego obtener las conclusiones.

# ÁRBOLES DE CLASIFICACIÓN


#arbol con metodo class
cfit<- rpart(clase ~ mean_radius + mean_texture +  mean_perimeter + mean_area + 
                     mean_smoothness + mean_compactness + mean_concavity + 
                     mean_concavepoints + mean_symmetry +  mean_fractaldim + 
               SE_radius + SE_texture +  SE_perimeter + SE_area + 
               SE_smoothness + SE_compactness + SE_concavity + 
               SE_concavepoints + SE_symmetry +  SE_fractaldim + 
               worst_radius + worst_texture +  worst_perimeter + worst_area + 
               worst_smoothness + worst_compactness + worst_concavity + 
               worst_concavepoints + worst_symmetry +  worst_fractaldim,
                data=datos, method='class')

print(cfit)
plot(cfit)
text(cfit)
library(car)
library(RcmdrMisc)
rcorr.adjust(wpbc, type = c("pearson", "spearman"), use = "complete.obs")

# 1º Dividimos el conjunto de entrenamiento y el de validación: no podemos dividir 
# nuestros datos en estos dos conjuntos ya que tenemos muy pocos datos (muy pocos
# positivos) por lo que realizaremos CV o Boostrap. Luego veremos cual es mejor.












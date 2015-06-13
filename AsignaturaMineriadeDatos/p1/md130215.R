# Iris data

data(iris)
iris2 <- iris
levels(iris2$Species) <- c('s','e','v')
# otherwise v*e*rsicolor are mixed ip with *v*irginica


## Exploramos los datos inicialmente


plot(iris[,1:4], col=as.numeric(iris$Species))

## Estimaci�n de la densidad condicional (estimadores kernel)

library(lattice)
densityplot( ~ Petal.Width+Petal.Length+Sepal.Length+Sepal.Width, groups=Species, data=iris2)



## Analisis Discriminante Lineal (lda).
## Estimaci�n de la matriz de mala clasificaci�n y del error de mala clasificacion aparente (peligro de optimismo)
require(MASS)
mod.lda <- lda(Species ~ ., data = iris2)
mod.lda

## representacion de las dos grupos de acuerdo con las dos primeras variables can�nica de clasificaci�n

plot(mod.lda)


## Estimaci�n de la matriz de mala clasificaci�n y del error de mala clasificacion aparente (peligro de optimismo)

predictLDA <- predict(mod.lda)
TAB <- table(iris2$Species, predictLDA$class)
TAB
mcrlda <- 1 - sum(diag(TAB))/sum(TAB)
mcrlda

## Representamos  de forma simultanea la frontera y los datos (utilizamos las variables can�nicas de discriminacion) 
## Uso de la funci�n partimat de la libreria klaR

library(klaR)

iris3<-cbind(iris2, predictLDA$x[,1:2])

partimat(Species ~ LD2+LD1, data=iris3, method="lda")

partimat(Species ~., data=iris2, method="lda", plot.matrix=TRUE)


## Uso de la libreria ggplot y representacion simultanea de datos y contornos de densidad (estimaci�n kernel)

library(ggplot2)

quickplot(iris3$LD1,iris3$LD2, colour=iris3$Species, geom=c("point","density2d"))



## Analisis Discriminante Cuadr�tico (qda).
## Estimaci�n de la matriz de mala clasificaci�n y del error de mala clasificacion aparente (peligro de optimismo)

mod.qda <- qda(Species ~ ., data = iris2)
mod.qda

## Estimaci�n de la matriz de mala clasificaci�n y del error de mala clasificacion aparente (peligro de optimismo)

predictQDA <- predict(mod.qda)
TAB <- table(iris2$Species, predictQDA$class)
TAB
mcrqda <- 1 - sum(diag(TAB))/sum(TAB)
mcrqda

## Representamos  de forma simultanea la frontera y los datos (utilizamos las variables can�nicas de discriminacion) 
## Uso de la funci�n partimat de la libreria klaR

library(klaR)



partimat(Species ~ ., data=iris2, method="qda", plot.matrix=TRUE)


## Naive Bayes (e1071). En v. métricas asume normalidad. Admite suavizado de laplace.

mod.NB <- naiveBayes(Species ~ ., data = iris2)
mod.NB

## Estimaci�n de la matriz de mala clasificaci�n y del error de mala clasificacion aparente (peligro de optimismo)


predictNB <- predict(mod.NB, newdata=iris2)
TAB <- table(iris2$Species, predictNB)
TAB
mcrNB <- 1 - sum(diag(TAB))/sum(TAB)
mcrNB

partimat(Species ~., data=iris2, method="naiveBayes", plot.matrix=TRUE)



## EStimacion por knn  (k=3 default)

mod.sknn <- sknn(Species ~ ., data = iris2, gamma=0)
mod.sknn

## Estimaci�n de la matriz de mala clasificaci�n y del error de mala clasificacion aparente (peligro de optimismo)


predictsknn <- predict(mod.sknn, newdata=iris2)
TAB <- table(iris2$Species, predictsknn$class)
TAB
mcrsknn <- 1 - sum(diag(TAB))/sum(TAB)
mcrsknn

partimat(Species ~., data=iris2, method="sknn", plot.matrix=TRUE, gamma=1.0)




## Seleccion de la variables

## Uso de la fucion stepclass de klaR. improvement=0.001 o m�s peque�o para que entren 2 o m�s variables
## Uso del error aparente.

# Lineal

stepclass(Species~., data=iris2, method="lda", improvement=0.001)

# Cuadr�tico

stepclass(Species~., data=iris2, method="qda", improvement=0.001)



### Seleccion de parámetros por "fuerza bruta" (grid search). 
### Funcion tune () (e1071)

## tune `knn' using a convenience function; this time with the
## conventional interface and bootstrap sampling:
x <- iris2[,-5]
y <- iris2[,5]
obj2 <- tune.knn(x,y, k = 1:15, tunecontrol = tune.control(sampling = "boot",nboot=500))
summary(obj2)
plot(obj2)

tune.knn(x,y, k = 12, tunecontrol = tune.control(sampling = "boot",nboot=500))






## funcion para resumir la capacidad clasificatoria de un m�todo
## admite prob. a priori

confusion <- function(actual, predicted, names=NULL, 
                      printit=TRUE, prior=NULL){
  if(is.null(names))names <- levels(actual)
  tab <- table(actual, predicted)
  acctab <- t(apply(tab, 1, function(x)x/sum(x)))
  dimnames(acctab) <- list(Actual=names,
                           "Predicted (cv)"=names) 
  if(is.null(prior)){
    relnum <- table(actual)
    prior <- relnum/sum(relnum)
    acc <- sum(tab[row(tab)==col(tab)])/sum(tab)
  } else
  {
    acc <- sum(prior*diag(acctab))
    names(prior) <- names
  }    
  if(printit)print(round(c("Overall accuracy"=acc, 
                           "Prior frequency"=prior),4))
  if(printit){    cat("\nConfusion matrix", "\n")
                  print(round(acctab,4))
  }
  invisible(acctab)
}


## Ejemplo con los datos iris

confusion(iris2$Species, predictLDA$class)

## Primer intento de obtener una estimaci�n con menor optimismo
## uso de CV=TRUE


# Lineal

confusion(iris2$Species, lda(Species ~., data=iris2,CV=TRUE)$class)

# Cuadratico


confusion(iris2$Species, qda(Species ~., data=iris2,CV=TRUE)$class)

table(iris2$Species, qda(Species ~., data=iris2,CV=TRUE)$class)






###################################
# Tarea: Repetir el an�lisis con los datos crabs de la libreria(MASS)
####################################


#lctura y creacion de 4 grupos

library(MASS)
data(crabs)
spsex <- 2 * as.numeric(crabs$sp) + as.numeric(crabs$sex)
spsex<-factor(spsex, levels=c(3,4,5,6), labels=c("b","B","o","O"))

plot(crabs[, 4:8], col = spsex)

densityplot( ~ FL+RW+CL+CW+BD, groups=spsex, data=crabs)


# ampliamos el conjunto de datos con las 3 primeras componentes principales

crabsPC <- predict(princomp(crabs[, 4:8]))
pairs(crabsPC[, 1:3], col = spsex)
crabs2<-cbind(crabs,spsex, crabsPC[,1:3])
names(crabs2)

mod.lda <- lda(spsex ~ Comp.1+Comp.2+Comp.3, data = crabs2)
mod.lda


# continuar .....




## Datos de los indios Pima
## Existe un grupo de entrenamiento y otro de validacion
## 2 grupos, numerosas variables (discretas y continuas)

library(MASS)

data(Pima.tr)

help(Pima.tr)



plot(Pima.tr[,1:7],col=Pima.tr[,8])

densityplot( ~ Pima.tr[,1:7], groups=Pima.tr$type, data=Pima.tr)





mod.lda <- lda(type ~ ., data = Pima.tr)
predictLDA <- predict(mod.lda, newdata = Pima.te)
TAB <- table(Pima.te$type, predictLDA$class)
TAB
mcrlda <- 1 - sum(diag(TAB))/sum(TAB)
mcrlda

confusion(Pima.te$type, predictLDA$class)



mod.qda <- qda(type ~ ., data = Pima.tr)
predictQDA <- predict(mod.qda, newdata = Pima.te)
TAB <- table(Pima.te$type, predictQDA$class)
TAB
mcrqda <- 1 - sum(diag(TAB))/sum(TAB)
mcrqda


confusion(Pima.te$type, predictQDA$class)


### MODELO de REGRESION LOG�STICA



pima.glm1<-glm(type~., binomial, data=Pima.tr)
summary(pima.glm1)

predictGLM<-predict(pima.glm1, newdata= Pima.te, type = "response")

TAB <- table(Pima.te$type, predictGLM > .5)
TAB
mcrglm1 <- 1 - sum(diag(TAB))/sum(TAB)
mcrglm1


# alternativamente

TAB <- table(Pima.te$type, round(predictGLM ))
TAB
mcrglm1 <- 1 - sum(diag(TAB))/sum(TAB)
mcrglm1



confusion(Pima.te$type, (predictGLM > .5))



## visualizacion de los efectos

termplot(pima.glm1, se=TRUE, rug=TRUE)

# Alternativa library(effecst)

require(effects)
plot(allEffects(pima.glm1))

## a�adir suavizadores para detectar efecto no lineal

termplot(pima.glm1, partial.resid=TRUE, se=TRUE, rug=TRUE, smooth=panel.smooth, span.smth=1/5)

## Suele ser necesario un n�mero grande de casos para
## poder visualizar un efecto no lineal de forma clara.


stepAIC(pima.glm1)

pima.glm2<-step(pima.glm1)

summary(pima.glm2)

## Ajustamos un modelo con posibles interacciones de orden 2.

stepAIC(pima.glm1, scope=c(upper=.~.^2, lower=.~1))

pima.glm3<-stepAIC(pima.glm1, scope=c(upper=.~.^2, lower=.~1))

summary(pima.glm3)

confusion(Pima.te$type ,round(predict(pima.glm3, newdata= Pima.te, type = "response")))

# visualización de interacciones

plot(allEffects(pima.glm3))


#notar la diferencia con 

stepAIC(glm(type~.^2,binomial, data=Pima.tr))





# evaluamos la capacidad predictiva del modelo 
# utilizamos una funci�n de la libreria boot

library(boot)
# esta funci�n cuenta la proporci�n de errores de predicci�n

cost <- function(r, pi=0) mean(abs(r-pi)>0.5)

# error de predicci�n con cross-validation K-fold
cv.glm(Pima.tr, pima.glm2, cost, K=50)


# si queremos un metodo de estimacion de tipo leave-one-out
# omitimos el valor de K (se puede hacer algo lento)

cv.glm(Pima.tr, pima.glm2, cost)

# en ambos casos, los errores estimados son mayores que
# el error aparente

fitpima.glm2<-predict(pima.glm2, type="response") 
table(Pima.tr$type,fitpima.glm2 > 0.5)






## ROC curve para lda y Pima

truepos <- numeric(19)
falsepos <- numeric(19)
p1 <- (1:19)/20
for(i in 1:19){
  p <- p1[i]
  Pima.CV1p <- lda(type ~ ., data=Pima.tr, CV=TRUE, prior=c(p, 1-p))
  confmat <- confusion(Pima.tr$type, Pima.CV1p$class, printit=FALSE)
  falsepos[i] <- confmat[1,2]
  truepos[i] <- confmat[2,2]
}

plot(truepos ~ falsepos, type="l", xlab="False positive rate",
     ylab="True positive rate (Sensitivity)")


for(i in 1:19){
  p <- p1[i]
  Pima.CV1p <- qda(type ~ ., data=Pima.tr, CV=TRUE, prior=c(p, 1-p))
  confmat <- confusion(Pima.tr$type, Pima.CV1p$class, printit=FALSE)
  falsepos[i] <- confmat[1,2]
  truepos[i] <- confmat[2,2]
}


lines(truepos ~ falsepos, type="l", col=2)

#ARBoles



library(rpart)
modeltree<-rpart(  type~. ,data=Pima.tr, method="class")
summary(modeltree)

plot(modeltree)
text(modeltree)
plotcp(modeltree)

mod.TREE<-rpart(type~. ,data=Pima.tr, cp=0.025)

# o bien
##  mod.TREE<-prune(modeltree, cp=0.025)

rpart(  type~. ,data=Pima.tr, method="class", parms=list( split='information'))



predictTREE <- predict(mod.TREE, newdata = Pima.te, type="class")
TAB <- table(Pima.te$type, predictTREE)
mcrtree <- 1 - sum(diag(TAB))/sum(TAB)
mcrtree



### Random Forest

library(randomForest)

Pima.rf <- randomForest(type~., data=Pima.tr, xtest=Pima.te[,-8],
                        ytest=Pima.te$type)
Pima.rf

plot(Pima.rf)

## ejermplo de uso de la libreria RWeka para una clasificacion


library(RWeka)

mod.J48 <- J48(type ~ ., data = Pima.tr)
## print and summary
mod.J48
summary(mod.J48) # calls evaluate_Weka_classifier()
TAB <- table(Pima.te$type, predict(mod.J48, type="class", newdata=Pima.te))
mcrJ48 <- 1 - sum(diag(TAB))/sum(TAB)
mcrJ48

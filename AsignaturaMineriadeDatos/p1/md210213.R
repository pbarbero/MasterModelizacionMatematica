



Random Forest








library(randomForest)

Pima.rf <- randomForest(type~., data=Pima.tr, xtest=Pima.te[,-8], ytest=Pima.te$type)

Pima.rf

mod.TREE<-rpart(type~. ,data=Pima.tr, cp=0.025)
predictTREE <- predict(mod.TREE, newdata = Pima.te, type="class")
TAB <- table(Pima.te$type, predictTREE)
mcrtree <- 1 - sum(diag(TAB))/sum(TAB)
mcrtree


##
##[1] 0.2439759
## 

confusion(Pima.te$type, predictTREE)






###############333
#  Bagging
#################

data(iris)
iris2 <- iris
levels(iris2$Species) <- c('s','e','v')
# otherwise v*e*rsicolor are mixed ip with *v*irginica

## rpart library should be loaded
library(rpart)

iris2.bagging <- bagging(Species~., data=iris2, mfinal=10)

irsi2.bagging

errorevol(iris2.bagging,Pima.tr)->evol.train
plot(evol.train$error, type="l", main="Adaboost error Vs number of trees",  col = "blue") 



##########
# BOOSTING
##########

cntrl<-rpart.control(maxdepth=1)

pima.adaboost <- boosting(type ~ ., data=Pima.tr, mfinal=100, control=cntrl)

#Error evolution along the iterations in training set 
errorevol(pima.adaboost,Pima.tr)->evol.train
plot(evol.train$error, type="l", main="Adaboost error Vs number of trees",  col = "blue") 


#comparing error evolution in training and test set
errorevol(pima.adaboost,Pima.te)->evol.test
plot(evol.test$error, type="l", ylim=c(0,1),  main="Adaboost error Vs number of trees",  
xlab="Iterations", ylab="Error", col = "red")
lines(evol.train$error, cex = .5 ,col="blue", lty=2)
legend("topright", c("test","train"), col = c("red", "blue"), lty=1:2)



### Cambiar los parámetros de control. P.ej.

cntrl<-rpart.control(maxdepth=3)


#############
# Rattle
##############

# Instalarlo

# Lectura de un conjunto de datos

# Clasificacion

# Evaluacion del clasificador
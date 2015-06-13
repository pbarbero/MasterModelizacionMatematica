################
# DATA READING #
################

library(MASS)
data(mcycle)
names(mcycle)
summary(mcycle)
mcycle$times
mcycle$accel
plot(mcycle)

mcycle2 <- mcycle


######################
# VALIDACION CRUZADA #
######################a

# gerardnico.com/wiki/r/cross_validation

#############
# REGRESION #
#############

===============
= PARAMETRICO =
===============

# Polinomial:
# ¡¡¡Elegir primero el grado, en este ejemplo es 20!!! 
## ajuste polinomial seleccionando el grado por CV
#    eep <- rep(NA,20)

#    for ( grado in 1:20 )  {
#          eep[grado] <- cv.glm( d, glm( y ~ poly(x, grado), data = d), K = 10)$delta[2]
#    }    

#   plot(eep,type = 'b')

==================
= NO PARAMETRICO =
==================


## Ejemplo de Función en R.
## Realiza una regresión de tipo k-nn
# como argumentos de entrada, los vectores xtrain, ytrain, el entero k (por defecto 1) y xtest como
# vector donde predecir

# la función devuelve un vector con los valores predichos en xtest

#knnreg = function(xtrain,ytrain,xtest,k=1)
#  {
#    n = length(xtrain)
#    N = length(xtest)
#    f.hat = rep(NA,N)
#    for(i in 1:N)
#      {
#        d = (xtest[i]-xtrain)^2
#        ind = c(1:n)[rank(d)<(k+1)]
#        f.hat[i] = mean(y[ind])
#      }
#    f.hat
#  }


#ajuste local constante (equiv al estimador de Nadaraya-Watson)

fit = loess(y~x,deg=0,span=0.4,data=d)
ess(formula = y ~ x, data = d, span = 0.4, degree = 0)

pred = predict(fit,d)
lines(hex,pred,lty=2)
summary(fit)

#Call:
#loess(formula = y ~ x, data = d, span = 0.4, degree = 0)
#
#Number of Observations: 100
#Equivalent Number of Parameters: 3.63
#Residual Standard Error: 0.3633
#Trace of smoother matrix: 4.23
#
#Control settings:
#  normalize:  TRUE
#  span     :  0.4
#  degree   :  0
#  family   :  gaussian
#  surface  :  interpolate    cell = 0.2


# SPLINES
library(splines)
library(ElemStatLearn)

## El paquete de funciones ElemStatLearn contiene un numero importante de funciones y la mayoría
## de los conjuntos de datos del libro de Hastie y Tibshirani


fit.spline<-smooth.spline(d$x,d$y)
lines(fit.spline, col=3)

library(mgcv) # biblioteca de GAM's

fit.gam<-gam(y~s(x), data=d)
plot(fit.gam)


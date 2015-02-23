
# Modelo simulado

###############

set.seed(353)
N = 100
x = sort(runif(N))
y = rnorm(N,sin(4*x),1/3)
plot(x,y,main="Regresión y Suavizadores")
x2 = seq(0,1,length=1000)
lines(x2,sin(4*x2))
d = data.frame(x=x,y=y)

#######################
# Modelado paramétrico
#######################


## función que construye la matriz del diseño asociado a polinomios trigonometricos

#
# Periodic model using sines and cosines
#


    cs <- function(t,harmonics=1) {
       # generate matrix with sin,cos pairs with period 1
       ret <- numeric(0)
       for ( i in harmonics) {
           ret <- cbind( ret, cos(2*pi*i*t), sin(2*pi*i*t))
       }
       if (missing(harmonics)) cnames <- c('c','s')
       else {
            cnames <- paste( c("c","s"), rep(harmonics, each = 2),sep=".")
       }
       colnames(ret) <- cnames
       ret
    }


# Ajustes con diferente grado de complejidad

    fit.cs= lm(y~cs(x),d)

    summary(fit.cs)   # can we drop a term??
    anova(fit.cs)     # here anova respects 'marginality'
    plot( d$x, d$y)
    dtest=data.frame(x=x2)
    lines( x2 , predict(fit.cs,dtest),col='red')

    lines( x2, predict( lm( y ~ cs(x,1:4),d),dtest),col='green')


    lines( x2, predict( lm( y ~ cs(x,1:8),d),dtest), col = 'black',lwd=1.5)


    lines( x2, predict( lm( y ~ cs(x,1:18),d),dtest), col = 'blue',lwd=1.5)

    lines( x2, predict( fit.cs <- lm( y ~ cs(x,1:50),d),dtest), col = 'purple',lwd=1.5)
    summary(fit.cs)

# seleccion del mejor modelo

## Using adjusted R2

    require(stats4)


    measures <- matrix(NA,nrow=50,ncol=8,
             ,dimnames=list(NULL,c("H","df","R.sq","Adj.R.sq","logLik",
             "-2logLik","AIC","BIC")))



    for ( maxh in 1:50) {
       fs <- summary( fit <- lm( y ~ cs(x,1:maxh), d))
       measures[maxh,] <- c(maxh, fit$df[1],100*fs$r.squared , 100*fs$adj.r.squared,logLik(fit),
                          -2*logLik(fit),  AIC(fit), BIC(fit))
    }

    round(measures,2)

## seleccion del número de harmonicos por CV

library(boot)


   eep <- rep(NA,20)
    for ( harmonic in 1:20 )  {
          eep[harmonic] <- cv.glm( d, glm( y ~ cs(x, 1:harmonic), data = d), K = 10)$delta[2]
    }
    
    plot(eep,type = 'b')


## minimum around 2   ????

# fit model using all data

    fit.cv <- lm( y ~ cs(x,1:2), d)
    plot( d$x, d$y)
    lines( d$x, predict(fit.cv))


## ajuste polinomial seleccionando el grado por CV


   eep <- rep(NA,20)
    for ( grado in 1:20 )  {
          eep[grado] <- cv.glm( d, glm( y ~ poly(x, grado), data = d), K = 10)$delta[2]
    }
    
    plot(eep,type = 'b')

# minimo en 7??

    fit.cv <- lm( y ~ poly(x,7), d)
    plot( d$x, d$y)
    lines( d$x, predict(fit.cv))

## Otras bases de funciones: ns (splines naturales)

  library(splines)
  eep <- rep(NA,20)
    for ( grado in 1:20 )  {
          eep[grado] <- cv.glm( d, glm( y ~ ns(x, grado), data = d), K = 10)$delta[2]
    }
    
    plot(eep,type = 'b')

   fit.cv <- lm( y ~ ns(x,3), d)
    plot( d$x, d$y)
    lines( d$x, predict(fit.cv))

#########################
# Modelado No paramétrico
#########################

## Ejemplo de Función en R.
## Realiza una regresión de tipo k-nn 
# como argumentos de entrada, los vectores xtrain, ytrain, el entero k (por defecto 1) y xtest como
# vector donde predecir

# la función devuelve un vector con los valores predichos en xtest

knnreg = function(xtrain,ytrain,xtest,k=1)
  {
    n = length(xtrain)
    N = length(xtest)
    f.hat = rep(NA,N)
    for(i in 1:N)
      {
        d = (xtest[i]-xtrain)^2
        ind = c(1:n)[rank(d)<(k+1)]
        f.hat[i] = mean(y[ind])
      }
    f.hat
  }
    
###############

# Modelos No-Paramétricos



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

# Ajuste por regresión mediante K-vecinos proximos, con k elegido en base al estimador de N-W

k=as.integer(100/4.23)
pred.knn = knnreg(x,y,x,k=k)
lines(x,pred.knn,lty=2)

# Ajuste localmente lineal (d=1)

fit.lp = loess(y~x,deg=1,span=0.5,data=d)
summary(fit.lp)
#Call:
#loess(formula = y ~ x, data = d, span = 0.5, degree = 1)
#
#Number of Observations: 100 
#Equivalent Number of Parameters: 3.67 
#Residual Standard Error: 0.3404 
#Trace of smoother matrix: 4.28 
#
#Control settings:
#  normalize:  TRUE 
#  span     :  0.5 
#  degree   :  1 
#  family   :  gaussian
#  surface  :  interpolate    cell = 0.2
pred.lp = predict(fit.lp,d)
lines(x,pred.lp,lty=3)
legend(0.2,-0.5,legend=c("True","N-W","knn18","loc-reg"),lty=1:4)

library(splines)
library(ElemStatLearn)

## El paquete de funciones ElemStatLearn contiene un numero importante de funciones y la mayoría
## de los conjuntos de datos del libro de Hastie y Tibshirani


fit.spline<-smooth.spline(d$x,d$y)
lines(fit.spline, col=3)

library(mgcv) # biblioteca de GAM's

fit.gam<-gam(y~s(x), data=d)
plot(fit.gam)


# Otra librería de uso frecuente en suavizados
library(KernSmooth)
help(locpoly)




#############################################################################3
# Propuesta: Analizar los datos motorcycle.  
# Utiliza loess, knnreg, smoothin.spline, gam, etc... 
#############################################################################
library(MASS)
data(mcycle)
plot(mcycle)
names(mcycle)

# ejemplo
#library(splines)
#lines(mcycle$times,lm(accel~ns(times,3), data=mcycle)$fit)


###################################################
# Selección automática de variables en regresion
######################################################

help(step)

library(ElemStatLearn) # En esta biblioteca se encuentra el conjunto de datos  prostate

help(prostate)

pairs(prostate[1:8])

# prostate = read.table("../../data/prostate.data")
prostate.lm = lm(lpsa~.-train,data=prostate)
step(prostate.lm)
#Start:  AIC=-58.32
#lpsa ~ lcavol + lweight + age + lbph + svi + lcp + gleason + 
#    pgg45
#
#          Df Sum of Sq     RSS     AIC
#- gleason  1     0.041  44.204 -60.231
#- pgg45    1     0.526  44.689 -59.174
#- lcp      1     0.674  44.837 -58.852
#<none>                  44.163 -58.322
#- age      1     1.550  45.713 -56.975
#- lbph     1     1.684  45.847 -56.693
#- lweight  1     3.586  47.749 -52.749
#- svi      1     4.935  49.099 -50.045
#- lcavol   1    22.372  66.535 -20.567
#
#Step:  AIC=-60.23
#lpsa ~ lcavol + lweight + age + lbph + svi + lcp + pgg45
#
#          Df Sum of Sq     RSS     AIC
#- lcp      1     0.662  44.867 -60.788
#<none>                  44.204 -60.231
#- pgg45    1     1.192  45.396 -59.650
#- age      1     1.517  45.721 -58.959
#- lbph     1     1.705  45.910 -58.559
#- lweight  1     3.546  47.751 -54.746
#- svi      1     4.898  49.103 -52.037
#- lcavol   1    23.504  67.708 -20.872
#
#Step:  AIC=-60.79
#lpsa ~ lcavol + lweight + age + lbph + svi + pgg45
#
#          Df Sum of Sq     RSS     AIC
#- pgg45    1     0.659  45.526 -61.374
#<none>                  44.867 -60.788
#- age      1     1.265  46.132 -60.092
#- lbph     1     1.647  46.513 -59.293
#- lweight  1     3.565  48.431 -55.373
#- svi      1     4.250  49.117 -54.009
#- lcavol   1    25.419  70.286 -19.248
#
#Step:  AIC=-61.37
#lpsa ~ lcavol + lweight + age + lbph + svi
#
#          Df Sum of Sq     RSS     AIC
#<none>                  45.526 -61.374
#- age      1     0.959  46.485 -61.352
#- lbph     1     1.857  47.382 -59.497
#- lweight  1     3.225  48.751 -56.735
#- svi      1     5.952  51.477 -51.456
#- lcavol   1    28.767  74.292 -15.870
#
#Call:
#lm(formula = lpsa ~ lcavol + lweight + age + lbph + svi, data = prostate)
#
#Coefficients:
#(Intercept)       lcavol      lweight          age         lbph          svi  
#0.95102      0.56561      0.42369     -0.01489      0.11184      0.72096

##
# Seleccion utilizando BIC en lugar de AIC
######

log(sum(prostate$train))


step(prostate.lm, k=log(sum(prostate$train)))

#Start:  AIC=-38.48
#lpsa ~ (lcavol + lweight + age + lbph + svi + lcp + gleason + 
#    pgg45 + train) - train
#
#          Df Sum of Sq     RSS     AIC
#- gleason  1     0.041  44.204 -42.594
#- pgg45    1     0.526  44.689 -41.536
#- lcp      1     0.674  44.837 -41.215
#- age      1     1.550  45.713 -39.337
#- lbph     1     1.684  45.847 -39.055
#<none>                  44.163 -38.479
#- lweight  1     3.586  47.749 -35.111
#- svi      1     4.935  49.099 -32.408
#- lcavol   1    22.372  66.535  -2.929
#
#Step:  AIC=-42.59
#lpsa ~ lcavol + lweight + age + lbph + svi + lcp + pgg45
#
#         Df Sum of Sq     RSS     AIC
#- lcp      1     0.662  44.867 -45.356
#- pgg45    1     1.192  45.396 -44.217
#- age      1     1.517  45.721 -43.526
#- lbph     1     1.705  45.910 -43.127
#<none>                  44.204 -42.594
#- lweight  1     3.546  47.751 -39.313
#- svi      1     4.898  49.103 -36.604
#- lcavol   1    23.504  67.708  -5.439
#
#Step:  AIC=-45.36
#lpsa ~ lcavol + lweight + age + lbph + svi + pgg45
#
#          Df Sum of Sq     RSS     AIC
#- pgg45    1     0.659  45.526 -48.146
#- age      1     1.265  46.132 -46.864
#- lbph     1     1.647  46.513 -46.064
#<none>                  44.867 -45.356
#- lweight  1     3.565  48.431 -42.145
#- svi      1     4.250  49.117 -40.781
#- lcavol   1    25.419  70.286  -6.020
#
#Step:  AIC=-48.15
#lpsa ~ lcavol + lweight + age + lbph + svi
#
#          Df Sum of Sq     RSS     AIC
#- age      1     0.959  46.485 -50.328
#- lbph     1     1.857  47.382 -48.473
#<none>                  45.526 -48.146
#- lweight  1     3.225  48.751 -45.712
#- svi      1     5.952  51.477 -40.433
#- lcavol   1    28.767  74.292  -4.847
#
#Step:  AIC=-50.33
#lpsa ~ lcavol + lweight + lbph + svi
#
#          Df Sum of Sq     RSS     AIC
#- lbph     1     1.300  47.785 -51.857
#<none>                  46.485 -50.328
#- lweight  1     2.801  49.286 -48.857
#- svi      1     5.806  52.291 -43.116
#- lcavol   1    27.830  74.315  -9.022
#
#Step:  AIC=-51.86
#lpsa ~ lcavol + lweight + svi
#
#          Df Sum of Sq     RSS     AIC
#<none>                  47.785 -51.857
#- svi      1     5.181  52.966 -46.076
#- lweight  1     5.892  53.677 -44.783
#- lcavol   1    28.045  75.830 -11.270
#
#Call:
#lm(formula = lpsa ~ lcavol + lweight + svi, data = prostate)
#
#Coefficients:
#(Intercept)       lcavol      lweight  
#    -0.2681       0.5516       0.5085  
#        svi  
#     0.6662  

##
# Acudir a la ayuda de prostate en ElemStatLearn

help(prostate)

# Realizamos una regresión ridge (o contraida)

pairs( prostate[,1:9], col="violet" )
train <- subset( prostate, train==TRUE )[,1:9]
test  <- subset( prostate, train=FALSE )[,1:9]


prostate.ridge <- simple.ridge( train[,1:8], train[,9], df=1:8 )
#
# coefficient traces:
#
matplot( prostate.ridge$df, t(prostate.ridge$beta), type="b", 
        col="blue", pch=17, ylab="coefficients" )

# Valor seleccionado para el parametro lambda
require(MASS)
lm.ridge(lpsa~., data=train)$kHKB


###
## Ajustamos el modelo mediante un GAM
##

prostate.gam<-gam(lpsa~s(lcavol)+s(lweight)+s(age)+s(lbph)+svi+s(lcp)+gleason+s(pgg45), data=train)

summary(prostate.gam)

plot(prostate.gam)


step(prostate.gam)

## Utilizamos la libreria leaps que 
## incluye una funcion regsubset de seleccion automatica
## de los nbest mejores modelos para cada uno de los tamaños
## de 1 a nvmax variables con el criterio BIC 


library(leaps)

# Si no se ha cargado anteriormente:
library(ElemStatLearn)



help(prostate)
names(prostate)
attach(prostate)


lm.prostate<-lm(lpsa~.-train,data=prostate)
step(lm.prostate)   # seleccion por AIC

step(lm.prostate, k=log(length(prostate$lpsa)))


lm.regsubset2 <- regsubsets(lpsa ~ .-train, data = prostate, nbest = 2, nvmax=8)
plot(lm.regsubset2)



## Ajustamos un modelo para predecir el lcavol

lm.regsubset2 <- regsubsets(lcavol ~ .-train, data = prostate, nbest = 2, nvmax=8)
plot(lm.regsubset2)


plot(lcavol~lpsa, data=prostate)

## permitimos interacción entre las variables

lm.regsubset3 <- regsubsets(lcavol ~ (.-train)^2, data = prostate, nbest = 2, nvmax=8)
plot(lm.regsubset3)

step(lm(lcavol~(.-train)^2,data=prostate), k=log(length(prostate$lpsa)))

lm.calvol.bic<-step(lm(lcavol~(.-train)^2,data=prostate), k=log(length(prostate$lpsa)))


plot(lcavol~lpsa, data=prostate)

## podemos analizar visualmente el efecto de la interacción de dos variables

coplot(lcavol ~ lpsa | lcp, data = prostate)
given.lcp <- co.intervals(prostate$lcp, number=2, overlap=.1)
coplot(lcavol ~ lpsa | lcp, data = prostate, given.v=given.lcp, rows=1)

## El paquete effects tb. facilita la representacion de efectos en 
## modelos lineales generales de forma cómoda
## Se carga automaticamente dentro de Rcmdr

library(effects)
plot(allEffects(lm.calvol.bic))

## Probad con 
### plot(allEffects(gam(lcavol~ ns(lpsa,2)*lcp, data=prostate)))



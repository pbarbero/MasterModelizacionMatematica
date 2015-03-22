# Iris data

data(iris)
iris2 <- iris
levels(iris2$Species) <- c('s','e','v')
# otherwise v*e*rsicolor are mixed ip with *v*irginica
library(klaR)
partimat(Species ~ ., data = iris2, method = "lda", mar=c(4,4,2,1))

mod.lda <- lda(Species ~ ., data = iris2)
summary(mod.lda)
predictLDA <- predict(mod.lda)
TAB <- table(iris2$Species, predictLDA$class)
TAB
mcrlda <- 1 - sum(diag(TAB))/sum(TAB)
mcrlda

      

#lctura y creacion de 4 grupos

library(MASS)
data(crabs)
spsex <- 2 * as.numeric(crabs$sp) + as.numeric(crabs$sex)
spsex<-factor(spsex, levels=c(3,4,5,6), labels=c("b","B","o","O"))
plot(crabs[, 4:8], col = spsex)

# ampliamos el conjunto de datos con las 3 primeras componentes principales

crabsPC <- predict(princomp(crabs[, 4:8]))
pairs(crabsPC[, 1:3], col = spsex)
crabs2<-cbind(crabs,spsex, crabsPC[,1:3])
names(crabs2)

mod.lda <- lda(spsex ~ Comp.1+Comp.2+Comp.3, data = crabs2)
mod.lda
predictLDA <- predict(mod.lda)
TAB <- table(crabs2$spsex, predictLDA$class)
TAB
mcrlda <- 1 - sum(diag(TAB))/sum(TAB)
mcrlda

mod.lda.cv <- lda(spsex ~ Comp.1+Comp.2+Comp.3, data = crabs2,CV=TRUE)
TAB <- table(crabs2$spsex, mod.lda.cv$class)
TAB
mcrlda <- 1 - sum(diag(TAB))/sum(TAB)
mcrlda
library(klaR)
partimat(spsex ~ Comp.1+Comp.2+Comp.3, plot.matrix=TRUE,data = crabs2, method = "lda")


library(rpart)
modeltree<-rpart(spsex~ Comp.1+Comp.2+Comp.3 ,data=crabs2)
summary(modeltree)

plot(modeltree)
text(modeltree)
plotcp(modeltree)
mod.TREE<-rpart(spsex~Comp.1+Comp.2+Comp.3 ,data=crabs2, cp=0.05)

predictTREE <- predict(mod.TREE, type="class")
TAB <- table(crabs2$spsex, predictTREE)
TAB
mcrtree <- 1 - sum(diag(TAB))/sum(TAB)
mcrtree
partimat(spsex ~ Comp.1+Comp.2+Comp.3, plot.matrix=TRUE,data = crabs2, method = "rpart",cp=0.05)

library(RWeka)

mod.J48 <- J48(spsex ~ Comp.1+Comp.2+Comp.3, data = crabs2)
## print and summary
mod.J48
summary(mod.J48) # calls evaluate_Weka_classifier()
TAB <- table(crabs2$spsex, predict(mod.J48, type="class"))
mcrJ48 <- 1 - sum(diag(TAB))/sum(TAB)
mcrJ48


if(require("party", quietly = TRUE)) plot(mod.J48)

mod.J48 <- J48(spsex ~ Comp.1+Comp.2+Comp.3, data = crabs2, control = Weka_control(R = TRUE))
mod.J48
summary(mod.J48)

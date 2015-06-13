library(TH.data)
library(ggplot2)
library(gridExtra)
data("wpbc")
View(wpbc)



#### histograma de las variables categoricas
p1 <- qplot(as.factor(status),data=wpbc,geom="bar",
            fill=I("black"), alpha=I(0.5),color=I("black"),
            main="Estado",xlab="",ylab="")
p2 <- qplot(as.factor(pnodes),data=wpbc,geom="bar",
            fill=I("blue"), alpha=I(0.5),color=I("black"),
            main="pnode",xlab="",
            ylab="")

grid.arrange(p1,p2,nrow=1)

# Script para la sesion del 28/2/14 de Introd. a la Mineria de Datos
# He tenido que terminar ejecutandolo en una versión 2.15 de R por
# problemas de la librería Rcpp en la versión 3.0.2 y similares.

## Necesitamo varios paquetes nuevos (tm: mineria de textos), 
## ('wordcloud': nubes de palabras)

## necesito un fichero (.txt) que contiene el texto a analizar. 
## He tomado algo aburrido como el debate entre Rubalcaba y Rajoy ()

## Instalar aquellas librerías que no estén todavía cargadas

require(tm)
require(wordcloud)
require(RColorBrewer)
require(SnowballC)

library(tm) # Framework for text mining.
library(SnowballC) # Provides wordStem() for stemming.
library(qdap) # Quantitative discourse analysis of transcripts.
library(qdapDictionaries)
library(dplyr) # Data preparation and pipes %>%.
library(RColorBrewer) # Generate palette of colours for plots.
library(ggplot2) # Plot word frequencies.
library(scales) # Include commas in numbers.
library(Rgraphviz) 


## Leer el fichero de texto, pero con cuidado de no leer palabras que no
## son semánticamente relevantes pero si muy frecuentes: de , por , la los,....
## El texto leído forma el Corpus


# leemos el archive con el texto

# creamos un corpus  



mensaje.corpus<-  Corpus(DirSource("txtmin/mensaje",encoding ="latin1"))

## En la versión 0.6 se debería  usar la siguiente sintaxis que asegura la herencia de atributos entre clases.


mensaje.corpus <- tm_map( mensaje.corpus,contet_transformer(removePunctuation)) #

getTransformations() # visualizamos las operaciones predefinidas

## es fácil definir nuevas operaciones adapatadas al usuario

# eliminamos puntuaciones y palabras específicas

mensaje.corpus <- tm_map( mensaje.corpus,contet_transformer(removePunctuation)) #

mensaje.corpus <- tm_map( mensaje.corpus, contet_transformer(tolower))  # todo a minusculas

mensaje.corpus <- tm_map( mensaje.corpus,function(x)removeWords(x,stopwords("spanish")))

stopwords("spanish")

mensaje.corpus <- tm_map( mensaje.corpus, contet_transformer(stemDocument))



# En la version 0.6 de tm conviene asegurar la conversion al tipo de dato adecuado tras las transformaciones
# realizadas por tm_map


mensaje.corpus<-Corpus(VectorSource(mensaje.corpus))


## En la versión 0.5-10 no era necesario el uso de contet_transformer

# eliminamos puntuaciones y palabras específicas

mensaje.corpus <- tm_map( mensaje.corpus, removePunctuation) #

mensaje.corpus <- tm_map( mensaje.corpus, tolower)  # todo a minusculas

mensaje.corpus <- tm_map( mensaje.corpus,function(x)removeWords(x,stopwords("spanish")))

stopwords("spanish")




mensaje.corpus <- tm_map(mensaje.corpus, stripWhitespace)


# Realizar o no stemming puede depender del problema concreto que se aborde.
# En este ejemplo no lo veo una mejora importante
#
# mensaje.corpus <- tm_map( mensaje.corpus, stemDocument, language = "spanish")
#



# inspeccion del corpus

inspect(mensaje.corpus)


# Construimos matrices de términos y documentos que podemos utilizar alternativamente al Corpus
# Calculamos la frecuencia con que aparecen las palabras
## usamos la matriz de documentos y terminos

tdm <- TermDocumentMatrix(mensaje.corpus)

dtm <- DocumentTermMatrix(mensaje.corpus)  ## traspuesta de la anterior

## Estadística de frecuencias básica

class(dtm)
dim(dtm)

freq<-colSums(as.matrix(dtm))


freq[tail(order(freq))]

findFreqTerms(dtm, lowfreq=5)

findAssocs(dtm, "ciudadanos", corlimit=0.000)

freq<-sort(colSums(as.matrix(dtm)), decreasing=TRUE)

freq

wf<- data.frame(word=names(freq), freq=freq)

head(wf)

library(ggplot2)
subset(wf, freq>500)
ggplot(aes(word, freq))
geom_bar(stat="identity")
theme(axis.text.x=element_text(angle=45,

# primera prueba

wordcloud(mensaje.corpus, scale=c(3,0.5), max.words=200, random.order=FALSE, rot.per=0.33, use.r.layout=FALSE, colors=brewer.pal(8, "Dark2"))


#Segunda forma

# Creamos  una matriz con la frecuencia de palabras

m <- as.matrix(tdm)

v <- sort(rowSums(m),decreasing=TRUE)

d <- data.frame(word = names(v),freq=v)

wordcloud(d$word,d$freq, min.freq=4)


nchar(as.character(d$word))

plot(density(nchar(as.character(d$word))))




# creamos un corpus

debate.corpus <- Corpus(DirSource("txtmin/"))

 inspect(debate.corpus)

            
# eliminamos puntuaciones y palabras específicas (esta sintaxis es adecuada para versiones anteriores)


debate.corpus <- tm_map( debate.corpus,removePunctuation) #

debate.corpus <- tm_map(debate.corpus, stripWhitespace)

debate.corpus <- tm_map( debate.corpus, tolower)  # todo a minusculas

debate.corpus <- tm_map( debate.corpus,function(x)removeWords(x,stopwords("spanish")))

debate.corpus <- tm_map(debate.corpus, stemDocument)


## En la versión 0.6 se debería  usar la siguiente sintaxis que asegura la herencia de atributos entre clases.


debate.corpus <- tm_map( debate.corpus,contet_transformer(removePunctuation)) #

getTransformations() # visualizamos las operaciones predefinidas

## es fácil definir nuevas operaciones adapatadas al usuario

# eliminamos puntuaciones y palabras específicas

mensaje.corpus <- tm_map( mensaje.corpus,contet_transformer(removePunctuation)) #

mensaje.corpus <- tm_map( mensaje.corpus, contet_transformer(tolower))  # todo a minusculas

mensaje.corpus <- tm_map( mensaje.corpus,function(x)removeWords(x,stopwords("spanish")))

stopwords("spanish")

mensaje.corpus <- tm_map( mensaje.corpus, stemDocument)

wordStem("Ni todos los españoles viven en España, ni todos los que viven en España son españoles", language="spanish")


# En la version 0.6 de tm conviene asegurar la conversion al tipo de dato adecuado tras las transformaciones
# realizadas por tm_map


debate.corpus<-Corpus(VectorSource(debate.corpus))








# hay palabras muy frecuentes que deberíamos eliminar

debate.corpus <- tm_map(debate.corpus, removeWords, c("usted", "ustedes", "señor","rajoy", "rubalcaba", "campo", "pérez"))





#Segunda forma

# Construimos matrices de términos que podemos utilizar alternativamente al Corpus
# Calculamos la frecuencia con que aparecen las palabras
# una matriz es la traspuesta de la otra

tdm <- TermDocumentMatrix(debate.corpus)

tdm <- TermDocumentMatrix(debate.corpus)

dtm <- DocumentTermMatrix(debate.corpus)




# Creamos  una matriz con la frecuencia de palabras

m <- as.matrix(tdm)

v <- sort(rowSums(m),decreasing=TRUE)

d <- data.frame(word = names(v),freq=v)

wordcloud(d$word,d$freq, min.freq=1)





# Tarea Abierta

# creamos con el texto del debate tres textos uno para cada interveniente 
# (comienzan en mayusculas)
#

x<-debate.corpus[[1]]

inspect(x)


# Idea: aprovechar la libreria qdap para "trocear" el debate.
# Ordenes para hacer algo parecido en un debate americano

library(qdap)
dat <- unlist(strsplit(x, "\\n"))

locs <- grep("STATEMENT OF ", dat)
nms <- sapply(strsplit(dat[locs], "STATEMENT OF |,"), "[", 2)
dat[locs] <- "SPLIT_HERE"
corp <- with(data.frame(person=nms, dialogue = 
                          Trim(unlist(strsplit(paste(dat[-1], collapse=" "), "SPLIT_HERE")))),
             df2tm_corpus(dialogue, person))

tm::inspect(corp)


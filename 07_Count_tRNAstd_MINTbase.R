library(readr)
library(dplyr)
library(Biobase) # Para construir y manipular el expressionSet
library(DESeq2)
library(sqldf)
library(ggplot2)
library(gridExtra)


#Laoding and ordering case data
SampleGr<-read_csv2("./Datos/Sample_groups.csv")
SampleGr <- SampleGr[order(SampleGr$sampleID), ]

Dataset1<-read.csv("./Datos/Dataset_Red.csv", header = TRUE)

# Rearrange and rename the matrix columns
dim(Dataset1)
str(Dataset1)
colnames(Dataset1)[5:44]<- SampleGr$sampleID
Seqs10<-Dataset1[2850:2855,c(1,4:8)]

Secuencias<-unique(Dataset1$Sequence)
length(Secuencias)

########################
# Annotations MINTbase #
########################

library(XML)
library(RCurl)
library(rlist)


## For RNA molecule:
## Read the URL. 
theurl <- getURL(paste("https://cm.jefferson.edu/MINTbase/InputController?v=v&g=GRCh37&e=1&search=submit&t=All&am=All&an=All&da=&tn=&fs=", seq, "&fn=&c=All&gs=&ge=", sep=""),.opts = list(ssl.verifypeer = FALSE) )
##  format table:
tables <- readHTMLTable(theurl)
## Keep table information
tables <- list.clean(tables, fun = is.null, recursive = FALSE)
str(tRNA_MINTb)

####### Tipo de tRNA #########

tRNA_MINTb<-data.frame(Type=character(),
            Sequence=character(), Length=numeric(), Only_tRNA_gene=character(),
            tRNA=character(), stringsAsFactors=FALSE)
for(i in 1:length(Secuencias)){
  seq<-Secuencias[i]
  tipo=character()
  molecula=character()
  solotRNA=character()
  Long=numeric()
  theurl2 <- getURL(paste("https://cm.jefferson.edu/MINTbase/InputController?v=f&g=GRCh37&e=1&search=submit&t=All&am=All&an=All&da=&tn=&fs=",seq,"&fn=&c=All&gs=&ge=", sep=""),.opts = list(ssl.verifypeer = FALSE) )
  tRFmintb <- readHTMLTable(theurl2)
  tRFmintb<-list.clean(tRFmintb, fun=is.null, recursive=FALSE)
  if(length(tRFmintb)==0){
    solotRNA= molecula= tipo="-"
    Long=NA
  } else if(!("Fragment Length" %in% names(tRFmintb$main_table))){
    tipo="tRF-1"
    molecula= solotRNA="-"
    Long=NA
  } else if(tRFmintb$main_table$`5'-half`!=""){
    tipo="5'-half"
    molecula=tRFmintb$main_table$`5'-half`
    Long<-as.numeric(tRFmintb$main_table$`Fragment Length`)
    solotRNA <- tRFmintb$main_table$`Exclusively within tRNA genes?`
  } else if(tRFmintb$main_table$`5'-tRF`!=""){
    tipo="5'-tRF"
    molecula=tRFmintb$main_table$`5'-tRF`
    Long<-as.numeric(tRFmintb$main_table$`Fragment Length`)
    solotRNA <- tRFmintb$main_table$`Exclusively within tRNA genes?`
  } else if (tRFmintb$main_table$`i-tRF`!=""){
    tipo="i-tRF"
    molecula=tRFmintb$main_table$`i-tRF`
    Long<-as.numeric(tRFmintb$main_table$`Fragment Length`)
    solotRNA <- tRFmintb$main_table$`Exclusively within tRNA genes?`
  } else if (tRFmintb$main_table$`3'-tRF`!=""){
    tipo="3'-tRF"
    molecula=tRFmintb$main_table$`3'-tRF`
    Long<-as.numeric(tRFmintb$main_table$`Fragment Length`)
    solotRNA <- tRFmintb$main_table$`Exclusively within tRNA genes?`
  } else if (tRFmintb$main_table$`3'-half`!=""){
    tipo="3'-half"
    molecula= tRFmintb$main_table$`3'-half`
    Long<-as.numeric(tRFmintb$main_table$`Fragment Length`)
    solotRNA <- tRFmintb$main_table$`Exclusively within tRNA genes?`
  } else {
    solotRNA= molecula = tipo="-"
    Long=NA
  }
  fila <- data.frame(Type=tipo, Sequence=seq, Length=Long, Only_tRNA_gene=solotRNA, tRNA=molecula)
  tRNA_MINTb <- rbind(tRNA_MINTb, fila)
}

table(tRNA_MINTb$Type)
write.csv(tRNA_MINTb, file = "./Datos/tRNAunique_fragments_MINTb.csv", row.names = FALSE)

tRNA_MINTb<- read.csv("./Datos/tRNAunique_fragments_MINTb.csv")

head(tRNA_MINTb)
tail(tRNA_MINTb)
sum(is.na(tRNA_MINTb$Length))
sum(tRNA_MINTb$tRNA=="-")
str(tRNA_MINTb)
NuclSeq<-nchar(tRNA_MINTb$Sequence)
tRNA_MINTb<-cbind(tRNA_MINTb,NuclSeq)
str(tRNA_MINTb)

table(tRNA_MINTb$Type)
NoData<-table(tRNA_MINTb$Type)[1]
NotRF<-table(tRNA_MINTb$Type)[7]
SumaNotRF<-NoData+NotRF

largo<-sqldf("SELECT * FROM tRNA_MINTb WHERE Length>35 AND Type LIKE '5''-tRF'")


## Grafs for proportions ##

# Long sequences
LongGr<-ggplot(tRNA_MINTb, aes(x = Length)) +
  geom_histogram(fill = "steelblue", color = "black", alpha = 0.75,bins=20) +
  labs(title = "Length of the sequences", x="Length (nt)", y = "Counts")+
  theme_bw()+ theme(plot.title=element_text(size=12),axis.title=element_text(size=10))

LongGr

NuclSGr<-ggplot(tRNA_MINTb, aes(x = NuclSeq)) +
  geom_histogram(fill = "salmon", color = "black", alpha = 0.75,bins=20) +
  labs(title = "Length of the sequences", x="Length (nt)", y = "Counts")+
  theme_bw()+ theme(plot.title=element_text(size=12),axis.title=element_text(size=10))

NuclSGr


Longitudes<-list(tRNAStudio=c(tRNA_MINTb$NuclSeq), MINTbase=c(tRNA_MINTb$Length))
head(Longitudes)

tStudio<-character()
Mintb<-character()
for (i in 1:length(tRNA_MINTb$NuclSeq)){
  tStudio<-c(tStudio,"tRNAStudio")
  Mintb<-c(Mintb,"MINTbase")
} 
Source<-c(tStudio,Mintb)
head(Source)
tail(Source)
length(Source)
Nucleot<-c(tRNA_MINTb$NuclSeq,tRNA_MINTb$Length)
class(Nucleot)
Longitudes<-data.frame(Source=Source,Nucleot=Nucleot)
dim(Longitudes)
str(Longitudes)

NtdsGr<-ggplot(Longitudes, aes(x = Nucleot , fill =Source)) + 
  geom_histogram(color="darkgray",alpha=0.5, position = "identity",bins=20)+
  scale_fill_manual(values = c("salmon","steelblue"))+
  labs(title = "MINTbase selecion of tRFs", x="Length (nt)", y = "Counts")+
  theme_bw()+ theme(plot.title=element_text(size=12),axis.title=element_text(size=10))+
  theme(legend.position = c(0.9,0.9), legend.text = element_text(size=10), legend.key.size = unit(0.5, 'cm'), legend.title =element_text(size=10))

NtdsGr


## Types of tRFs

TipostRNA<-as.data.frame(table(tRNA_MINTb$Type))
names(TipostRNA)<-c("Type", "Freq")
porcent<-round(TipostRNA$Freq/sum(TipostRNA$Freq)*100, digits=1)
TipostRNA <- data.frame(TipostRNA, porcent)
TipostRNA <- TipostRNA %>% mutate(porcent = paste0(porcent, ' %'))

tiposGr<-ggplot(TipostRNA,aes(x="",y=Freq, fill=Type))+
  geom_bar(stat = "identity",color="white")+
  geom_text(aes(label=porcent), position=position_stack(vjust=0.5),color="black",size=4)+
  scale_fill_manual(values=c("lightblue","salmon","cornflowerblue","darkgoldenrod2","cadetblue","darkgray","azure2"))+
  coord_polar(theta="y")+theme_void()+ labs(title="Proportion of tRF types")

onlytRNA<-as.data.frame(table(tRNA_MINTb$Only_tRNA_gene))
names(onlytRNA)<-c("Exclusive", "Freq")
porcOtRNA<-round(onlytRNA$Freq/sum(onlytRNA$Freq)*100, digits=1)
onlytRNA <- data.frame(onlytRNA, porcOtRNA)
onlytRNA <- onlytRNA %>% mutate(porcOtRNA = paste0(porcOtRNA, ' %'))

ExclGr<-ggplot(onlytRNA,aes(x="",y=Freq, fill=Exclusive))+
  geom_bar(stat = "identity",color="white")+
  geom_text(aes(label=porcOtRNA), position=position_stack(vjust=0.5),color="black",size=4)+
  scale_fill_manual(values=c("lightblue","salmon","steelblue"))+
  coord_polar(theta="y")+theme_void()+ labs(title="Exclusive form tRNA genes")

### parental tRNA species

tRNAsps<-as.data.frame(tRNA_MINTb$tRNA)
names(tRNAsps)<-"tRNA"
tRNAsps$tRNA<-as.character(tRNAsps$tRNA)
tRNAsps2 <- tRNAsps %>% filter(!grepl("-", tRNA))

SpcGrq<-ggplot(tRNAsps2, aes(x = tRNA)) +
  geom_bar(fill = "steelblue", color = "white", alpha = 0.75) +
  labs(title = "Identified tRNA species", x = NULL, y = "Frequency (number)") +
  theme(axis.text = element_text(size = 5)) +
  coord_flip() + theme_bw()


Patron <- matrix(c(1,2), ncol = 2, byrow = TRUE)
grid.arrange(tiposGr,ExclGr, layout_matrix=Patron)
SpcGrq

table(tRNA_MINTb$Length)


Largos5half<-sqldf("SELECT * FROM tRNA_MINTb WHERE Length > 20 AND Type LIKE '%5''-half%'")
Cortos5half<-sqldf("SELECT * FROM tRNA_MINTb WHERE Length < 20 AND Type LIKE '%5''-half%'")
Largos3half<-sqldf("SELECT * FROM tRNA_MINTb WHERE Length > 20 AND Type LIKE '%3''-half%'")
Cortos3half<-sqldf("SELECT * FROM tRNA_MINTb WHERE Length < 20 AND Type LIKE '%3''-half%'")

which(tRNA_MINTb$Type=="3'-tRF")
Long3tRF<-sqldf("SELECT * FROM tRNA_MINTb WHERE Type LIKE '%3''-tRF%'")
which(Long3tRF$Length==37)
Long3tRF[c(4,9,10),]


LongGr<-ggplot(tRNA_MINTb, aes(x = Length)) +
  geom_histogram(fill = "steelblue", color = "black", alpha = 0.75,bins=20) +
  labs(title = "Length of the sequences", x="Length (nt)", y = "Counts")+
  theme_bw()+ theme(plot.title=element_text(size=12),axis.title=element_text(size=10))

LongGr

Todos5prima<-sqldf("SELECT * FROM tRNA_MINTb WHERE Type LIKE '%5''-half%' OR Type LIKE '%5''-tRF%'")
todos5pGr<-ggplot(Todos5prima, aes(x = Length , fill =Type)) + 
  geom_histogram(color="darkgray",alpha = 0.6, position = "identity",bins=20)+
  scale_fill_manual(values = c("steelblue","darkgoldenrod2"))+
  labs(title = "5' fragments", x="Length (nt)", y = "Counts")+
  theme_bw()+ theme(plot.title=element_text(size=12),axis.title=element_text(size=10))+
  theme(legend.position = c(0.85,0.8), legend.text = element_text(size=7), legend.key.size = unit(0.3, 'cm'), legend.title =element_text(size=7))

todos5pGr

#Para el 3 prima
Todos3prima<-sqldf("SELECT * FROM tRNA_MINTb WHERE Type LIKE '%3''-half%' OR Type LIKE '%3''-tRF%'")
todos3pGr<-ggplot(Todos3prima, aes(x = Length , fill =Type)) + 
  geom_histogram(color="darkgray",alpha = 0.6, position = "identity",bins=20)+
  scale_fill_manual(values = c("cadetblue", "salmon"))+
  labs(title = "3' fragments", x="Length (nt)", y = "Counts")+
  theme_bw()+ theme(plot.title=element_text(size=12),axis.title=element_text(size=10))+
  theme(legend.position = c(0.85,0.8), legend.text = element_text(size=7), legend.key.size = unit(0.3, 'cm'), legend.title =element_text(size=7))

todos3pGr

#Para el intermedio
TodosI<-sqldf("SELECT * FROM tRNA_MINTb WHERE Type LIKE '%i-tRF%'")
todosIGr<-ggplot(TodosI, aes(x = Length , fill =Type))+
  geom_histogram(color="darkgray",alpha = 0.6, position = "identity",bins=20)+
  scale_fill_manual(values = "gray")+
  labs(title = "i-tRF fragments", x="Length (nt)", y = "Counts")+
  theme_bw()+ theme(plot.title=element_text(size=12),axis.title=element_text(size=10))+
  theme(legend.position = c(0.85,0.8), legend.text = element_text(size=7), legend.key.size = unit(0.3, 'cm'), legend.title =element_text(size=7))

todosIGr

grid.arrange(LongGr,todos5pGr,todos3pGr, todosIGr,ncol=2, nrow=2)



##### tRF type and aligment  #####

# Dataframe using I() to allow lists
AlignMINTb<-data.frame(Type=character(), SeqID=character(), Sequence=character(),tRNA=as.character(), tRNA_name=I(data.frame()),gtRNAdb_name=I(data.frame()),num_tRFs=I(data.frame()), stringsAsFactors=FALSE)
for(i in 1:nrow(SeqsPr)){
     seq<-SeqsPr[i,3]
     # tRNA molecule
     tipo=as.character()
     molecula=as.character()
     theurl2 <- getURL(paste("https://cm.jefferson.edu/MINTbase/InputController?v=f&g=GRCh37&e=1&search=submit&t=All&am=All&an=All&da=&tn=&fs=",seq,"&fn=&c=All&gs=&ge=", sep=""),.opts = list(ssl.verifypeer = FALSE) )
     tRFmintb <- readHTMLTable(theurl2)
     tRFmintb<-list.clean(tRFmintb, fun=is.null, recursive=FALSE)
     if(length(tRFmintb)==0){
       molecula= tipo="-"
     } else if(!("Fragment Length" %in% names(tRFmintb$main_table))){
       tipo="No_tRF"
       molecula="-"
     } else if(tRFmintb$main_table$`5'-half`!=""){
       tipo="5'-half"
       molecula=tRFmintb$main_table$`5'-half`
     } else if(tRFmintb$main_table$`5'-tRF`!=""){
       tipo="5'-tRF"
       molecula=tRFmintb$main_table$`5'-tRF`
     } else if (tRFmintb$main_table$`i-tRF`!=""){
       tipo="i-tRF"
       molecula=tRFmintb$main_table$`i-tRF`
     } else if (tRFmintb$main_table$`3'-tRF`!=""){
       tipo="3'-tRF"
       molecula=tRFmintb$main_table$`3'-tRF`
     } else if (tRFmintb$main_table$`3'-half`!=""){
       tipo="3'-half'"
       molecula= tRFmintb$main_table$`3'-half`
     } else molecula= tipo="-"
     # tRNA alignment
     theurl <- getURL(paste("https://cm.jefferson.edu/MINTbase/InputController?v=v&g=GRCh37&e=1&search=submit&t=All&am=All&an=All&da=&tn=&fs=", seq, "&fn=&c=All&gs=&ge=", sep=""),.opts = list(ssl.verifypeer = FALSE))
     ResMintb <- readHTMLTable(theurl)
     ResMintb <- list.clean(ResMintb, fun=is.null, recursive=FALSE)
     if(length(ResMintb)==0){
       nueva_fila <- data.frame(Type=tipo, SeqID=SeqsPr[i,1], Sequence=seq, tRNA=molecula, tRNA_name=NA, gtRNAdb_name=NA, num_tRFs=NA)
     } else if(!("MINTbase tRNA name" %in% names(ResMintb$main_table))){
         nueva_fila <- data.frame(Type=tipo, SeqID=SeqsPr[i,1], Sequence=seq, tRNA=molecula, tRNA_name="No_tRF_name", gtRNAdb_name="No_tRF",num_tRFs=NA)
         } else {
         nueva_fila <- data.frame(Type=tipo, SeqID=SeqsPr[i,1], Sequence=seq, tRNA=molecula, tRNA_name=I(ResMintb$main_table$`MINTbase tRNA name`), gtRNAdb_name=I(ResMintb$main_table$`gtRNAdb name`), num_tRFs=I(ResMintb$main_table$`Number of tRFs returned`))
         }
     AlignMINTb <- rbind(AlignMINTb, nueva_fila)
}

table(AlignMINTb$Type)
table(AlignMINTb$tRNA)
tail(AlignMINTb)
str(AlignMINTb)
AlignMINTb[25,]

which(AlignMINTb$tRNA_name=="No_tRF_name")
AlignMINTb[c(1:20,410:425),]

# Guardar el dataframe como archivo CSV
write.csv(AlignMINTb, file = "./Datos/AlignMINTb.csv", row.names = FALSE)


##### Loading the dataframe with the MINTbase annotations ####
library(ggVennDiagram)
library(ggplot2)
library(VennDiagram) # No he usado este

AnnotMINTb<-read_csv("./Datos/AlignMINTb.csv")
sum(AnnotMINTb$gtRNAdb_name=="no_tRF")

tStudio_MINTb<-list(tRNAstudio=Dataset_Red$tRNA, MINTbase=AnnotMINTb$gtRNAdb_name)

SoloMINTb<- setdiff(tStudio_MINTb$MINTbase, tStudio_MINTb$tRNAstudio)
SolotRNAst<- setdiff(tStudio_MINTb$tRNAstudio, tStudio_MINTb$MINTbase)


dim(tStudio_MINTb)
ggVennDiagram(tStudio_MINTb, label_alpha = 0, category.names = c("tStd", "MTb"), set_color=c("gray","darkgray"))+
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+ theme(legend.position = "none")

length(table(AnnotMINTb$gtRNAdb_name))
length(table(Dataset_Red$tRNA))

table(Dataset_Red$tRNA)
table(AnnotMINTb$gtRNAdb_name)

# Those that aremnot in MINTbase
Errores<-sqldf("SELECT SeqID, Sequence FROM AnnotMINTb WHERE gtRNAdb_name=='no_tRF'")
SolMINT<-sqldf("SELECT SeqID, Sequence FROM AnnotMINTb WHERE gtRNAdb_name=='no_tRF'")
Faltan<-subset(AnnotMINTb, gtRNAdb_name=='no_tRF', select=c(SeqID,Sequence))

dim(Faltan)
dim(Errores)
table(Faltan$SeqID)

hola<-subset(AnnotMINTb, SeqID=='Seq_106', select=c(SeqID,Sequence, gtRNAdb_name))
hola
adios<-subset(Dataset_Red, SeqID=='Seq_106', select=c(SeqID,Sequence, tRNA))
adios
sum(is.na(AnnotMINTb$gtRNAdb_name))



class(AlignMINTb)
summary(AlignMINTb)
str(AlignMINTb)
AlignMINTb$tRNA_name <- as.character(AlignMINTb$tRNA_name)
AlignMINTb$gtRNAdb_name <- as.character(AlignMINTb$gtRNAdb_name)
AlignMINTb$num_tRFs <- as.character(AlignMINTb$num_tRFs)


#Keep the results
write.csv(AlignMINTb, file = "./Datos/AlignMINTb.csv", row.names = FALSE)



###### FUNCION DE BÚSQUEDA DE MOLÉCULAS EN MINTbase ######
#Funciona con una columna de un dataframe o de una matriz.

Molecules_MINTb <- function(SeqsPr) {
   tRNA_MINTb<-data.frame(Sequence=character(), Type=character(),
        Only_tRNA_gene=character(), tRNA=character(), stringsAsFactors=FALSE)
   
   for(i in 1:nrow(SeqsPr)){
      seq<-SeqsPr[i,]
      tipo=character()
      molecula=character()
      solotRNA=character()
      theurl2 <- getURL(paste("https://cm.jefferson.edu/MINTbase/InputController?v=f&g=GRCh37&e=1&search=submit&t=All&am=All&an=All&da=&tn=&fs=",seq,"&fn=&c=All&gs=&ge=", sep=""),.opts = list(ssl.verifypeer = FALSE) )
      tRFmintb <- readHTMLTable(theurl2)
      tRFmintb<-list.clean(tRFmintb, fun=is.null, recursive=FALSE)
      if(length(tRFmintb)==0){
        solotRNA= molecula= tipo="-"
      } else if(!("Fragment Length" %in% names(tRFmintb$main_table))){
        tipo="No_tRF"
        molecula= solotRNA="-"
      } else if(tRFmintb$main_table$`5'-half`!=""){
        tipo="5'-half"
        molecula=tRFmintb$main_table$`5'-half`
        solotRNA <- tRFmintb$main_table$`Exclusively within tRNA genes?`
      } else if(tRFmintb$main_table$`5'-tRF`!=""){
        tipo="5'-tRF"
        molecula=tRFmintb$main_table$`5'-tRF`
        solotRNA <- tRFmintb$main_table$`Exclusively within tRNA genes?`
      } else if (tRFmintb$main_table$`i-tRF`!=""){
        tipo="i-tRF"
        molecula=tRFmintb$main_table$`i-tRF`
        solotRNA <- tRFmintb$main_table$`Exclusively within tRNA genes?`
      } else if (tRFmintb$main_table$`3'-tRF`!=""){
        tipo="3'-tRF"
        molecula=tRFmintb$main_table$`3'-tRF`
        solotRNA <- tRFmintb$main_table$`Exclusively within tRNA genes?`
      } else if (tRFmintb$main_table$`3'-half`!=""){
        tipo="3'-half"
        molecula= tRFmintb$main_table$`3'-half`
        solotRNA <- tRFmintb$main_table$`Exclusively within tRNA genes?`
      } else {
        solotRNA= molecula = tipo="-"
    }
    fila <- data.frame(Sequence=seq, Type=tipo, Only_tRNA_gene=solotRNA, tRNA=molecula)
    tRNA_MINTb <- rbind(tRNA_MINTb, fila)
   }
   return(tRNA_MINTb)
}   

Hola<- Molecules_MINTb(Secuencias)

Secuencias<-c(Dataset_Red[160:165,3])
class(Secuencias)





###### FUNCION DE ALINEAMIENTO DE MOLÉCULAS EN MINTbase ######
#Funciona con una columna de un dataframe o de una matriz.

# En el dataframe uso I() para permitir añadir listas
Alignment_MINTb<- function(SeqsPr) {
   AlignMINTb<-data.frame(SeqID=character(), Sequence=character(),tRNA_name=I(data.frame()),gtRNAdb_name=I(data.frame()),num_tRFs=I(data.frame()), stringsAsFactors=FALSE)
   for(i in 1:nrow(SeqsPr)){
     seq<-SeqsPr[i,3]
     # tRNA alignment
     theurl <- getURL(paste("https://cm.jefferson.edu/MINTbase/InputController?v=v&g=GRCh37&e=1&search=submit&t=All&am=All&an=All&da=&tn=&fs=", seq, "&fn=&c=All&gs=&ge=", sep=""),.opts = list(ssl.verifypeer = FALSE))
     ResMintb <- readHTMLTable(theurl)
     ResMintb <- list.clean(ResMintb, fun=is.null, recursive=FALSE)
     if(length(ResMintb)==0){
       nueva_fila <- data.frame(SeqID=SeqsPr[i,1], Sequence=seq, tRNA_name=NA, gtRNAdb_name=NA, num_tRFs=NA)
     } else if(!("MINTbase tRNA name" %in% names(ResMintb$main_table))){
       nueva_fila <- data.frame(SeqID=SeqsPr[i,1], Sequence=seq, tRNA_name="No_tRF_name", gtRNAdb_name="No_tRF",num_tRFs=NA)
     } else {
       nueva_fila <- data.frame(SeqID=SeqsPr[i,1], Sequence=seq, tRNA_name=I(ResMintb$main_table$`MINTbase tRNA name`), gtRNAdb_name=I(ResMintb$main_table$`gtRNAdb name`), num_tRFs=I(ResMintb$main_table$`Number of tRFs returned`))
     }
     AlignMINTb <- rbind(AlignMINTb, nueva_fila)
   }
   return(AlignMINTb)
}

Secuencias<-as.data.frame(Dataset_Red[160:165,])
class(Secuencias)
Prueba<- Alignment_MINTb(Secuencias)





##### Tipo de tRNA y Alineamiento #####

# En el dataframe uso I() para permitir añadir listas
AlignMINTb<-data.frame(Type=character(), SeqID=character(), Sequence=character(),tRNA=as.character(), tRNA_name=I(data.frame()),gtRNAdb_name=I(data.frame()),num_tRFs=I(data.frame()), stringsAsFactors=FALSE)
for(i in 1:nrow(SeqsPr)){
  seq<-SeqsPr[i,3]
  # tRNA molecule
  tipo=as.character()
  molecula=as.character()
  theurl2 <- getURL(paste("https://cm.jefferson.edu/MINTbase/InputController?v=f&g=GRCh37&e=1&search=submit&t=All&am=All&an=All&da=&tn=&fs=",seq,"&fn=&c=All&gs=&ge=", sep=""),.opts = list(ssl.verifypeer = FALSE) )
  tRFmintb <- readHTMLTable(theurl2)
  tRFmintb<-list.clean(tRFmintb, fun=is.null, recursive=FALSE)
  if(length(tRFmintb)==0){
    molecula= tipo="-"
  } else if(!("Fragment Length" %in% names(tRFmintb$main_table))){
    tipo="No_tRF"
    molecula="-"
  } else if(tRFmintb$main_table$`5'-half`!=""){
    tipo="5'-half"
    molecula=tRFmintb$main_table$`5'-half`
  } else if(tRFmintb$main_table$`5'-tRF`!=""){
    tipo="5'-tRF"
    molecula=tRFmintb$main_table$`5'-tRF`
  } else if (tRFmintb$main_table$`i-tRF`!=""){
    tipo="i-tRF"
    molecula=tRFmintb$main_table$`i-tRF`
  } else if (tRFmintb$main_table$`3'-tRF`!=""){
    tipo="3'-tRF"
    molecula=tRFmintb$main_table$`3'-tRF`
  } else if (tRFmintb$main_table$`3'-half`!=""){
    tipo="3'-half'"
    molecula= tRFmintb$main_table$`3'-half`
  } else molecula= tipo="-"
  # tRNA alignment
  theurl <- getURL(paste("https://cm.jefferson.edu/MINTbase/InputController?v=v&g=GRCh37&e=1&search=submit&t=All&am=All&an=All&da=&tn=&fs=", seq, "&fn=&c=All&gs=&ge=", sep=""),.opts = list(ssl.verifypeer = FALSE))
  ResMintb <- readHTMLTable(theurl)
  ResMintb <- list.clean(ResMintb, fun=is.null, recursive=FALSE)
  if(length(ResMintb)==0){
    nueva_fila <- data.frame(Type=tipo, SeqID=SeqsPr[i,1], Sequence=seq, tRNA=molecula, tRNA_name=NA, gtRNAdb_name=NA, num_tRFs=NA)
  } else if(!("MINTbase tRNA name" %in% names(ResMintb$main_table))){
    nueva_fila <- data.frame(Type=tipo, SeqID=SeqsPr[i,1], Sequence=seq, tRNA=molecula, tRNA_name="No_tRF_name", gtRNAdb_name="No_tRF",num_tRFs=NA)
  } else {
    nueva_fila <- data.frame(Type=tipo, SeqID=SeqsPr[i,1], Sequence=seq, tRNA=molecula, tRNA_name=I(ResMintb$main_table$`MINTbase tRNA name`), gtRNAdb_name=I(ResMintb$main_table$`gtRNAdb name`), num_tRFs=I(ResMintb$main_table$`Number of tRFs returned`))
  }
  AlignMINTb <- rbind(AlignMINTb, nueva_fila)
}


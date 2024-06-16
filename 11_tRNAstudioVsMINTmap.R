library(readr)
library(dplyr)
library(sqldf)
library(ggplot2)
library(gridExtra)
library(ggvenn)


## Ful datasets
tStudio_MINTmap<-read.csv("Datos/tRNAunique_fragments_MINTb.csv", header=TRUE)
MINTdata<-read.csv("./Datos/MINTmapData_Red.csv")
dim(MINTdata)

ListaSeqs<-list(tRNAstudio=tStudio_MINTmap$Sequence,MINTmap=MINTdata$Sequence)

SoloMINT<- setdiff(ListaSeqs$MINTmap, ListaSeqs$tRNAstudio)
NsoloMINT<-length(SoloMINT)
NsoloMINT*100/length(ListaSeqs$MINTmap)

SolotRNA<-setdiff(ListaSeqs$tRNAstudio, ListaSeqs$MINTmap)
NsolotRNA<-length(SolotRNA)
NsolotRNA*100/length(ListaSeqs$tRNAstudio)


#Venn's Diagrams
jpeg("./Grafics/Preliminares/VennDiagram.jpeg", units="cm", width=11, height=8, res=300)
ggvenn(ListaSeqs, show_percentage=F, stroke_size = 0.5, text_size=5, set_name_size=4)
dev.off()

# Anotations tRFs de MINTmap
AnotacMINTmap<-read.csv("./Datos/MINTmap_Anotation_seqs.csv")
AnotacMINTmap$type<-sapply(strsplit(AnotacMINTmap$tRF_type, "\t"), `[`, 1)
table(AnotacMINTmap$Type)
AnotacMINTmap$Type[7]
AnotMINTmapUnic<-unique(AnotacMINTmap[,c(3,6)])
table(AnotMINTmapUnic$type)
# Dividir la segunda columna en elementos individuales
TodostRF <- unlist(strsplit(AnotMINTmapUnic$type, ", "))
table(TodostRF)

TiposMINT<-as.data.frame(table(TodostRF))
names(TiposMINT)<-c("Tipus", "Freq")
percent<-round(TiposMINT$Freq/sum(TiposMINT$Freq)*100, digits=1)
TiposMINT <- data.frame(TiposMINT, percent)
TiposMINT <- TiposMINT %>% mutate(percent = paste0(percent, ' %'))


colortRF<-c("orchid","cornflowerblue","darkgoldenrod2","cadetblue","salmon", "darkgray","lightblue")

jpeg("./Grafics/Preliminares/TypesMINT.jpeg", units="cm", width=11.5, height=10, res=300)
ggplot(TiposMINT,aes(x="",y=Freq, fill=Tipus))+
  geom_bar(stat = "identity",color="white")+
  geom_text(aes(label=percent), position=position_stack(vjust=0.5),color="black",size=4)+
  scale_fill_manual(values=colortRF)+
  coord_polar(theta="y")+theme_void()+ labs(title="MINTmap")
dev.off()


#Anotations of tRFs in tRNAstudio

tRNA_Anot<- read.csv("./Datos/tRNAunique_fragments_MINTb.csv")
str(tRNA_Anot)
tRNA_Anot$Type <- gsub("tRF-1", "i-tRF", tRNA_Anot$Type)

TipostRNAst<-as.data.frame(table(tRNA_Anot$Type))
names(TipostRNAst)<-c("Tipus", "Freq")
porcent<-round(TipostRNAst$Freq/sum(TipostRNAst$Freq)*100, digits=1)
TipostRNAst <- data.frame(TipostRNAst, porcent)
TipostRNAst <- TipostRNAst %>% mutate(porcent = paste0(porcent, ' %'))

jpeg("./Grafics/Preliminares/TypestRNAstd.jpeg", units="cm", width=11.5, height=10, res=300)
ggplot(TipostRNAst,aes(x="",y=Freq, fill=Tipus))+
   geom_bar(stat = "identity",color="white")+
   geom_text(aes(label=porcent), position=position_stack(vjust=0.5),color="black",size=4)+
   scale_fill_manual(values=c("azure2","salmon","cornflowerblue","darkgoldenrod2","cadetblue","darkgray"))+
   coord_polar(theta="y")+theme_void()+ labs(title="tRNAstudio")
dev.off()

grid.arrange(grStd,grM, ncol=2)



##############################
#### Comparison of PLSDAs ####
##############################

# tRNAstudio data
Todo_tRNAStd<-read.csv("./Datos/all_sequence_tRNA_level_ok.csv", sep="\t")
tRNADataRed<-read.table("./Datos/tRNAstd_Dataset_uniq.txt", header=TRUE)
IDsPLSDA_std<-read.table("./Datos/FinalVIPS_PLSDA.txt")
SeqPLSDA_std <- data.frame(subset(tRNADataRed$Sequence, tRNADataRed$SeqID %in% IDsPLSDA_std$x))
colnames(SeqPLSDA_std)<-"Sequences"

# MINTmap data
MINTDataRed<-read.csv("./Datos/MINTmapData_Red.csv", header = TRUE, sep=",")
tRNAfullData<-read.csv("./Datos/MINTmap_Anotation_seqs.csv", header = TRUE)
IDsPLSDA_MINT<-read.table("./Datos/FinalVIPS_PLSDA_mint.txt")
SeqPLSDA_mint <- data.frame(subset(MINTDataRed$Sequence, MINTDataRed$SeqID %in% IDsPLSDA_MINT$x))
colnames(SeqPLSDA_mint)<-"Sequences"

ListaPLSDA<-list(tRNAstd=SeqPLSDA_std$Sequences, MINTmap=SeqPLSDA_mint$Sequences)

#Those that are equal:

jpeg("./Grafics/PLSDA/Venn_SeqsComuns.jpeg", units="cm", width = 12, height=10,   res=300)
ggvenn(ListaPLSDA, show_percentage=F, stroke_size = 0.5, text_size=5, set_name_size=4)
dev.off()

ComunesPLSDA<-intersect(ListaPLSDA$tRNAstd, ListaPLSDA$MINTmap)
write.table(ComunesPLSDA,"./Datos/SeqsComunes_PLSDA.txt", quote = FALSE)

#Parental tRNAs
tRNAs_com<-data.frame(subset(Todo_tRNAStd$tRNA, Todo_tRNAStd$Sequence %in% ComunesPLSDA))
colnames(tRNAs_com)<-"Comuns"

tRNAs_com$Comuns<-factor(tRNAs_com$Comuns)
TablatRNAs<-table(tRNAs_com$Comuns)
dim(table(tRNAs_com$Comuns))

tRNAspecies<-sort(unique(tRNAs_com$Comuns))
length(tRNAspecies)
freqs<-table(tRNAs_com$Comuns)
OrdFreqs<-freqs[order(names(freqs))]
tRNA_freq.df<-data.frame(tRNA=rep(NA,length(tRNAspecies)), freqs=rep(NA,length(tRNAspecies)))
tRNA_freq.df$tRNA<-tRNAspecies
tRNA_freq.df$freqs<-OrdFreqs
dim(tRNA_freq.df)
str(tRNA_freq.df)

#Lolliplot  
jpeg("./Grafics/PLSDA/Lolliplot_tRNAComuns.jpeg", units="cm", width = 12, height=11,   res=300)
ggplot(data= tRNA_freq.df, aes(x=tRNA, y=freqs))+
  geom_segment(aes(x=reorder(tRNA,freqs), xend=tRNA, y=0, yend=freqs))+
  xlab("") + ylab("Nombre") + labs(title="tRNAs de les seqüències comuns")+
  coord_flip() +
  geom_point(size=3, color=("#a0522d"), fill=("#a1887f"), shape=21, stroke=1)+
  theme_bw()
dev.off()


#Those that anotate in MINTbase

tRNAfullData$Type<-sapply(strsplit(tRNAfullData$tRF_type, "\t"), `[`, 1)
table(tRNAfullData$Type)
tRNAfullData$Type[7]
str(tRNAfullData)
AnotMINTmapUnic<-unique(tRNAfullData[,c(3,6)])

table(AnotMINTmapUnic$Type)
TodostRF <- unlist(strsplit(AnotMINTmapUnic$Type, ", "))

Seqs<- rep(AnotMINTmapUnic$tRF_seq, times = sapply(strsplit(AnotMINTmapUnic$Type, ", "), length))
# DataFrame
MINT_tRFs <- data.frame(Seqs=Seqs,Type=TodostRF)
table(TodostRF)
str(MINT_tRFs)
table(MINT_tRFs$Type)

# freqs tRFs VIPs 
Comm_tRFs<-subset(MINT_tRFs$Type, MINT_tRFs$Seqs %in% ComunesPLSDA)
Comm_tRFs<-data.frame(Comm_tRFs)
colnames(Comm_tRFs)<-"Tipus"
dim(Comm_tRFs)
str(Comm_tRFs)
table(Comm_tRFs)

colortRF<-c("cornflowerblue","cadetblue", "darkgray")
jpeg("./Grafics/PLSDA/Histograma_tRFComuns.jpeg", units="cm", width = 12, height=11,   res=300)
ggplot(data = Comm_tRFs, aes(x = "", fill = Tipus)) +
  geom_bar(width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "tRFs de les seqüències comuns", x = "", y = "Nombre") +
  scale_fill_manual(values = colortRF) +
  theme_void()+
  geom_text(aes(label = ..count..), stat = "count", position = position_stack(vjust = 0.5))
dev.off()

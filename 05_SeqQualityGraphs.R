library(readr)
library(ggplot2)
library(ggfortify)
library(gridExtra)


####### READS ORIGINALES #######

PerSeq_OR<- read_csv2("./Datos/Perseq_Original.csv")
summary(PerSeq_OR)
str(PerSeq_OR)

PhredScore_OR<- read_csv2("./Datos/PhredSeq_Original.csv")
PhredScore_OR<-PhredScore_OR[1:29,]
summary(PhredScore_OR)
str(PhredScore_OR)

Adapters_OR<- read_csv2("./Datos/Adapters3_original.csv")
Adapters_OR<-Adapters_OR[1:40,]
summary(Adapters_OR)
str(Adapters_OR)

# Quality per secuencia
jpeg("./Grafics/Preliminares/QualitatXbase1.jpeg", units="cm", width=11, height=10, res=300)
matplot(PerSeq_OR[,2:41], type="l", ylim=c(0,40), main = "Qualitat per base",
        xlab = "Posició (pb)", ylab = "Puntuació")
abline(h=c(20,28), col =c("red","orange"), lty=c(2,2))
dev.off()  

# Phred Score
jpeg("./Grafics/Preliminares/PhredScore1.jpeg", units="cm", width=11, height=10, res=300)
matplot(x=PhredScore_OR$Q_Score, y=PhredScore_OR$`01BLACK`, ylim=c(0,43000000), type = "l", main = "Puntuació Phred",
        xlab = "Puntuació mitjana", ylab = "Nombre de lectures")
matlines(x=PhredScore_OR$Q_Score, y=PhredScore_OR[,3:41], type = "l", lty = 1)     
dev.off()   

# Adapters
jpeg("./Grafics/Preliminares/Adaptadores1.jpeg", units="cm", width=11, height=10, res=300)
matplot(x=Adapters_OR[,2:41], type="l", ylim=c(0,100), main = "Adaptador 3' d'ARN petit",
        xlab = "Posició (pb)", ylab = "Seqüències (%)")
dev.off()

#######  SMALL RNA DATA #######

PerSeqQ_sRNA<- read_csv2("./Datos/PerseqQ_sRNA.csv")
summary(PerSeqQ_sRNA)
str(PerSeqQ_sRNA)

PhredScore_sRNA<- read_csv2("./Datos/PhredScore_sRNA.csv")
PhredScore_sRNA<-PhredScore_sRNA[1:25,]
summary(PhredScore_sRNA)
str(PhredScore_sRNA)

Adapters_sRNA<- read_csv2("./Datos/Adapters_sRNA.csv")
Adapters_sRNA<-Adapters_sRNA[1:32,]
summary(Adapters_sRNA)
str(Adapters_sRNA)

# Quality per seqüence
jpeg("./Grafics/Preliminares/QualitatXbase2.jpeg", units="cm", width=11, height=10, res=300)
matplot(PerSeqQ_sRNA[,2:41], type="l", ylim=c(0,40), main = "Qualitat per base",
        xlab = "Posició (pb)", ylab = "Puntuació")
abline(h=c(20,28), col =c("red","orange"), lty=c(2,2))
dev.off()

# Phred Score
jpeg("./Grafics/Preliminares/PhredScore2.jpeg", units="cm", width=11, height=10, res=300)
matplot(x=PhredScore_sRNA$Q_Score, y=PhredScore_sRNA$`01BLACK`, ylim=c(0,23000000), type = "l", main = "Puntuació Phred",
        xlab = "Puntuació mitjana", ylab = "Nombre de lectures")
matlines(x=PhredScore_sRNA$Q_Score, y=PhredScore_sRNA[,3:41], type = "l", lty = 1)  
dev.off()

# Adaptaers
jpeg("./Grafics/Preliminares/Adaptadores2.jpeg", units="cm", width=11, height=10, res=300)
matplot(x=Adapters_sRNA[,2:41], type="l", ylim=c(0,100), main = "Adaptador 3' d'ARN petit",
        xlab = "Posició (pb)", ylab = "Seqüències (%)")
dev.off()



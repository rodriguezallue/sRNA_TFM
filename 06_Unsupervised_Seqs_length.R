library(readr)
library(grid)
library(ggplot2)
library(ggfortify)
library(gridExtra)

# Loading data
# Loading 2 data sets: Read lengths and Groups
LengthFrq<- read_csv2("./Datos/Length_reads.csv")
SampleGr<-read_csv2("./Datos/Sample_groups.csv")
str(SampleGr, give.attr=FALSE)
SampleGr$MutGroup<- gsub("Healthy", "Ctrl", SampleGr$MutGroup)

SampleGr$Batch<-as.factor(SampleGr$Batch)
SampleGr$Group<-factor(SampleGr$MutGroup, levels=c("C9orf72", "Ctrl","MAPT","Pick","TDP43"), labels=c("C9orf72", "Ctrl","MAPT","Pick","TDP43"))
SampleGr$Group <- relevel(SampleGr$Group, "Ctrl")
SampleGr$Gender<-factor(SampleGr$GenderB, levels=c("F", "M","unknown"), labels=c("Famella", "Mascle","Desconegut"))

Longitud<-c(17:43)
LengT<-as.data.frame(t(LengthFrq[-1]))
colnames(LengT) <- LengthFrq$seq_length 
sampleID<- names(LengthFrq[-1])
LengT<-cbind(sampleID,LengT)
DadesSample<-SampleGr[,c("sampleID","Grup","Genere","AgeDeath","Batch")]

LengthDades<- merge(DadesSample, LengT, by = "sampleID", all = TRUE)
SamGrupo<-table(LengthDades$Group)
SamBatch<-table(LengthDades$Batch)
SamSex<-table(LengthDades$Gender)

jpeg("./Grafics/Preliminares/Demografia.jpeg", units="cm", width=11, height=10, res=300)
ggplot(data=DadesSample, aes(x = Group, y= AgeDeath)) + geom_point(aes(col=Gender), size=3)+
  labs(title= "Demografia de la cohort")+ ylab("Edat de mort (anys)") + xlab("Grup")+
  theme_bw()+theme(legend.title= element_blank())+ scale_color_manual(values = c("#f0b27a", "#ca6f1e","#808080"))
dev.off()

# Select subsets to make graphs

Ctrl = subset(LengthDades, Group=="Ctrl", drop=FALSE)
Ctrl_t<-cbind(Longitud,t(Ctrl[6:32]),rowMeans(t(Ctrl[6:32])))

C9orf72 =subset(LengthDades, Group=="C9orf72", drop=FALSE)
C9orf_t<-cbind(Longitud,t(C9orf72[6:32]), rowMeans(t(C9orf72[6:32])))

TDP43 =subset(LengthDades, Group=="TDP43", drop=FALSE)
TDP43_t<-cbind(Longitud,t(TDP43[6:32]), rowMeans(t(TDP43[6:32])))

Pick =subset(LengthDades, Group=="Pick", drop=FALSE)
Pick_t<-cbind(Longitud,t(Pick[6:32]), rowMeans(t(Pick[6:32])))

MAPT =subset(LengthDades, Group=="MAPT", drop=FALSE)
MAPT_t<-cbind(Longitud,t(MAPT[6:32]), rowMeans(t(MAPT[6:32])))


par(mfrow=c(2,3))

jpeg("./Grafics/Preliminares/LongControl.jpeg", units="cm", width=11, height=8, res=300)
matplot(x=Ctrl_t[,1], y=Ctrl_t[,2:8], type="l", col=2, ylim=c(0,0.5), main = "Grup Control",
        xlab = "Longitud (nt)", ylab = "Freqüència")
matlines(x=Ctrl_t[,1],Ctrl_t[,9], type = "l",col = 1, lwd = 2, lty = 1)
dev.off()

jpeg("./Grafics/Preliminares/LongC9orf72.jpeg", units="cm", width=11, height=8, res=300)
matplot(x=C9orf_t[,1], y=C9orf_t[,2:12], type="l", col=2, ylim=c(0,0.5), main = "Grup C9orf72",
        xlab = "Longitud (nt)", ylab = "Freqüència")
matlines(x=C9orf_t[,1],C9orf_t[,13], type = "l",col = 1, lwd = 2, lty = 1)
dev.off()

jpeg("./Grafics/Preliminares/LongTDP43.jpeg", units="cm", width=11, height=8, res=300)
matplot(x=TDP43_t[,1], y=TDP43_t[,2:10], type="l", col=2, ylim=c(0,0.5), main = "Grup TDP43",
        xlab = "Longitud (nt)", ylab = "Freqüència")
matlines(x=TDP43_t[,1],TDP43_t[,11], type = "l",col = 1, lwd = 2, lty = 1)
dev.off()

jpeg("./Grafics/Preliminares/LongPick.jpeg", units="cm", width=11, height=8, res=300)
matplot(x=Pick_t[,1], y=Pick_t[,2:8], type="l", col=2, ylim=c(0,0.5), main = "Grup Pick",
        xlab = "Longitud (nt)", ylab = "Freqüència")
matlines(x=Pick_t[,1],Pick_t[,9], type = "l",col = 1, lwd = 2, lty = 1)
dev.off()

jpeg("./Grafics/Preliminares/LongMAPT.jpeg", units="cm", width=11, height=8, res=300)
matplot(x=MAPT_t[,1], y=MAPT_t[,2:7], type="l", col=2, ylim=c(0,0.5), main = "Grup MAPT",
        xlab = "Longitud (nt)", ylab = "Freqüència")
matlines(x=MAPT_t[,1],MAPT_t[,8], type = "l",col = 1, lwd = 2, lty = 1)
dev.off()

par(mfrow=c(1,1))


# Ordenadas por batch
batchA = subset(LengthDades, Batch=="a", drop=FALSE)
batchA_t<-cbind(Longitud,t(batchA[6:32]),rowMeans(t(batchA[6:32])))
batchB = subset(LengthDades, Batch=="b", drop=FALSE)
batchB_t<-cbind(Longitud,t(batchB[6:32]),rowMeans(t(batchB[6:32])))
batchC = subset(LengthDades, Batch=="c", drop=FALSE)
batchC_t<-cbind(Longitud,t(batchC[6:32]),rowMeans(t(batchC[6:32])))
batchD = subset(LengthDades, Batch=="d", drop=FALSE)
batchD_t<-cbind(Longitud,t(batchD[6:32]),rowMeans(t(batchD[6:32])))

par(mfrow=c(2,3))
matplot(x=batchA_t[,1], y=batchA_t[,2:28], type="l", col=2, ylim=c(0,0.5), main = "Batch A",
        xlab = "Size (nt)", ylab = "Frequency")
matlines(x=batchA_t[,1],batchA_t[,29], type = "l",col = 1, lwd = 2, lty = 1)

matplot(x=batchB_t[,1], y=batchB_t[,2:5], type="l", col=2, ylim=c(0,0.5), main = "Batch B",
        xlab = "Size (nt)", ylab = "Frequency")
matlines(x=batchB_t[,1],batchB_t[,6], type = "l",col = 1, lwd = 2, lty = 1)

matplot(x=batchC_t[,1], y=batchC_t[,2:8], type="l", col=2, ylim=c(0,0.5), main = "Batch C",
        xlab = "Size (nt)", ylab = "Frequency")
matlines(x=batchC_t[,1],batchC_t[,9], type = "l",col = 1, lwd = 2, lty = 1)

matplot(x=batchD_t[,1], y=batchD_t[,2:3], type="l", col=2, ylim=c(0,0.5), main = "Batch D",
        xlab = "Size (nt)", ylab = "Frequency")
matlines(x=batchD_t[,1],batchD_t[,4], type = "l",col = 1, lwd = 2, lty = 1)

par(mfrow=c(1,1))

#Samples per batch
colores<-c("cadetblue","darkgoldenrod2","MediumPurple","salmon","RosyBrown")
jpeg("./Grafics/Preliminares/MostresXtanda.jpeg", units="cm", width=11, height=10, res=300)
ggplot(data=DadesSample, aes(x = Batch, fill = Group)) + geom_bar()+
  ylab("Nombre de mostres") + xlab("Tanda")+labs(title= "Mostres per tanda") + theme_bw()+
  scale_fill_manual(values=colores)+theme(legend.title= element_blank())
dev.off()

# Graphs of freq length per group (means)
ctrlgr<- ggplot(data = LengthFrq, aes(seq_length, Ctrl_t)) +   
    geom_point() + geom_line()+ theme_bw()+
    xlab("Length (nt)")+ ylab("Relative frequency")+ ggtitle("Controls")
c9orfgr<- ggplot(data = LengthFrq, aes(seq_length, C9orf_t)) +   
    geom_point() + geom_line()+ theme_bw()+
    xlab("Length (nt)")+ ylab("Relative frequency")+ ggtitle("C9orf72")
tdp43gr<-ggplot(data = LengthFrq, aes(seq_length, TDP43_t)) +   
      geom_point() + geom_line()+theme_bw()+
      xlab("Length (nt)")+ ylab("Relative frequency")+ ggtitle("TDP43")
pickgr<-ggplot(data = LengthFrq, aes(seq_length, Pick_t)) +   
      geom_point() + geom_line()+theme_bw()+
      xlab("Length (nt)")+ ylab("Relative frequency")+ ggtitle("Pick")
maptgr<- ggplot(data = LengthFrq, aes(seq_length, MAPT_t)) +   
      geom_point() + geom_line()+theme_bw()+
      xlab("Length (nt)")+ ylab("Relative frequency")+ ggtitle("MAPT")
grid.arrange(ctrlgr, c9orfgr, maptgr, tdp43gr, pickgr, nrow=2)

# Graphs of freq length per Batch
grafA<- ggplot(data = LengthFrq, aes(seq_length, batchA_t)) +   
    geom_point() + geom_line()+ 
    xlab("Length (nt)")+ ylab("Relative frequency") + ggtitle("Batch A")
grafB<-ggplot(data = LengthFrq, aes(seq_length, batchB_t)) +   
    geom_point() + geom_line()+
    xlab("Length (nt)")+ ylab("Relative frequency") + ggtitle("Batch B")
grafC<- ggplot(data = LengthFrq, aes(seq_length, batchC_t)) +   
    geom_point() + geom_line()+
    xlab("Length (nt)")+ ylab("Relative frequency") + ggtitle("Batch C")
grafD<- ggplot(data = LengthFrq, aes(seq_length, batchD_t)) +   
    geom_point() + geom_line()+
    xlab("Length (nt)")+ ylab("Relative frequency") + ggtitle("Batch D")
grid.arrange(grafA, grafB, grafC, grafD, nrow=2)

# samples per batch
Distrgr<- ggplot(data=SampleGr, aes(x = Batch, fill = MutGroup)) + geom_bar()+
  ylab("Number of samples") + labs("Group")

grid.arrange(grafA,grafB,grafC,grafD, Distrgr, ncol = 4, 
             layout_matrix = rbind(c(1,1,2,2), c(3,3,4,4), c(NA,5,5,5)))

# Arrange data for PCA

MatriuLeng <- as.matrix(LengT[,2:28],rownames=1)
#Perform the PCA
pcaLeng <- prcomp(MatriuLeng, scale = TRUE)
summary(pcaLeng)

jpeg("./Grafics/Preliminares/PCAtanda.jpeg", units="cm", width=11, height=10, res=300)
autoplot(pcaLeng, data = DadesSample, main="PCA per tanda", colour="Batch", size=3) +
  theme_bw()+ xlim(-0.3, 0.4)+ylim(-0.45,0.5)+
  theme(legend.position = c(0.1,0.85), legend.text = element_text(size=10),
          legend.key.size = unit(0.3, 'cm'), legend.title = element_blank())

jpeg("./Grafics/Preliminares/PCAsex.jpeg", units="cm", width=11, height=10, res=300)
autoplot(pcaLeng, data = DadesSample, main="PCA per gènere", colour="Gender", size=3)+
  theme_bw()+ xlim(-0.3, 0.4)+ylim(-0.45,0.5)+
  scale_color_manual(values = c("#f0b27a", "#ca6f1e","#808080"))+
  theme(legend.position = c(0.15,0.9), legend.text = element_text(size=10),
        legend.key.size = unit(0.3, 'cm'), legend.title = element_blank())
dev.off() 

jpeg("./Grafics/Preliminares/PCAmutacio.jpeg", units="cm", width=11, height=10, res=300)
autoplot(pcaLeng, data = DadesSample, main="PCA per grup de diagnòstic", colour="Group", size=3)+
  theme_bw()+ xlim(-0.3, 0.4)+ylim(-0.45,0.5)+
  scale_color_manual(values=colores)+
  theme(legend.position = c(0.15,0.85), legend.text = element_text(size=10),
        legend.key.size = unit(0.3, 'cm'), legend.title = element_blank())
dev.off()

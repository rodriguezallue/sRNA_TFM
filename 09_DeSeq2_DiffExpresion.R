library(readr)
library(DESeq2)
library(sqldf)
library(ggplot2)
library(factoextra)
library(ggfortify)
library(gridExtra)
library(RColorBrewer)
library(cluster)
library(ComplexHeatmap)
library(scatterplot3d)

#Loading and ordering case data
SampleGr<-read_csv2("./Datos/Sample_groups.csv")
SampleGr <- SampleGr[order(SampleGr$sampleID), ]
SampleGr$MutGroup<- gsub("Healthy", "Ctrl", SampleGr$MutGroup)
SampleGr$Group<-factor(SampleGr$MutGroup, levels=c("C9orf72", "Ctrl","MAPT","Pick","TDP43"), labels=c("C9orf72", "Ctrl","MAPT","Pick","TDP43"))
SampleGr$Group <- relevel(SampleGr$Group, "Ctrl")

#Loading ComBat-modified count data
CB_contajes<-read.csv("./Datos/MINT_CBatConts.csv", header = TRUE)
colnames(CB_contajes)[2:41]<-SampleGr$sampleID
rownames(CB_contajes)<-CB_contajes$X
CB_contajes<-CB_contajes[,2:41]

str(CB_contajes[1:8])
is.data.frame(CB_contajes)

#El set de datos original
tRNADataRed<-read.csv("./Datos/MINTmapData_Red.csv", header = TRUE)

#Boxplot de Raw data
boxplot(log10(CB_contajes), range = 0, las = 2, main=" ComBat-transformed counts")
class(CB_contajes)


#######################################
# DESeq analysis of CB_contajes uniq  #
#######################################

all(rownames(SampleGr$sampleID) == colnames(CB_contajes))

SampleGr$MutGroup<-as.factor(SampleGr$MutGroup)
class(SampleGr$MutGroup)
MutGr <- relevel(SampleGr$MutGroup, ref="Ctrl")
DsDataFin <- DESeqDataSetFromMatrix(countData = CB_contajes, colData = SampleGr,
              design= ~ GenderB + AgeDeath + MutGroup)

DsDataFin$MutGroup <- relevel(DsDataFin$MutGroup, ref="Ctrl")
rownames(DsDataFin) <- rownames(CB_contajes)

class(DsDataFin)
dim(DsDataFin)
head(DsDataFin)
ddsAnalisis<-DESeq(DsDataFin)

#Dispersion plot
jpeg("./Grafics/Unsupervised/DispPlot_MINT.jpeg", units="cm", width=12, height=10, res=300)
plotDispEsts(ddsAnalisis, main ="Gràfic de dispersió", xlab="Mitjana dels comptatges normalitzats",ylab="Dispersió" )
dev.off()

ResDiffExp<-results(ddsAnalisis)
ResDiffExpOrd_p <- ResDiffExp[order(ResDiffExp$pvalue),]
head(ResDiffExpOrd_p)

# Estimate the number of sequences with Differential Expression
sum(ResDiffExp$padj < 0.1, na.rm=TRUE)
sum(ResDiffExp$padj < 0.2, na.rm=TRUE)
sum(ResDiffExp$pvalue < 0.1, na.rm=TRUE)
sum(ResDiffExp$pvalue < 0.2, na.rm=TRUE)


##### Select all data with pval < 0.1 #######
resDExp <- as.data.frame(ResDiffExp)
resDExp$SeqID<-rownames(resDExp)
dim(resDExp)
str(resDExp)
Seleccion01<-sqldf("SELECT* FROM resDExp WHERE pvalue < 0.1 ORDER BY pvalue")
dim(Seleccion01)
str(Seleccion01)

ResDiffExpOrd_padj <- ResDiffExp[order(ResDiffExp$padj),]
head(ResDiffExpOrd_padj)

# Exploring Results
# Make general MA plot with labelled points with p-value < 0.1
umbral_p <- 0.1
# Filter p< threshold
pvalues <- ResDiffExp$pvalue
ind_sign <- pvalues < umbral_p
jpeg("./Grafics/DiffExpre/Signif_MINT.jpeg", units="cm", width=12, height=10, res=300)
plotMA(ResDiffExp, alpha = 0.1, ylim = c(-5, 5), main = "Sequences with p.value < 0.1")
points(ResDiffExp$baseMean[ind_sign], ResDiffExp$log2FoldChange[ind_sign], col = "blue", pch = 16, cex = 0.5)
dev.off() 

resultsNames(ddsAnalisis)

###### Results per Grups Vs Ctrl ########

colores<-c("cadetblue","darkgoldenrod2","MediumPurple","salmon","RosyBrown")

## TDP43 ##

res_TDP43<-results(ddsAnalisis, contrast=c("MutGroup","TDP43", "Ctrl"))
resTDP43_Ord_padj <- res_TDP43[order(res_TDP43$padj),]
head(resTDP43_Ord_padj)
length(which(res_TDP43$padj<0.1))

#Volcano plot
ResTDP43VPlot<-data.frame(res_TDP43)
head(ResTDP43VPlot)
ResTDP43VPlot$signif<-"NO"
ResTDP43VPlot$signif[ResTDP43VPlot$padj<0.1]<-"YES"
table(ResTDP43VPlot$signif)

jpeg("./Grafics/DiffExpre/Volcano_TDP43_MINT.jpeg", units="cm", width=12, height=10, res=300)
ggplot(data=ResTDP43VPlot, aes(x=log2FoldChange, y=-log10(padj), col=signif))+ 
   geom_point()+ labs(title="TDP43 vs control", x= "log2 (FCh)", y="-log10(p ajustada)") +
   xlim(-6,6.5) + ylim(0,3)+
   geom_hline(yintercept=-log10(0.1), linetype="dashed", col = "black")+ 
   geom_vline(xintercept = c(-1.4, 1.4), linetype = "dashed", col = "black")+
   scale_color_manual (values = c("NO" = "darkgray","YES" = "RosyBrown"))+ 
   theme_bw()+ theme(legend.position ="none")
dev.off() 

## MAPT ##

res_MAPT<-results(ddsAnalisis, contrast=c("MutGroup","MAPT","Ctrl"))
res_MAPT_Ord_padj <- res_MAPT[order(res_MAPT$padj),]
head(res_MAPT_Ord_padj)
length(which(res_MAPT$padj<0.1))

#Volcano plot
ResMAPTVPlot<-data.frame(res_MAPT)
head(ResMAPTVPlot)
ResMAPTVPlot$signif<-"NO"
ResMAPTVPlot$signif[ResMAPTVPlot$padj<0.1]<-"YES"
table(ResMAPTVPlot$signif)

jpeg("./Grafics/DiffExpre/Volcano_MAPT_MINT.jpeg", units="cm", width=12, height=10, res=300)
ggplot(data=ResMAPTVPlot, aes(x=log2FoldChange, y=-log10(padj), col=signif))+ 
   geom_point()+ labs(title="MAPT vs control", x= "log2 (FCh)", y="-log10(p ajustada)") +
   xlim(-6,6.5) + ylim(0,3)+
   geom_hline(yintercept=-log10(0.1), linetype="dashed", col = "black")+ 
   geom_vline(xintercept = c(-1.4, 1.4), linetype = "dashed", col = "black")+
   scale_color_manual (values = c("NO" = "darkgray","YES" = "MediumPurple"))+ 
   theme_bw()+ theme(legend.position ="none")
dev.off()

## PICK ##

res_Pick<-results(ddsAnalisis, contrast=c("MutGroup","Pick","Ctrl"))
res_Pick_Ord_padj <- res_Pick[order(res_Pick$padj),]
head(res_Pick_Ord_padj)
length(which(res_Pick$padj<0.1))

# Volcano_plot
ResPickVPlot<-data.frame(res_Pick)
head(ResPickVPlot)
ResPickVPlot$signif<-"NO"
ResPickVPlot$signif[ResPickVPlot$padj<0.1]<-"YES"
table(ResPickVPlot$signif)

jpeg("./Grafics/DiffExpre/Volcano_Pick_MINT.jpeg", units="cm", width=12, height=10, res=300)
ggplot(data=ResPickVPlot, aes(x=log2FoldChange, y=-log10(padj), col=signif))+ 
     geom_point()+ labs(title="Pick vs control", x= "log2 (FCh)", y="-log10(p ajustada)") +
     xlim(-6,6.5) + ylim(0,3)+
     geom_hline(yintercept=-log10(0.1), linetype="dashed", col = "black")+ 
     geom_vline(xintercept = c(-1.4, 1.4), linetype = "dashed", col = "black")+
     scale_color_manual (values = c("NO" = "darkgray","YES" = "salmon"))+ 
     theme_bw()+ theme(legend.position ="none")
dev.off() 

## C9orf72 ##

res_C9orf<-results(ddsAnalisis, contrast=c("MutGroup","C9orf72","Ctrl"))
res_C9orf_Ord_padj <- res_C9orf[order(res_C9orf$padj),]
head(res_C9orf_Ord_padj)
length(which(res_C9orf$padj<0.1))

# Volcano_plot
ResC9orfVPlot<-data.frame(res_C9orf)
head(ResC9orfVPlot)
ResC9orfVPlot$signif<-"NO"
ResC9orfVPlot$signif[ResC9orfVPlot$padj<0.1]<-"YES"
table(ResC9orfVPlot$signif)

jpeg("./Grafics/DiffExpre/Volcano_C9orf_MINT.jpeg", units="cm", width=12, height=10, res=300)
ggplot(data=ResC9orfVPlot, aes(x=log2FoldChange, y=-log10(padj), col=signif))+ 
  geom_point()+ labs(title="C9orf72 vs Ctrl", x= "log2 (FCh)", y="-log10 (p ajustada)") +
  xlim(-6,6.5) + ylim(0,3)+
  geom_hline(yintercept=-log10(0.1), linetype="dashed", col = "black")+ 
  geom_vline(xintercept = c(-1.4, 1.4), linetype = "dashed", col = "black")+
  scale_color_manual (values = c("NO" = "darkgray","YES" = "darkgoldenrod2"))+ 
  theme_bw()+ theme(legend.position ="none")
dev.off() 

## Summary of results ###
summary(res_C9orf)
summary(res_MAPT)
summary(res_Pick)
summary(res_TDP43)


### NON SUPERVISED ANALYSIS WITH SIGNIFICANT SEQUENCESS 

colorsPick<-c("cadetblue","gray","gray","salmon","gray")


#####  Seqs with padj <0.1 from Pick  and make PCA ####

res_Pick.df<-as.data.frame((res_Pick))
dim(res_Pick.df)[1]
#SeqIDs significativos (P-adj < 0.1)
PickSigIDs <- rownames(res_Pick.df)[which(res_Pick.df$padj<0.1)]
#Secuencias significativas (P-adj < 0.1)
SeqsPickSign <- tRNADataRed$Sequence [tRNADataRed$SeqID %in% PickSigIDs]


#Keep them in files
write.table(PickSigIDs,"./Datos/Diff_exp/SeqID_MINT_Pick01.txt",quote=FALSE)
write.table(SeqsPickSign,"./Datos/Diff_exp/Seqs_MINT_Pick01.txt",quote=FALSE)

Cont_Pick05 <- CB_contajes[rownames(CB_contajes) %in% PickSigIDs, 1:40]
dim(Cont_Pick05)
# rlog Normalization
rlgCont_Pick<-rlog(as.matrix(Cont_Pick05))


# Plotting PCAs to compare results by groups
Pick_PCA<-prcomp(t(rlgCont_Pick), scale=TRUE)
fviz_eig(Pick_PCA, addlabels = TRUE)


pca_signPick <- data.frame(SampleGr$Group, Pick_PCA$x)
colnames(pca_signPick)[1]<-"Group"

pca_signPick$Group<-factor(pca_signPick$Group, levels=c("C9orf72", "Ctrl","MAPT","Pick","TDP43"),
                         labels=c("C9orf72", "Ctrl","MAPT","Pick","TDP43"))
pca_signPick$Group <- relevel(pca_signPick$Group, "Ctrl")

jpeg("./Grafics/DiffExpre/PCApick_MINT.jpeg", units="cm", width=11, height=10, res=300)
autoplot(Pick_PCA, data = SampleGr, main="Pick vs control - adj. p < 0.1", colour="Group", size=3) +
  theme_bw()+ xlim(-0.35, 0.4)+ylim(-0.4,0.6)+
  scale_color_manual(values=colorsPick)+
  theme(legend.position = c(0.15,0.87), legend.text = element_text(size=10),
        legend.key.size = unit(0.3, 'cm'), legend.title = element_blank())
dev.off()

jpeg("./Grafics/DiffExpre/PCA3Dpick_MINT.jpeg", units="cm", width=11, height=10, res=300)
scatterplot3d(pca_signPick[2:4], pch = 16, highlight.3d = FALSE, type = "h", main = "3 PC", angle = 50, color = colorsPick[pca_signPick$Group])
dev.off()
scatterplot3d(pca_signPick[2:4], pch = 16, highlight.3d = FALSE, type = "h", main = "3D Scatterplot - angle 45", angle = 75, color = colorsPick[pca_signPick$Group])

scatterplot3d(pca_signPick[2:4], pch = 16, highlight.3d = FALSE, type = "h", main = "3D Scatterplot - angle 50", angle = 50, color = colorsPick[pca_signPick$Group])




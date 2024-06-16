library(readr)
library(dplyr)
library(Biobase) # Para construir y manipular el expressionSet
library(DESeq2)
library(ggplot2)
library(apeglm)
library(gridExtra)
library(ggfortify)
library(gg3D) #graficas 3D
library(cowplot) # grids graficas 3D
library(factoextra)
library(sva) # Analiza variables subrogadas y el efecto batch
library(sqldf)


#Laoding and ordering case data
SampleGr<-read_csv2("./Datos/Sample_groups.csv")
SampleGr <- SampleGr[order(SampleGr$sampleID), ]

Dataset1<-read.table("./Datos/all_sequence_tRNA_level_ok.csv", sep = "\t", header = TRUE)

# Rearrange and rename the matrix columns
dim(Dataset1)
str(Dataset1)
Orden<- c(1,2,3,4,17,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,
        35,36,37,38,39,40,41,42,43,5,6,7,8,9,10,11,12,13,14,15,16,18)

Dataset1 <- Dataset1[,Orden]
colnames(Dataset1)[(ncol(Dataset1) - 39):ncol(Dataset1)] <- SampleGr$sampleID

#Add a column of our sequence IDs with dplyr, to identify the sequences 
Dataset1 <- Dataset1 %>% mutate(SeqID = paste("Seq_", row_number(), sep = "")) %>%
    select(SeqID, everything())
dim(Dataset1)
str(Dataset1)
Seqs10<-Dataset1[32850:32855,c(1,4:8)]

head(Dataset1)
tail(Dataset1)

# Generate the DesSeqDataSet 
MutGroup <- factor(SampleGr$MutGroup)
Batch <- factor(SampleGr$Batch)
MatriuDataSet <- DESeqDataSetFromMatrix(countData = Dataset1[,5:44], colData = SampleGr,
          design= ~ Batch + GenderB + AgeDeath + MutGroup)
rownames(MatriuDataSet) <- Dataset1$SeqID

# Remove sequences with too low counts (50% column/groups with less than 5 counts)
smallestGroupSize <- 4
keepLarge <- rowSums(counts(MatriuDataSet) >= 5) >= smallestGroupSize
MatriuDataRed <- MatriuDataSet[keepLarge,]
Seqs2<-dim(MatriuDataRed)[1]

class(MatriuDataRed)

# Reduced dataframe with the selected counts.
filasMatR <- rownames(counts(MatriuDataRed))
DadesSel<-which(Dataset1$SeqID %in% filasMatR)
Dataset_Red<- Dataset1[DadesSel,]
dim(Dataset_Red)

if(length((filasMatR)) == length(DadesSel)){print("TRUE")}         #It is not the same length!!!!!!!!

table(Dataset1$Sequence)[1:15]
table(filasMatR)[1:10]
table(Dataset_Red$Sequence)[1:10]


class(Dataset1)
class(MatriuDataRed)

#Have a look to the selected data.
MatriuSel <- counts(MatriuDataRed)
RedSeqs10<-MatriuSel[856:860,1:5]
class(MatriuSel)

# Make the DESeq model of the reduced dataset
ModRed <- DESeq(MatriuDataRed, )
# Graphs with the data distribution 
boxplot(log2(assays(ModRed)[["cooks"]]), range=0, las=2,
        main="Distribution of counts by sample")


# Custom Transformation
ModRedSF <- estimateSizeFactors(ModRed)
# shifted log of normalized counts
se <- SummarizedExperiment(log2(counts(ModRedSF, normalized=TRUE)+1), colData=colData(ModRedSF))
# the call to DESeqTransform() is needed to trigger our plotPCA method.
plotPCA(DESeqTransform(se), intgroup ="MutGroup")

# Other different normalizations
vst_ModRed <- vst(ModRed)
rlog_ModRed <- rlog(ModRed)

boxplot(log10(counts(vst_ModRed)), range = 0, las = 2)

plotDispEsts(rlog_ModRed)
vst_matriz<-as.data.frame(log10(assay(vst_ModRed)))
boxplot(vst_matriz, range = 0, las = 2)
boxplot(log10(assays(vst_ModRed)[["cooks"]]), range=0, las=2)

class(vst_ModRed)
Mat_rlog <- assay(rlog_ModRed)

# Plotting PCAs to compare results by groups
ModPCA<-plotPCA(DESeqTransform(se), intgroup ="MutGroup")+ ggtitle("DESeq2 Nomalization")+
  xlim (-50, 75)+ ylim(-40,50)+ geom_point(size=2)
vstPCA<-plotPCA(vst_ModRed, intgroup = "MutGroup") + ggtitle("VST Nomalization")+ 
  xlim (-50, 75)+ ylim(-40,50)+ geom_point(size=2)
rlogPCA<-plotPCA(rlog_ModRed, intgroup = "MutGroup")+ ggtitle("rlog Normalization")+
  xlim (-50, 75)+ ylim(-40,50)+ geom_point(size=2)
grid.arrange(ModPCA, vstPCA, rlogPCA, ncol=1, nrow=3)

# Plotting PCAs to compare results by batch
ModPCAb<-plotPCA(DESeqTransform(se), intgroup ="Batch")+ ggtitle("DESeq2 Nomalization")+
  xlim (-50, 75)+ ylim(-40,50)+ geom_point(size=2)
vstPCAb<-plotPCA(vst_ModRed, intgroup = "Batch") + ggtitle("VST Nomalization")+ 
  xlim (-50, 75)+ ylim(-40,50)+ geom_point(size=2)
rlogPCAb<-plotPCA(rlog_ModRed, intgroup = "Batch")+ ggtitle("rlog Normalization")+
  xlim (-50, 75)+ ylim(-40,50)+ geom_point(size=2)
grid.arrange(ModPCAb, vstPCAb, rlogPCAb, ncol=1, nrow=3)

grid_arr <- grid.arrange(ModPCA, vstPCA, rlogPCA, ncol=1, nrow=3)

# Ajustar los márgenes y tamaños de los gráficos
grid_arr$widths <- unit(c(11, "cm"), "null")  
grid_arr$heights <- unit(c(8.5, 8.5, 8.5), "cm")
print(grid_arr)


###############################
#  ANALYSIS OF THE BACH EFECT #
###############################

# Extrae la matriz de diseño del objeto DESeqDataSet
diseny_Datos <- MatriuDataRed@colData
disenyMat <- model.matrix(~ Batch + GenderB + AgeDeath + MutGroup, data = diseny_Datos)
class(disenyMat)
Mod0 = model.matrix(~Batch, data=diseny_Datos)

###### Analizo el batch effect usando comBat ######

CB_contajes<-ComBat_seq(counts(MatriuDataRed), batch=Batch, group=NULL)

pVals_CB = f.pvalue(CB_contajes,disenyMat, Mod0)
qVals_CB = p.adjust(pVals_CB,method="BH")
# % de genes que son diferentes de 0: Son el 25%
sum(qVals_CB<0.05)/length(qVals_CB)

class(CB_contajes)
tail(CB_contajes)
min(CB_contajes)

write.table(CB_contajes, file = "./Datos/CB_contajes.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

#Boxplot de Raw data
boxplot(log10(CB_contajes), range = 0, las = 2, main=" ComBat-transformed counts")

# Normalizaciones
rlogCB_cont <- rlog(CB_contajes)
boxplot(rlogCB_cont, range = 0, las = 2, main=" ComBat-transformed counts, rlog normalization")

vstCB_cont <- vst(CB_contajes)
boxplot(vstCB_cont, range = 0, las = 2, main=" ComBat-transformed counts, vst normalization")


# Plotting PCAs to compare results by groups
CB_contPCA<-prcomp(t(CB_contajes), scale=FALSE)
rlogCB_PCA<-prcomp(t(rlogCB_cont), scale=FALSE)

autoplot(CB_contPCA, data = SampleGr, colour = 'Batch', label=FALSE)+
  ggtitle("ComBat-transformed data")
autoplot(rlogCB_PCA, data = SampleGr, colour = 'Batch', label = FALSE)+
  ggtitle("ComBat-transformed data after rlog normalization")


class(rlogCB_PCA)
plot(rlogCB_PCA)
fviz_eig(rlogCB_PCA, addlabels = TRUE)


# Grafics 3D amb diferents paràmetres

# Gràfics PCA per phenotype
pca_rlogCB <- data.frame(SampleGr$MutGroup, rlogCB_PCA$x)

names(pca_rlogCB) = c("MutGroup", "PC1", "PC2", "PC3")

angle1p <- ggplot(pca_rlogCB, aes(x=PC1, y=PC2, z=PC3, geom="blank", colour = MutGroup)) +
  theme_void()+
  axes_3D()+
  stat_3D()

angle2p <- ggplot(pca_rlogCB, aes(x=PC1, y=PC2, z=PC3, geom="blank", colour = MutGroup)) +
  theme_void()+
  axes_3D(theta=100)+
  stat_3D(theta=100)

angle3p <- ggplot(pca_rlogCB, aes(x=PC1, y=PC2, z=PC3, geom="blank", colour = MutGroup)) +
  theme_void()+
  axes_3D(theta=170)+
  stat_3D(theta=170)

angle4p <- ggplot(pca_rlogCB, aes(x=PC1, y=PC2, z=PC3, geom="blank", colour = MutGroup)) +
  theme_void()+
  axes_3D(phi=20)+
  stat_3D(phi=20)

plot_grid(angle1p,angle2p,angle3p,angle4p, ncol=2)


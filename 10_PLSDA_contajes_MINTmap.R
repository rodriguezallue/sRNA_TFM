library(readr)
library(mixOmics)
library(matrixStats)
library(ggplot2)
library(ggfortify)
library(gridExtra)
library(factoextra)
library(RColorBrewer)
library(scatterplot3d)
library(sqldf)
library(cluster)
library(ComplexHeatmap)
library(dendextend)
library(circlize)
library(ggVennDiagram)


#Loading and ordering case data
SampleGr<-read_csv2("./Datos/Sample_groups.csv")
SampleGr <- SampleGr[order(SampleGr$sampleID), ]
SampleGr$MutGroup<- gsub("Healthy", "Ctrl", SampleGr$MutGroup)
SampleGr$Group<-factor(SampleGr$MutGroup, levels=c("C9orf72", "Ctrl","MAPT","Pick","TDP43"), labels=c("C9orf72", "Ctrl","MAPT","Pick","TDP43"))
SampleGr$Group <- relevel(SampleGr$Group, "Ctrl")

#Loading count data
rlogConts<-read.table("./Datos/rlogMINTconts.txt", header = TRUE)
colnames(rlogConts)<-SampleGr$sampleID
str(rlogConts[1:8])
is.data.frame(rlogConts)

# PERFORM PLS-DA 

# Rearrange data into a matrix
Xconts <- t(as.matrix(rlogConts,rownames=1))
Y <- SampleGr$Group

# Execute the model
Counts.plsda <-plsda(Xconts,Y, ncomp=15, scale=TRUE, near.zero.var = TRUE)

# Plotting the results
colores<-c("cadetblue","darkgoldenrod2","MediumPurple","salmon","RosyBrown")

#2D plot
jpeg("./Grafics/PLSDA/PLSDA1_MINT.jpeg", units="cm", width=12, height=11, res=300)
plotIndiv(Counts.plsda , comp = 1:2, 
          group = SampleGr$MutGroup, ind.names = FALSE,  # colour points by class
          ellipse = FALSE,
          col=colores, X.label='Comp 1 (19%)', Y.label='Comp 2 (8%)',
          legend = TRUE, title = 'PLS-DA with confidence ellipses', size.title = 1, style = 'ggplot2')
dev.off()

#3D plot
plotIndiv(Counts.plsda , comp = 1:3, 
          group = SampleGr$MutGroup, ind.names = FALSE,    # color points by class
          ellipse = TRUE, col=colores, #X.label = 'Comp 1 (19%)', Y.label = 'Comp 2 (8%)',Z.label='Comp 3 (5%)',
          legend = TRUE, title='PLS-DA with confidence ellipses', style="3d")

# PREDICTION MAP use the max.dist measure to form decision boundaries between classes based on PLS-DA
background = background.predict(Counts.plsda, comp.predicted=2, dist = "max.dist")
#tiff("PredictionMap.tiff", units="cm", width = 12, height=11,   res=400)
plotIndiv(Counts.plsda, comp = 1:2,
          group = SampleGr$MutGroup, ind.names = FALSE, col=colores,
          background = background, X.label = 'Comp 1 (19%)', Y.label = 'Comp 2 (8%)',
          legend = TRUE, title = "PLS-DA prediction background", size.title = 1, style="ggplot2")
#dev.off()

# Cross-validation with perf() to calculate the AUC and the BER of the model
set.seed(43210)
perf.plsda <- perf(Counts.plsda, validation="Mfold", folds=4, nrepeat=50, near.zero.var = TRUE, progressBar = FALSE, auc=TRUE) 
set.seed(NULL)

plot(perf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

# The optimal number of components
perf.plsda$choice.ncomp
perf.plsda$choice.ncomp[2,2]

# The error.rate and the AUC of the last optimal component
perf.plsda$error.rate$BER[[3]]
mean(perf.plsda$auc$comp1[,1])

perf.plsda$error.rate
perf.plsda$auc

# ESTIMATION OF THE VARIABLE IMPORTANCE IN THE PROJECTION (VIP) AND SELECTION OF VIPs >= 1.2
vip.Counts<-as.data.frame(vip(Counts.plsda))
vip_medio<- rowMeans(vip.Counts)
SeqID<-rownames(rlogConts)
VIPs_Cont<-cbind(SeqID,rlogConts,vip_medio)
#write.csv(VIPs_lip, file = "./Dades/VIPs_lipids.csv", row.names = FALSE, col.names=FALSE)
VIPs <-sqldf("SELECT * FROM VIPs_Cont WHERE vip_medio >= 1.2 ORDER BY SeqID")
dim(VIPs)
summary(VIPs)
VIPs[c("SeqID","vip_medio")]
ListaVIPs<-VIPs$SeqID

#Select VIPs counts for a second PLS-DA
dfCont_vip<-VIPs_Cont[VIPs_Cont$SeqID %in% ListaVIPs,]
dim(dfCont_vip)
str(dfCont_vip)

#Second PLS-DA
# Rearrange data into a matrix
Xvip <- t(as.matrix(dfCont_vip[,2:41], rownames=1))
dim(Xvip)
Xvip[1:5,1:5]

# Execute the model
VIP.plsda <-plsda(Xvip,Y, ncomp=15, scale=TRUE, near.zero.var = TRUE)

#2D plot
jpeg("./Grafics/PLSDA/PLSDAvips_MINT.jpeg", units="cm", width=12, height=11, res=300)
plotIndiv(VIP.plsda , comp = 1:2, 
          group = SampleGr$MutGroup, ind.names = FALSE,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          col=colores, X.label='Comp 1 (16%)', Y.label='Comp 2 (10%)',
          legend = TRUE, title = 'PLS-DA with confidence ellipses', size.title = 1, style = 'ggplot2')
dev.off()

#3D plot
plotIndiv(VIP.plsda , comp = 1:3, 
          group = SampleGr$MutGroup, ind.names = FALSE,    # color points by class
          ellipse = TRUE, col=colores, X.label = 'Comp 1 (16%)', Y.label = 'Comp 2 (10%)',Z.label='Comp 3 (7%)',
          legend = TRUE, title='PLS-DA with confidence ellipses', style="3d")

# PREDICTION MAP use the max.dist measure to form decision boundaries between classes based on PLS-DA
background2 = background.predict(VIP.plsda, comp.predicted=2, dist = "max.dist")
jpeg("./Grafics/PLSDA/PredMap_VIPS_MINT.jpeg", units="cm", width=12, height=11, res=300)
plotIndiv(VIP.plsda, comp = 1:2,
          group = SampleGr$MutGroup, ind.names = FALSE, col=colores,
          background = background2, X.label = 'Comp 1 (13%)', Y.label = 'Comp 2 (12%)',
          legend = TRUE, title = "PLS-DA prediction background", size.title = 1, style="ggplot2")
dev.off()

# Cross-validation with perf() to calculate the AUC and the BER of the model
set.seed(43210)
perf2.plsda <- perf(VIP.plsda, validation="Mfold", folds=4, nrepeat=150, near.zero.var = TRUE, progressBar = FALSE, auc=TRUE) 
set.seed(NULL)

plot(perf2.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

# The optimal number of components
perf2.plsda$choice.ncomp
NumComp<-perf2.plsda$choice.ncomp[2,1]

# The error.rate and the AUC of the last optimal component
perf2.plsda$error.rate$BER[6,]
perf2.plsda$error.rate$BER[[perf2.plsda$choice.ncomp[2,1]]]
perf2.plsda$auc$comp6

mean(perf2.plsda$auc$comp6[,1])
max(perf2.plsda$auc$comp6[,1])
perf2.plsda$error.rate
mean(perf2.plsda$auc$comp6[,1])
perf2.plsda$auc[7][[1]][[1]]

############ ITERATIONS TO GET THE P VALUE OF THE MODEL ###########

set.seed(43210)
# Define vectors and counters to collect data
AUC<-c()    
ERp<-c() 
numAUC<-0
numERp<-0
SeqID<-rownames(rlogConts)

# Permutation system
for(i in 1:1){
   # Vector 'Group' with 40 elements in 5 random groups
   Yval<-sample(1:5, size = length(SampleGr$MutGroup), replace = TRUE)
   #Run the model
   Iter_Counts.plsda <- plsda(Xconts, Yval, scale=TRUE, ncomp=15, near.zero.var=TRUE )
   #Collect VIPS > 1.2
   Iter.vips<-as.data.frame(vip(Iter_Counts.plsda))
   Iter.vip_medio<- rowMeans(Iter.vips)
   Iter.VIPs_Count<-cbind(SeqID,rlogConts,Iter.vip_medio)
   Iter.VIPs12 <-sqldf("SELECT * FROM Iter.VIPs_Count WHERE Iter.vip_medio >= 1.2 ORDER BY SeqID")
   dim(Iter.VIPs12)
   summary(Iter.VIPs12)
   Iter.VIPs12[c("SeqID","vip_medio")]
   Iter.ListaVIPs<-Iter.VIPs12$SeqID
   #Select counts from VIPS for second PLS-DA
   dfCont_vipI<-data.frame()
   dfCont_vipI<-Iter.VIPs_Count[Iter.VIPs_Count$SeqID %in% Iter.ListaVIPs,]
   dim(dfCont_vipI)
   str(dfCont_vipI)
   
   XvipI <- t(as.matrix(dfCont_vipI[,2:41], rownames=1))
   dim(XvipI)
   Xvip[1:5,1:5]
   
   Iter_VIPs.plsda <- plsda(XvipI, Yval, scale=TRUE, ncomp = 15 )
  
   perf.Iter_VIPS <- perf(Iter_VIPs.plsda.spls, validation="Mfold",fold =5, nrepeat=150,
              progressBar = FALSE, auc = TRUE)
  # Select the optimal number of Components (max.dist, BER)
   NumCompI<-perf.Iter_VIPS$choice.ncomp[2,1]                      
   # Vectors with the BER and the AUC of the optimal components
   ERp<-c(ERp, perf.Iter_VIPS$error.rate$BER[[NumCompI]])
   AUC<-c(AUC, perf.Iter_VIPS$auc[NumCompI][[1]][[1]])
   # counting the number of analysis with better BER and AUC values
   if (perf.Iter_VIPS$error.rate$BER[[NumCompI]] <= 0.2206022) {numERp=numERp+1}  
   if (perf.Iter_VIPS$auc[NumCompI][[1]][[1]] >= 0.97152) {numAUC=numAUC+1}
}
set.seed(NULL)

# Plots of the iterations statistics.
AUC_df<-as.data.frame(AUC)
plot_AUC<- ggplot(AUC_df, aes(AUC)) + 
  geom_histogram (aes(y=..density..), colour="black", fill="white") + labs(x="AUC") +
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=0.975364), color="blue", linetype="dashed")

ERp_df<-as.data.frame(ERp)
plot_ERp<-ggplot(ERp_df, aes(ERp)) + 
  geom_histogram (aes(y=..density..), colour="black", fill="white") + labs(x="Error rate") +
  geom_density(alpha=.2, fill="#5499C7") +
  geom_vline(aes(xintercept=0.0777778), color="Red", linetype="dashed")

grid.arrange(plot_AUC, plot_ERp, ncol=2, nrow =1)


####SELECCTION OF SEQUENCES (VIP >1.1) #########

# SECOND ESTIMATION OF THE VIP AND SELECTION OF VIPs >= 1.1
vip.plsda2<-as.data.frame(vip(VIP.plsda))
mean_vip<- rowMeans(vip.plsda2)
min(mean_vip)
max(mean_vip)
VIPs_fin<-data.frame()
str(dfCont_vip)
VIPs_fin<-cbind(dfCont_vip[,-42],mean_vip)
str(VIPs_fin)
dim(dfCont_vip)
dim(VIPs_fin)
#write.csv(VIPs_lip, file = "./Dades/VIPs_lipids.csv", row.names = FALSE, col.names=FALSE)
VIPs_fin <-sqldf("SELECT * FROM VIPs_fin WHERE mean_vip >= 1.1 ORDER BY SeqID")
dim(VIPs_fin)
str(VIPs_fin)
summary(VIPs_fin)
FinalSeqs<-VIPs_fin[c("SeqID","mean_vip")]
FinListVIPs<-VIPs_fin$SeqID

#Keep SeqIDs from VIPs
write.table(FinListVIPs,"./Datos/FinalVIPs_PLSDA_mint.txt")

VIPsGraf <-sqldf("SELECT * FROM VIPs_fin WHERE mean_vip >= 1.2 ORDER BY SeqID")

#Lolliplot of VIPS 
#tiff("VIP_Analysis.tiff", units="cm", width = 12, height=11,   res=400)
ggplot(data= VIPsGraf, aes(x=SeqID, y=mean_vip))+
  geom_segment(aes(x=reorder(SeqID,mean_vip), xend=SeqID, y=1.1, yend=mean_vip))+
  xlab("Sequence ID") + ylab("VIP score") + labs(title="Sequences with the highest VIP scores")+
  coord_flip() +
  geom_point(size=3, color=("#a0522d"), fill=("#a1887f"), shape=21, stroke=1)+
  theme_bw()
dev.off()

# Generate the PCA with the final VIPs
MatVIPf<- as.matrix(VIPs_fin[,2:41], rownames=1)
pcaVIPs <- prcomp(t(MatVIPf), scale = TRUE)
pcaVIPs$x
summary(pcaVIPs)
#Optimal number of simensions (10)
fviz_eig(pcaVIPs, ncp = 40, addlabels = TRUE, main = "Optimal number of dimensions")

# Graphic of PCA results
pca_VIPsdf <- data.frame(SampleGr$MutGroup, pcaVIPs$x)
colnames(pca_VIPsdf)[1]<-"Group"
pca_VIPsdf$Group

pca_VIPsdf$Group<-factor(pca_VIPsdf$Group, levels=c("C9orf72", "Ctrl","MAPT","Pick","TDP43"),
                         labels=c("C9orf72", "Ctrl","MAPT","Pick","TDP43"))
pca_VIPsdf$Group <- relevel(pca_VIPsdf$Group, "Ctrl")


jpeg("./Grafics/PLSDA/PCA_VIPfin_MINT.jpeg", units="cm", width=12, height=11, res=300)
autoplot(pcaVIPs, data = SampleGr, main="PCA amb seqüènces VIP > 1.1", colour="Group",size=2)+
  scale_colour_manual(values=colores)+
  theme(legend.position="bottom") + theme_bw()
dev.off()

#3D graph of PCA results
colores<-c("darkgoldenrod2","cadetblue","MediumPurple","salmon","RosyBrown")

colores2 <- c("RosyBrown","darkgoldenrod2","salmon","MediumPurple","cadetblue") 

scatterplot3d(pca_VIPsdf[2:4], pch = 16, highlight.3d = FALSE, type = "h", main = "3D Scatterplot - angle 75", angle = 75, color = colores[pca_VIPsdf$Group])
legend("topright", legend = unique(pca_VIPsdf$Group), col = colores2, pch = 16, title = "Group")

scatterplot3d(pca_VIPsdf[2:4], pch = 16, highlight.3d = FALSE, type = "h", main = "3D Scatterplot - angle 45", angle = 45, color = colores[pca_VIPsdf$Group])
legend("topright", legend = unique(pca_VIPsdf$Group), col = colores2, pch = 16, title = "Group")

#Annotations of colums and lipids of the Heatmap
# Data standarization
Zscore<- function(x) {
  Z=(x-rowMeans(x))/rowSds(as.matrix(x), useNames = FALSE)
  return(Z)}

MatVIPz<- Zscore(MatVIPf)
rownames(MatVIPz)<-VIPs_fin$SeqID

ColorHM=colorRampPalette(c("navy", "white", "firebrick3"))(50)
FiveCols=list(Group=c("TDP43"="RosyBrown","Ctrl"="cadetblue","Pick"="salmon","MAPT"="MediumPurple","C9orf72"="darkgoldenrod2"))
column_ha<-HeatmapAnnotation(Group=SampleGr$MutGroup, col =FiveCols)

# Distancess Pearson
Heatmap(MatVIPz, name="Z-score", col=ColorHM, show_row_names=FALSE, show_column_names=FALSE, top_annotation = column_ha,
        show_row_dend = FALSE)
# Clustering Agnes
Heatmap(MatVIPz, name="Z-score", col=ColorHM, show_row_names=FALSE, show_column_names=FALSE, top_annotation = column_ha,
        show_row_dend = FALSE, cluster_columns = agnes(t(MatVIPz)))

##### Columns Clustering Diana  #######
Heatmap(MatVIPz, name="Z-score", col=ColorHM, show_row_names=FALSE, show_column_names=FALSE, top_annotation = column_ha,
        show_row_dend = FALSE, cluster_columns = diana(t(MatVIPz)))

##### rows Diana and Columns Clustering Agnes  #######
jpeg("./Grafics/PLSDA/HeatmapVIPs_MINT.jpeg", units="cm", width=14, height=18, res=300)
Heatmap(MatVIPz, name="Z-score", col=ColorHM, show_row_names=FALSE, show_column_names=FALSE,
        top_annotation = column_ha, width = unit(10, "cm"), height = unit(15, "cm"),
        show_row_dend = FALSE, cluster_rows = diana(MatVIPz),cluster_columns = agnes(t(MatVIPz)))
dev.off()





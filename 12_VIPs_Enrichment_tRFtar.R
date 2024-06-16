library(ggplot2)



##########################
#### DATOS DE MINTmap ####
##########################

VIPsmint_GO.BP<-read.csv("./Datos/Enrichment/MINTb_VIPs_GO_BP.csv")
VIPsmint_KEGG<-read.csv("./Datos/Enrichment/MINTb_VIPs_GO_KEGG.csv")

str(VIPsmint_GO.BP)
colnames(VIPsmint_GO.BP)[colnames(VIPsmint_GO.BP) == "Term.ID"] <- "ID"
colnames(VIPsmint_GO.BP)[1] <-"ID"
colnames(VIPsmint_GO.BP)[2] <-"Description"

##### Biological process ######
#Gene ratios
VIPsmint_GO.BP$Gene.List <- strsplit(as.character(VIPsmint_GO.BP$Gene.List), ", ")
# Unlist + unique
all_genesMnt <- unique(unlist(VIPsmint_GO.BP$Gene.List))
length(all_genesMnt)

totalGenesMnt<-length(all_genesMnt)

VIPsmint_GO.BP$Counts<-sapply(VIPsmint_GO.BP$Gene.List, length)
VIPsmint_GO.BP$GeneRatio<-sapply(VIPsmint_GO.BP$Gene.List, length)/totalGenesMnt

top15Mnt_GO.BP <- VIPsmint_GO.BP %>% arrange(Q.value) %>% head(15)


jpeg("./Grafics/Enrichment/MINT_GO_BP.jpeg", units="cm", width=23, height=12, res=300)
ggplot(top15Mnt_GO.BP, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Counts, color = Q.value)) +
  geom_point() + labs(title = "Processos biològics enriquits",x ="Proporció de gens", y = "")+
  scale_color_continuous(low = "red", high = "blue", name = "p ajustada") +
  theme(axis.text.y = element_text(size = 8), legend.position = "right")+
  theme_bw()
dev.off()

##### KEGG ######

str(VIPsmint_KEGG)
colnames(VIPsmint_KEGG)[colnames(VIPs_KEGG) == "Term.ID"] <- "ID"
colnames(VIPsmint_KEGG)[1] <-"ID"
colnames(VIPsmint_KEGG)[2] <-"Pathway"

# Gene ratios
VIPsmint_KEGG$Gene.List <- strsplit(as.character(VIPsmint_KEGG$Gene.List), ", ")
# Unlist + unique
Mnt_allgenesKg <- unique(unlist(VIPsmint_KEGG$Gene.List))


MntGenesKg<-length(Mnt_allgenesKg)
VIPsmint_KEGG$Counts<-sapply(VIPsmint_KEGG$Gene.List, length)
VIPsmint_KEGG$GeneRatio<-sapply(VIPsmint_KEGG$Gene.List, length)/MntGenesKg

length(which(VIPsmint_KEGG$Q.value <0.05))


Mntop15_Kegg <- VIPsmint_KEGG %>% arrange(Q.value) %>% head(15)

jpeg("./Grafics/Enrichment/MINT_GO_KEGG.jpeg", units="cm", width=20, height=12, res=300)
ggplot(Mntop15_Kegg, aes(x = reorder(Pathway,  Counts), y = Counts, fill = Q.value)) +
  geom_bar(stat = "identity") + labs(title = "Vies KEGG enriquides", x = "",y = "Nombre de gens")+
  scale_fill_gradient(low = "red", high = "blue", name = "p ajustada") +
  coord_flip() + theme_bw()
dev.off()  


#### MINTbase in STRING ####

##BIOL PROCESS NODE 1 ##

String1.BP<-read.table("./Datos/Enrichment/enrichProcess_MainNode95.txt", sep ="\t", header=FALSE)
colnames(String1.BP)<-c("ID","Description","Counts","BackgroundCount","Strength","Q.value","Prot_ID","Prot_Name")
str(String1.BP)


String_1top15 <- String1.BP %>% arrange(Q.value) %>% head(15)

jpeg("./Grafics/Enrichment/String1_GO_BP.jpeg", units="cm", width=20, height=12, res=300)
ggplot(String_1top15, aes(x = reorder(Description,  Counts), y = Counts, fill = Q.value)) +
  geom_bar(stat = "identity") + labs(title = "Processos biològics enriquits - node 1", x = "",y = "Nombre de proteïnes")+
  scale_fill_gradient(low = "red", high = "blue", name = "p ajustada") +
  coord_flip() + theme_bw()
dev.off()

## Mol ecular Function  NODE 1 ##

String1.MF<-read.table("./Datos/Enrichment/enrichFunction_MainNode95.txt", sep ="\t", header=FALSE)
colnames(String1.MF)<-c("ID","Function","Counts","BackgroundCount","Strength","Q.value","Prot_ID","Prot_Name")
str(String1.MF)

String_1topMF <- String1.MF %>% arrange(Q.value) %>% head(10)

jpeg("./Grafics/Enrichment/String1_GO_MF.jpeg", units="cm", width=20, height=12, res=300)
ggplot(String_1topMF, aes(x = reorder(Function,  Counts), y = Counts, fill = Q.value)) +
  geom_bar(stat = "identity") + labs(title = "Funció molecular - node 1", x = "",y = "Nombre de proteïnes")+
  scale_fill_gradient(low = "red", high = "blue", name = "p ajustada") +
  coord_flip() + theme_bw() + theme(aspect.ratio = 0.65)
dev.off()


## BIOLOGIC process NODE 2 ##

String2.BP<-read.table("./Datos/Enrichment/enrichProcess_SecondNode95.txt", sep ="\t", header=FALSE)
colnames(String2.BP)<-c("ID","Description","Counts","BackgroundCount","Strength","Q.value","Prot_ID","Prot_Name")
str(String2.BP)

String_2top <- String2.BP %>% arrange(Q.value) %>% head(15)

jpeg("./Grafics/Enrichment/String2_GO_BP.jpeg", units="cm", width=20, height=12, res=300)
ggplot(String_2top, aes(x = reorder(Description,  Counts), y = Counts, fill = Q.value)) +
  geom_bar(stat = "identity") + labs(title = "Processos biològics enriquits - node 2", x = "",y = "Nombre de proteïnes")+
  scale_fill_gradient(low = "red", high = "blue", name = "p ajustada") +
  coord_flip() + theme_bw() + theme(aspect.ratio = 0.75)
dev.off()

## molecular function NODE 2 ##

String2.MF<-read.table("./Datos/Enrichment/enrichFunction_SecondNode95.txt", sep ="\t", header=FALSE)
colnames(String2.MF)<-c("ID","Function","Counts","BackgroundCount","Strength","Q.value","Prot_ID","Prot_Name")
str(String2.MF)

String_2topMF <- String2.MF %>% arrange(Q.value) %>% head(15)

jpeg("./Grafics/Enrichment/String2_GO_MF.jpeg", units="cm", width=20, height=12, res=300)
ggplot(String_2topMF, aes(x = reorder(Function,  Counts), y = Counts, fill = Q.value)) +
  geom_bar(stat = "identity") + labs(title = "Funció molecular - node 2", x = "",y = "Nombre de proteïnes")+
  scale_fill_gradient(low = "red", high = "blue", name = "p ajustada") +
  coord_flip() + theme_bw()+ theme(aspect.ratio = 0.6)
dev.off()


## BIOL process NODE 3 ##

String3.BP<-read.table("./Datos/Enrichment/enrichProcess_ThirdNode95.txt", sep ="\t", header=FALSE)
colnames(String3.BP)<-c("ID","Description","Counts","BackgroundCount","Strength","Q.value","Prot_ID","Prot_Name")
str(String3.BP)

String3top <- String3.BP %>% arrange(Q.value) %>% head(15)

jpeg("./Grafics/Enrichment/String3_GO_BP.jpeg", units="cm", width=20, height=12, res=300)
ggplot(String3top, aes(x = reorder(Description,  Counts), y = Counts, fill = Q.value)) +
  geom_bar(stat = "identity") + labs(title = "Processos biològics enriquits - node 3", x = "",y = "Nombre de proteïnes")+
  scale_fill_gradient(low = "red", high = "blue", name = "p ajustada") +
  coord_flip() + theme_bw()+ theme(aspect.ratio = 0.65)
dev.off()

## Molecular function NODE 3 ##

String3.MF<-read.table("./Datos/Enrichment/enrichFunction_ThirdNode95.txt", sep ="\t", header=FALSE)
colnames(String3.MF)<-c("ID","Function","Counts","BackgroundCount","Strength","Q.value","Prot_ID","Prot_Name")
str(String3.MF)

String3topMF <- String3.MF %>% arrange(Q.value) %>% head(15)

jpeg("./Grafics/Enrichment/String3_GO_MF.jpeg", units="cm", width=20, height=12, res=300)
ggplot(String3topMF, aes(x = reorder(Function,  Counts), y = Counts, fill = Q.value)) +
  geom_bar(stat = "identity") + labs(title = "Funció molecular - node 3", x = "",y = "Nombre de proteïnes")+
  scale_fill_gradient(low = "red", high = "blue", name = "p ajustada") +
  coord_flip() + theme_bw()+ theme(aspect.ratio = 0.1) 
dev.off()




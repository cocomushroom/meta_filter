#### This file create FPKM table for nutrient-related genes. 
### Input file is generated using BOWTIE2 tabular file with prior text editing.
### The goal is to identify top ranking expressed genes in different layers. We plan to create a large file with manually curated annotation
## for all possible interesting genes. Using this mega annotation, different gene sets will be investigated.
## ADDITIONAL WORK NEED TESTING : Normalize data using original sequencing depth (probably more suitable for heatmap drawing). 
##This is the total reads: read_ori<-c(28125340,43088278,34631176,41315262,61686300,32507584,52059190,42425280,47060558)

### Add genes detected by ordinary DEG test (DEseq2) and compare if the results differ.
### At the end of the file there is DEseq2 analysis



#unwant_blast<- read.delim("/Users/kc178/Documents/nutrient/check_fungi_uni/unwant_15")
#unwant_trino <- read.delim("/Users/kc178/Documents/nutrient/clean_up/BOWTIE2_filter/trinotate_based_filter_Trinityname_0115")
#unwant_blast_unlist <- as.vector(unwant_blast)
#unwant_trino_unlist <- as.vector(unwant_trino)
#combine_unwant<-c(unwant_blast_unlist,unwant_trino_unlist)
#unwant<-unique(combine_unwant)
#unwant <- as.vector(unwant)


filter<-read.delim("/Users/kc178/Documents/nutrient/clean_up/BOWTIE2_filter/combine_filter_0115",header=F)
unwant<-unique(filter)

nrow(unwant)

#test <- c(unwant_blast,unwant_trino)





### Filter out unwanted genes##
#setwd("/Users/kc178/Documents/nutrient/check_fungi_uni")
#phy_e6<-read.delim("/Users/kc178/Documents/nutrient/check_fungi_uni/filtered.name", header=FALSE)
#all_name <- read.delim("/Users/kc178/Documents/nutrient/check_fungi_uni/blast_format6",header=FALSE)
del<-table[,1] %in% unwant[,1]
non_del<-table[!del,]
nrow(non_del)
### double check the removed number ##
all<-dim(table)
removed<-dim(unwant)
keep<-dim(non_del)
supposed_keep<-all-removed ## The number doesn't exactly match, the "unwant" list might be generated prior to this batch of sequences


## output transmembrane annotation ##
del_sq <-sq[,2] %in% unwant[,1]
sq_clean <- sq[!del_sq,]
sq_transport<-sq_clean[,2] %in% ANPC_transport_name
write.table(sq_clean[sq_transport,], row.names=FALSE, col.names=FALSE, quote=FALSE, "sq_transport_ACPN.txt")

nrow(sq_clean[sq_transport,2])

transport_row<-non_del[,1] %in% ANPC_transport_name
transporter_name <- non_del[transport_row,]


######## Start here ####
setwd("/Users/kc178/Documents/nutrient/clean_up/BOWTIE2_filter/")
table<-read.delim("/Users/kc178/Documents/nutrient/clean_up/BOWTIE2_filter//TMBTMBTMB_length",header=TRUE)
sq<-read.delim("/Users/kc178/Documents/nutrient/clean_up/BOWTIE2_filter//report.sqlite",header=TRUE)

#### Append nutrient type information (P, N, C, A)

A_name <- read.delim("/Users/kc178/Documents/nutrient/clean_up/BOWTIE2_filter/transport_name/AA_name",header=F)
N_name <-read.delim("/Users/kc178/Documents/nutrient/clean_up/BOWTIE2_filter/transport_name/nitrogen_name",header=F)
P_name <- read.delim("/Users/kc178/Documents/nutrient/clean_up/BOWTIE2_filter/transport_name/pho_name",header=F)
C_name <- read.delim("/Users/kc178/Documents/nutrient/clean_up/BOWTIE2_filter/transport_name/sugar_name",header=F)


A_name<-unique(A_name)
N_name<-unique(N_name)
P_name<-unique(P_name)
C_name<-unique(C_name)

for (i in 1:nrow(A_name)) {
  A_name[i,2] <- "AA"
}
for (i in 1:nrow(N_name)) {
  N_name[i,2] <- "N"
}
for (i in 1:nrow(P_name)) {
  P_name[i,2] <- "P"
}
for (i in 1:nrow(C_name)) {
  C_name[i,2] <- "C"
}


ANPC_transport_name_cat <-rbind(A_name,N_name,P_name,C_name)
colnames(ANPC_transport_name_cat)[1]<-"ID"
colnames(ANPC_transport_name_cat)[2]<-"cat"

ANPC_transport_name<-unique(c(as.vector(A_name[,1]),as.vector(P_name[,1]), as.vector(N_name[,1]), as.vector(C_name[,1])))

ANPC_bowtie_row <- table[,1] %in% ANPC_transport_name
ANPC_bowtie<-table[ANPC_bowtie_row,]


## To look at genes related to transport (transmembrane transport, transmembrane transporter)
#transport <- read.delim("/Users/kc178/Documents/nutrient/clean_up/BOWTIE2_filter/trinotate_based_transport", header=F)

#ANPC_transport_name<-c(as.vector(A_name[,1]),as.vector(P_name[,1]), as.vector(N_name[,1]), as.vector(C_name[,1]))





## Just for double check if something is duplicated
##dup<-duplicated(ANPC_transport_name)
##ANPC_transport_name[dup]

### prepare table format##
table_num<-ANPC_bowtie[2:11] ## this variable is used for the DEseq2 too
row.names(table_num)<-ANPC_bowtie[,1]
head(table_num,n=10)
nrow(table_num)



##FPKM count##Please note it's better to use the word "FPKM" here since these are paired end reads
## Normalized by total read counts
name<-list()
for (p in 1:9) {
  name[[p]]<-paste(colnames(table_num[p]),"RPKM", sep="_")
  table_num[[name[[p]]]]<-(table_num[,p]*10^9)/sum(table_num[,p])/table_num$length
}

dim(table_num)
#################
####Order combine layers####


OTU_raw <- table_num
OTU<-OTU_raw[,11:19]

colIndex<- vector(mode="list")
colIndex$b<- grep("b.", names(OTU))
colIndex$m<- grep("m.", names(OTU))
colIndex$t<- grep("t.", names(OTU))


OTU3<-t( apply(OTU, 1, function(x){
  c(mean(x[colIndex$b]/3),
    mean(x[colIndex$m]/3),
    mean(x[colIndex$t]/3) )
})  )
colnames(OTU3)<- c("B","M","T")


order_t <- order(OTU3[,1], decreasing=T)
order_m<-order(OTU3[,2], decreasing=T)
order_b<-order(OTU3[,3], decreasing=T)

ID_t100 <- row.names(OTU3[(order_t[1:100]),])
ID_m100 <- row.names(OTU3[(order_m[1:100]),])
ID_b100 <- row.names(OTU3[(order_b[1:100]),])


#####find annotation for the transcripts with high ranking in combined set####

t100_cat<-ANPC_transport_name_cat[ANPC_transport_name_cat$ID %in% ID_t100,2]
m100_cat<-ANPC_transport_name_cat[ANPC_transport_name_cat$ID %in% ID_m100,2]
b100_cat<-ANPC_transport_name_cat[ANPC_transport_name_cat$ID %in% ID_b100,2]

#tmb_bind<-cbind(t100_cat,m100_cat, b100_cat)
#table(tmb_bind)

table(t100_cat)
table(m100_cat)
table(b100_cat)

#anno_b100<-sq[,2] %in% ID_b100

#t100_function <-sq[anno_t100,]
#m100_function <-sq[anno_m100,]
#b100_function <-sq[anno_b100,]

#write.table(t100_function, row.names=FALSE, col.names=FALSE, quote=FALSE, "t100_unwant_remove.txt")
#write.table(m100_function, row.names=FALSE, col.names=FALSE, quote=FALSE, "m100_unwant_remove.txt")
#write.table(b100_function, row.names=FALSE, col.names=FALSE, quote=FALSE, "b100_unwant_remove.txt")
##### 
########



###################

########Order by per sample#########
t1_order<- order(table_num$t1_RPKM, decreasing=TRUE)
t1_100_ID<-row.names(table_num[t1_order[1:100],])

t2_order<- order(table_num$t2_RPKM, decreasing=TRUE)
t2_100_ID<-row.names(table_num[t2_order[1:100],])

t3_order<- order(table_num$t3_RPKM, decreasing=TRUE)
t3_100_ID<-row.names(table_num[t3_order[1:100],])

#t12<-intersect(t1_order[1:100],t2_order[1:100])
#t123<-intersect(t12,t3_order[1:100])

#t12_u<-union(t1_order[1:100],t2_order[1:100])
#t123_u<-union(t12_u,t3_order[1:100])

#write.table(t1_100_ID, row.names=FALSE, col.names=FALSE, quote=FALSE, "t1.txt")

#ID_t123_union<-row.names(table_num[t123_u,])
#write.table(ID_t123_union, row.names = FALSE, col.names=FALSE, quote=FALSE, "t123_union.txt")


m1_order<- order(table_num$m1_RPKM, decreasing=TRUE)
m1_100_ID<-row.names(table_num[m1_order[1:100],])

m2_order<- order(table_num$m2_RPKM, decreasing=TRUE)
m2_100_ID<-row.names(table_num[m2_order[1:100],])

m3_order<- order(table_num$m3_RPKM, decreasing=TRUE)
m3_100_ID<-row.names(table_num[m3_order[1:100],])



b1_order<- order(table_num$b1_RPKM, decreasing=TRUE)
b1_100_ID<-row.names(table_num[b1_order[1:100],])

b2_order<- order(table_num$b2_RPKM, decreasing=TRUE)
b2_100_ID<-row.names(table_num[b2_order[1:100],])

b3_order<- order(table_num$b3_RPKM, decreasing=TRUE)
b3_100_ID<-row.names(table_num[b3_order[1:100],])

#####  Summarize SeqID names for manual annotation####

tmb_123tmb_100_ID <- c(ID_t100,ID_m100,ID_b100,t1_100_ID,t2_100_ID, t3_100_ID, m1_100_ID, m2_100_ID, m3_100_ID, b1_100_ID, b2_100_ID, b3_100_ID)
length(tmb_123tmb_100_ID)
uniq_tmb_123tmb_100_ID<-unique(tmb_123tmb_100_ID)
length(uniq_tmb_123tmb_100_ID)

anno_uniq_tmb_123tmb_100_ID<-sq[,2] %in% uniq_tmb_123tmb_100_ID

uniq_tmb_123tmb_100_ID_function <-sq[anno_uniq_tmb_123tmb_100_ID,]

write.table(uniq_tmb_123tmb_100_ID_function, row.names=FALSE, col.names=FALSE, quote=FALSE, "uniq_tmb_123tmb_100_ID_function_transport0111.txt")




#######################
####### DEseq2 ########
#######################

library("DESeq2", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/resources/library")
library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
library("gdata", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")


countC <- table_num[,1:9]
type <- read.delim("condition_type", header=TRUE)
a <- type[,2:3]
row.names(a) <- as.character(type[,1])
as.character(type[,1])
row.names(a) <- as.character(type[,1])

print(type)
dim(type)
ddsC <- DESeqDataSetFromMatrix(countData = countC,
                               colData = a,
                               design = ~ condition)
ddsC <- ddsC[, c(1,4,7,2,5,8,3,6,9)]
ddsC <- DESeq(ddsC)

C_resTB<- results(ddsC)
C_resTB
summary(C_resTB)
#change which condition to compare with
ddsC$condition <- relevel(ddsC$condition, "middle")
# C_resTM gives the C_results of top vs. middle group
C_resTM <- results(ddsC)
C_resTM
summary(C_resTM)
# C_resMB gives the C_results of middle vs. bottom group
# sometimes run "top" first and change it to bottom will work 
ddsC$condition <- relevel(ddsC$condition, "top")
ddsC$condition <- relevel(ddsC$condition, "bottom")
C_resMB <- results(ddsC)
C_resMB
summary(C_resMB)
ddsC$condition <- relevel(ddsC$condition, "top")

library("RColorBrewer")
library("gplots")
hmcol <- colorRampPalette(brewer.pal(9, "YlOrBr"))(100)
blue_yellow <- colorRampPalette(brewer.pal(9,"YlGnBu"))(100)

plotMA(C_resTB, ylim=c(-15,15), main="DESeq2")
rld <- rlog(ddsC)
vsd <- varianceStabilizingTransformation(ddsC)
rlogMat <- assay(rld)
vstMat <- assay(vsd)

#clustering heatmap
distsRL <- dist(t(assay(vsd)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(ddsC), paste(condition, type, sep=" : "))
hc <- hclust(distsRL)
heatmap.2(mat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col = rev(blue_yellow), margin=c(10, 10))


## 
C_resTBP <- C_resTB$padj
C_resTBP_T <- C_resTBP < 0.05
#c: "which" can give the index of whatever are true
#list of genes with p< certain thC_reshold  in Top vs. Bottom comparison
C_resTBP_T_index  <- which(C_resTBP_T)
length(C_resTBP_T_index)

C_resTMP <- C_resTM$padj
C_resTMP_T <- C_resTMP < 0.05
#c: "which" can give the index of whatever are true
#list of genes with p< certain thC_reshold  in Top vs. Bottom comparison
C_resTMP_T_index  <- which(C_resTMP_T)
length(C_resTMP_T_index)

C_resMBP <- C_resMB$padj
C_resMBP_T <- C_resMBP < 0.05
#c: "which" can give the index of whatever are true
#list of genes with p< certain thC_reshold  in Top vs. Bottom comparison
C_resMBP_T_index  <- which(C_resMBP_T)
length(C_resMBP_T_index)


#common listc
common <- abs(C_resTB$log2FoldChange) < log2(2)
commonnames <- C_resTB@rownames[common]

commonP <- C_resTB$padj > 0.05
commonP[is.na(commonP)] <- F
commonPnames <- C_resTB@rownames[commonP]

common2P <- common & commonP
common2Pnames <- C_resTB@rownames[common2P]

##check up- down- regulated ones

sigTB <- C_resTB$padj < 0.05
TBup <- C_resTB$log2FoldChange > log2(4) ## fold change >1=> logfoldchange >0
sigTBup <- sigTB & TBup
sigTBup[is.na(sigTBup)] <- F## this truns NA to false so only the true ones will be printed out

TBnamesup <- C_resTB@rownames[sigTBup]
TBup[is.na(TBup)]<-F
TBlogup<-as.vector(C_resTB@rownames[TBup])
write.csv(TBlogup,row.names=FALSE, "TBlogup.csv")

gup <- genus[TBnamesup]
write.csv(TBnamesup, row.names = FALSE, "TBnamesup.csv")
length(TBnamesup)

TBdown <- C_resTB$log2FoldChange < -log2(4)
sigTBdown <- sigTB & TBdown
sigTBdown[is.na(sigTBdown)] <- F 
TBdown[is.na(TBdown)]<-F
TBlogdown<-as.vector(C_resTB@rownames[TBdown])
TBnamesdown <- C_resTB@rownames[sigTBdown]
length(TBnamesdown)

sigTBup_cat<-ANPC_transport_name_cat$ID %in% TBnamesup
sigTBdown_cat<-ANPC_transport_name_cat$ID %in% TBnamesdown

heatmap.2(assay(vsd)[sigTBup,], labRow=ANPC_transport_name_cat$cat[sigTBup_cat], main = c("top layer", "up regulated genes"), col = blue_yellow, Rowv = TRUE, Colv = FALSE, scale="none", 
          dendrogram="row", trace="none", margin=c(3, 8), cexCol = 1, cexRow=0.9, srtCol=70, keysize=1.5, density.info="none")

heatmap.2(assay(vsd)[sigTBdown,], labRow=ANPC_transport_name_cat$cat[sigTBdown_cat], main = c("bottom layer", "up regulated genes"), col = blue_yellow, Rowv = TRUE, Colv = FALSE, scale="none", 
          dendrogram="row", trace="none", margin=c(3, 8), cexCol = 1, cexRow=0.9, srtCol=70, keysize=1.5, density.info="none")

table(ANPC_transport_name_cat$cat[sigTBup_cat])
table(ANPC_transport_name_cat$cat[sigTBdown_cat])

#test
png("LR3_MA_2P.png", width=1200, height=1000)
par(font.axis=3)#make label italic#
par(cex.main=3) #change main title size
heatmap.2(assay(vsd)[common2P,], labRow=genus[common2Pnames], 
          main = c("taxa equally presented: LR3-MA"), 
          col = blue_yellow, Rowv = TRUE, Colv = FALSE, scale="none", 
          dendrogram="row", trace="none", margin=c(6, 24.2), 
          cexCol = 3.2, cexRow=2.9, srtCol=70, keysize=0.75, 
          density.info="none", key.par=list(cex.lab=2.2, cex.axis=2.2, cex.main=2.5),
          RowSideColors=unlist(pcol3[common2Pnames],
                               
          )  )
dev.off()







#PCA plot

dataPCA <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=TRUE)
data<-dataPCA
percentVar <- round(100 * attr(data, "percentVar"))
data$condition<- factor(data$condition, levels=c("top",    "middle", "bottom") )

ggplot(data, aes(PC1, PC2, color=condition, shape=type)) +
  #scale_fill_manual(values=c("greenyellow","chocolate4","darkolivegreen")) +
  #scale_fill_manual(values=plotColor) +
  #scale_color_hue(l=40, c=35) +
  geom_point(colour="black",size=4.5)+
  scale_color_manual(values=col2hex(c("greenyellow","darkolivegreen","chocolate4"))) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  ggtitle("nutri_variance")




### Failed , work on this later###

##### Obtain the DEseq2-normalized counts (note some debates, also these probably didn't get normalized by bases)####
#####
dds <- estimateSizeFactors(ddsC)
DE_count<-counts(dds, normalized=TRUE)
head(DE_count,n=10)
unlist(DE_count)
DE_count$t1
attributes(DE_count)$b1
t1_DE_order<- order(DE_count$t1, decreasing=TRUE)
t1_150_DE_ID<-row.names(DE_count[t1_DE_order[1:150],])

t2_DE_order<- order(DE_count$t2, decreasing=TRUE)
t2_150_DE_ID<-row.names(DE_count[t2_DE_order[1:150],])

t3_DE_order<- order(DE_count$t3, decreasing=TRUE)
t3_150_DE_ID<-row.names(DE_count[t3_DE_order[1:150],])

m1_DE_order<- order(DE_count$m1, decreasing=TRUE)
m1_150_DE_ID<-row.names(DE_count[m1_DE_order[1:150],])

m2_DE_order<- order(DE_count$m2, decreasing=TRUE)
m2_150_DE_ID<-row.names(DE_count[m2_DE_order[1:150],])

m3_DE_order<- order(DE_count$m3, decreasing=TRUE)
m3_150_DE_ID<-row.names(DE_count[m3_DE_order[1:150],])


b1_DE_order<- order(DE_count$b1, decreasing=TRUE)
b1_150_DE_ID<-row.names(DE_count[b1_DE_order[1:150],])

b2_DE_order<- order(DE_count$b2, decreasing=TRUE)
b2_150_DE_ID<-row.names(DE_count[b2_DE_order[1:150],])

b3_DE_order<- order(DE_count$b3, decreasing=TRUE)
b3_150_DE_ID<-row.names(DE_count[b3_DE_order[1:150],])


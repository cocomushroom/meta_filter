#### This file create FPKM table for nutrient-related genes. 
### Input file is generated using BOWTIE2 tabular file with prior text editing.
### The goal is to identify top ranking expressed genes in different layers. We plan to create a large file with manually curated annotation
## for all possible interesting genes. Using this mega annotation, different gene sets will be investigated.
## ADDITIONAL WORK NEED TESTING : Normalize data using original sequencing depth (probably more suitable for heatmap drawing). 
##This is the total reads: read_ori<-c(28125340,43088278,34631176,41315262,61686300,32507584,52059190,42425280,47060558)

### Add genes detected by ordinary DEG test (DEseq2) and compare if the results differ.
### At the end of the file there is DEseq2 analysis


setwd("/Users/kc178/Documents/nutrient/clean_up/BOWTIE2_filter/")
table<-read.delim("/Users/kc178/Documents/nutrient/clean_up/BOWTIE2_filter//TMBTMBTMB_length",header=TRUE)

#unwant_blast<- read.delim("/Users/kc178/Documents/nutrient/check_fungi_uni/unwant_15")
#unwant_trino <- read.delim("/Users/kc178/Documents/nutrient/clean_up/BOWTIE2_filter/trinotate_based_filter_Trinityname")

sq<-read.delim("/Users/kc178/Documents/nutrient/clean_up/BOWTIE2_filter//report.sqlite",header=TRUE)

filter<-read.delim("/Users/kc178/Documents/nutrient/clean_up/BOWTIE2_filter/combine_filter_0111",header=F)
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


## To look at genes related to transport (transmembrane transport, transmembrane transporter)
transport <- read.delim("/Users/kc178/Documents/nutrient/clean_up/BOWTIE2_filter/trinotate_based_transmembrane", header=F)
transport_row<-non_del[,1] %in% transport[,1]
transporter_name <- non_del[transport_row,]



### prepare table format##
table_num<-transport_name[2:11] ## this variable is used for the DEseq2 too
row.names(table_num)<-non_del[,1]
head(table_num,n=10)
nrow(table_num)



##FPKM count##Please note it's better to use the word "FPKM" here since these are paired end reads
## Normalized by total read counts
name<-list()
for (p in 1:9) {
  name[[p]]<-paste(colnames(table_num[p]),"RPKM", sep="_")
  table_num[[name[[p]]]]<-(table_num[,p]*10^9)/sum(table_num[,p])/table_num$length
}

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

ID_t150 <- row.names(OTU3[(order_t[1:150]),])
ID_m150 <- row.names(OTU3[(order_m[1:150]),])
ID_b150 <- row.names(OTU3[(order_b[1:150]),])


#####find annotation for the transcripts with high ranking in combined set####

anno_t150<-sq[,2] %in% ID_t150
anno_m150<-sq[,2] %in% ID_m150
anno_b150<-sq[,2] %in% ID_b150

t150_function <-sq[anno_t150,]
m150_function <-sq[anno_m150,]
b150_function <-sq[anno_b150,]

write.table(t150_function, row.names=FALSE, col.names=FALSE, quote=FALSE, "t150_unwant_remove.txt")
write.table(m150_function, row.names=FALSE, col.names=FALSE, quote=FALSE, "m150_unwant_remove.txt")
write.table(b150_function, row.names=FALSE, col.names=FALSE, quote=FALSE, "b150_unwant_remove.txt")
##### 

########



###################

########Order by per sample#########
t1_order<- order(table_num$t1_RPKM, decreasing=TRUE)
t1_150_ID<-row.names(table_num[t1_order[1:150],])

t2_order<- order(table_num$t2_RPKM, decreasing=TRUE)
t2_150_ID<-row.names(table_num[t2_order[1:150],])

t3_order<- order(table_num$t3_RPKM, decreasing=TRUE)
t3_150_ID<-row.names(table_num[t3_order[1:150],])

#t12<-intersect(t1_order[1:150],t2_order[1:150])
#t123<-intersect(t12,t3_order[1:150])

#t12_u<-union(t1_order[1:150],t2_order[1:150])
#t123_u<-union(t12_u,t3_order[1:150])

#write.table(t1_150_ID, row.names=FALSE, col.names=FALSE, quote=FALSE, "t1.txt")

#ID_t123_union<-row.names(table_num[t123_u,])
#write.table(ID_t123_union, row.names = FALSE, col.names=FALSE, quote=FALSE, "t123_union.txt")


m1_order<- order(table_num$m1_RPKM, decreasing=TRUE)
m1_150_ID<-row.names(table_num[m1_order[1:150],])

m2_order<- order(table_num$m2_RPKM, decreasing=TRUE)
m2_150_ID<-row.names(table_num[m2_order[1:150],])

m3_order<- order(table_num$m3_RPKM, decreasing=TRUE)
m3_150_ID<-row.names(table_num[m3_order[1:150],])



b1_order<- order(table_num$b1_RPKM, decreasing=TRUE)
b1_150_ID<-row.names(table_num[b1_order[1:150],])

b2_order<- order(table_num$b2_RPKM, decreasing=TRUE)
b2_150_ID<-row.names(table_num[b2_order[1:150],])

b3_order<- order(table_num$b3_RPKM, decreasing=TRUE)
b3_150_ID<-row.names(table_num[b3_order[1:150],])

#####  Summarize SeqID names for manual annotation####

tmb_123tmb_150_ID <- c(ID_t150,ID_m150,ID_b150,t1_150_ID,t2_150_ID, t3_150_ID, m1_150_ID, m2_150_ID, m3_150_ID, b1_150_ID, b2_150_ID, b3_150_ID)
length(tmb_123tmb_150_ID)
uniq_tmb_123tmb_150_ID<-unique(tmb_123tmb_150_ID)
length(uniq_tmb_123tmb_150_ID)

anno_uniq_tmb_123tmb_150_ID<-sq[,2] %in% uniq_tmb_123tmb_150_ID

uniq_tmb_123tmb_150_ID_function <-sq[anno_uniq_tmb_123tmb_150_ID,]

write.table(uniq_tmb_123tmb_150_ID_function, row.names=FALSE, col.names=FALSE, quote=FALSE, "uniq_tmb_123tmb_150_ID_function_transport0111.txt")



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


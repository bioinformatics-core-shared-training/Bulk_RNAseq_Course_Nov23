#!/usr/bin/Rscript
#Description: Day3 for bulk RNA-seq course: annotation and visualisation
#load libraries 
library(AnnotationHub)
library(AnnotationDbi)
library(ensembldb)
library(DESeq2)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(ComplexHeatmap)
library(circlize)

# Annotation ----
## load data ----
ddsObj.interaction <- readRDS("RObjects/DESeqDataSet.interaction.rds")
results.interaction.11 <- readRDS("RObjects/DESeqResults.interaction_d11.rds")
results.interaction.33 <- readRDS("RObjects/DESeqResults.interaction_d33.rds")

### Annotationhub
ah<-AnnotationHub()
ah

## Select the data that we need
MouseEnsDb<-query(ah,c("EnsDb", "Mus musculu","102"))
length(MouseEnsDb)
#MouseEnsDb<-query(ah,c("EnsDb", "Mus musculu"))[1]
MouseEnsDb<-query(ah,c("EnsDb", "Mus musculu","102"))[[1]]

## Extract genes
annotations<-genes(MouseEnsDb, return.type='data.frame')
colnames(annotations)
View(annotations)
annot<-annotations %>%
  dplyr::select(gene_id,gene_name,entrezid) %>%
  dplyr::filter(gene_id %in% rownames(results.interaction.11))
head(annot)
length(annot$entrezid)
length(unique(annot$entrezid))
sum(is.na(annot$entrezid))

## query to an already made annotation file 
ensemblAnnot<-readRDS("RObjects/Ensembl_annotations.rds")
head(ensemblAnnot)

## Annotate results ----
annot.interaction.11<-as.data.frame(results.interaction.11) %>%
  rownames_to_column("GeneID") %>%
  left_join(ensemblAnnot,"GeneID") %>%
  dplyr::rename(logFC=log2FoldChange, FDR=padj)
annot.interaction.11

#if want to check if the gene of interest is there
annot.interaction.11[annot.interaction.11$Symbol=='Sox17',]

#Save the results
write_tsv(annot.interaction.11, "results/Interaction.11_Results_Annotated.txt")

# Visulization ----
## sanity check 
hist(annot.interaction.11$pvalue)

## shrink values 
ddsShrink.11<-lfcShrink(ddsObj.interaction,res=results.interaction.11, type="ashr")

shrinkTab.11<-as.data.frame(ddsShrink.11) %>%
  rownames_to_column("GeneID") %>%
  left_join(ensemblAnnot, "GeneID") %>%
  dplyr::rename(logFC=log2FoldChange, FDR=padj)

head(shrinkTab.11)


## MA plot ----
par(mfrow=c(1,2))
plotMA(results.interaction.11, alpha=0.05, main="raw")
plotMA(ddsShrink.11, alpha=0.05, main="shrinked")

## Volcano plots 
volcano.Tab.11<-shrinkTab.11 %>%
  mutate(`-log10(pvalue)`=-log10(pvalue))

ggplot(volcano.Tab.11, aes(x=logFC, y=`-log10(pvalue)`)) +geom_point(aes(col=FDR<0.05),size=1) +
  geom_text(data=~top_n(.x, 1, wt=-FDR), aes(label=Symbol)) +
  labs(x="log2(fold change)",y ="-log10(p-value)", colour="FDR <5%", 
       tilte="Infected vs Uninfected (day 11)")

ggplot(volcano.Tab.11, aes(x=logFC, y=`-log10(pvalue)`)) +geom_point(aes(col=FDR<0.05),size=1) +
  geom_text(data=~top_n(.x, 4, wt=-FDR), aes(label=Symbol)) +
  labs(x="log2(fold change)",y ="-log10(p-value)", colour="FDR <5%", 
       tilte="Infected vs Uninfected (day 11)")

ggplot(volcano.Tab.11, aes(x=logFC, y=`-log10(pvalue)`)) +geom_point(aes(col=FDR<0.05),size=1) +
  geom_text_repel(data=~top_n(.x, 4, wt=-FDR), aes(label=Symbol)) +
  labs(x="log2(fold change)",y ="-log10(p-value)", colour="FDR <5%", 
       tilte="Infected vs Uninfected (day 11)")


# Exercise 1 ----
#Shrink the results for the day 33 contrast 
ddsShrink.33 <- lfcShrink(ddsObj.interaction, 
                          res = results.interaction.33,
                          type = "ashr")

volcano.Tab.33<-as.data.frame(ddsShrink.33) %>% 
  mutate(`-log10(pvalue)`=-log10(pvalue))
head(volcano.Tab.33)

shrinkTab.33 <- as.data.frame(volcano.Tab.33) %>%
  rownames_to_column("GeneID") %>% 
  left_join(ensemblAnnot, "GeneID") %>%
  dplyr::rename(logFC=log2FoldChange, FDR=padj)

head(shrinkTab.33)

#Create the vocalno plot 
ggplot(shrinkTab.33, aes(x=logFC, y=`-log10(pvalue)`)) + geom_point(aes(col = FDR < 0.05), size=1) +
  geom_text_repel(data = ~top_n(.x, 4, wt = -FDR), aes(label = Symbol)) 

#compare d11 vs d33
#par(mfrow=c(1,2))
vol1<-ggplot(volcano.Tab.11, aes(x=logFC, y=`-log10(pvalue)`)) + geom_point(aes(col = FDR < 0.05), size=1) +
  geom_text_repel(data = ~top_n(.x, 4, wt = -FDR), aes(label = Symbol)) + ggtitle('Day 11')
vol2<-ggplot(shrinkTab.33, aes(x=logFC, y=`-log10(pvalue)`)) + geom_point(aes(col = FDR < 0.05), size=1) +
  geom_text_repel(data = ~top_n(.x, 4, wt = -FDR), aes(label = Symbol)) + ggtitle('Day 33')
vol1+vol2


# Exercise 2 ----
plotMA(ddsShrink.11, alpha=0.05, main="shrinked") 
#transform the data first to have the same format
maTab.33<-shrinkTab.33 %>%
  mutate(M=log2(baseMean))
ggplot(maTab.33, aes(x=M, y=logFC)) + geom_point((aes(col=FDR<0.05)), size =1) +
  scale_y_continuous(limits=c(-4,4)) + theme_bw()
ggplot(maTab.33, aes(x=M, y=logFC)) + geom_point((aes(col=FDR<0.05)), size =1) +
  scale_y_continuous(limits=c(-4,4), oob=scales::squish) + theme_bw()

## Strip Charts ----
# Extract ensembl id using the gene name
geneID<-dplyr::filter(shrinkTab.11, Symbol=="Il10ra") %>% pull(GeneID)
plotCounts(ddsObj.interaction, 
           gene = geneID, 
           intgroup = c("TimePoint", "Status", "Replicate"),
           returnData = T) %>% 
  ggplot(aes(x=Status, y=log2(count))) +
  geom_point(aes(fill=Replicate), shape=21, size=2) +
  facet_wrap(~TimePoint) +
  expand_limits(y=0) +
  labs(title = "Normalised counts - Interleukin 10 receptor, alpha")

# Exercise 3 ----
## Stripcharts with another gene Jchain
geneID <- dplyr::filter(shrinkTab.11, Symbol=="Jchain") %>% pull(GeneID)
plotCounts(ddsObj.interaction, 
           gene = geneID, 
           intgroup = c("TimePoint", "Status", "Replicate"),
           returnData = T) %>% 
  ggplot(aes(x=Status, y=log2(count))) +
  geom_point(aes(fill=Replicate), shape=21, size=2) +
  facet_wrap(~TimePoint) +
  expand_limits(y=0) +
  labs(title = "Normalised counts for Jchain")



# Venn Diagram ----
vennDat <- tibble(Geneid=rownames(results.interaction.11)) %>% 
  mutate(Upregulated_11 = results.interaction.11$padj < 0.05 & 
           !is.na(results.interaction.11$padj) & 
           results.interaction.11$log2FoldChange > 0) %>% 
  mutate(Downregulated_11 = results.interaction.11$padj < 0.05 & 
           !is.na(results.interaction.11$padj) & 
           results.interaction.11$log2FoldChange < 0) %>%
  mutate(Upregulated_33 = results.interaction.33$padj < 0.05 & 
           !is.na(results.interaction.33$padj) & 
           results.interaction.33$log2FoldChange > 0) %>%
  mutate(Downregulated_33 = results.interaction.33$padj < 0.05 & 
           !is.na(results.interaction.33$padj) & 
           results.interaction.33$log2FoldChange < 0) 

head(vennDat)

## plot Venn diagram
ggvenn(vennDat,set_name_size=4)


# Heatmap ----
#filtering 
#get the top genes
sigGenes<- shrinkTab.11 %>%
  top_n(300, wt= -FDR) %>%
  pull("GeneID")
head(sigGenes)

plotDat<-vst(ddsObj.interaction)[sigGenes,] %>%
  assay()
head(plotDat)
Heatmap(plotDat)

#z-score transformation
z.mat<-t(scale(t(plotDat), center=TRUE, scale=TRUE))

#change colour with colour palette
myPalette <- c("royalblue3", "ivory", "orangered3")
myRamp <- colorRamp2(c(-2, 0, 2), myPalette)

Heatmap(z.mat, name="z-score",col=myRamp,
        show_row_names=FALSE)

#split the heatmap into clusters and add annotation 
ha1= HeatmapAnnotation(df=colData(ddsObj.interaction)[,c("Status","TimePoint")])
Heatmap(z.mat, name = "z-score",
        col = myRamp,
        show_row_names = FALSE,  top_annotation = ha1) 

Heatmap(z.mat, name = "z-score",
        col = myRamp,            
        show_row_name = FALSE,
        split=3,
        rect_gp = gpar(col = "lightgrey", lwd=0.3),
        top_annotation = ha1)

#how to change the colours of the bars at the top
ha1 = HeatmapAnnotation(df = colData(ddsObj.interaction)[,c("Status", "TimePoint")], 
                        col = list(Status = c("Uninfected" = "darkgreen", 
                                              "Infected" = "palegreen"), 
                                   TimePoint = c("d11" = "lightblue", 
                                                 "d33" = "darkblue")))

Heatmap(z.mat, name = "z-score",
        col = myRamp,            
        show_row_name = FALSE,
        split=3,
        rect_gp = gpar(col = "lightgrey", lwd=0.3),
        top_annotation = ha1)


saveRDS(annot.interaction.11, file="results/Annotated_Results.d11.rds")
saveRDS(shrinkTab.11, file="results/Shrunk_Results.d11.rds")
saveRDS(annot.interaction.33, file="results/Annotated_Results.d33.rds")
saveRDS(shrinkTab.33, file="results/Shrunk_Results.d33.rds")

# pre-processed course files
saveRDS(annot.interaction.11, file="RObjects/Annotated_Results.d11.rds")
saveRDS(shrinkTab.11, file="RObjects/Shrunk_Results.d11.rds")
# saveRDS(annot.interaction.33, file="RObjects/Annotated_Results.d33.rds")
saveRDS(shrinkTab.33, file="RObjects/Shrunk_Results.d33.rds")


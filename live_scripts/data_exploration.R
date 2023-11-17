library(tximport)
library(DESeq2)
library(tidyverse)

# read in some sample information

sampleinfo <- read_tsv("data/samplesheet.tsv", col_types = c("cccc"))
arrange(sampleinfo, Status, TimePoint, Replicate)

files <- file.path("salmon", sampleinfo$SampleName, "quant.sf")
files <- set_names(files, sampleinfo$SampleName)
tx2gene <- read_tsv("references/tx2gene.tsv")

# import salmon files
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
str(txi)
head(txi$counts)
saveRDS(txi, "salmon_outputs/txi.rds")

# Exercise 1

tpm <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
str(tpm)
head(tpm$counts)

# prepare raw counts

rawCounts <- round(txi$counts, 0)
dim(rawCounts)

keep <- rowSums(rawCounts) > 5
table(keep, useNA = "always")

filtCounts <- rawCounts[keep,]
dim(filtCounts)


# Start exploration

summary(filtCounts)

boxplot(filtCounts, main= "Raw Counts", las = 2)
plot(rowMeans(filtCounts), rowSds(filtCounts),
     main="Raw Counts sd vs mean",
     xlim=c(0,10000),
     ylim=c(0,5000))


logcounts <- log2(filtCounts + 1)

# make a colour vector
statusCols <- case_when(sampleinfo$Status=="Infected" ~ "red", 
                        sampleinfo$Status=="Uninfected" ~ "orange")

# Check distributions of samples using boxplots
boxplot(logcounts,
        xlab="",
        ylab="Log2(Counts)",
        las=2,
        col=statusCols,
        main="Log2(Counts)")
# Let's add a blue horizontal line that corresponds to the median
abline(h=median(logcounts), col="blue")

plot(rowMeans(logcounts), rowSds(logcounts),
     main="Log2 Counts sd vs mean")

# vst

vst_counts <- vst(filtCounts)
boxplot(vst_counts,
        xlab="",
        ylab="VSTcounts",
        las = 2,
        col = statusCols)
abline(h=median(vst_counts), col="blue")

plot(rowMeans(vst_counts), rowSds(vst_counts),
     main="VST Counts sd vs mean")

# PCA

library(ggfortify)

rlogcounts <- rlog(filtCounts)

pcDat <- prcomp(t(rlogcounts))
autoplot(pcDat)

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5)

# Exercise 3

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         x=2,
         y=3,
         size = 5)

# labels

library(ggrepel)
autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5) +
  geom_text_repel(aes(x=PC1, y=PC2, label=SampleName, box.padding=0.8))

sampleinfo <- mutate(sampleinfo, Status= case_when(
  SampleName=="SRR7657882" ~ "Uninfected",
  SampleName=="SRR7657873" ~ "Infected", 
  TRUE ~ Status))

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5) +
  geom_text_repel(aes(x=PC1, y=PC2, label=SampleName, box.padding=0.8))

# Hclust

library(ggdendro)
 
hclDat <- t(rlogcounts) %>%
  dist(method = "euclidean") %>%
  hclust()

ggdendrogram(hclDat, rotate = TRUE)  

hclDat2 <- hclDat
hclDat2$labels <- str_c(sampleinfo$Status, ":", sampleinfo$TimePoint)
ggdendrogram(hclDat2, rotate = TRUE)




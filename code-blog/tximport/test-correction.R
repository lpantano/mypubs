library(tximportData)
dir <- system.file("extdata", package="tximportData")
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
names(files) <- paste0("sample",1:6)

tx2gene <- read.csv(file.path(dir, "tx2gene.gencode.v27.csv"))

txi.no <- tximport(files, type="salmon", tx2gene=tx2gene, countsFromAbundance = "no",txOut = T)

txi.stpm <- tximport(files, type="salmon", tx2gene=tx2gene, countsFromAbundance = "scaledTPM", txOut = T)

txi.lstpm <- tximport(files, type="salmon", tx2gene=tx2gene, countsFromAbundance = "lengthScaledTPM", txOut = T)

head(txi.no$abundance)
head(txi.no$counts) # different
colSums(txi.no$abundance) # 1e+06
colSums(txi.no$counts) #

head(txi.stpm$abundance)
head(txi.stpm$counts) # different
colSums(txi.stpm$counts) # tpm scaled to original lib size

head(txi.lstpm$abundance)
head(txi.lstpm$counts) # different
colSums(txi.lstpm$counts) # tpm scaled to original lib size
library(magrittr)
dds.no = DESeqDataSetFromTximport(txi.no, samples, desig = ~1) %>% 
    varianceStabilizingTransformation()

dds.stpm = DESeqDataSetFromTximport(txi.stpm, samples, desig = ~1) %>% 
    varianceStabilizingTransformation()

dds.lstpm = DESeqDataSetFromTximport(txi.lstpm, samples, desig = ~1) %>% 
    varianceStabilizingTransformation()

head(assay(dds.no))
head(assay(dds.stpm))
head(assay(dds.lstpm))
head(txi.no$length)

cor(rowMeans(assay(dds.no)), rowMeans(txi.no$length), method = "spearman")
plot(rowMeans(assay(dds.no)), log2(rowMeans(txi.no$length)))
cor(rowMeans(assay(dds.stpm)), rowMeans(txi.stpm$length), method = "spearman")
plot(rowMeans(assay(dds.stpm)), log2(rowMeans(txi.stpm$length)))
cor(rowMeans(assay(dds.lstpm)), rowMeans(txi.lstpm$length), method = "spearman")
plot(rowMeans(assay(dds.lstpm)), log2(rowMeans(txi.lstpm$length)))

cor(rowMeans(assay(dds.no)), rowMeans(assay(dds.stpm)))
cor(rowMeans(assay(dds.no)), rowMeans(assay(dds.lstpm)))



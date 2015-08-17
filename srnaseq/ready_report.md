---
output:
  knitrBootstrap::bootstrap_document:
    theme: readable
    highlight: zenburn
    theme.chooser: TRUE
    highlight.chooser: TRUE
  html_document:
    toc: true
    highlight: zenburn
---

last update Fri Aug 14 15:20:39 2015 by @lopantano



For replication [go here](http://seqcluster.readthedocs.org/example_pipeline.html)


```r
library(knitr)

library(ggplot2)
library(reshape)
library(DESeq2)
```

```
## Loading required package: S4Vectors
## Loading required package: stats4
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following objects are masked from 'package:dplyr':
## 
##     combine, intersect, setdiff, union
## 
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
##     Position, rank, rbind, Reduce, rep.int, rownames, sapply,
##     setdiff, sort, table, tapply, union, unique, unlist, unsplit
## 
## 
## Attaching package: 'S4Vectors'
## 
## The following object is masked from 'package:dplyr':
## 
##     rename
## 
## The following object is masked from 'package:reshape':
## 
##     rename
## 
## Loading required package: IRanges
```

```
## Warning: replacing previous import by 'BiocGenerics::mad' when loading
## 'IRanges'
```

```
## Warning: replacing previous import by 'BiocGenerics::IQR' when loading
## 'IRanges'
```

```
## Warning: multiple methods tables found for 'mad'
```

```
## Warning: multiple methods tables found for 'IQR'
```

```
## 
## Attaching package: 'IRanges'
## 
## The following objects are masked from 'package:dplyr':
## 
##     collapse, desc, slice
## 
## The following object is masked from 'package:reshape':
## 
##     expand
## 
## Loading required package: GenomicRanges
## Loading required package: GenomeInfoDb
```

```
## Warning: replacing previous import by 'IRanges::mad' when loading
## 'GenomeInfoDb'
```

```
## Warning: replacing previous import by 'IRanges::IQR' when loading
## 'GenomeInfoDb'
```

```
## Warning: replacing previous import by 'IRanges::mad' when loading
## 'GenomicRanges'
```

```
## Warning: replacing previous import by 'IRanges::IQR' when loading
## 'GenomicRanges'
```

```
## Warning: replacing previous import by 'IRanges::mad' when loading 'XVector'
```

```
## Warning: replacing previous import by 'IRanges::IQR' when loading 'XVector'
```

```
## Loading required package: SummarizedExperiment
## Loading required package: Biobase
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Warning: replacing previous import by 'IRanges::mad' when loading
## 'SummarizedExperiment'
```

```
## Warning: replacing previous import by 'IRanges::IQR' when loading
## 'SummarizedExperiment'
```

```
## Loading required package: Rcpp
## Loading required package: RcppArmadillo
```

```
## Warning: replacing previous import by 'IRanges::mad' when loading 'DESeq2'
```

```
## Warning: replacing previous import by 'IRanges::IQR' when loading 'DESeq2'
```

```
## Warning: replacing previous import by 'IRanges::mad' when loading
## 'AnnotationDbi'
```

```
## Warning: replacing previous import by 'IRanges::IQR' when loading
## 'AnnotationDbi'
```

```r
library(genefilter)
```

```
## 
## Attaching package: 'genefilter'
## 
## The following object is masked from 'package:base':
## 
##     anyNA
```

```r
library(CHBUtils)
library(gtools)
library(gridExtra)
library(devtools)
library(dplyr)
library(isomiRs)
```

```
## Loading required package: DiscriMiner
```

```
## Warning: replacing previous import by 'IRanges::mad' when loading 'isomiRs'
```

```
## Warning: replacing previous import by 'IRanges::IQR' when loading 'isomiRs'
```

```r
knitr::opts_chunk$set(tidy=TRUE, highlight=TRUE, dev="png", fig.width=6,fig.heigh=6,
               cache=FALSE, highlight=TRUE, autodep=TRUE, warning=FALSE, error=FALSE,
               message=FALSE, prompt=TRUE, comment='', fig.cap='', bootstrap.show.code=FALSE)

root_file = paste0(root_path,"/report/")

condition = "condition"
```



# Overview

mirRQC project

[paper](http://www.nature.com/nmeth/journal/v11/n8/full/nmeth.3014.html)

samples overview:

>> Universal Human miRNA reference RNA (Agilent Technologies, #750700), human brain total RNA (Life Technologies, #AM6050), human liver total RNA (Life Technologies, #AM7960) and MS2-phage RNA (Roche, #10165948001) were diluted to a platform-specific concentration. RNA integrity and purity were evaluated using the Experion automated gel electrophoresis system (Bio-Rad) and Nanodrop spectrophotometer. All RNA samples were of high quality (miRQC A: RNA quality index (RQI, scale from 0 to 10) = 9.0; miRQC B: RQI = 8.7; human liver RNA: RQI = 9.2) and high purity (data not shown). RNA was isolated from serum prepared from three healthy donors using the miRNeasy mini kit (Qiagen) according to the manufacturer's instructions, and RNA samples were pooled. Informed consent was obtained from all donors (Ghent University Ethical Committee). Different kits for isolation of serum RNA are available; addressing their impact was outside the scope of this work. Synthetic miRNA templates for let-7a-5p, let-7b-5p, let-7c, let-7d-5p, miR-302a-3p, miR-302b-3p, miR-302c-3p, miR-302d-3p, miR-133a and miR-10a-5p were synthesized by Integrated DNA Technologies and 5′ phosphorylated. Synthetic let-7 and miR-302 miRNAs were spiked into MS2-phage RNA and total human liver RNA, respectively, at 5 × 106 copies/μg RNA. These samples do not contain endogenous miR-302 or let-7 miRNAs, which allowed unbiased analysis of cross-reactivity between the individual miR-302 and let-7 miRNAs measured by the platform and the different miR-302 and let-7 synthetic templates in a complex RNA background. Synthetic miRNA templates for miR-10a-5p, let-7a-5p, miR-302a-3p and miR-133a were spiked in human serum RNA at 6 × 103 copies per microliter of serum RNA or at 5-times higher, 2-times higher, 2-times lower and 5-times lower concentrations, respectively. All vendors received 10 μl of each serum RNA sample.

[samples](http://www.nature.com.ezp-prod1.hul.harvard.edu/nmeth/journal/v11/n8/images/nmeth.3014-F1.jpg)


```r
> # setwd(root_path)
> files = read.table(file.path(root_path, "summary_re.csv"), sep = ",", header = T, 
+     colClasses = "character")
> 
> samples = files[, "sample_id"]
> 
> names_stats = files[, "size_stats"]
> names(names_stats) = samples
> 
> groups = files[, "group"]
> names(groups) = samples
> 
> summarydata = data.frame(row.names = samples, samples = samples, group = groups)
> design <- data.frame(row.names = files$sample_id, condition = files$group)
```

# Exploratory analysis

In this section we will see exploratory figures about quality of the data, 
reads with adapter, reads mapped to miRNAs and reads mapped to other small RNAs. 

## size distribution

After adapter removal, we can plot the size distribution of the small RNAs. In a normal
small RNA sample, we should see a peak at 22/23 and maybe another at 26 or 31 depending on the biological background.


```r
> tab = data.frame()
> for (sample in samples) {
+     d = read.table(file.path(root_path, names_stats[sample]), sep = " ")
+     tab = rbind(tab, d %>% mutate(sample = sample, group = groups[sample]))
+ }
> 
> 
> reads_adapter = tab %>% group_by(sample, group) %>% summarise(total = sum(V2))
> ggplot(reads_adapter, aes(x = sample, y = total, fill = group)) + geom_bar(stat = "identity", 
+     position = "dodge") + ggtitle("total number of reads with adapter") + theme(axis.text.x = element_text(angle = 90, 
+     vjust = 0.5, hjust = 1))
```

![](figure/adapter-1.png) 

```r
> ggplot(tab, aes(x = V1, y = V2, group = sample)) + geom_bar(stat = "identity", 
+     position = "dodge") + facet_wrap(~group, ncol = 2) + ggtitle("size distribution") + 
+     ylab("abundance") + xlab("size") + theme(axis.text.x = element_text(angle = 90, 
+     vjust = 0.5, hjust = 1))
```

![](figure/adapter-2.png) 

## miRNA

### total miRNA expression annotated with mirbase


```r
> mi_files = file.path(root_path, files[, "miraligner"])
> row.names(design) = samples
> 
> obj <- IsomirDataSeqFromFiles(files = mi_files, design = design, header = T, 
+     cov = 1)
```



```r
> ggplot(data.frame(sample = colnames(counts(obj)), total = colSums(counts(obj)))) + 
+     geom_bar(aes(x = sample, y = total), stat = "identity") + theme(axis.text.x = element_text(angle = 90, 
+     vjust = 0.5, hjust = 1))
```

![](figure/mirna-mirbase-1.png) 

```r
> mirna_step <- as.data.frame(colSums(counts(obj)))
```

### Distribution of mirna expression


```r
> ggplot(melt(counts(obj))) + geom_boxplot(aes(x = X2, y = value)) + scale_y_log10() + 
+     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
```

![](figure/depth-1.png) 

### cumulative distribution of miRNAs


```r
> cs <- as.data.frame(apply(counts(obj), 2, cumsum))
> cs$pos <- 1:nrow(cs)
> 
> ggplot((melt(cs, id.vars = "pos"))) + geom_point(aes(x = pos, y = value, color = variable)) + 
+     scale_y_log10()
```

![](figure/cum-1.png) 

## Others small RNA

The data was analyzed with [seqcluster](http://seqcluster.readthedocs.org/)

This tools used all reads, uniquely mapped and multi-mapped reads. The first
step is to cluster sequences in all locations they overlap. The second step is to 
create meta-clusters: is the unit that merge all clusters that share the 
same sequences. This way the output are meta-clusters, common sequences that could
come from different region of the genome.


### genome covered


```r
> cov_stats <- read.table(file.path(root_path, "align", "seqs_clean_rmlw.bam_cov.tsv"), 
+     sep = "\t", check.names = F)
> 
> kable(cov_stats[cov_stats$V1 == "genome", ] %>% select(coverage = V2, ratio_genome = V5))
```



|    | coverage| ratio_genome|
|:---|--------:|------------:|
|129 |        0|    0.9997890|
|130 |        1|    0.0002112|

The normal value for data with strong small RNA signal is: 0.0002

### classification

Number of reads in the data after each step:

* raw: initial reads
* cluster: after cluster detection
* multimap: after meta-cluster detection using all hits


```r
> reads_stats <- read.table(file.path(root_path, "seqcluster", "cluster", "read_stats.tsv"), 
+     sep = "\t", check.names = F)
> ggplot(reads_stats, aes(x = V2, y = V1, fill = V3)) + geom_bar(stat = "identity", 
+     position = "dodge") + labs(list(x = "samples", y = "reads")) + scale_fill_brewer("steps", 
+     palette = "Set1") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
+     hjust = 1))
```

![](figure/reads-track-1.png) 



```r
> clus <- read.table(file.path(root_path, files[1, "clusters"]), header = T, sep = "\t", 
+     row.names = 1, check.names = FALSE)
> ann <- clus[, 2]
> toomany <- clus[, 1]
> clus_ma <- clus[, 3:ncol(clus)]
> clus_ma = clus_ma[, row.names(design)]
```


Check complex meta-clusters: This kind of events happen when there are small RNA over the whole genome, and all
repetitive small RNAs map to thousands of places and sharing many sequences in many positions.
If any meta-cluster is > 40% of the total data, maybe it is worth to add some filters
like: minimum number of counts `-e` or `--min--shared` in `seqcluster prepare`



```r
> clus_ma_norm = sweep(clus_ma, 2, colSums(clus_ma), "/")
> head(clus_ma_norm[toomany > 0, ])
```

```
 [1] miRQC_A                                        
 [2] miRQC_A_repeat                                 
 [3] miRQC_B                                        
 [4] miRQC_B_repeat                                 
 [5] miRQC_C                                        
 [6] miRQC_C_repeat                                 
 [7] miRQC_D                                        
 [8] miRQC_D_repeat                                 
 [9] liver_miR-302a                                 
[10] liver_miR-302b                                 
[11] liver_miR-302c                                 
[12] liver_miR-302d                                 
[13] MS2_let-7a-5p                                  
[14] MS2_let-7b-5p                                  
[15] MS2_let-7c                                     
[16] MS2_let-7d-5p                                  
[17] serum_spiked_miRs_constant_concentration       
[18] serum_spiked_miRs_constant_concentration_repeat
[19] serum_spiked_miRs_variable_concentration       
[20] serum_spiked_miRs_variable_concentration_repeat
<0 rows> (or 0-length row.names)
```

Until here is an example of the `Rmd` template that the user can get from `seqcluster-helper` and render directly in `R` / `Rstudio`.

Contribution by class:


```r
> rRNA <- colSums(clus_ma[grepl("rRNA", ann) & grepl("miRNA", ann) == F, ])
> miRNA <- colSums(clus_ma[grepl("miRNA", ann), ])
> tRNA <- colSums(clus_ma[grepl("tRNA", ann) & grepl("rRNA", ann) == F & grepl("repeat", 
+     ann) == F & grepl("miRNA", ann) == F, ])
> piRNA <- colSums(clus_ma[grepl("piRNA", ann) & grepl("rRNA", ann) == F & grepl("miRNA", 
+     ann) == F, ])
> rmsk <- colSums(clus_ma[grepl("repeat", ann) & grepl("rRNA", ann) == F & grepl("miRNA", 
+     ann) == F & grepl("piRNA", ann) == F, ])
> total <- colSums(clus_ma)
> 
> dd <- data.frame(samples = names(rRNA), rRNA = rRNA, miRNA = miRNA, tRNA = tRNA, 
+     piRNA = piRNA, rmsk = rmsk, total = total)
> 
> ggplot(melt(dd)) + geom_bar(aes(x = samples, y = value, fill = variable), stat = "identity", 
+     position = "dodge") + scale_fill_brewer(palette = "Set1") + theme(axis.text.x = element_text(angle = 90, 
+     vjust = 0.5, hjust = 1))
```

![](figure/cluster_type-1.png) 

```r
> dd_norm = dd
> dd_norm[, 2:5] = sweep(dd[, 2:5], 1, dd[, 6], "/")
> ggplot(melt(dd_norm[, 1:5])) + geom_bar(aes(x = samples, y = value, fill = variable), 
+     stat = "identity", position = "dodge") + scale_fill_brewer(palette = "Set1") + 
+     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
+     labs(list(title = "relative proportion of small RNAs", y = "% reads"))
```

![](figure/cluster_type-2.png) 


# Comparison



```r
> library(DESeq2)
> library(DEGreport)
> library(vsn)
```


```r
> filter_handle <- function(res) {
+     res_nona <- res[!is.na(res$padj), ]
+     keep <- res_nona$padj < 0.1
+     res_nona[keep, ]
+ }
> 
> handle_deseq2 = function(dds, summarydata, column) {
+     all_combs = combn(levels(summarydata[, column]), 2, simplify = FALSE)
+     all_results = list()
+     contrast_strings = list()
+     for (comb in all_combs) {
+         contrast_string = paste(comb, collapse = "_vs_")
+         contrast = c(column, comb)
+         res = results(dds, contrast = contrast)
+         res = res[order(res$padj), ]
+         all_results = c(all_results, res)
+         contrast_strings = c(contrast_strings, contrast_string)
+     }
+     names(all_results) = contrast_strings
+     return(all_results)
+ }
> 
> plot_MA = function(res) {
+     for (i in seq(length(res))) {
+         plotMA(res[[i]])
+         title(paste("MA plot for contrast", names(res)[i]))
+     }
+ }
> 
> plot_volcano = function(res) {
+     for (i in seq(length(res))) {
+         stats = as.data.frame(res[[i]][, c(2, 6)])
+         p = volcano_density_plot(stats, title = names(res)[i], lfc.cutoff = 1)
+         print(p)
+     }
+ }
> 
> do_de = function(raw, summarydata, condition) {
+     dss = DESeqDataSetFromMatrix(countData = raw[rowMeans(raw) > 3, ], colData = summarydata, 
+         design = ~condition)
+     dss = DESeq(dss)
+     plotDispEsts(dss)
+     dss
+ }
> 
> do_norm = function(dss, root_path, prefix) {
+     rlog_ma = assay(rlog(dss))
+     count_ma = counts(dss, normalized = TRUE)
+     raw = counts(dss, normalized = FALSE)
+     fn_log = paste0(root_file, prefix, "log_matrix.txt")
+     write.table(rlog_ma, fn_log, sep = "\t")
+     fn_count = paste0(root_file, prefix, "count_matrix.txt")
+     write.table(count_ma, fn_count, sep = "\t")
+     fn_raw = paste0(root_file, prefix, "raw_matrix.txt")
+     write.table(count_ma, fn_raw, sep = "\t")
+ }
```

## mirna


```r
> mirna_c = counts(obj)
```


### abundance detection of miRQC samples

There are 4 samples:

* A: universal human RNA sample
* B: human brain sample
* C: 0.25 * A + 0.75 * B
* D: 0.25 * B + 0.75 * A

If A > B then A > C > D > B

If A < B then A < C < D < B


Note that C/D samples are swapped in the paper and in the GEO web. Text from the paper:

>> These samples (termed miRQC A–D) consist of 100% Universal Human miRNA Reference RNA (UHmiRR; A), 100% human brain RNA (HBR; B) and two titrations thereof (C = 0.75A + 0.25B and D = 0.25A + 0.75B). 

while in the GEO:

>> Source name 	miRQC C
Organism 	Homo sapiens
Characteristics 	biomaterial: Mixture of 25% UHRR and 75% HuBr Total RNA

>> Source name 	miRQC D
Organism 	Homo sapiens
Characteristics 	biomaterial: Mixture of 75% UHRR and 25% HuBr Total RNA




```r
> keep = grepl("miRQC", colnames(mirna_c))
> keep_d = grepl("miRQC", rownames(design))
> library(edgeR)
> dge = DGEList(mirna_c[, keep])
> dge = calcNormFactors(dge, method = "upperquartile")
> mirna_n = cpm(dge, normalized.lib.sizes = T, log = T)
> mirna_n = mirna_n[rowMedians(mirna_n) > 5, ]
```


```r
> library(GGally)
> GGally::ggpairs(mirna_n, axisLabels = "internal")
```

![](figure/mirqc-cor-1.png) 


```r
> mirqca = rowMeans(mirna_n[, 1:2])
> mirqcb = rowMeans(mirna_n[, 3:4])
> mirqcc = rowMeans(mirna_n[, 5:6])
> mirqcd = rowMeans(mirna_n[, 7:8])
> 
> top_a = mirqca > mirqcb
> conc_a = top_a & mirqca > mirqcd & mirqcd > mirqcc
> top_b = mirqca < mirqcb
> conc_b = top_b & mirqcb > mirqcc & mirqcc > mirqcd
```

miRNAs which mirqca > mirqcd > mrqcc are 111 out of 111
miRNAs which mirqcb > mirqcc > mrqcd are 174 out of 181


ratio expression summary of A/D


```r
> kable(t(as.matrix(summary(mirqca[top_a] - mirqcd[top_a]))))
```



|   Min.| 1st Qu.| Median|   Mean| 3rd Qu.|   Max.|
|------:|-------:|------:|------:|-------:|------:|
| 0.2856|  0.4861| 0.5816| 0.5789|  0.6748| 0.9794|

the logFC is 0.5 that is a FC of log2(1/0.75)=0.5 that should be D according to description


ratio expression summary of A/C


```r
> kable(t(as.matrix(summary(mirqca[top_a] - mirqcc[top_a]))))
```



|   Min.| 1st Qu.| Median|  Mean| 3rd Qu.|  Max.|
|------:|-------:|------:|-----:|-------:|-----:|
| 0.3806|   1.071|  2.007| 1.667|   2.156| 2.891|

the logFC is 1.6 that is a FC of log2(1/0.25)=2 that should be C according to description

Same happens when comparing B vs D 


```r
> kable(t(as.matrix(summary(mirqcb[top_b] - mirqcd[top_b]))))
```



|   Min.| 1st Qu.| Median|  Mean| 3rd Qu.|  Max.|
|------:|-------:|------:|-----:|-------:|-----:|
| 0.3703|  0.9992|  1.389| 1.409|   1.822| 2.544|

and B vs C


```r
> kable(t(as.matrix(summary(mirqcb[top_b] - mirqcc[top_b]))))
```



|  Min.| 1st Qu.| Median|   Mean| 3rd Qu.|  Max.|
|-----:|-------:|------:|------:|-------:|-----:|
| 0.268|  0.6431| 0.7505| 0.7431|  0.8395| 1.044|

according to this: 

miRQC_C = 0.75 * miRQC_B + 0.25 * miRQC_A

miRQC_D = 0.75 * miRQC_A + 0.25 * miRQC_B

that is the same that is in the GEO data set description file.


### specificity

>> we spiked in 8 synthetic miRNAs from two miRNA families into human liver RNA (miR-302a/b/c/d) or MS2-phage RNA (let-7a/b/c/d)

We should only see those miRNAs in those samples and not in anywhere else.


```r
> kable(isoSelect(obj, mirna = "hsa-let-7a-5p")[, grepl("MS", rownames(design))])
```



|                                                     | MS2_let-7a-5p| MS2_let-7b-5p| MS2_let-7c| MS2_let-7d-5p|
|:----------------------------------------------------|-------------:|-------------:|----------:|-------------:|
|hsa-let-7a-5p.iso.t5:0.seed:0.t3:0.ad:u-A.mm:0       |             0|             0|          0|             0|
|hsa-let-7a-5p.iso.t5:0.seed:0.t3:d-T.ad:0.mm:0       |             0|             0|          0|             0|
|hsa-let-7a-5p.iso.t5:0.seed:0.t3:u-GTT.ad:u-ATT.mm:0 |            14|             0|          0|             0|
|hsa-let-7a-5p.iso.t5:0.seed:0.t3:u-T.ad:0.mm:0       |             0|             0|          0|             0|
|hsa-let-7a-5p.ref.t5:0.seed:0.t3:0.ad:0.mm:0         |           453|             0|          0|             0|

```r
> kable(isoSelect(obj, mirna = "hsa-let-7b-5p")[, grepl("MS", rownames(design))])
```



|                                                  | MS2_let-7a-5p| MS2_let-7b-5p| MS2_let-7c| MS2_let-7d-5p|
|:-------------------------------------------------|-------------:|-------------:|----------:|-------------:|
|hsa-let-7b-5p.iso.t5:0.seed:0.t3:0.ad:u-A.mm:0    |             0|             0|          0|             0|
|hsa-let-7b-5p.iso.t5:0.seed:0.t3:0.ad:u-AT.mm:0   |             0|             0|          0|             0|
|hsa-let-7b-5p.iso.t5:0.seed:0.t3:d-T.ad:0.mm:0    |             0|             0|          0|             0|
|hsa-let-7b-5p.iso.t5:0.seed:0.t3:d-T.ad:u-T.mm:0  |             0|             0|          0|             0|
|hsa-let-7b-5p.iso.t5:0.seed:0.t3:u-T.ad:0.mm:0    |             0|             0|          0|             0|
|hsa-let-7b-5p.iso.t5:0.seed:0.t3:u-T.ad:u-A.mm:0  |             0|             0|          0|             0|
|hsa-let-7b-5p.iso.t5:0.seed:0.t3:u-T.ad:u-AT.mm:0 |             0|             0|          0|             0|
|hsa-let-7b-5p.iso.t5:0.seed:0.t3:u-TT.ad:0.mm:0   |             0|             0|          0|             0|
|hsa-let-7b-5p.ref.t5:0.seed:0.t3:0.ad:0.mm:0      |             0|           753|          0|             0|

```r
> kable(isoSelect(obj, mirna = "hsa-let-7c-5p")[, grepl("MS", rownames(design))])
```



|                                               | MS2_let-7a-5p| MS2_let-7b-5p| MS2_let-7c| MS2_let-7d-5p|
|:----------------------------------------------|-------------:|-------------:|----------:|-------------:|
|hsa-let-7c-5p.iso.t5:0.seed:0.t3:0.ad:u-A.mm:0 |             0|             0|          0|             0|
|hsa-let-7c-5p.iso.t5:0.seed:0.t3:d-T.ad:0.mm:0 |             0|             0|          0|             0|
|hsa-let-7c-5p.iso.t5:0.seed:0.t3:u-T.ad:0.mm:0 |             0|             0|          0|             0|
|hsa-let-7c-5p.ref.t5:0.seed:0.t3:0.ad:0.mm:0   |             0|             0|       1150|             0|

```r
> kable(isoSelect(obj, mirna = "hsa-let-7d-5p")[, grepl("MS", rownames(design))])
```



|                                                | MS2_let-7a-5p| MS2_let-7b-5p| MS2_let-7c| MS2_let-7d-5p|
|:-----------------------------------------------|-------------:|-------------:|----------:|-------------:|
|hsa-let-7d-5p.iso.t5:0.seed:0.t3:0.ad:0.mm:1TA  |             0|             0|          0|             0|
|hsa-let-7d-5p.iso.t5:0.seed:0.t3:0.ad:u-A.mm:0  |             0|             0|          0|             0|
|hsa-let-7d-5p.iso.t5:0.seed:0.t3:d-T.ad:0.mm:0  |             0|             0|          0|             0|
|hsa-let-7d-5p.iso.t5:0.seed:0.t3:u-T.ad:0.mm:0  |             0|             0|          0|             0|
|hsa-let-7d-5p.iso.t5:0.seed:0.t3:u-TT.ad:0.mm:0 |             0|             0|          0|             0|
|hsa-let-7d-5p.ref.t5:0.seed:0.t3:0.ad:0.mm:0    |             0|             0|          0|           262|


```r
> kable(isoSelect(obj, mirna = "hsa-miR-302a-3p")[, grepl("liver", rownames(design))])
```



|                                                   | liver_miR-302a| liver_miR-302b| liver_miR-302c| liver_miR-302d|
|:--------------------------------------------------|--------------:|--------------:|--------------:|--------------:|
|hsa-miR-302a-3p.iso.t5:0.seed:0.t3:0.ad:u-A.mm:0   |              0|              0|              0|              0|
|hsa-miR-302a-3p.iso.t5:0.seed:0.t3:u-A.ad:0.mm:0   |              0|              0|              0|              0|
|hsa-miR-302a-3p.iso.t5:0.seed:0.t3:u-A.ad:u-G.mm:0 |              0|              0|              0|              0|
|hsa-miR-302a-3p.iso.t5:0.seed:0.t3:u-A.ad:u-T.mm:0 |              0|              0|              0|              0|
|hsa-miR-302a-3p.iso.t5:0.seed:0.t3:u-GA.ad:0.mm:0  |              0|              0|              0|              0|
|hsa-miR-302a-3p.iso.t5:d-T.seed:0.t3:0.ad:0.mm:0   |              0|              0|              0|              0|
|hsa-miR-302a-3p.ref.t5:0.seed:0.t3:0.ad:0.mm:0     |             41|              0|              0|              0|

```r
> kable(isoSelect(obj, mirna = "hsa-miR-302b-3p")[, grepl("liver", rownames(design))])
```



|                                                   | liver_miR-302a| liver_miR-302b| liver_miR-302c| liver_miR-302d|
|:--------------------------------------------------|--------------:|--------------:|--------------:|--------------:|
|hsa-miR-302b-3p.iso.t5:0.seed:0.t3:0.ad:u-A.mm:0   |              0|              0|              0|              0|
|hsa-miR-302b-3p.iso.t5:0.seed:0.t3:0.ad:u-T.mm:0   |              0|              0|              0|              0|
|hsa-miR-302b-3p.iso.t5:0.seed:0.t3:u-AG.ad:0.mm:0  |              0|              0|              0|              0|
|hsa-miR-302b-3p.iso.t5:0.seed:0.t3:u-G.ad:0.mm:0   |              0|              0|              0|              0|
|hsa-miR-302b-3p.iso.t5:0.seed:0.t3:u-G.ad:u-A.mm:0 |              0|              0|              0|              0|
|hsa-miR-302b-3p.iso.t5:0.seed:0.t3:u-G.ad:u-T.mm:0 |              0|              0|              0|              0|
|hsa-miR-302b-3p.iso.t5:0.seed:0.t3:u-TAG.ad:0.mm:0 |              0|              0|              0|              0|
|hsa-miR-302b-3p.iso.t5:d-T.seed:0.t3:0.ad:u-T.mm:0 |              0|              0|              0|              0|
|hsa-miR-302b-3p.ref.t5:0.seed:0.t3:0.ad:0.mm:0     |              0|             72|              0|              0|

```r
> kable(isoSelect(obj, mirna = "hsa-miR-302c-3p")[, grepl("liver", rownames(design))])
```



|                                                    | liver_miR-302a| liver_miR-302b| liver_miR-302c| liver_miR-302d|
|:---------------------------------------------------|--------------:|--------------:|--------------:|--------------:|
|hsa-miR-302c-3p.iso.t5:0.seed:0.t3:0.ad:0.mm:18GC   |              0|              0|              0|              0|
|hsa-miR-302c-3p.iso.t5:0.seed:0.t3:0.ad:u-T.mm:0    |              0|              0|              0|              0|
|hsa-miR-302c-3p.iso.t5:0.seed:0.t3:u-GG.ad:0.mm:0   |              0|              0|              0|              0|
|hsa-miR-302c-3p.iso.t5:0.seed:0.t3:u-TGG.ad:0.mm:0  |              0|              0|              0|              0|
|hsa-miR-302c-3p.iso.t5:d-T.seed:0.t3:0.ad:0.mm:0    |              0|              0|              0|              0|
|hsa-miR-302c-3p.iso.t5:d-T.seed:0.t3:0.ad:u-C.mm:0  |              0|              0|              0|              0|
|hsa-miR-302c-3p.iso.t5:d-T.seed:0.t3:0.ad:u-T.mm:0  |              0|              0|              0|              0|
|hsa-miR-302c-3p.iso.t5:d-T.seed:0.t3:0.ad:u-TT.mm:0 |              0|              0|              0|              0|
|hsa-miR-302c-3p.iso.t5:d-TA.seed:0.t3:0.ad:u-T.mm:0 |              0|              0|              0|              0|
|hsa-miR-302c-3p.ref.t5:0.seed:0.t3:0.ad:0.mm:0      |              0|              0|             85|              0|

```r
> kable(isoSelect(obj, mirna = "hsa-miR-302d-3p")[, grepl("liver", rownames(design))])
```



|                                                   | liver_miR-302a| liver_miR-302b| liver_miR-302c| liver_miR-302d|
|:--------------------------------------------------|--------------:|--------------:|--------------:|--------------:|
|hsa-miR-302d-3p.iso.t5:0.seed:0.t3:0.ad:u-A.mm:0   |              0|              0|              0|              0|
|hsa-miR-302d-3p.iso.t5:0.seed:0.t3:0.ad:u-T.mm:0   |              0|              0|              0|              0|
|hsa-miR-302d-3p.iso.t5:0.seed:0.t3:u-GT.ad:0.mm:0  |              0|              0|              0|              0|
|hsa-miR-302d-3p.iso.t5:0.seed:0.t3:u-T.ad:0.mm:0   |              0|              0|              0|              0|
|hsa-miR-302d-3p.iso.t5:0.seed:0.t3:u-T.ad:u-A.mm:0 |              0|              0|              0|              0|
|hsa-miR-302d-3p.iso.t5:0.seed:0.t3:u-T.ad:u-G.mm:0 |              0|              0|              0|              0|
|hsa-miR-302d-3p.iso.t5:0.seed:0.t3:u-TGT.ad:0.mm:0 |              0|              0|              0|              0|
|hsa-miR-302d-3p.iso.t5:d-T.seed:0.t3:0.ad:0.mm:0   |              0|              0|              0|              0|
|hsa-miR-302d-3p.ref.t5:0.seed:0.t3:0.ad:0.mm:0     |              0|              0|              0|            178|


According to the text they saw cross-mapping between these miRNAs, but here we are seeing perfect annotation.

[Figure-e](http://www.nature.com.ezp-prod1.hul.harvard.edu/nmeth/journal/v11/n8/images/nmeth.3014-F4.jpg)

Probably the method was causing that, since I am not seeing this with seqbuster.

### accuracy


```r
> target = c("hsa-miR-10a-5p", "hsa-let-7a-5p", "hsa-miR-302a-3p", "hsa-miR-133a-3p")
> keep = grepl("serum", colnames(mirna_c))
> keep_d = grepl("serum", rownames(design))
> dge = DGEList(mirna_c[, keep])
> dge = calcNormFactors(dge, method = "TMM")
> serum_n = cpm(dge, normalized.lib.sizes = T, log = T)
> 
> serum_n[target, ]
```

```
                serum_spiked_miRs_constant_concentration
hsa-miR-10a-5p                                 17.569177
hsa-let-7a-5p                                  11.664742
hsa-miR-302a-3p                                 5.037919
hsa-miR-133a-3p                                 9.042439
                serum_spiked_miRs_constant_concentration_repeat
hsa-miR-10a-5p                                        17.766612
hsa-let-7a-5p                                         11.793600
hsa-miR-302a-3p                                        5.037919
hsa-miR-133a-3p                                        5.037919
                serum_spiked_miRs_variable_concentration
hsa-miR-10a-5p                                 19.903714
hsa-let-7a-5p                                  12.667415
hsa-miR-302a-3p                                 5.037919
hsa-miR-133a-3p                                 5.037919
                serum_spiked_miRs_variable_concentration_repeat
hsa-miR-10a-5p                                        19.619057
hsa-let-7a-5p                                         13.105691
hsa-miR-302a-3p                                        5.037919
hsa-miR-133a-3p                                        5.037919
```

For miR10 and let-7a the changes are clear although not proportional.

I cannot detect sequences for miR-302 and miR-133 (< 5 reads). I checked directly in the raw data and these sequences are not there.

I think they weren't captured by the sequencing at all. The text shows the same for miR-302, and not 
detection of changes in the miR-133 ( I guess for the same reason).

## isomiRs

As an example of some figures you can do with this package.  (read more [here](http://lpantano.github.io/isomiRs)). 

There is one figure per type of isomiR. 
The y-axes shows the percentage of unique sequences with that change.
The x-axes shows the percentage of abundance with that change.


```r
> obj = isoPlot(obj)
```

![](figure/isomir-example-1.png) 

```r
> obj = isoPlot(obj, type = "iso3")
```

![](figure/isomir-example-2.png) 

```r
> obj = isoPlot(obj, type = "add")
```

![](figure/isomir-example-3.png) 

```r
> obj = isoPlot(obj, type = "subs")
```

![](figure/isomir-example-4.png) 

Seems there are some nt-changes for serum and MS2 samples at position 13/14, 
and at position 9-11 for miRQC and liver samples. These are just
exploratory figures and could lead to a deep analysis of differential
expression of isomiRs.

## clusters
 
The same logic was applied to clusters detection.


```r
> keep = grepl("miRQC", colnames(clus_ma))
> keep_d = grepl("miRQC", rownames(design))
> dge = DGEList(clus_ma[, keep])
> dge = calcNormFactors(dge, method = "upperquartile")
> clus_ma_norm = cpm(dge, normalized.lib.sizes = T, log = T)
> clus_ma_norm = clus_ma_norm[rowMedians(clus_ma_norm) > 5, ]
```

### matrix correlation among samples


```r
> GGally::ggpairs(clus_ma_norm, axisLabels = "internal")
```

![](figure/cluster-cor-1.png) 

### abundance detection


```r
> mirqca = rowMeans(clus_ma_norm[, 1:2])
> mirqcb = rowMeans(clus_ma_norm[, 3:4])
> mirqcc = rowMeans(clus_ma_norm[, 5:6])
> mirqcd = rowMeans(clus_ma_norm[, 7:8])
> 
> top_a = mirqca > mirqcb
> conc_a = top_a & mirqca > mirqcd & mirqcd > mirqcc
> 
> top_b = mirqca < mirqcb
> conc_b = top_b & mirqcb > mirqcc & mirqcc > mirqcd
```

clusters which mirqca > mirqcd > mrqcc are 139 out of 147
clusters which mirqcb > mirqcc > mrqcd are 222 out of 230

ratio expression summary of A/D


```r
> kable(t(as.matrix(summary(mirqca[top_a] - mirqcd[top_a]))))
```



|    Min.| 1st Qu.| Median|   Mean| 3rd Qu.|   Max.|
|-------:|-------:|------:|------:|-------:|------:|
| -0.1564|   0.208| 0.3339| 0.3144|  0.4182| 0.7166|

the logFC is 0.5 that is a FC of log2(1/0.75)=0.5 that should be D


ratio expression summary of A/C


```r
> kable(t(as.matrix(summary(mirqca[top_a] - mirqcc[top_a]))))
```



|    Min.| 1st Qu.| Median| Mean| 3rd Qu.|  Max.|
|-------:|-------:|------:|----:|-------:|-----:|
| -0.4284|  0.5452|  1.588| 1.28|   1.937| 2.526|

the logFC is 1.6 that is a FC of log2(1/0.25)=2 that should be C

Same happens when comparing B vs D 


```r
> kable(t(as.matrix(summary(mirqcb[top_b] - mirqcd[top_b]))))
```



|    Min.| 1st Qu.| Median|  Mean| 3rd Qu.|  Max.|
|-------:|-------:|------:|-----:|-------:|-----:|
| 0.06527|  0.7847|   1.15| 1.218|   1.606| 3.599|

and B vs C


```r
> kable(t(as.matrix(summary(mirqcb[top_b] - mirqcc[top_b]))))
```



|    Min.| 1st Qu.| Median|   Mean| 3rd Qu.| Max.|
|-------:|-------:|------:|------:|-------:|----:|
| -0.7351|  0.3412| 0.4821| 0.4784|  0.5845| 2.35|

The exactly same thing than we saw for miRNA analysis and in concordance
with the samples description file.

# Conclusions

We can conclude that the miRNAs and clusters quantification is accurate. The mapping annotation for miRNAs is perfect with those examples, and the difference detection is good for 2 miRNAs, and bad for the other two miRNAs due to a lack of reads support for those miRNAs.

In general, `seqbuster/seqcluster` show good accuracy for abundance detection and mapping accuracy for miRNAs.

# Thanks

special thanks to the author of that papers to make data available. I encourage to use this data for any tool that analyze small RNA data.

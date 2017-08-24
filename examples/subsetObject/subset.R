library(DESeq2)
library(SummarizedExperiment)
library(dplyr)

dummyClass <- setClass("dummyClass", contains = "SummarizedExperiment")

dds <- makeExampleDESeqDataSet(n = 10000, m = 100)
ddsSmall <- makeExampleDESeqDataSet(n = 100, m = 5)


dummyBig <- new("dummyClass",
                SummarizedExperiment(assays = SimpleList(counts(dds))),
                colData = colData(dds),
                metadata = list(dds = dds))
save(dummyBig, file = "big.rda")

dummySmall <- new("dummyClass",
                SummarizedExperiment(assays = SimpleList(counts(ddsSmall))),
                colData = colData(ddsSmall),
                metadata = list(dds = ddsSmall))
save(dummySmall, file = "micro.rda")


setMethod("[", "dummyClass",
          function(x, i, j, ..., drop) {
              se <- SummarizedExperiment(assays = SimpleList(assays(x)[[1]]),
                                         colData = colData(x))
              tmp <- se[i, j]
              set <- DESeqDataSetFromMatrix(countData = assays(tmp)[[1]],
                                            colData = colData(tmp),
                                            design = ~1)
              new("dummyClass",
                  SummarizedExperiment(assays = SimpleList(counts(set))),
                  colData = colData(set),
                  metadata = list(dds = set))
})

dummySilly <- dummyBig[1:100, 1:5]
save(dummySilly, file = "microSilly.rda")


.intSubset <- function(x){
    DESeqDataSetFromMatrix(countData = assays(x)[[1]],
                           colData = colData(x),
                           design = ~1)
}
setMethod("[", "dummyClass",
          function(x, i, j, ..., drop) {
              se <- SummarizedExperiment(assays = SimpleList(assays(x)[[1]]),
                                         colData = colData(x))
              tmp <- se[i, j]
              set <- .intSubset(tmp)
              new("dummyClass",
                  SummarizedExperiment(assays = SimpleList(counts(set))),
                  colData = colData(set),
                  metadata = list(dds = set))
})


dummySmart <- dummyBig[1:100, 1:5]
save(dummySmart, file = "microSmart.rda")

rdaSize <- sapply(list.files(".", pattern = "rda"),
                  function(x) file.size(x) / 1024 )
objSize <- sapply(list(dummyBig, dummySmall, dummySilly, dummySmart),
                  function(x) format(object.size(x), units = "Kb"))

data.frame(file = list.files(".", pattern = "rda"),
           rdaSize = rdaSize,
           objSize = objSize) %>% knitr::kable()

lapply(list.files(".", pattern = "rda"), function(fn){
    unlink(fn)
})

sessionInfo()

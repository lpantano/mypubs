load("fa_isonetwork.rda")
library(isomiRs)
library(DESeq2)
library(SummarizedExperiment)

mi_cold = fa_cold[colnames(fa_mirna),, drop = F]
mi_cold$day = droplevels(mi_cold$day)
mi_dse = DESeqDataSetFromMatrix(round(2^fa_mirna), mi_cold, design = ~ day)
mi_dse = DESeq(mi_dse)
mi_res = results(mi_dse)
mi_res = mi_res[!is.na(mi_res$padj),]
mi_top = row.names(mi_res[mi_res$padj < 0.05,])
mi_rse = SummarizedExperiment(assays = SimpleList(norm=fa_mirna),
                              colData = mi_cold,
                              metadata = list(sign=mi_top))

m_cold = fa_cold[colnames(fa_mrna),, drop = F]
m_cold$day = droplevels(m_cold$day)
m_dse = DESeqDataSetFromMatrix(round(2^fa_mrna), m_cold, design = ~ day)
m_dse = DESeq(m_dse)
m_res = results(m_dse)
m_res = m_res[!is.na(m_res$padj),]
m_top = row.names(m_res[m_res$padj < 0.05,])
m_cold_fixed = m_cold[m_cold$day  %in%  mi_cold$day,, drop = F]
m_cold_fixed$day = droplevels(m_cold_fixed$day)
m_rse = SummarizedExperiment(assays=SimpleList(norm=fa_mrna[,m_cold$day  %in%  mi_cold$day]),
                              colData = m_cold_fixed,
                              metadata=list(sign=m_top))

library(org.Mm.eg.db)
library(clusterProfiler)
ego <- enrichGO(m_top,
                org.Mm.eg.db,
                "ENSEMBL",
                ont = "BP",
                universe = rownames(m_res))

library(targetscan.Mm.eg.db)
m2t = mirna2targetscan(mi_top, species = "mmu", org = org.Mm.eg.db, keytype = "ENSEMBL")

mirna_targets = findTargets(mi_rse, m_rse, m2t[,c("ENSEMBL", "mir")], summarize = "day", min_cor = -0.7)

data <- isoNetwork(mi_rse, m_rse, min_fc = 0.3,
                   summarize = "day", target = mirna_targets,
                   enrich = ego)
library(cowplot)
isoPlotNet(data, minGenes = 3) + ggsave("network.pdf", width = 12, height = 9)


library(multiMiR)
multimir_results <- get_multimir(org     = 'mmu',
                                 mirna   = mi_top,
                                 table   = 'validated',
                                 summary = TRUE)

library(magrittr)
m2t_multimir = slot(multimir_results, "data")[,c("target_ensembl", "mature_mirna_id")] %>%  dplyr::filter(target_ensembl != "") %>% dplyr::distinct()

mirna_targets = findTargets(mi_rse, m_rse, m2t_multimir, summarize = "day", min_cor = -0.7)
data <- isoNetwork(mi_rse, m_rse, min_fc = 0.3,
                   summarize = "day", target = mirna_targets,
                   enrich = ego) # this error for not having enough genes


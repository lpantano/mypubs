---
title: "miRNA annotation"
author: "Lorena Pantano"
date: "10/02/2015"
output:
  knitrBootstrap::bootstrap_document:
    theme: readable
    highlight: zenburn
    theme.chooser: TRUE
    highlight.chooser: TRUE
  html_document:
    highlight: zenburn

---



# Methods
 * isomiRs are simulated with the script inside `seqcsluter`. see `mirannotation.sh`.
 * To see reproducibility, please check: `mirannotation.sh` script. `installer.sh` helps with installation of some specific tools that are less likely you have installed.
 * miRBase 21 is the version used.
 * scoring hits to precursor: in case the tool gives a score, the best score will be used. If not the first hit is the one used. If score is the same, the first hit is used. Tools that can be scored are: `bowtie, bowtie2, chimira_blast, GEM, microzer, miraligner, miraligner-python, novoaling, razer3, STAR`
 * only `miraligner*` and `srnabench `gives miRNA annotation, so these tools should have an advantage since they are parsing the hits to get the best annotation. Anyway, I am trying to get the best of all of them using the precursor name instead of the miRNA. In the last section, since the rest tools are only alignerts, I call miRNA the firs three names in the precursor (i.e hsa-let-7a in hsa-let-7a-1/2/3).
 
*Note* Future work is the use of the [python miraligner](http://seqcluster.readthedocs.org/mirna_annotation.html#miraligner-inside-seqcluster) version to annotate miRNA/isomiR from BAM files with precursor hits coming from any tool that gives a BAM file as output.



# Mapped
Proportion of mapped and no-mapped sequences
![](figure/mapped-mir-1.png) 


# Size effect
How size affects the alignments
![](figure/size-mir-1.png) 

# Isomirs effect
How changes in the mature miRNA affect the alignment:

* red: addition + mutation
* blue: addition
* green: mutation
* purple: just trimming precursor at different position

![](figure/iso-mir-1.png) 


# Specificity at precursor level
How many were assigned to the correct miRNA using the precursor name as the true positive.
Red would be "not correct" and blue "correct". This is only considering mapped sequences.
![](figure/sp-precursor-1.png) 


# Specificity at miRNA level
Same logic than before but with the miRNA names in case the tool gives the miRNA names, if not, only the three first field in the name are used as annotation (i.e hsa-let-7a-1 will ignore any character beyond 7a). This will increase the number of TP, since many miRNAs has multiple precursors being the same mature miRNA at the end.
![](figure/sp-mir-1.png) 


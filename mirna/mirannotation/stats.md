---
title: "miRNA annotation"
author: "Lorena Pantano"
date: "18/02/2015"
output:
  knitrBootstrap::bootstrap_document:
    theme: readable
    highlight: zenburn
    theme.chooser: TRUE
    highlight.chooser: TRUE
  html_document:
    highlight: zenburn

---






# Mapped
Proportion of mapped and no-mapped sequences
${image?fileName=mapped%2Dmir%2D1%2Epng&align=none&scale=100}


# Size effect
How size affects the alignments
${image?fileName=size%2Dmir%2D1%2Epng&align=none&scale=100}

# Isomirs effect
How changes affect the alignment
${image?fileName=iso%2Dmir%2D1%2Epng&align=none&scale=100}


# Specificity at precursor level
How many were assigned to only the correct miRNA. Normally miRNA can map to other miRNAs that are from the same family/precursor. "Red" will be sequences mapping to multiple miRNAs.
${image?fileName=sp%2Dprecursor%2D1%2Epng&align=none&scale=100}


# Specificity at miRNA level

${image?fileName=sp%2Dmir%2D1%2Epng&align=none&scale=100}


library(tidyverse)
library(ggplot2)
library(qqman)

df = data.table::fread("out/EOPC_assoc.pancreas_cancer.glm.logistic.hybrid") %>% dplyr::filter(TEST == "ADD")

# save manhattan plot as PDF
pdf("manhattanPlotCandidateEOPCvariants.pdf",width=7,height=5)

# manhattan plot and label top candidate
qqman::manhattan(df,chr = "#CHROM",bp="POS",snp="ID",p="P",annotatePval = 0.01,main="Candidate variants associated with EOPC risk") 
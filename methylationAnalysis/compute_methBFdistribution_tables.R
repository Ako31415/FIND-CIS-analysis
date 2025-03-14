#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(dplyr)


geno <- args[1]
meth <- args[2]
binsize <- as.numeric(as.character(args[3]))

inputfile <- paste0(geno, ".WW.q255.41bp.", meth, ".noNP.GT1.justMPs.tsv")

NAM_CG <- read.delim(inputfile, header=FALSE, na.strings="n.a.")
colnames(NAM_CG) <- c("B73_chr", "B73_pos", "B73_allele", "NAM_allele", "NAM_ID", "GT", "AMP", "BF", "B73_meth", "NAM_meth", "meth_diff")

NAM_CG <- subset(NAM_CG, BF!=1 & BF!=0)

NAM_CG$sign_meth_diff <- ifelse((NAM_CG$B73_meth < 0.1 & NAM_CG$NAM_meth > 0.7), "hypo/hyper", 
                                ifelse((NAM_CG$NAM_meth < 0.1 & NAM_CG$B73_meth > 0.7), "hyper/hypo", "EM"))


readyData <- NAM_CG %>%
  mutate(binCounts = cut(BF, breaks = seq(0, 1, by = binsize))) %>%
  group_by(binCounts,sign_meth_diff) %>%
  count() %>%
  group_by(binCounts) %>%
  mutate(Prop = n/(sum(n)))

readyData$NAMline <- geno

fout <- paste0(geno, ".WW.q255.41bp.", meth, ".noNP.GT1.justMPs.dmProportionsPerBFBin.", binsize, "bins.csv")
write.csv(readyData, file = fout)



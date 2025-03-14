#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(reshape)


methylContext <- args[1]

fin <- paste0("ALL_lines.WW.q255.41bp.", methylContext, ".noNP.GT1.justMPs.dmProportionsPerBFBin.0.025bins.csv")

alllines <- read.delim(fin, header=FALSE, sep=",")
colnames(alllines) <- c("runningnumber", "binCounts", "sign_meth_diff", "n", "Prop", "NAMline")


fout1 <- paste0("ALL_lines.WW.q255.41bp.", methylContext, ".noNP.GT1.justMPs.dmProportionsPerBFBin.0.025bins.noLinesPlot.svg")
svg(fout1, width=5.8, height=5.32)
ggplot(alllines, aes(x=binCounts, y=Prop, color=sign_meth_diff, fill=sign_meth_diff, group=sign_meth_diff)) +
  geom_boxplot(outlier.shape = NA, width=1, lwd=1, aes(group = NULL)) +
  theme_classic() +
  scale_color_manual(values = c("#66ccff", "#ff0000", "#000000")) +
  scale_fill_manual(values = c("#66ccff", "#ff0000", "#000000")) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=5), legend.position="none")
dev.off()




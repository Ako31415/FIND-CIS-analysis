library("ggplot2")

dat <- read.delim("ALL_Lines.41bpVS11bp.CG.hyperBindAFPChanges.tsv", sep='\t', header=FALSE)

colnames(dat) <- c("NAM_line", "total_AFP_count", "hyperbindCount_41bp", "to_hyper_count", "to_hypo_count", "to_NoDM_count", "hyperbind_percent_41bp", "to_hyper_percent", "to_hypo_percent", "to_NoDM_percent")



summary(dat$"hyperbind_percent_41bp")
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.009602 0.014854 0.018765 0.017704 0.020317 0.024768 


summary(dat$"to_hyper_percent")
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.5116  0.5443  0.5616  0.5609  0.5777  0.6275 

summary(dat$"to_NoDM_percent")
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.3725  0.4223  0.4384  0.4391  0.4557  0.4884 


summary(dat$"to_hypo_percent")
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#      0       0       0       0       0       0 


# Plot

dat$xlab <- "F1s"

fout1 <- "plot_hyperbindAFPs_boxplot.pdf"
  
pdf(fout1, width=4, height=8)
print(ggplot(dat, aes(y = hyperbind_percent_41bp*100, x = xlab, fill = xlab, color=xlab)) +
          geom_boxplot(outlier.shape = NA, lwd=0.3) +
          geom_point(position = position_jitterdodge(jitter.width = 0.3), pch = 20, color="black", size = 2, alpha = 0.7) +
          #geom_jitter(width = 0.1, size = 0.4, alpha = 0.9, size = 15) +
          labs(x="", y = "AMPs with binding bias to hypermethylated allele [%]", fill="") +
          theme_classic() +
          theme(legend.position = "none") +
          scale_fill_manual(values = c("lightgrey")) +
          scale_color_manual(values=c("black")) +
          scale_y_continuous(expand = c(0,0), limits = c(0,100)))
dev.off()






  fout1 <- "plot_hyperbindToNoDmAFPs_boxplot.pdf"
  
  pdf(fout1, width=4, height=8)
  print(ggplot(dat, aes(y = to_NoDM_percent*100, x = xlab, fill = xlab, color=xlab)) +
          geom_boxplot(outlier.shape = NA, lwd=0.6) +
          geom_point(position = position_jitterdodge(jitter.width = 0.3), pch = 20, color="black", size = 2, alpha = 0.7) +
          #geom_jitter(width = 0.1, size = 0.4, alpha = 0.9, size = 15) +
          labs(x="", y = "AMPs with affinity to hypermethylated allele with no methylation difference directly at AMP site [%]", fill="") +
          theme_classic() +
          theme(legend.position = "none") +
          scale_fill_manual(values = c("lightgrey")) +
          scale_color_manual(values=c("black")) +
          scale_y_continuous(expand = c(0,0), limits = c(0,100)))

  dev.off()



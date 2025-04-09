#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(ggplot2)


fin <- paste0("ALL_Lines.WW.q255.41bp.CG.noNP.GT1.justMPs.bed")
all_SNPs <- read.delim(fin, header=FALSE, na.strings="n.a.")
colnames(all_SNPs) <- c("B73_chr", "B73_pos", "B73_allele", "NAM_allele", "NAM_ID", "GT", "AMP", "BF", "B73_meth", "NAM_meth", "meth_diff")

# make n.r. in BF column to n.a.
all_SNPs[all_SNPs$BF=="n.r.", ] <- NA

all_SNPs$BF <- as.numeric(as.character(all_SNPs$BF))

all_SNPs$B73_ID <- paste0(all_SNPs$B73_chr, ":", all_SNPs$B73_pos)

all_SNPs$NAM_line <- sapply(strsplit(as.character(all_SNPs$NAM_ID), "-"), "[", 1)


d <- subset(all_SNPs, BF!=1 & BF!=0 & GT=="1/1")


rm(all_SNPs)
gc()

# not diff. methylated is either both alleles <10% or >70% methylated:
d$diff_meth <- ifelse(((d$B73_meth<0.1 & d$NAM_meth>0.7) | (d$NAM_meth<0.1 & d$B73_meth>0.7)), "diff_meth", 
                  ifelse(((d$B73_meth<0.1 & d$NAM_meth<0.1) | (d$B73_meth>0.7 & d$NAM_meth>0.7)), "no_diff_meth", "unclear_meth"))


d$dmAMP <- ifelse((d$AMP=="AMP" & d$diff_meth=="diff_meth"), "dmAMP", 
                  ifelse((d$AMP=="AMP" & d$diff_meth=="no_diff_meth"), "noDM_AMP", 
                         ifelse((d$AMP=="MP" & d$diff_meth=="diff_meth"), "dmMP",
                                ifelse((d$AMP=="MP" & d$diff_meth=="no_diff_meth"), "noDM_MP",
                                       ifelse((d$AMP=="AMP" & d$diff_meth=="unclear_meth"), "ucm_AMP",
                                              ifelse((d$AMP=="MP" & d$diff_meth=="unclear_meth"), "ucm_MP", "mistake"))))))



# remove line's sites that have unclear methylation
data <- d[d$dmAMP!="ucm_AMP", ]
data <- data[data$dmAMP!="ucm_MP", ]


######## Try 2 dm 2 nodm 1 amp as requirements for the sites
######## as in the number of inbred lines that are dm/noDm/amp at the specified site
######## only take sites that meet the above requirement
######## (dm=differentially methylated; noDm=not differentially methylated; amp=allele-specifically bound MOA-footprint)


# get IDs of sites

dm_groups <- data %>% group_by(B73_ID) %>% count(diff_meth)

dm2IDs <- unique(dm_groups[(dm_groups$diff_meth=="diff_meth" & dm_groups$n>1),]$B73_ID)

nodm2IDs <- unique(dm_groups[(dm_groups$diff_meth=="no_diff_meth" & dm_groups$n>1),]$B73_ID)

dm2nodm2IDs <- dm2IDs[dm2IDs %in% nodm2IDs]


# get sites
dm2nodm2 <- data[data$B73_ID %in% dm2nodm2IDs, ]

dm2nodm2AMPIDs <- unique(dm2nodm2[dm2nodm2$AMP=="AMP",]$B73_ID)

dm2nodm2AMP <- dm2nodm2[dm2nodm2$B73_ID %in% dm2nodm2AMPIDs, ]



dm2nodm2AMP_counts <- dm2nodm2AMP %>% group_by(B73_ID,diff_meth) %>% count(AMP)

dm2nodm2AMP_props <- dm2nodm2AMP_counts %>% group_by(B73_ID,diff_meth) %>% mutate(Prop=n/sum(n))


write.csv(dm2nodm2AMP_props, file = "dm2nodm2AMP_props.csv")



# Plot:
plot_dm2nodm2AMP_props <- dm2nodm2AMP_props[dm2nodm2AMP_props$AMP=="AMP", ]



  fout1 <- "plot_dm2nodm2AMP_props.violin.pdf"
  
  pdf(fout1, width=4, height=8)
  print(ggplot(plot_dm2nodm2AMP_props, aes(x = diff_meth, y = Prop)) +
          geom_violin(lwd=1.3, aes(fill = diff_meth, color = diff_meth)) +
          geom_boxplot(width = 0.1, fill = c("red", "#66ccff"), color = c("black", "#005580"), outlier.shape = NA, lwd=1, alpha = 0.8) +
          # geom_point(position = position_jitterdodge(jitter.width = 0.3), pch = 20, color="black", size = 2, alpha = 0.7) +
          # geom_jitter(width = 0.1, size = 0.4, alpha = 0.9, size = 15) +
          labs(x="", y = "Proportion of AMPs per site [%]", fill="") +
          theme_classic() +
          theme(legend.position = "none") +
          scale_fill_manual(values = c("#ff8080", "#b3e7ff")) +
          scale_color_manual(values=c("#999999", "#80d4ff")) +
          scale_y_continuous(expand = c(0,0), limits = c(0,1)))

  dev.off()



  fout1 <- "plot_dm2nodm2AMP_props.violin_beeswarm.pdf"
  
  pdf(fout1, width=4, height=8)
  print(ggplot(plot_dm2nodm2AMP_props, aes(x = diff_meth, y = Prop)) +
          geom_violin(lwd=1.3, aes(fill = diff_meth, color = diff_meth)) +
          geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA, lwd=1, alpha = 0.8) +
          geom_point(aes(fill = diff_meth), position = position_jitterdodge(jitter.width = 0.3), pch = 20, color="black", size = 2, alpha = 0.7) +
          geom_jitter(width = 0.1, size = 0.4, alpha = 0.9, size = 15) +
          labs(x="", y = "Proportion of AMPs per site [%]", fill="") +
          theme_classic() +
          theme(legend.position = "none") +
          scale_fill_manual(values = c("red", "#66ccff")) +
          scale_color_manual(values=c("black", "#0077b3")) +
          scale_y_continuous(expand = c(0,0), limits = c(0,1)))

  dev.off()



  fout1 <- "plot_dm2nodm2AMP_props.density.pdf"
  
  pdf(fout1, width=8, height=4)
  print(ggplot(plot_dm2nodm2AMP_props, aes(x = Prop, color = diff_meth)) +
          geom_density() +
          # geom_point(position = position_jitterdodge(jitter.width = 0.3), pch = 20, color="black", size = 2, alpha = 0.7) +
          # geom_jitter(width = 0.1, size = 0.4, alpha = 0.9, size = 15) +
          labs(x="Proportion of AMPs per site [%]") +
          theme_classic() +
          # theme(legend.position = "none") +
          scale_color_manual(values = c("red", "#66ccff")) +
          scale_y_continuous(expand = c(0,0)))

  dev.off()
  




# plot density of all combinations

dm2nodm2AMP_props$label <- paste0("proportion_of_", dm2nodm2AMP_props$AMP, "_out_of_all_", dm2nodm2AMP_props$diff_meth, "_F1s_per_site")



  fout1 <- "plot_dm2nodm2AMP_props.density.combinations.pdf"
  
  pdf(fout1, width=7, height=3.5)
  # print(ggplot(dm2nodm2AMP_props, aes(x = Prop, color = diff_meth, linetype=AMP)) +
  print(ggplot(dm2nodm2AMP_props, aes(x = Prop, color = label)) +
          geom_density( linewidth=1.8 ) +
          labs(x="Proportion per site [%]") +
          theme_classic() +
          theme(legend.position = "none") +
          scale_color_manual(values = c("red", "#66ccff", "black", "#00cc44")) +
          scale_y_continuous(expand = c(0,0)))

  dev.off()



# Solid red is (dmAMP F1s/ allDM F1s), dotted red is (dmNonAMP F1s/ allDM F1s), solid blue is (noDmAMP F1s/ allNoDM F1s), dotted blue is (noDmNonAMP F1s/ allNoDM F1s)

# red is (dmAMP F1s/ allDM F1s), black is (dmNonAMP F1s/ allDM F1s), blue is (noDmAMP F1s/ allNoDM F1s), green is (noDmNonAMP F1s/ allNoDM F1s)






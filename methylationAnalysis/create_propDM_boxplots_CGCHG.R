#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(dplyr)


all_lines <- data.frame()

# List of all the NAM parent names
hybrids <- c("A188", "B97", "CML247", "CML103", "CML277", "CML322", "CML333", "CML69", "HP301", "IL14H", "Ki11", "Ki3", "Ky21", "M162W", "M37W", "Mo17", "Mo18W", "Ms71", "NC358", "Oh43", "Oh7b", "P39", "Tx303", "W22")



  
for(aline in hybrids) {
    
    
  geno <- aline

  NAMlineCG <- read.delim(paste0(geno, ".WW.q255.41bp.CG.noNP.GT1.justMPs.tsv"), header=FALSE, na.strings="n.a.")
  colnames(NAMlineCG) <-  c("B73_chr", "B73_pos", "B73_allele", "NAM_allele", "NAM_ID", "GT", "AMP", "BF", "B73_CG", "NAM_CG", "CG_diff")


  NAMlineCHG <- read.delim(paste0(geno, ".WW.q255.41bp.CHG.noNP.GT1.justMPs.tsv"), header=FALSE, na.strings="n.a.")
  colnames(NAMlineCHG) <-  c("B73_chr", "B73_pos", "B73_allele", "NAM_allele", "NAM_ID", "GT", "AMP", "BF", "B73_CHG", "NAM_CHG", "CHG_diff")
      
  NAMline <- merge(NAMlineCG, NAMlineCHG, sort=TRUE)

  
    
  # make n.r. in BF column to n.a.
  NAMline[NAMline$BF=="n.r.", ] <- NA
  NAMline$BF <- as.numeric(as.character(NAMline$BF))
    
  NAMline$diff_meth <- ifelse((NAMline$B73_CG < 0.1 & NAMline$NAM_CG > 0.7), "DM", 
                                ifelse((NAMline$NAM_CG < 0.1 & NAMline$B73_CG > 0.7), "DM", 
                                    ifelse((NAMline$B73_CHG < 0.1 & NAMline$NAM_CHG > 0.7), "DM",
                                        ifelse((NAMline$NAM_CHG < 0.1 & NAMline$B73_CHG > 0.7), "DM", "EM"))))
    
    
  d <- subset(NAMline, !(BF==0.0 | BF==1.0) & GT=="1/1") %>%
    select(NAM_ID, BF, AMP, diff_meth) 
    
  d$strict_AMP <- ifelse(d$BF <= 0.15 | d$BF >= 0.85, "strAMP", "no_strAMP")
    
    
  all_counts <- d %>%
    count(diff_meth) %>%
    mutate(Prop = n/(sum(n)))
    
    
  AMP_counts <- subset(d, AMP=="AMP") %>%
    count(diff_meth) %>%
    mutate(Prop = n/(sum(n)))
    
  strAMP_counts <- subset(d, strict_AMP=="strAMP") %>%
    count(diff_meth) %>%
    mutate(Prop = n/(sum(n)))
    
    
    
  counts <- rbind.data.frame(
    cbind.data.frame(all_counts, cat = "all MPs", line = geno),
    cbind.data.frame(AMP_counts, cat = "Ball AMPs", line = geno),  # all AMPs (Ball AMPs written to have a simple cheat for ordering the bars´)
    cbind.data.frame(strAMP_counts, cat = "str. AMPs", line = geno)
  )
    
  all_lines <- rbind.data.frame(all_lines, counts)
    
    
}
  
  
  
fout <- paste0("ALL_lines.WW.q255.41bp.CGandCHG.noNP.GT1.justMPs.count_DMs.data.txt")
write.csv(all_lines, file = fout)
  
dms <- subset(all_lines, all_lines$diff_meth=="DM")[ , -1]
all_lines <- data.frame()
  

fout1 <- paste0("ALL_lines.WW.q255.41bp.CGandCHG.noNP.GT1.justMPs.proportions.boxplot.plot.pdf")
pdf(fout1, width=4, height=8)
print(ggplot(dms, aes(x = cat, y = Prop, fill = cat)) +
        geom_boxplot(outlier.shape = NA, lwd=1.2) +
        geom_point(position = position_jitterdodge(jitter.width = 0.3), pch = 20, color="black", size = 3) +
        # geom_jitter(width = 0.1, size = 0.4, alpha = 0.9, size = 15) +
        labs(x="", y = paste0("Number of CG and/or CHG diff. methylated MPs [%]"), fill="") +
        theme_classic() +
        theme(legend.position = "none") +
        scale_x_discrete(labels=c("all MPs" = "all MPs", "Ball AMPs" = "all AMPs", "str. AMPs" = "str. AMPs")) +
        scale_fill_manual(values = c("grey11", "lightgrey", "white")) +
        scale_y_continuous(expand = c(0,0), limits = c(0,1)))
# scale_y_continuous(expand = c(0,0), labels = scales::percent))
# scale_y_continuous(expand = c(0,0), labels = scales::percent, limits = c(0,1)))
dev.off()
  
  
fout <- paste0("ALL_lines.WW.q255.41bp.CGandCHG.noNP.GT1.justMPs.StatTests.txt")
sink(fout, append = TRUE)
#sink(fout, append = TRUE, type = "message")
 
print("CGandCHG")
## Test assumptions for ANOVA
  
print("Test normality of residuals:")

# Normality of the residuals:
res_aov <- aov(Prop ~ cat, data = dms)

# Check normality statistically
print(shapiro.test(res_aov$residuals))
 
print("_____________")

print("Perform Kruskal test:")
## Perform Kruskal-Wallis test
print(kruskal.test(Prop ~ cat, data = dms))

print("_____________")

print("Perform pairwise wilcox test:")
# Pairwise test with multiple testing correction:
print(pairwise.wilcox.test(dms$Prop, dms$cat, p.adjust.method = "BH"))

sink()
sink(type = "message")

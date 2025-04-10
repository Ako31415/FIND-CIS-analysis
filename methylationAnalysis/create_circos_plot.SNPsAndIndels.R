library(circlize)
library(stringr)
library(dplyr)

circos.clear()

####### Load data #######


### GWAS data ###
# get cleaned data from all tissues
GWAS_atlas <- read.delim("GWAS_SNPs_from_GWAS_Atlas_database.ALL_tissues.cleaned.bedGraph", header=FALSE)
colnames(GWAS_atlas) <- c("chr", "start", "stop", "pval")


# To see better the y-axis is set to a max of 40 All values >40 are set to 49 so that they are not drawn over the plot region
GWAS40 <- GWAS_atlas
GWAS40$pval <- ifelse(GWAS40$pval>40, 40, GWAS40$pval)



### MOA bQTL data ###
bQTLs <- read.delim("bQTL_WW_SNPsAndINDEL.jbrowse.tidy.bed", header=FALSE)
colnames(bQTLs) <- c("chr", "start", "stop", "value")

# To see something the yscale is set to 45. Since some would be drawn over the bounds of the plots, every value greater 45 is set to 45 instead
bQTLs45 <- bQTLs
bQTLs45$value <- ifelse(bQTLs45$value>45, 45, bQTLs45$value)

# dim(bQTLs[bQTLs$value>45,])



## Selective sweeps ###
# Selective sweeps
Sel_sweeps <- read.delim("remapped_GSE145586_Table_S2_Selective_sweep.cleaned.bedGraph", header=FALSE)
colnames(Sel_sweeps) <- c("chr", "start", "stop", "val")


# To see better the yscale is set to 250. Since some would be drawn over the bounds of the plots, every value greater 250 is set to 250 instead
sweeps250 <- Sel_sweeps
sweeps250$val <- ifelse(sweeps250$val>=250, 250, sweeps250$val)


### Set up the plot ###

# Create Matrix for Chomosome lenghts (B73 v5)
chromlenghts <- matrix(c(0, 308452471, 
                         0, 243675191,
                         0, 238017767,
                         0, 250330460,
                         0, 226353449,
                         0, 181357234,
                         0, 185808916,
                         0, 182411202,
                         0, 163004744,
                         0, 152435371), ncol = 2, byrow = T)





# Draw plot
circos.clear()
col_text <- "black"
circos.par("track.height"=0.8, gap.degree=2, cell.padding=c(0, 0, 0, 0))
circos.initialize(factors=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10"),
                  xlim=chromlenghts)


# Labels for genes
geneLabels <- data.frame(matrix(c("chr1", 80823669, 80823670, "TRE1", "chr1", 177940469, 177940470, "BIF2", "chr1", 198205028, 198205029, "SO", "chr2", 15685108, 15685109, "UPA2", "chr2", 28118441, 28118442, "TIP3d", "chr3", 4730702, 4730703, "TINY", "chr3", 161185978, 161185979, "MADS69", "chr3", 198724865, 198724866, "SBT11", "chr4", 192254463, 192254464, "PGM1", "chr5", 130302870, 130302871, "SWEET4c", "chr8", 126680771, 126680773, "ZCN8", "chr8", 135590787, 135590788, "VGT1", "chr9", 120909923, 120909924, "CCT9", "chr9", 141030560, 141030561, "ZmPHYB2"), ncol = 4, byrow = T))
colnames(geneLabels) <- c("chr", "start", "stop", "label")

geneLabels$start <- as.numeric(as.character(geneLabels$start))
geneLabels$stop <- as.numeric(as.character(geneLabels$stop))

circos.genomicLabels(geneLabels, labels.column = 4, side = "outside", col = "red", line_col = "red", cex = 1, line_lwd = 2, font = 3, connection_height = mm_h(4))



## bar plot bQTL pvalues
circos.genomicTrackPlotRegion(bQTLs45, panel.fun = function(region, value, ...){
  circos.genomicLines(region, value, type = "h", col="#000080ff",...)}
  , bg.col = "grey97", bg.border = T,track.height = 0.15
)



## genes

genes <- data.frame(matrix(c("chr1", 80823669, 80823670, 40, "chr1", 177940469, 177940470, 40, "chr1", 198205028, 198205029, 40, "chr2", 15685108, 15685109, 40, "chr2", 28118441, 28118442, 40, "chr3", 4730702, 4730703, 40, "chr3", 161185978, 161185979, 40, "chr3", 198724865, 198724866, 40, "chr4", 192254463, 192254464, 40, "chr5", 130302870, 130302871, 40, "chr8", 126680771, 126680773, 40, "chr8", 135590787, 135590788, 40, "chr9", 120909923, 120909924, 40, "chr9", 141030560, 141030561, 40), ncol = 4, byrow = T))
colnames(genes) <- c("chr", "start", "stop", "val")

genes$start <- as.numeric(as.character(genes$start))
genes$stop <- as.numeric(as.character(genes$stop))
genes$val <- as.numeric(as.character(genes$val))


circos.genomicTrack(genes, track.index = get.current.track.index(), ylim = c(0, 45), numeric.column = 4, 
                    panel.fun = function(region, value, ...) {
                      # numeric.column is automatically passed to `circos.genomicPoints()`
                      circos.genomicPoints(region, value, col = "red", pch = 19, cex = 0.7, ...)
                    }, bg.border=F)




# bar plot GWAS pvalues
circos.genomicTrackPlotRegion(GWAS40, ylim=c(0, 40), panel.fun = function(region, value, ...){
  circos.genomicLines(region, value, type = "h", col="#07bc3dff",...)}
  , bg.col = "#d9e8f0", bg.border = NA,track.height = 0.15
)


# bar plot with selective sweeps
circos.genomicTrackPlotRegion(sweeps250, ylim=c(0, 250), panel.fun = function(region, value, ...){
  circos.genomicLines(region, value, type = "h", col="#fdcc00ff",...)}
  , bg.col = "grey97", bg.border = T,track.height = 0.15
)


# chromosomes
circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), substring(chr, 4, 5), cex=1.5, col="black", 
              facing="bending.inside", niceFacing=TRUE, font = 2)
}, bg.col="#d9e8f0", bg.border=F, track.height=0.05)

text(0,0, "GWAS", cex = 1.3, col="#07bc3dff")
text(0,0.06, "bQTL", cex = 1.3, col="#000080ff")
text(0,-0.06, "Selection", cex = 1.3, col="#fdcc00ff")








# sessionInfo()
# R version 4.4.3 (2025-02-28)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 20.04.6 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0 
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=de_DE.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Berlin
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] dplyr_1.1.4     stringr_1.5.1   circlize_0.4.16
# 
# loaded via a namespace (and not attached):
#   [1] Matrix_1.7-3        compiler_4.4.3      tidyselect_1.2.1    Rcpp_1.0.13-1       DHARMa_0.4.7       
# [6] splines_4.4.3       boot_1.3-31         yaml_2.3.8          fastmap_1.2.0       lattice_0.22-5     
# [11] R6_2.5.1            generics_0.1.3      shape_1.4.6.1       knitr_1.46          MASS_7.3-64        
# [16] tibble_3.2.1        nloptr_2.1.1        minqa_1.2.8         pillar_1.9.0        rlang_1.1.3        
# [21] utf8_1.2.4          stringi_1.8.4       xfun_0.44           GlobalOptions_0.1.2 cli_3.6.2          
# [26] magrittr_2.0.3      digest_0.6.35       grid_4.4.3          rstudioapi_0.16.0   lme4_1.1-35.5      
# [31] lifecycle_1.0.4     nlme_3.1-167        vctrs_0.6.5         evaluate_0.23       glue_1.7.0         
# [36] fansi_1.0.6         colorspace_2.1-0    rmarkdown_2.27      pkgconfig_2.0.3     tools_4.4.3        
# [41] htmltools_0.5.8.1  





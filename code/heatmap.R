
library(argparser)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)
argv <- arg_parser('')
argv <- add_argument(argv, "--geneID", help="gene ID")
argv <- add_argument(argv, "--normalized_dataset", help="fpkm.xls")
argv <- add_argument(argv, "--out_dir", help="the output directory")
argv <- parse_args(argv)

geneID <- argv$geneID
fpkm_group <-argv$fpkm_group
out_dir <- argv$out_dir

############################File preparation######################################################################


fpkm_file <- read.delim(fpkm_group, header=TRUE)
geneID_file <- read.delim(geneID, header=TRUE)
rownames(fpkm_file) <- fpkm_file[,1]
print(geneID_file[,1])
unionion_fpkm <- subset(fpkm_file, rownames(fpkm_file) %in% geneID_file[,1])

#####################comparison###################################################################################

rownames(unionion_fpkm) <- unionion_fpkm[,1]
unionion_fpkm <- unionion_fpkm[, -1]
unionion_fpkm <- log10(unionion_fpkm + 1)
if(length(unionion_fpkm[,1])<= 50){
  showname <- TRUE
}else{
  showname <- FALSE
}

if(length(unionion_fpkm[1,])<=40){
  cell_widths <- 32
  cell_width <- 10
}else{
  cell_widths <- floor(600/length(unionion_fpkm[1,]))
  cell_width <- floor(400/length(unionion_fpkm[1,]))
}

num <- length(unionion_fpkm[,1])
if(dim(unionion_fpkm)[2]==2){
  scale_row_col = "column"
}else{
  scale_row_col = "row"
} 



png(paste(out_dir, "/heatCluster.png", sep=""), type="cairo-png", res=1200, width = 10, height = 8, units = 'in')



pheatmap(unionion_fpkm,
color=colorRampPalette(rev(c("red", "light yellow")))(500),
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         legend = T,
         show_rownames = showname,
         cellwidth = cell_widths,
         fontsize = cell_width,
         main = "48 BCR clonotypes spike-in UMI count Heatmap")
dev.off()
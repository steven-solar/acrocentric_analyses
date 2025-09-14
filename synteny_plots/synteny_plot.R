library(SVbyEye)
library(dplyr)
library(tibble)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

paf.file <- args[1]
output <- args[2]
method <- args[3]
sats <- args[4]

paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")

annot <- sats 
annot.df <- read.table(annot, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
annot.gr <- GenomicRanges::makeGRangesFromDataFrame(annot.df, keep.extra.columns = TRUE)

if (grepl("hsa2a", paf.file)) {
  seqnames.order <- c("4.mPonAbe1.chr12_hap1_hsa2a", "5.mPonPyg2.chr12_hap1_hsa2a", "3.mGorGor1.chr12_pat_hsa2a", "1.mPanTro3.chr12_hap1_hsa2a", "2.mPanPan1.chr12_mat_hsa2a", "0.chm13.chr2")
} else if (grepl("hsa2b", paf.file)) {
  seqnames.order <- c("4.mPonAbe1.chr11_hap1_hsa2b", "5.mPonPyg2.chr11_hap1_hsa2b", "3.mGorGor1.chr11_mat_hsa2b", "1.mPanTro3.chr13_hap1_hsa2b", "2.mPanPan1.chr13_mat_hsa2b", "0.chm13.chr2")
} else if (grepl("hsa9", paf.file)) {
  seqnames.order <- c("4.mPonAbe1.chr13_hap1_hsa9", "5.mPonPyg2.chr13_hap1_hsa9", "3.mGorGor1.chr13_pat_hsa9", "0.chm13.chr9", "1.mPanTro3.chr11_hap1_hsa9.rev", "2.mPanPan1.chr11_mat_hsa9.rev")
} else if (grepl("hsa13", paf.file)) {
  seqnames.order <- c("4.mPonAbe1.chr14_hap1_hsa13", "5.mPonPyg2.chr14_hap1_hsa13", "3.mGorGor1.chr14_pat_hsa13", "1.mPanTro3.chr14_hap1_hsa13", "2.mPanPan1.chr14_pat_hsa13", "0.chm13.chr13")
} else if (grepl("hsa14", paf.file)) {
  seqnames.order <- c("4.mPonAbe1.chr15_hap1_hsa14", "5.mPonPyg2.chr15_hap1_hsa14", "1.mPanTro3.chr15_hap1_hsa14", "2.mPanPan1.chr15_mat_hsa14", "0.chm13.chr14", "3.mGorGor1.chr15_pat_hsa14")
} else if (grepl("hsa15", paf.file)) {
  seqnames.order <- c("4.mPonAbe1.chr16_hap1_hsa15", "5.mPonPyg2.chr16_hap1_hsa15", "3.mGorGor1.chr16_pat_hsa15", "0.chm13.chr15", "1.mPanTro3.chr16_hap1_hsa15", "2.mPanPan1.chr16_pat_hsa15")
} else if (grepl("hsa18", paf.file)) {
  seqnames.order <- c("4.mPonAbe1.chr17_hap1_hsa18", "5.mPonPyg2.chr17_hap1_hsa18", "3.mGorGor1.chr17_mat_hsa18", "1.mPanTro3.chr17_hap1_hsa18", "2.mPanPan1.chr17_mat_hsa18", "0.chm13.chr18")
} else if (grepl("hsa21", paf.file)) {
  seqnames.order <- c("4.mPonAbe1.chr22_hap1_hsa21", "5.mPonPyg2.chr22_hap1_hsa21", "3.mGorGor1.chr22_mat_hsa21", "1.mPanTro3.chr22_hap1_hsa21", "2.mPanPan1.chr22_pat_hsa21", "0.chm13.chr21")
} else if (grepl("hsa22", paf.file)) {
  seqnames.order <- c("4.mPonAbe1.chr23_hap1_hsa22", "5.mPonPyg2.chr23_hap1_hsa22", "3.mGorGor1.chr23_mat_hsa22", "1.mPanTro3.chr23_hap1_hsa22", "2.mPanPan1.chr23_pat_hsa22", "0.chm13.chr22")
} else if (grepl("hsaY", paf.file)) {
  seqnames.order <- c("4.mPonAbe1.chrY_hap2_hsaY", "5.mPonPyg2.chrY_hap2_hsaY", "3.mGorGor1.chrY_pat_hsaY", "1.mPanTro3.chrY_hap2_hsaY.rev", "2.mPanPan1.chrY_pat_hsaY.rev", "0.chm13.chrY")
} 

pdf(output, width = 7, height = 5)

if (method == "plotMiro") {
	plt <- plotMiro(paf.table = paf.table, color.by = "direction")
} else if (method == "plotAVA") {
	plt <- plotAVA(paf.table = paf.table, color.by = "direction", seqnames.order = seqnames.order)
} else {
	stop("Invalid method")
}

combined_colors <- c(annot.df$color)
unique_colors <- unique(combined_colors)
colors <- setNames(paste0("#", unique_colors), unique_colors)
plt <- addAnnotation(ggplot.obj = plt, annot.gr = annot.gr, shape = "rectangle", coordinate.space = "self", y.label.id = "ID", annotation.level = 0, fill.by = "color", color.palette = colors)
plt <- plt + theme(legend.position = "none")
plt <- plt + labs(x = "Genomic Position (Mbp)")
plt <- plt + scale_x_continuous(labels = scales::label_number(scale = 1e-6))
print(plt)
dev.off()


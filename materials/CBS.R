# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
if (!requireNamespace("DNAcopy", quietly = TRUE))
    BiocManager::install("DNAcopy")
library(DNAcopy)

# setwd("E:/Genomics - MOD 2 [Computational Human Genomics]/chg/project/codes")
folder = "../resources/outputs"
cn           <- read.table(file.path(folder, "SCNA.copynumber.called"), header = TRUE)

pdf("../resources/statistics/SegPlot.pdf")
print(cn)

# plot(cn$raw_ratio,          pch = ".", ylim = c(-2.5, 2.5))
# plot(cn$adjusted_log_ratio, pch = ".", ylim = c(-2.5, 2.5))

# CNA.object   <- CNA(
#     genomdat  = cn$adjusted_log_ratio, 
#     chrom     = cn$chrom,
#     maploc    = cn$chr_start,
#     data.type = 'logratio'
# )

# CNA.smoothed <- smooth.CNA(CNA.object)
# segs         <- segment(
#     CNA.smoothed,
#     alpha       = 0.001,
#     min.width   = 2,
#     undo.splits = "sdundo", #undoes splits that are not at least this many SDs apart.
#     undo.SD     = 3,
#     verbose     = 1
# )

# plot(segs, plot.type = "w", xlab = "chromosomal position", ylab = "Log2 ratio")

dev.off()

# segs2        <- segs$output
# write.table(
#     segs2,
#     file      = file.path(folder, "SCNA.copynumber.called.seg"),
#     row.names = FALSE,
#     col.names = TRUE,
#     quote     = FALSE,
#     sep       = "\t"
# )

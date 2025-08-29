library(data.table)
if (!requireNamespace("CLONETv2", quietly = TRUE))
    BiocManager::install("CLONETv2")
library(CLONETv2)
if (!requireNamespace("TPES", quietly = TRUE))
    BiocManager::install("TPES")
library(TPES)
library(ggplot2)

folder = "../data/output_2"

setwd("E:/Genomics - MOD 2 [Computational Human Genomics]/chg/project/code")

normal = fread(file.path(folder, "Control.ASEReadCounter.csv"),data.table=F)
normal$af = normal$altCount/normal$totalCount
tumor = fread(file.path(folder, "Tumor.ASEReadCounter.csv"),data.table=F)
tumor$af = tumor$altCount/tumor$totalCount

pileup.normal = normal[,c(1,2,4,5,14,8)]
colnames(pileup.normal) = c("chr","pos","ref","alt","af","cov")

pileup.tumor = tumor[,c(1,2,4,5,14,8)]
colnames(pileup.tumor) = c("chr","pos","ref","alt","af","cov")

seg.tb <- fread(file.path(folder, "SCNA.copynumber.called.seg"), data.table=F)

setwd("E:/Genomics - MOD 2 [Computational Human Genomics]/chg/project/data/clonet/")

bt <- compute_beta_table(seg.tb, pileup.tumor, pileup.normal)

## Compute ploidy table with default parameters
pl.table <- compute_ploidy(bt, beta_limit_for_neutral_reads = 0.8)

adm.table <- compute_dna_admixture(beta_table = bt, ploidy_table = pl.table)

allele_specific_cna_table <- compute_allele_specific_scna_table(beta_table = bt,
                                                                ploidy_table = pl.table, 
                                                                admixture_table = adm.table)


check.plot <- check_ploidy_and_admixture(beta_table = bt, ploidy_table = pl.table,
                                         admixture_table = adm.table)
# plot(check.plot)
print(check.plot)

p<-ggplot(data=allele_specific_cna_table,aes(x=cnA,y=cnB,col=log2.corr)) + coord_fixed()  + geom_point()
# plot(p)

# TPES
snv.reads = fread("../output_2/somatic.purity.pm",data.table=F)
snv.reads = snv.reads[which(snv.reads$somatic_status=="Somatic"),]
snv.reads = snv.reads[,c("chrom","position","position","tumor_reads1","tumor_reads2")]
colnames(snv.reads) = c("chr","start","end","ref.count","alt.count")
snv.reads$sample = "Sample.1"

TPES_purity(ID = "Sample.1", SEGfile = seg.tb,
            SNVsReadCountsFile = snv.reads,
            ploidy = pl.table,
            RMB = 0.47, maxAF = 0.6, minCov = 10, minAltReads = 10, minSNVs = 1)

TPES_report(ID = "Sample.1", SEGfile = seg.tb,
            SNVsReadCountsFile = snv.reads,
            ploidy = pl.table,
            RMB = 0.47, maxAF = 0.6, minCov = 10, minAltReads = 10, minSNVs = 1)

# #some examples from real samples (data included in the library)
# TPES_report(ID = "TCGA-HT-8564", SEGfile = TCGA_HT_8564_seg,
#      SNVsReadCountsFile = TCGA_HT_8564_maf, ploidy = TCGA_HT_8564_ploidy,
#      RMB = 0.47, maxAF = 0.55, minCov = 10, minAltReads = 5, minSNVs = 10)
# 
# TPES_purity(ID = "TCGA-HT-8564", SEGfile = TCGA_HT_8564_seg,
#      SNVsReadCountsFile = TCGA_HT_8564_maf, ploidy = TCGA_HT_8564_ploidy,
#      RMB = 0.47, maxAF = 0.55, minCov = 10, minAltReads = 5, minSNVs = 10)



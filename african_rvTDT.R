#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# command line: Rscript --no-save --no-restore african_rvTDT.R /dcl01/beaty/data/agu/african_dataset/african_chr8_127to129_phased.vcf.vcf hg38 /dcl01/beaty/data/agu/african_dataset/african_chr08_v6.fam /dcl01/beaty/data/mtaub/targetedSeq/code/rv-tdt/rvTDT /dcl01/beaty/data/agu/african_dataset/african_chr8_127to129rare.var.vcf fp.african127to129.rvTDT.results.txt /dcl01/beaty/data/agu/african_dataset/african_127to129_rvTDTManhattan_test.pdf RV-TDT_AfricanResults_for_window_size_=_25_SNPs_24_SNP_overlap Position_hg38 -log10p-value_at_center_of_window

library(VariantAnnotation)
library(trio)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(dplyr)
library(ggplot2)

# filepath for phased vcf file; ex: "/home/directory/file"; note: make sure you add a .vcf to your fp.phased.vcf from the fp.phased.vcf in the haplotypePhasing file
fp.phased.vcf = as.character(args[1])

# reference genome; ex: "hg38"
refgenome = as.character(args[2])

# filepath to the PED file; ex: "/home/directory/file"
fp.ped = as.character(args[3])

# filepath for rvTDT function; ex: "/home/directory/file"
fp.RV_TDT = as.character(args[4])

# filepath for vcf that only consists of rare variants and monomorphic SNPs; ex: "/home/directory/file" 
fp.rare.vcf = as.character(args[5])

# .txt file name for rvTDT results; ex: "filepath.rvTDT.results.txt"
fp.rvtdt.results.txt = as.character(args[6])

# filepath for PDF of rvTDT Manhattan plot; ex: "/home/directory/file"
rvTDTManhattan.pdf = as.character(args[7])

# rvTDT Manhattan plot title; ex: "Plot Title"
rvTDTManhattan.title = as.character(args[8])

# rvTDT Manhattan plot x-label; ex: "Position (hg38)"
rvTDTManhattan.xlab = as.character(args[9])

# rvTDT Manhattan plot y-label; ex: "-log10p at Center of Window"
rvTDTManhattan.ylab = as.character(args[10])

vcf_rvtdt <- readVcf(fp.phased.vcf, refgenome)
table(geno(vcf_rvtdt)$GT)

maf.matrix <- geno(vcf_rvtdt)$GT
maf.matrix

maf.matrix[maf.matrix == "0|0"]<- 0
maf.matrix[maf.matrix == "1|0"]<- 1
maf.matrix[maf.matrix == "0|1"]<- 1
maf.matrix[maf.matrix == "1|1"]<- 2

maf.matrix <- matrix(as.numeric(maf.matrix), nrow = nrow(maf.matrix))

maf.matrix <- 0.5*rowMeans(maf.matrix)

# Identify SNPs with MAF < 0.01 and remove monomorphic SNPS
rare.var.snps <- which(maf.matrix < 0.01 & maf.matrix != 0)
#take the first 1000 of the rare variants as a test
#rare.var.snps.1000 = rare.var.snps[1:1000]

# Filter VCF to rare variants only
vcf.rare <- vcf_rvtdt[rare.var.snps, ]
vcf.rare

# Write out VCF
writeVcf(vcf.rare, fp.rare.vcf)
print("hello")
# rvTDT
RV_TDT.results<-rvtrio::RV_TDT(vcf = vcf.rare, ped = read.table(fp.ped, header = TRUE),
                               fp.RV_TDT, window.size = 100, window.type = "M")

# Visualization 
n.windows <- nrow(RV_TDT.results)

# Convert to long format
RV_TDT.results.long <- RV_TDT.results %>%
  tidyr::gather(key = test, value = pval,
                CMC.Analytical,BRV.Haplo,CMC.Haplo,VT.BRV.Haplo,VT.CMC.Haplo,WSS.Haplo)

# Bonferroni Correction 
bonferroni.sig.level <- -log10(0.05/n.windows)

# Plot 
write.table(RV_TDT.results, file = fp.rvtdt.results.txt, sep = "        ", quote = FALSE, row.names = FALSE) 

pdf(rvTDTManhattan.pdf, height = 8, width = 16)

ggplot() +
  geom_line(data = RV_TDT.results.long, aes(group=test, color = test,
                                            x = mid.window.pos, y = -log10(pval)))+
  geom_hline(yintercept=bonferroni.sig.level, linetype=2, color = "red", size=2) +
  labs(title = rvTDTManhattan.title,
       x = rvTDTManhattan.xlab, y = rvTDTManhattan.ylab) +
  guides(color=guide_legend("RV-TDT test type")) +
  scale_linetype_manual(name = "Bonferroni-corrected significance", values = 2,
                        guide = guide_legend(override.aes = list(color = c("red"))))

dev.off()



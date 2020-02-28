# filepath for rare variants vcf 
fp.rare.vcf = 
  
# .txt file name for rvTDT results
fp.rvtdt.results.txt = 

# filepath for PDF of rvTDT Manhattan plot
rvTDTManhattan.pdf = 

# rvTDT Manhattan plot title
rvTDTManhattan.title = 
  
# rvTDT Manhattan plot x-label
rvTDTManhattan.xlab = 

# rvTDT Manhattan plot y-label
rvTDTManhattan.ylab = 

vcf_rvtdt <- readVcf(fp.phased.vcf, refgenome)

maf.matrix <- geno(vcf_rvtdt)$GT

maf.matrix[maf.matrix == "0|0"]<- 0
maf.matrix[maf.matrix == "1|0"]<- 1
maf.matrix[maf.matrix == "0|1"]<- 1
maf.matrix[maf.matrix == "1|1"]<- 2

maf.matrix <- matrix(as.numeric(maf.matrix), nrow = nrow(maf.matrix))

# Identify SNPs with MAF < 0.01 and remove monomorphic SNPS
rare.var.snps <- which(maf.matrix < 0.01 & maf.matrix != 0)

# Filter VCF to rare variants only
vcf.rare <- vcf[rare.var.snps, ]

# Write out VCF
writeVcf(vcf.rare, fp.rare.vcf)

# rvTDT
RV_TDT.results<-rvtrio::RV_TDT(vcf=vcf.rare, ped = ped,
                               filepath.RV_TDT = "/dcl01/beaty/data/mtaub/targetedSeq/code/rv-tdt/rvTDT")

# Visualization 
n.windows<-nrow(RV_TDT.results)

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
       x = rvTDTManhattan.xlab, y = rvTDTManhattan.ylab) s+
  guides(color=guide_legend("RV-TDT test type")) +
  scale_linetype_manual(name = "Bonferroni-corrected significance", values = 2,
                        guide = guide_legend(override.aes = list(color = c("red"))))


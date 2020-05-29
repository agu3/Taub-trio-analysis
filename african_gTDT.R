#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# command line: Rscript --no-save --no-restore african_gTDT.R 127978136 129978136 8 hg38 /dcl01/beaty/data/agu/african_dataset/african_chr08_v6.vcf /dcl01/beaty/data/agu/african_dataset/african_chr08_v6.fam /dcl01/beaty/data/agu/african_chr8_allfiltered_127to129.recode.vcf gTDTAfricanResults_127to129.txt /dcl01/beaty/data/agu/African_gTDTManhattan_127to129_test.pdf /dcl01/beaty/data/agu/african_dataset/african_chr08_v6.vcf.bgz 

library(VariantAnnotation)
library(trio)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(dplyr)
library(ggplot2)

# start of the sequence you want to analyze; enter as vector with one element; ex: c(127978136)
start = as.numeric(args[1])

# end of the sequence you want to analyze; enter as vector with one element; ex: c(12998136)
end = as.numeric(args[2])

# chromosome number; ex: "chr8" or "8" to match vcf file 
chr = as.character(args[3])

# reference genome; ex: "hg38"
refgenome = as.character(args[4])

# filepath to the vcf file of the entire chromosome of interest; ex:        "/home/directory/file"
fp.vcf = as.character(args[5])

# filepath to the bgzipped vcf file of the entire chromosome of interest
#fp.zipped = paste0(fp.vcf, ".bgz") if you already have a bgzipped vcf file
#fp.zipped = bgzip(fp.vcf, as.character(args[10])) to create the bgzip if you don't have the file already
fp.zipped = as.character(args[10])

# filepath to the tabix indexed file of the bgzipped vcf file of the entire chromosome of interest
fp.indexed.tabix = paste0(fp.zipped, ".tbi") # if you have the tbi file already 
#fp.indexed.tabix = indexTabix(fp.zipped, format = "vcf") if you don't have the tbi file already

# filepath to the PED file; ex: "/home/directory/file"
fp.ped = as.character(args[6])

# filepath for the vcf file that has been subsetted to the region of interest; ex: "/home/directory/file"
fp.filtered.vcf = as.character(args[7])

# filepath for the results of the gTDT: ex: "/home/directory/file"
gTDT.results.fp = as.character(args[8])

# filepath for pdf of Manhattan plot generated from the gTDT: ex: "/home/directory/file"
gTDT.Manhattan = as.character(args[9])

fp.tabix = TabixFile(fp.zipped, fp.indexed.tabix)

# Read in subset of the chromosome and write out region of interest to a vcf file 

rng = GRanges(seqnames = chr,
              ranges = IRanges(
                start, end
              )
)

vcf.rng <- readVcf(fp.tabix, refgenome, param=rng)

writeVcf(vcf.rng,fp.filtered.vcf)

ped <- read.table(fp.ped, header=TRUE)

vcf <- readVcf(fp.filtered.vcf, refgenome)

sum(colnames(geno(vcf)$GT) %in% ped$pid)
colnames(geno(vcf)$GT)[!colnames(geno(vcf)$GT)%in% ped$pid]
sum(gsub("\\_2","",colnames(geno(vcf)$GT)) %in% ped$pid)
colnames(vcf) <- gsub("\\_2","",colnames(vcf))

# Perform gTDT and create dataframe of genotypic TDT p-values and snp_names

trio.geno <- trio::vcf2geno(vcf, ped, na.string="./.")
gTDT.results <- trio::colTDT(trio.geno, model = c("additive"))

gTDT.df <- data.frame(
  names(gTDT.results$pval),
  gTDT.results$stat, 
  gTDT.results$pval
)

colnames(gTDT.df)<-c("snp","stat","pval")

gTDT.df <- gTDT.df %>% 
  arrange(pval) %>% 
  mutate(neglogp = -log(pval)) 

# Visualization - extract position and SNP IDs from the filtered vcf and Manhattan Plot 

svp <- ScanVcfParam(fixed="NA", info="NA", geno="NA")
vcf.svp <- readVcf(fp.vcf, refgenome, svp)

snp <- names(vcf.svp) 
pos <- start(rowRanges(vcf.svp)) 
snp.pos.df <- data.frame(snp,pos)

gTDT.df <- left_join(gTDT.df,snp.pos.df)

write.table(gTDT.df, file = gTDT.results.fp, sep = "        ", quote = FALSE, row.names = FALSE)

pdf(gTDT.Manhattan, height = 8, width = 16)

ggplot(gTDT.df)+
  geom_point(aes(x=pos,y=neglogp)) +
  labs(title="Genotypic TDT results",
       x="SNP position", y = "-log10p")

dev.off()

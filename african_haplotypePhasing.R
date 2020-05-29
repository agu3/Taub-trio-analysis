#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# command line: Rscript --no-save --no-restore african_haplotypePhasing.R 127978136 129978136 8 hg38 /dcl01/beaty/data/agu/african_dataset/african_chr08_v6.vcf /dcl01/beaty/data/agu/african_dataset/african_chr08_v6.fam /dcl01/beaty/data/agu/filepath.beagle4 1000 /dcl01/beaty/data/agu/african_dataset/african_chr8_127to129_tophase.vcf /dcl01/beaty/data/agu/african_dataset/african_chr8_127to129_phased.vcf /dcl01/beaty/data/agu/african_dataset/african_chr08_v6.vcf.bgz /dcl01/beaty/data/agu/african_dataset/african_chr08_v6.vcf.bgz.tbi

library(VariantAnnotation)
library(trio)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(dplyr)
library(ggplot2)

devtools::install_github("lindagai/rvtrio")
library(rvtrio)

# start of the sequence you want to analyze; ex: 127978136
start = as.numeric(args[1])

# end of the sequence you want to analyze; ex: 129978136
end = as.numeric(args[2])

# chromosome number; ex: "chr8" or "8" depending on if you used plink
chr = as.character(args[3])

# reference genome; ex: "hg38"
refgenome = as.character(args[4])

# filepath to the vcf file of the entire chromosome of interest; ex:        "/home/directory/file"
fp.vcf = as.character(args[5])

# filepath to the bgzipped vcf file of the entire chromosome of interest
#fp.zipped = paste0(fp.vcf, ".bgz")
fp.zipped = as.character(args[11])

# filepath to the tabix indexed file of the bgzipped vcf file of the entire chromosome of interest
#fp.indexed.tabix = paste0(fp.zipped, ".tbi")
fp.indexed.tabix = as.character(args[12])

fp.tabix <- TabixFile(fp.zipped, fp.indexed.tabix)

# filepath to the PED file; ex: "/home/directory/file"
fp.ped = as.character(args[6])

# filepath for beagle; ex: "/home/directory/file"
filepath.beagle4 = as.character(args[7])

# specify buffer window size; ex: 1000
buffer.window.size = as.numeric(args[8])

# filepath of vcf for subsetted genome region (including buffer) that you would like to phase; ex: "/home/directory/file"
fp.vcf.to.phase = as.character(args[9])

# filepath for phased vcf file; ex: "/home/directory/file"
fp.phased.vcf = as.character(args[10])


dl.beagle4<-paste0("wget -O ",filepath.beagle4," https://faculty.washington.edu/browning/beagle/beagle.r1399.jar")
system(dl.beagle4)

rng_rvtdt <- GRanges(seqnames = chr,
                     ranges=IRanges(
                       start = start - buffer.window.size ,
                       end = end + buffer.window.size)
)

vcf.rng_rvtdt = readVcf(fp.tabix, refgenome, param = rng_rvtdt)

writeVcf(vcf.rng_rvtdt, fp.vcf.to.phase)

phase.command = paste0("java -Xmx10000m -jar ", filepath.beagle4,
                       " gt=",fp.vcf.to.phase,
                       " ped=",fp.ped,
                       " out=",fp.phased.vcf)
phase.command

system(phase.command)
system(paste0("gunzip ", fp.phased.vcf))

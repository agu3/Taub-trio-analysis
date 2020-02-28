library(VariantAnnotation)
library(trio)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(dplyr)
library(ggplot2)
library(rvTDT)

# start of the sequence you want to analyze; enter as vector with one element
start = 

# end of the sequence you want to analyze; enter as vector with one element
end = 

# chromosome number 
chr = 

# reference genome 
refgenome = 
  
# filepath to the vcf file of the entire chromosome of interest
fp.vcf = 

# filepath to the bgzipped vcf file of the entire chromosome of interest
fp.zipped = paste0(fp.vcf, ".bgz")

# filepath to the tabix indexed file of the bgzipped vcf file of the entire chromosome of interest
fp.indexed.tabix = paste0(fp.zipped, ".tbi")

# filepath to the PED file 
fp.ped = 

# filepath for the vcf file that has been subsetted to the region of interest
fp.filtered.vcf = 
  
# filepath for the results of the gTDT
gTDT.results.fp = 

# filepath for the Manhattan plot generated from the gTDT
gTDT.Manhattan = 


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
  

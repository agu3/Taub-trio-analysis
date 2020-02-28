devtools::install_github("lindagai/rvtrio")
library(rvtrio)

# filepath for beagle 
filepath.beagle4 = 

# specify buffer window size
buffer.window.size = 
  
# filepath of vcf for subsetted genome region (including buffer) that you would like to phase
fp.vcf.to.phase = 
  
# filepath for phased vcf file 
fp.phased.vcf = 

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

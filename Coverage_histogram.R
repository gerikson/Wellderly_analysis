
Total_coverage <- read.csv("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/MedianCovereage_wellderly/total_coverage.csv", header=FALSE)
cov <- Total_coverage$V2
Coverage <- cov[cov < 200]

tiff("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/MedianCovereage_wellderly/cov.tiff", width=2000, height=2000, res=300, compression="lzw")
hist(Coverage)
dev.off()
q()

whiteOnly_wg_and_exome <- read.delim("~/Desktop/Coverage/CG_coverage/whiteOnly_wg_and_exome.txt", header=FALSE)
hist(whiteOnly_wg_and_exome$V3, main = paste("Whole genome coverage"), xlab="Coverage")
hist(whiteOnly_wg_and_exome$V4, main = paste("Exome coverage), xlab="Coverage")


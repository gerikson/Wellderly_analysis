#in R
dat <- read.table("final_vcf_allChrom_snps_AF0.01.pheno",F)
well <- dat[dat$V3 == 2,1:2]
inova <- dat[dat$V3 == 1,1:2]
write.table(well, "ids-wellderly.ids", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(inova, "ids-inova.ids", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")


/gpfs/home/nwineing/plink --bfile final_vcf_allChrom_snps_AF0.01 --keep ids-wellderly.ids --freq --out final_vcf_allChrom_snps_AF0.01-wellderly
/gpfs/home/nwineing/plink --bfile final_vcf_allChrom_snps_AF0.01 --keep ids-inova.ids --freq --out final_vcf_allChrom_snps_AF0.01-inova


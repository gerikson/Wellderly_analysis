install.packages("qqman")
library(qqman)

#Manhattan plot
cat test-results.assoc.logistic | awk '{print $2,"\t",$1,"\t",$3,"\t",$9}' >gwasResults.txt
#Extract the hw variants e-5 and AF in either cohort to be >0.05 and <0.95
python remove_hwe.py

#in R
results <- read.table("gwasResults.txt",T)
#Extract the ones with AF>0.05 (leave the LD variants)
#snp_names <- read.table("final_vcf_allChrom_snps_AF0.01-0.05_MAF.bim.short", F)

#snp_names <- read.table("final_association.txt")

pruned_results <- results[results$SNP %in% snp_names$V1,]
f1_results <- pruned_results[!is.na(pruned_results$P),]

#Extract wellderly AF and inova AF <0.01
well_freq <-read.table("Wellderly_AF.frq",T)
good_well_freq <- well_freq[well_freq$MAF>=0.05,]
ino_freq <-read.table("Inova_AF.frq",T)
good_ino_freq <- ino_freq[ino_freq$MAF>=0.05,]
f2_results <- f1_results[f1_results$SNP %in% good_well_freq$SNP]
f2_results <- f1_results[f1_results$SNP %in% good_well_freq$SNP]

max_results <- final_results[final_results$P<0.00001,]
max_results_plus <- snp_names[snp_names$V1 %in% max_results$SNP,]
write.table(max_results_plus,"smallest_pvalues.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

tiff("manhattan.tiff", width=2000, height=2000, res=300, compression="lzw")
manhattan(final_results)
dev.off()
results <- read.table("final_vcf_allChrom_snps_AF0.01-results.assoc.logistic.short", T)
snp_names <- read.table("final_vcf_allChrom_snps_AF0.01-0.05_MAF.bim.short", F)
pruned_results <- results[results$SNP %in% snp_names$V1,]
final_results <- pruned_results[!is.na(pruned_results$P),]
max_results <- final_results[final_results$P<0.000001,]
write.table(max_results,"smallest_pvalues.txt")
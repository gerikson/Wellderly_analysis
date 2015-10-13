results <- read.table("final_vcf_allChrom_snps_AF0.01-results.assoc.logistic.short", T)
snp_names <- read.table("final_vcf_allChrom_snps_AF0.01-0.05_MAF.bim.short", F)
pruned_results <- results[results$SNP %in% snp_names$V1,]
final_results <- pruned_results[!is.na(pruned_results$P),]

# create QQ plot (same R session) #
obs_P <- -log10(final_results$P[order(final_results$P, decreasing=FALSE)])
exp_P <- -log10(c(1:length(obs_P))/length(obs_P))

library(ggplot2)
p <- ggplot() +
    geom_point(aes(exp_P, obs_P)) +
    geom_line(aes(c(0,max(exp_P)), c(0,max(exp_P))), color="red") +
    ylim(c(0,max(exp_P)))

tiff("test-qq_plot.tiff", width=2000, height=2000, res=300, compression="lzw")
p
dev.off()
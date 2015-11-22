#results <- read.table("final_association.txt", T)
results <- read.table("Wellderly_final_association.txt", T)
#filter by LD
snps_ld <- read.table("test-MAF_LD_pruned.bim")
#keep only the snps in ld
pruned_results <- results[results$SNP %in% snps_ld$V2,]
final_results <- pruned_results[!is.na(pruned_results$P),]
#> dim(final_results)
#[1] 518609      9
obs_P <- -log10(final_results$P[order(final_results$P, decreasing=FALSE)])
exp_P <- -log10(c(1:length(obs_P))/length(obs_P))

library(ggplot2)
p <- ggplot() +
    geom_point(aes(exp_P, obs_P), size = 0.8) +
    geom_line(aes(c(0,max(exp_P)), c(0,max(exp_P))), color="red") +
    ylim(c(0,max(exp_P))) +
    xlab("Expected (-logP)")+ 
    ylab("Observed (-logP)") +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))
tiff("Figure 2B - QQ.tiff", width=3, height=3, res=300, units = 'in', compression="lzw")
p
dev.off()


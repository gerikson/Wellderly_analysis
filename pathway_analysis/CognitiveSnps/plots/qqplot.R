pvals <- read.table("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/cognitive_snps/Cognitive_snps_P_Vals.txt", header=T)

observed <- sort(pvals$PVAL)
lobs <- -(log10(observed))

expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))



pdf("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/cognitive_snps/qqplot_CognitiveSNPs.pdf", width=6, height=6)
plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
dev.off()
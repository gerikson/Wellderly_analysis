results <- read.table("Wellderly_final_association.txt", T)
snps_ld <- read.table("test-MAF_LD_pruned.bim")
pruned_results <- results[results$SNP %in% snps_ld$V2,]
final_results <- pruned_results[!is.na(pruned_results$P),]
head(final_results)
p_values <- final_results$P
data <- qchisq(p_values, 1, lower.tail = FALSE)
data <- sort(data)
ppoi <- ppoints(data)
ppoi <- sort(qchisq(ppoi, df = 1, lower.tail = FALSE))
out <- list()
s <- summary(lm(data ~ 0 + ppoi))$coeff
out$estimate <- s[1, 1]
out$se <- s[1, 2]
out$estimate <- median(data, na.rm = TRUE)/qchisq(0.5, 1)
out$estimate
[1] 0.9718326


#See if p-value changes based on the low inflation factor
smallest_pval <-min(small_pvalues$V4)
smallest_pval
[1] 6.081e-07
chisq <- qchisq(smallest_pval)
div <- chisq/0.9718326
res_pVal <- pchisq(div,1)
res_pVal
[1] 6.168496e-07
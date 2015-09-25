# create plink format file #
/gpfs/home/nwineing/plink --vcf wellderly_inova.VQHIGH.0.95white.chr19.vcf.gz --double-id --vcf-half-call m --make-bed --out test

# name markers according to position #
awk '{print $1, $1"-"$4"-"$5"-"$6, $3, $4, $5, $6}' < test.bim > test.temp
rm test.bim
mv test.temp test.bim

## ANY OTHER FILTERS (biallelic) ##
## ASSSOCIATION TESTING ON WHATEVER FILES CREATED AFTER FILTERS HERE ##

##### PCA #####

# filter by MAF #
/gpfs/home/nwineing/plink --bfile test --maf 0.05 --make-bed --out test-0.05_MAF

# filter by LD #
/gpfs/home/nwineing/plink --bfile test-0.05_MAF --indep-pairwise 50 5 0.5 --out test-LD
/gpfs/home/nwineing/plink --bfile test-0.05_MAF --extract test-LD.prune.in --make-bed --out test-MAF_LD_pruned

# calculate PCs (first 10) #
/gpfs/home/nwineing/plink --bfile test-MAF_LD_pruned --pca 10 --out test-PCA

### plot PCs in R ###
library(ggplot2)
dat <- read.table("test-PCA.eigenvec", F)
# first 519 are wellderly, remaining are inova #
COLOR <- c(rep("red", times=519), rep("blue", times=dim(dat)[1]-519))
p <- ggplot() +
    geom_point(aes(dat$V3, dat$V4), col=COLOR) +
    xlab("pc1") +
    ylab("pc2")
tiff("pca.tiff", width=2000, height=2000, res=300, compression="lzw")
p
dev.off()
q()


##### ASSOCIATION TESTING #####

# create phenotype file in R #
# wellderly = "affected" #
dat <- read.table("test-MAF_LD_pruned.fam", F)
pheno <- c(rep(2, times=519), rep(1, times=dim(dat)[1]-519))
out <- data.frame(dat$V1, dat$V2, pheno)
write.table(out, "test.pheno", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
q()

# perform association testing (logistic) using pcs as covariate #
plink --bfile test --logistic --pheno test.pheno --covar test-PCA.eigenvec --allow-no-sex --hide-covar --out test-results


##### POST-ASSOCIATION PRUNING (for QQ plot) #####
##### THIS SCRIPT ASSUMES test-0.05_MAF IS THE FINAL FILTERED DATASET SET #####
##### OTHER FILTERS WILL BE USED (and a different dataset) IN THE FINAL ANALYSIS #####
##### filter by LD for QQ plot, but not for results #####

# reduce asssociation testing output #
awk '{ print $2, $9}' < test-results.assoc.logistic > test-results.assoc.logistic.short

# reduce marker names output #
awk '{ print $2}' < test-0.05_MAF.bim > test-0.05_MAF.bim.short         # this is the filtered marker list

# prune in R #
results <- read.table("test-results.assoc.logistic.short", T)
snp_names <- read.table("test-0.05_MAF.bim.short", F)
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









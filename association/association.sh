i# create plink format file #
/gpfs/home/nwineing/plink --vcf /gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_allFilters_36kmer_snpsOnly_AF0.01/final_vcf_nokmer_snps_AF0.01.noRelated.vcf.gz --double-id --vcf-half-call m --make-bed --out test

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

#Hardy-Weinberg 1e-4
/gpfs/home/nwineing/plink --bfile final_vcf_allChrom_snps_AF0.01-0.05_MAF --hwe 0.0001 --write-snplist

#Hardy-Weinberg for case-control
gpfs/home/nwineing/plink --bfile final_vcf_allChrom_snps_AF0.01-0.05_MAF --pheno final_vcf_allChrom_snps_AF0.01.pheno --allow-no-sex –hardy

# calculate PCs (first 10) #
/gpfs/home/nwineing/plink --bfile test-MAF_LD_pruned --pca 10 --out test-PCA

### plot PCs in R ###
library(ggplot2)
dat <- read.table("final_vcf_allChrom_snps_AF0.01-PCA.eigenvec", F)
# first 519 are wellderly, remaining are inova #
COLOR <- c(rep("blue", times=511), rep("red", times=dim(dat)[1]-511))
p <- ggplot() +
    geom_point(aes(dat$V3, dat$V4), col=COLOR) +
    xlab("pc1") +
    ylab("pc2")
tiff("pca.tiff", width=2000, height=2000, res=300, compression="lzw")
p
dev.off()
q()

COLOR <- c(rep("bue", times=dim(dat)[1]-511), rep("red", times=511))
p <- ggplot() +
    geom_point(aes(dat$V3, dat$V4), col=COLOR) +
    xlab("pc1") +
    ylab("pc2")
tiff("pca.tiff", width=2000, height=2000, res=300, compression="lzw")
p
dev.off()
q()

### do IBD

plink --file mydata --genome --min 0.125  (http://pngu.mgh.harvard.edu/~purcell/plink/ibdibs.shtml)

##### ASSOCIATION TESTING #####

# create phenotype file in R #
# wellderly = "affected" #
dat <- read.table("final_vcf_allChrom_snps_AF0.01-MAF_LD_pruned.fam", F)
pheno <- c(rep(2, times=511), rep(1, times=dim(dat)[1]-511))
out <- data.frame(dat$V1, dat$V2, pheno)
write.table(out, "final_vcf_allChrom_snps_AF0.01.pheno", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
q()

# perform association testing (logistic) using pcs as covariate #
plink --bfile final_vcf_allChrom_snps_AF0.01 --logistic --pheno final_vcf_allChrom_snps_AF0.01.pheno --covar final_vcf_allChrom_snps_AF0.01-PCA.eigenvec --allow-no-sex --hide-covar --out final_vcf_allChrom_snps_AF0.01-results

##### POST-ASSOCIATION PRUNING (for QQ plot) #####
##### THIS SCRIPT ASSUMES test-0.05_MAF IS THE FINAL FILTERED DATASET SET #####
##### OTHER FILTERS WILL BE USED (and a different dataset) IN THE FINAL ANALYSIS #####
##### filter by LD for QQ plot, but not for results #####

# reduce asssociation testing output #
awk '{ print $2, $9}' < final_vcf_allChrom_snps_AF0.01-results.assoc.logistic > final_vcf_allChrom_snps_AF0.01-results.assoc.logistic.short

# reduce marker names output #
awk '{ print $2}' < final_vcf_allChrom_snps_AF0.01-0.05_MAF.bim > final_vcf_allChrom_snps_AF0.01-0.05_MAF.bim.short         # this is the filtered marker list

# prune in R #
#results <- read.table("final_vcf_allChrom_snps_AF0.01-results.assoc.logistic.short", T)
#snp_names <- read.table("final_vcf_allChrom_snps_AF0.01-0.05_MAF.bim.short", F)
#Keep only the snps that passed Hardy-Weinberg 1e-5
#snp_names <- read.table("plink.snplist", F)



pruned_results <- results[results$SNP %in% snp_names$V1,]
final_results <- pruned_results[!is.na(pruned_results$P),]

#Extract smallest p-values (<1e-5)
max_results <- final_results[final_results$P<0.00001,]
write.table(max_results,"smallest_pvalues.txt")

# create QQ plot (same R session) #
obs_P <- -log10(final_results$P[order(final_results$P, decreasing=FALSE)])
exp_P <- -log10(c(1:length(obs_P))/length(obs_P))

library(ggplot2)
p <- ggplot() +
    geom_point(aes(exp_P, obs_P)) +
    geom_line(aes(c(0,max(exp_P)), c(0,max(exp_P))), color="red") +
    ylim(c(0,max(exp_P)))

tiff("qq_plot.tiff", width=2000, height=2000, res=300, compression="lzw")
p
dev.off()

#Manhattan plot
cat final_vcf_allChrom_snps_AF0.01-results.assoc.logistic | awk '{print $2,"\t",$1,"\t",$3,"\t",$9}' >gwasResults.txt
#in R
#results <- read.table("gwasResults.txt",T)
library(qqman)
#Read in the snps with HWE excluded
python remove_hwe.py
results <- read.table("final_association.txt",T)
#Extract the ones with AF>0.01
snp_names <- read.table("final_vcf_allChrom_snps_AF0.01-0.05_MAF.bim.short", F)
pruned_results <- results[results$SNP %in% snp_names$V1,]
final_results <- pruned_results[!is.na(pruned_results$P),]
max_results <- final_results[final_results$P<0.00001,]
write.table(max_results,"smallest_pvalues.txt")

tiff("manhattan.tiff", width=2000, height=2000, res=300, compression="lzw")
manhattan(final_results)
dev.off()




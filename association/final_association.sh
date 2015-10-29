#Final association steps
#From /gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/final_assoc/

#Set VQLOW as missing after all filters
#Generate plink files
/gpfs/home/nwineing/plink --vcf /gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_final_allFilters_noVQLOW_snpsOnly_AF0.01/final_vcf_allChrom_snps_AF0.01.vcf.gz --double-id --vcf-half-call m --make-bed --out test

awk '{print $1, $1"-"$4"-"$5"-"$6, $3, $4, $5, $6}' < test.bim > test.temp

rm test.bim
mv test.temp test.bim

#See if there are any individuals with more then 10% missing rate
/gpfs/home/nwineing/plink --bfile test --mind 0.1 --make-bed --out cleaned


#Extract AF > 0.05 and missing rate < 0.1
/gpfs/home/nwineing/plink --bfile test --maf 0.05 --geno 0.1 --make-bed --out test-0.05_MAF
#After 0.95% white, no related
#1197 people: 511 wellderly, 686 inova
#6997503 start snps only AF>0.01 after all filters
#3844648 test-0.05_MAF.bim -->without missing
#3511 variants removed due to missing genotype data (--geno).
#3152416 variants removed due to MAF threshold(s) (--maf/--max-maf).
#3841576 variants and 1197 people pass filters and QC.
#0 people removed due to missing genotype data (--mind).
#3747785 final_association.txt

#Hardy-Weinberg for case-control, this generate a file will ALL hwe
/gpfs/home/nwineing/plink --bfile test-0.05_MAF -hardy --allow-no-sex --pheno final_vcf_allChrom_snps_AF0.01.pheno

#Extracting the ones that are <1e-4 in either wellderlt, inova or all


#Extract the variants in LD
/gpfs/home/nwineing/plink --bfile test-0.05_MAF --indep-pairwise 50 5 0.5 --out test-LD
#Prune the variants in LD
/gpfs/home/nwineing/plink --bfile test-0.05_MAF --extract test-LD.prune.in --make-bed --out test-MAF_LD_pruned

# calculate PCs (first 10) #
/gpfs/home/nwineing/plink --bfile test-MAF_LD_pruned --pca 10 --out test-PCA

### plot PCs in R ###
library(ggplot2)
dat <- read.table("test-PCA.eigenvec", F)

#First 511 are wellderly in red up top, inova in blue behind
COLOR <- c(rep("blue", times=dim(dat)[1]-511), rep("red", times=511))
p <- ggplot() +
    geom_point(aes(dat$V3, dat$V4), col=COLOR) +
    xlab("pc1") +
    ylab("pc2")
p
dev.off()
q()

# create phenotype file in R #
# wellderly = "affected" #
dat <- read.table("final_vcf_allChrom_snps_AF0.01-MAF_LD_pruned.fam", F)
pheno <- c(rep(2, times=511), rep(1, times=dim(dat)[1]-511))
out <- data.frame(dat$V1, dat$V2, pheno)
write.table(out, "final_vcf_allChrom_snps_AF0.01.pheno", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
q()

#Run the association
plink --bfile test-0.05_MAF --logistic hide-covar --pheno final_vcf_allChrom_snps_AF0.01.pheno --covar test-PCA.eigenvec --allow-no-sex --out test-results

#Reduce the association result
awk '{ print $2, $9}' < test-results.assoc.logistic > test-results.assoc.logistic.short
# reduce marker names output #
awk '{ print $2}' < test-0.05_MAF.bim > test-0.05_MAF.bim.short  


#Extract AF by population
#in R
dat <- read.table("final_vcf_allChrom_snps_AF0.01.pheno",F)
well <- dat[dat$V3 == 2,1:2]
inova <- dat[dat$V3 == 1,1:2]
write.table(well, "ids-wellderly.ids", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(inova, "ids-inova.ids", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

/gpfs/home/nwineing/plink --bfile test-0.05_MAF --keep ids-wellderly.ids --freq --out Wellderly_AF
/gpfs/home/nwineing/plink --bfile test-0.05_MAF --keep ids-inova.ids --freq --out Inova_AF

#Remove the variants <0.05 in either inova or wellderly and the hwe <1e-4 in either cohort or both
#Remove hwe <1e-4 in either cohort or both
python ~/scripts/Wellderly_scripts/GitHub/association/remove_hwe.py

#Manhattan plot
install.packages("qqman")
library(qqman)

#in R
#Load in final association with everything filtered out
snp_names <- read.table("final_association.txt",T)
#Extract N/As
final_results <- snp_names[!is.na(snp_names$P),]

#Save to file smallest p-values
max_results <- final_results[final_results$P<0.00001,]
write.table(max_results,"smallest_pvalues.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

tiff("manhattan.tiff", width=4000, height=2000, res=300, compression="lzw")
#manhattan(final_results, col = c("red1", "green2","blue2","cyan","blueviolet","yellow","snow3","gray0"))
manhattan(final_results, col = c("red1","blue2","darkorange1","blueviolet","green2","gray0"))
dev.off()

#Make QQ plot
#in R

results <- read.table("final_association.txt", T)
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
    geom_point(aes(exp_P, obs_P)) +
    geom_line(aes(c(0,max(exp_P)), c(0,max(exp_P))), color="red") +
    ylim(c(0,max(exp_P)))

tiff("qq_plot.tiff", width=2000, height=2000, res=300, compression="lzw")
p
dev.off()



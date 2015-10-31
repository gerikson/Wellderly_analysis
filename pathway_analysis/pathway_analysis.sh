#Extract/combine transcripts by Gene add 100kb to begin end, split 
#genes by chrom
#From pathway_analysis folder
python split_UCSC_geneByChrom.py

#Generate bim file with the snps that passed ALL filters
#in R from /gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/final_assoc

#snp_names <- read.table("test-0.05_MAF.bim")
#assoc <- read.table("final_association.txt",T)
#pruned_snps <- snp_names[snp_names$V2 %in% assoc$SNP,]
#write.table(pruned_snps, "final.bim")

#Cleaning bim file remove first column, header and '"'
#cat final.bim | awk '{print $2,"\t",$3,"\t",$4,"\t",$5,"\t",$6,"\t",$7}' >final_fixed.bim

#Extract the variants in bim file not found in the association file,
#will use plink to remove them

assoc <- read.table("final_association.txt",T)
snp_names <- read.table("test-0.05_MAF.bim")
pruned_snps <- snp_names[!(snp_names$V2 %in% assoc$SNP),]
write.table(pruned_snps,"snps_to_be_removed.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
cat snps_to_be_removed.txt | awk '{print $2}' >remove_snps.txt

#Make new plink files
/gpfs/home/nwineing/plink --bfile test-0.05_MAF --exclude remove_snps.txt --make-bed
#--exclude: 3799647 variants remaining.

#Add genes to the bim file
python add_gene_to_bim.py

#Run the association with corrected file (make sure no problems)
plink --bfile plink --logistic hide-covar --pheno plink.pheno --covar plink-PCA.eigenvec --allow-no-sex --out results

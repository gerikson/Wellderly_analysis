PATHWAY ANALYSIS
1) Extracted the kgXref genes and mapped them to the UCSC known genes using the mysql database:
mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -N -A -e 'select kgXref.kgID, kgXref.geneSymbol,knownGene.name,knownGene.chrom,knownGene.txStart,knownGene.txEnd from kgXref, knownGene where knownGene.name=kgXref.kgID' >genes_positions.txt

2) From kgXref file kept only the genes that have valid entries in the 5th column and mapped it to the UCSC Know Genes:
Unique genes in kgXref file:
zcat kgXref.txt.gz | awk -F "\t" '{print $5}' | sort -u | wc -l
28517

3)Removed the genes that have same names, are on the same chr but don’t overlap even when adding 100kb to begin and end

4) Assigned each snp to the closest gene (within 100kb of begin or end of the snp)

5) If snp is inside multiple genes that overlapped assigned the gene that has the begin closest to the snp
6) Generated 10k *.pheno files with the casa/control shuffled
7) Ran association (with covariates) on all of the files
8) For each simulation file extract the lowest p-values by gene
9) Per each gene calculated a New_p-value that is the number of simulations that have lower p-value then the real association + 1 divided by 10001
10) Extracted -log10(New_p_value)
11) If there are still same gene name on different chroms, kept only the gene with the lowest p-value
12) Uploaded in GSEA the *.rnk file that has gene_name + “\t” + [-log10(New_p_value)]
13) Run GSEA on a Pre-Ranked gene list with gene set databases: Hallmarks, Reactome, GO: biological processes




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

#Start 10k jobs to generate random *.pheno and run association with simulated Pheno files
python generate_pheno_files.py

#Restart the jobs that failed
#From directory check file sizes
find -size 348M >jobs.succesfully.complete
wc -l jobs.succesfully.completed
10000 jobs.succesfully.completed
python restart_failed_jobs.py

#Extract gene names and index of the variants
# in plink files where the gene is located
python Extract_gene_from_plink.py >no_assigned_gene.txt
sort -k3n,3 gene_names_coordinates_plink.txt >gene_names_coordinates_plink.sorted.txt

3799647 plink.withGene.bim
881747 no_assigned_gene.txt

#Extract smallest p-value by gene:
python Extract_smallest_p_value_byGene.py

#Real data smallest p-values by gene
/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/backup_plink_files/REAL_DATA_smallest_p_value_per_gene.txt


#Start all simulation jobs, extract minimum p-value
python Create_jobs_extract_min_perGene_simulation.py

#Combine simulations
python combine_simulations.py

#Remove duplicate genes
python Remove_duplicate_genes.py

#Run GSEA 
Run GSEA on a Pre-Ranked gene list with gene set databases: Hallmarks, Reactome, GO: biological processes

#SECOND TRY!!!! ONE GENE PER SNP
#Create list of valid genes, verify that the gene is present in the
#kgXref file 5th column:
python ./Reasign_genes/Create_list_of_valid_genes.py

#ADD genes to bim file:
python add_gene_to_bim.py

#Reasign genes only 1 gene per snp (for the overlapping genes add the gene
#whose begin is closest to the gene)
python ./Reasign_genes/genes_reasignment.py

#Sort files
/gpfs/home/gerikson/wellderly/resources/
sort -k 1n,1 -k 3n,3 gene_names_coordinates_plink_final.txt > gene_names_coordinates_plink_final.sorted.txt

#Still some issues, for the overlapping genes the snps gets the gene that starts first
#Fixed that issue and we have
#22883 gene_names_coordinates_plink_final.txt
#vs before 
#21047 gene_names_coordinates_plink_final.txt_v1
#When we were assigning each snp to ALL genes we had
#26228 gene_names_coordinates_plink.sorted.txt

#Extract minimum p-value per simulation
python ./Reassign_genes/Create_jobs_extract_min_perGene.py

#Extract smallest p-value real data
python ./Reasign_genes/Extract_smallest_pvalue_real_data.py

#combine simulations
python ./Reasign_genes/combine_simulations.py
22883 FINAL_RESULTS_p_values

#Remove duplicate genes
python ./Reasign_genes/Remo.py
22845 FINAL_FINAL_RESULTS_negLog10

#Filter out the variants in LD
python ./No_LD/Gene_reasignment_noLD.py

python ./No_LD/Create_jobs_extract_min_perGene.py



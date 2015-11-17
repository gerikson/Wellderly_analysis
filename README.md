# Wellderly_analysis

* Extract Wellderly genotypes only, remove variants that are not found in the wellderly, transform to vcf

python Create_job_extract_wellderly_vcf.py

* Parse to vcf all of the data

python Create_jobs_parse_genomeComb.py

* Sort wellderly vcf

python Create_jobs_sort_vcf.py

* Extract the varinats with at least one VQHIGH in white individuals

python remove_vqlow.py

*Extract variants that are clustered in >0.1 wellderly or inova

python create_jobs_remove_clustered.py

*Extract the variants with AF >0.01
python create_jobs_0.01AF.py  ---> DID'T work with vcftools, will have to do it manually

*Extract the repeats, homopolymers, etc
python Create_jobs_extractRepeats_etc.py

*Count the number of VQHIGH passed filters by AF
python create_jobs_count_totalVQHIGH_byAF.py

*Count rest of the filters (in the Count_filters folder)

*Remove variants with >10% missing in either wellderly or inova
python Create_jobs_extract_missing.py

*Remove variants with coverage <10 or >100
python  create_jobs_remove_coverage.py

* Extract snp position based on rsID
python ./snps_of_interest/Extract_position_of_snp.py

*Extract the snps of interest
python Exract_snps_of_interest.py


*Extracting separatelly snps and delins with AF>0.01
python Extract_snpsOnly_AFmoreThen0.01.py

*Concatenating the vcf file by chrom into a final one:
vcf-concat vcf_snps_AF0.01.chr1.vcf.gz vcf_snps_AF0.01.chr2.vcf.gz vcf_snps_AF0.01.chr3.vcf.gz vcf_snps_AF0.01.chr4.vcf.gz vcf_snps_AF0.01.chr5.vcf.gz vcf_snps_AF0.01.chr6.vcf.gz vcf_snps_AF0.01.chr7.vcf.gz vcf_snps_AF0.01.chr8.vcf.gz vcf_snps_AF0.01.chr9.vcf.gz vcf_snps_AF0.01.chr10.vcf.gz vcf_snps_AF0.01.chr11.vcf.gz vcf_snps_AF0.01.chr12.vcf.gz vcf_snps_AF0.01.chr13.vcf.gz vcf_snps_AF0.01.chr14.vcf.gz vcf_snps_AF0.01.chr15.vcf.gz vcf_snps_AF0.01.chr16.vcf.gz vcf_snps_AF0.01.chr17.vcf.gz vcf_snps_AF0.01.chr18.vcf.gz vcf_snps_AF0.01.chr19.vcf.gz vcf_snps_AF0.01.chr20.vcf.gz vcf_snps_AF0.01.chr21.vcf.gz vcf_snps_AF0.01.chr22.vcf.gz | gzip -c >final_vcf_allChrom_snps_AF0.01.vcf.gz


*Run the first step of the association on all data:
python create_job_association_final.py

*Second step:
python run_association.py

*Association shows biggest p-values in repeat regions, removing all of them --> skiped this for final analysis:
python ./association/create_jobs_remove_ALLrepeats_association.py

*Extracted the related individuals <-- This doesn't work
/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_AFmore0.01> vcftools --gzvcf final_vcf_allChrom_snps_AF0.01.vcf.gz --remove eliminate_individuals.txt --out final_vcf_allChroms_snps_AD0.01_noRelated.vcf.gz

*Extract related better version
python Excluded_related.py

*Concatenate files
python concatenate.py

*Extracted inova coverage

*Concatenating v1 (still need to concatenate chr1-chr4)
vcf-concat final_vcf_nokmer_snps_AF0.01.noRelated.chr5.vcf.gz final_vcf_nokmer_snps_AF0.01.noRelated.chr6.vcf.gz final_vcf_nokmer_snps_AF0.01.noRelated.chr7.vcf.gz final_vcf_nokmer_snps_AF0.01.noRelated.chr8.vcf.gz final_vcf_nokmer_snps_AF0.01.noRelated.chr9.vcf.gz final_vcf_nokmer_snps_AF0.01.noRelated.chr10.vcf.gz final_vcf_nokmer_snps_AF0.01.noRelated.chr11.vcf.gz final_vcf_nokmer_snps_AF0.01.noRelated.chr12.vcf.gz final_vcf_nokmer_snps_AF0.01.noRelated.chr13.vcf.gz final_vcf_nokmer_snps_AF0.01.noRelated.chr14.vcf.gz final_vcf_nokmer_snps_AF0.01.noRelated.chr16.vcf.gz final_vcf_nokmer_snps_AF0.01.noRelated.chr17.vcf.gz final_vcf_nokmer_snps_AF0.01.noRelated.chr18.vcf.gz final_vcf_nokmer_snps_AF0.01.noRelated.chr19.vcf.gz final_vcf_nokmer_snps_AF0.01.noRelated.chr20.vcf.gz final_vcf_nokmer_snps_AF0.01.noRelated.chr21.vcf.gz final_vcf_nokmer_snps_AF0.01.noRelated.chr22.vcf.gz | gzip -c >final_vcf_nokmer_snps_AF0.01.noRelated.temp.vcf.gz


*Extract snps of interest from vcf file
python Extract_snps_of_interest.vcf.py

*Extract coveage by individual snps of interest:
python Extract_coverage_SnpSOfInterest_wellderly.py

*Calculate AF and median coverage (from median coverage file) SNPs of interest, missing geno and VQHIGH for both wellderly and inova:
python Calculate_AF_welldVsInova.py

* Match p-values from association based on location
final_association.sh

*Extract rare variants
python Create_jobs_extract_rareVariants.py

*Extract AF, p-value (from association) snps of interest 
python Extract_AF_p-values.py

ASSOCIATION
./association/final_association.sh


PATHWAY ANALYSIS
FOR pathway analysis, extract genes/positions (in ~/wellderly/resources)
mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -N -A -e 'select kgXref.kgID, kgXref.geneSymbol,knownGene.name,knownGene.chrom,knownGene.txStart,knownGene.txEnd from kgXref, knownGene where knownGene.name=kgXref.kgID' >genes_positions.txt

* Extract genes with 100kb interval (from pathway analysis folder)
python split_UCSC_geneByChrom.py

*Add gene to bim file (not needed in the end)
python add_gene_to_bim.py

*Generate 10k *.pheno files and run the simulation
python generate_pheno_files.py


*For pathway analysis read ./pathway_analysis/pathway_analysis.sh

*Extract all filters:
python Create_jobs_apply_all_filters.py

FOR Rare variants
*Extracting ALL clustered variants:
python ./Rare_variant_analysis/Create_jobs_remove_ALL_clustered.py




*Extracting snp position for cognitive snps
python ./pathway_analysis/CognitiveSnps/Extract_snp_position.py

#Extract cognitive snps p-values from all 10k simulations
python ./pathway_analysis/CognitiveSnps/Create_jobs_extract_sim_pvalues.py


*Count 36mers:
python ./Count_filters/Create_jobs_count_36mers.py








* Apply the filters: missing/uncertain genotype > 10 perc in either wellderly or inova, covereage <10 or >100, whites only (testing 0.85 white and 0.95)

python Extract_white_filter.py

* Run the plink analysis, maf > 0.01, in LD

python data_analysis.py

* Add filters to the association file (reapeat, homopolymer, segDup, Microsat)

python Add_filters.py

* Test different filters to see which one work better
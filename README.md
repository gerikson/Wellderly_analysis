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












* Apply the filters: missing/uncertain genotype > 10 perc in either wellderly or inova, covereage <10 or >100, whites only (testing 0.85 white and 0.95)

python Extract_white_filter.py

* Run the plink analysis, maf > 0.01, in LD

python data_analysis.py

* Add filters to the association file (reapeat, homopolymer, segDup, Microsat)

python Add_filters.py

* Test different filters to see which one work better
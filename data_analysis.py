"""
Data Analysis Wellderly

author: gerikson
date: July 28th

* fixed the padding issue from the first version of mcompar_to_vcf
* implemented all of the fixed from fix_vcf.py
    1. From the wellderly data excludes the variants that have no alternate alleles
    2. Excluded the other alles if not found in the wellderly datased
    3. Gets rid of the '@' allele
    4. Sorts output file
* combinning the dels at the same position

"""

import os, sys, gzip, datetime 

WRITE_BUFFER_SIZE = 100

#Work directory for test:
#/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_results_v2/Plink/July29
def run_plink(chrom):
    input_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_withCorrectGenotypes/chr22.filtered.vcf.gz"
    #input_filename= "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_withCorrectGenotypes/"+chrom+".backup.vcf.gz"
    output_path="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/plink_data_analysis/"

    #Transform from vcf file to plink
    os.system("plink --vcf " + input_filename + " --double-id --vcf-half-call m --out "+ output_path+"data_in_plink_format/"+chrom)
    
    #Sort file
    os.system("plink --bfile " + output_path+"data_in_plink_format/"+chrom +" --make-bed --out " + output_path+"bed_files/"+chrom + ".sorted")

    #USE THIS:
    #os.system(plink --bfile chr22.sorted --maf 0.05  --pheno /gpfs/group/stsi/data/projects/wellderly/GenomeComb/plink_data_analysis/wellderly.pheno --assoc --allow-no-sex --out)
    '''
    filter Maximum minor allele frequency 0.01
    filter by LD:the command above that specifies 50 5 0.5 would 
    a) consider a window of 50 SNPs, 
    b) calculate LD between each pair of SNPs in the window, 
    b) remove one of a pair of SNPs if the LD is greater than 0.5, 
    c) shift the window 5 SNPs forward and repeat the procedure.
    '''
    os.system("plink --bfile "+ output_path + "bed_files/"+chrom + ".sorted" + " --max-maf 0.01 --indep-pairwise 50 5 0.5 --out "+ output_path + "pruned_data/" + chrom + ".pruneddata")
    
    #Extract the pruned variants
    os.system("plink --bfile "+ output_path + "bed_files/"+chrom + ".sorted" + " --extract " + output_path + "pruned_data/" + chrom + ".pruneddata.in" + " --make-bed --out " + output_path + "final_pruned_data/" + chrom+".pruneddata.final")
    
    #Run the association
    os.system("plink --bfile "+ output_path + "final_pruned_data/" + chrom+".pruneddata.final" + " --pheno /gpfs/group/stsi/data/projects/wellderly/GenomeComb/plink_data_analysis/wellderly.pheno --assoc --allow-no-sex --out "+ output_path + "association/" + chrom)

    #Add the filters to the associated file



    #plink --bfile chr22.sorted --extract plink.prune.in --make-bed --out chr22.pruneddata

    #This removes 90% of variants
    #plink --bfile chr22.pruneddata --maf 0.01 --make-bed --out chr22.pruneddata.AF

    #Doing the association
    #plink --bfile chr22.pruneddata.AF --pheno /gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_results_v2/Plink/July29/welldery.pheno --assoc --allow-no-sex --out chr22

    #Principal component
    #plink --bfile chr22.pruneddata.AF -pca --out chr22.pca

def main(chrom):
    run_plink(chrom)

if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)


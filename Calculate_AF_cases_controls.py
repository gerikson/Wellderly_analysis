"""
Filter by HWE

"""
import os, sys, gzip, datetime


def main(chrom):

    workingdir = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/Plink_files/"
    input_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/wellderly_all_filters_withHWE"+str(chrom)+".vcf.gz"
    plink_command="/gpfs/home/nwineing/plink --vcf "+input_filename+" --double-id --vcf-half-call m --make-bed --out "+workingdir+chrom
    os.system(plink_command)
    
    make_snpID="awk '{print $1, $4"+'"-"'+"$5, $3, $4, $5, $6}' < "+workingdir+chrom+".bim > "+workingdir+chrom+".temp"
    os.system(make_snpID)
    
    rm_snp = "rm "+workingdir+chrom+".bim"
    os.system(rm_snp)
    os.system("mv "+workingdir+chrom+".temp "+workingdir+chrom+".bim")
    
    pheno_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/final_assoc/final_vcf_allChrom_snps_AF0.01.pheno"
    plink_command_wellderly = "/gpfs/home/nwineing/plink --bfile "+workingdir+chrom+" --filter-cases --freq --pheno "+pheno_file+" --allow-no-sex --out "+ workingdir+chrom+".wellderly"
    os.system(plink_command_wellderly)
    plink_command_inova = "/gpfs/home/nwineing/plink --bfile "+workingdir+chrom+" --filter-controls --freq --pheno "+pheno_file+" --allow-no-sex --out "+ workingdir+chrom+".inova"
    os.system(plink_command_inova)    


    
if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)
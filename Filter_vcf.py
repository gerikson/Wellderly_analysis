"""
Filter vcf file by:


Whites only Missing/Uncertain Genotypes in >10% of Wellderly or >10% of Inova Genome calls. (this is the standard threshold used in PLINK for GWAS) 

Median coverage <10 or >100 in Wellderly or Inova Genomes.

Variants Clustered in >10% of Wellderly or >10% of Inova Genome calls. (the GenomeComb output has a column that indicates whether the variant is clustered or not) 


"""

def main(chrom):
	infile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_withCorrectGenotypes"


if __name__ == '__main__':

    print "Python Version: " + sys.version
    args = parse_args()
    import time
    start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    print 'This took {0} seconds'.format(end - start)
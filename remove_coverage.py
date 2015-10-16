"""
Remove the variants that are bellow 10 coverage or >110

"""
import os, sys, gzip, datetime


def main(chrom):
    ch = chrom[3:]
    print "chrom is: " + ch
    
    #V1 for the wellderly filter
    #input_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_clustered_repeats_homopoly_etc_missing/wellderly_inova.VQHIGH.0.95white.nocluster.repeats.etc.missing."+str(chrom)+".vcf.gz"
    #output_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_clustered_repeats_homopoly_etc_missing_cov/v1_wellderly_inova.VQHIGH.0.95white.nocluster.repeats.etc.missing.cov."+str(chrom)+".vcf.gz"
    #counter_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_clustered_repeats_homopoly_etc_missing_cov/v1_counter_wellderly_inova.VQHIGH.0.95white.nocluster.repeats.etc.missing.cov.txt"
    #coverage_file="/gpfs/group/torkamani/bhuvan/wellderly/coverage/CoverageInfo/MediansCompiled/medians_chrm_"+str(ch)+".csv"
    
    #V2 for the inova filter
    input_filename= "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_allFilters_36kmer_snpsOnly_AF0.01/vcf_nokmer_snps_AF0.01.noRelated." + chrom + ".vcf.gz"
    output_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_allFilters_36kmer_snpsOnly_AF0.01/final_vcf_nokmer_snps_AF0.01.noRelated." + chrom + ".vcf.gz"
    counter_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_allFilters_36kmer_snpsOnly_AF0.01/counter.txt"
    coverage_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/inova_median_coverage/MediansCompiled/medians_chrm_"+str(ch)+".csv"


    o = gzip.open(output_filename, 'w')
    #Create a dictionary of the covereage
    cov_file = open(coverage_file)
    covereage = {}
    print "Start creating dictionary"
    for i in cov_file:
        ln = i.split(",")
        covereage[ln[0]] = ln[1]
    print "Dictionary created"


    c = open(counter_file, "a")
    print 'Calculating...'

    counter = 0
    good_lines = 0
    block = ""
    f = gzip.open(input_filename)

    for line in f:

        if line[:1] == "#":
        	print "header"
        	o.write(line)

        else:

            counter += 1

            if counter%100 == 0:
                o.write(block)
                block = ""
                o.flush()

            if counter%10000 == 0:
                print datetime.datetime.now().time()
                print "total lines " + str(counter)
                print "Good lines " + str(good_lines)
                sys.stdout.flush()

            tp_line = line.strip().split()
            pos = tp_line[1] + "_" + str(len(tp_line[3]))

            #Extract the coverage from the dictionary
            cov = "0.0"
            try:
                cov = covereage[pos]
            except:
                print "coverage not found"
                print tp_line[:5]


            if float(cov) < 10.0 or float(cov) > 110.0:
                continue
            else:
                block = block+line.strip()+"\n"
                good_lines += 1
        	

    o.write(block)

    print "Total lines " + str(counter)
    print "Total good lines " + str(good_lines)

    c.write(chrom + "\t" + str(counter) + "\t" + str(good_lines)+ "\n")
    c.close()
    f.close()
    o.close()


if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)
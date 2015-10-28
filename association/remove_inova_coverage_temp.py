"""
October 27th
This is a temp scripts
Checking if I removed inova covereage correctly from the association file

"""
import os, sys, gzip, datetime


def main():
    ch = "1"
    print "chrom is: " + ch

    input_filename= "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/final_assoc/test-0.05_MAF.bim"
    output_filename= "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/final_assoc/test-0.05_MAF_correct.bim"
    
    counter_file= "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/final_assoc/inova_cov.txt"
    coverage_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/inova_median_coverage/MediansCompiled/medians_chrm_"+str(ch)+".csv"


    o = open(output_filename, 'w')
    #Create a dictionary of the covereage
    cov_file = open(coverage_file)
    covereage = {}
    print "Start creating dictionary"
    for i in cov_file:
        ln = i.split(",")
        covereage[ln[0]] = ln[1]
    print "Dictionary created"
    cov_file.close()

    c = open(counter_file, "a")
    print 'Calculating...'

    counter = 0
    good_lines = 0
    block = ""
    f = open(input_filename)

    for line in f:

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
        if tp_line[0] != ch:
                ch = tp_line[0]
                print "New chrom: " + ch
                coverage_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/inova_median_coverage/MediansCompiled/medians_chrm_"+str(ch)+".csv"
                #Create a dictionary of the covereage
                cov_file = open(coverage_file)
                covereage = {}
                print "Start creating dictionary"
                for i in cov_file:
                    ln = i.split(",")
                    covereage[ln[0]] = ln[1]
                print "Dictionary created"
                cov_file.close()

        pos = tp_line[3] + "_1"

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
    main()
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)
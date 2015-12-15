"""
Count final dataset after all filters

"""

import os, sys, gzip, datetime 

def main(chrom):


	#AF_file = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/FINAL_AF/final_"+chrom+".txt.gz")
	#Before AD filter
	AF_file = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/AF_noAD_filter/final_"+chrom+".txt.gz")
	counter = 0
	counter_unique_no_pass = 0
	good_lines = 0
	
	well_001 = 0
	well_005 = 0
	well_005_plus = 0

	inova_001 = 0
	inova_005 = 0
	inova_005_plus = 0
	part_count = 0
	for line in AF_file:
		counter += 1
		tp_line = line.strip().split("\t")


		counter_unique_no_pass += 1
	
		well_AF = float(tp_line[5])
		inova_AF = float(tp_line[6])

		if well_AF > 0.0:
			if well_AF < 0.01:
			    well_001 += 1
			elif well_AF < 0.05:
			    well_005 += 1
			else:
			    well_005_plus += 1

		if inova_AF > 0.0:
			if inova_AF < 0.01:
			    inova_001 += 1
			elif inova_AF < 0.05:
			    inova_005 += 1
			else:
			    inova_005_plus +=1

	AF_file.close()


	print "Total lines"
	print str(counter)
	print "New Filtered out lines"
	print str(counter_unique_no_pass)
	
	counter_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_count_filters/final/FINAL_counter_after_ALL_filters_exceptAD.txt"
	counterf = open(counter_file, 'a')
	
	AF_counter = chrom + "\t" + str(counter_unique_no_pass) + "\t" + str(well_001) + "\t" + str(well_005) + "\t" + str(well_005_plus) + \
				"\t" + str(inova_001) + "\t" + str(inova_005) + "\t" + str(inova_005_plus) + "\n"
	
	print AF_counter
	counterf.write(AF_counter)
	counterf.close()
	print "DONE!"




if __name__ == '__main__':

    print "Python Version: " + sys.version
    #args = parse_args()
    #import time
    #start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    #end = time.time()
    #print 'This took {0} seconds'.format(end - start)

"""
Count after all filters except AD

"""

import os, sys, gzip, datetime 

def main(chrom):

	
	infile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/sanity_check_wellderly_all_filters_withHWE.noZeroAF"+chrom+".vcf.gz"
	AF_file = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/AF/final_"+chrom+".txt.gz")

	counter_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_count_filters/final/All_filters_except_AD.txt"
	
	counterf = open(counter_file, 'a')


	i = gzip.open(infile)


	#Index the AF file
	AF_dictionary = {}
	for line in AF_file:
		af_array = []
		tp_line = line.strip().split("\t")
		dict_key = tp_line[1]+"_"+tp_line[3]+"_"+tp_line[4]
		af_array.append(float(tp_line[5]))
		af_array.append(float(tp_line[6]))
		AF_dictionary[dict_key] = af_array

	counter = 0
	well_001 = 0
	well_005 = 0
	well_005_plus = 0

	inova_001 = 0
	inova_005 = 0
	inova_005_plus = 0
	for line in i:
		counter = counter + 1
		if line[:1] == "#":
			#o.write(line)
			continue
		else:
			tp_line = line.split("\t")
			dict_key = tp_line[1] + "_" + tp_line[3] + "_" + tp_line[4]
			
			#Won't do try, all variants should have a AF
			af_array = AF_dictionary[dict_key]
			well_AF = float(af_array[0])
			inova_AF = float(af_array[1])

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


		if counter%10000 ==0:
			print "Total lines"
			print str(counter)
			sys.stdout.flush()

	file_buffer = ""
	print "Total lines"
	print str(counter)
	AF_counter = chrom + "\t" + str(counter) + "\t" + str(well_001) + "\t" + str(well_005) + "\t" + str(well_005_plus) + \
				"\t" + str(inova_001) + "\t" + str(inova_005) + "\t" + str(inova_005_plus) + "\n"
	
	print AF_counter
	counterf.write(AF_counter)
	counterf.close()
	i.close()




if __name__ == '__main__':

    print "Python Version: " + sys.version
    #args = parse_args()
    #import time
    #start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    #end = time.time()
    #print 'This took {0} seconds'.format(end - start)

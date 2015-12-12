"""
Count total before AD filter

"""

import os, sys, gzip, datetime 

def main(chrom):


	infile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/sanity_check_wellderly_all_filters_withHWE.noZeroAF"+chrom+".vcf.gz"
	AF_file = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_noRelated_AF/AF_"+chrom+".gz")

	counter_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_count_filters/final/cluster_filters_counter.txt"
	
	counterf = open(counter_file, 'a')

	cl = gzip.open(cluster)
	i = gzip.open(infile)
	#o = gzip.open(outputfile, 'w')
	
	#counterf.write("\t\t\tWellderly\t\t\tInova\n")
	#counterf.write("AF\tTotal Var\t Filtered Var\t<0.01\t0.01-0.05\t>0.05\t<0.01\t0.01-0.05\t>0.05\n")


	cluster_dict = dict()
	#PUT the cluster file in a dictionary
	for line in cl:
		if "#" in line:
			continue
		#keep only the entries with yes	
		elif "YES" in line:
			tp_line = line.split("\t")
			key = tp_line[1] + "_" + tp_line[3] + "_" + tp_line[4]
			cluster_dict[key] = "YES"
	cl.close()

	#Index the AF file
	AF_dictionary = {}
	for line in AF_file:
		af_array = []
		tp_line = line.strip().split("\t")
		dict_key = tp_line[1]+"_"+tp_line[3]+"_"+tp_line[4]
		af_array.append(float(tp_line[5]))
		af_array.append(float(tp_line[6]))
		AF_dictionary[dict_key] = af_array

	#file_buffer = ""
	counter = 0
	counter_unique_no_pass = 0
	good_lines = 0
	
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
			
			try:
				if cluster_dict[key]:

					counter_unique_no_pass += 1

					try:
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

					except:
						continue
			except:
				good_lines = good_lines + 1

		if counter%10000 ==0:
			print "Total lines"
			print str(counter)
			print "No pass lines"
			print str(counter_unique_no_pass)
			sys.stdout.flush()

	#o.write(file_buffer)
	file_buffer = ""
	#o.flush()
	print "Total lines"
	print str(counter)
	print "New Filtered out lines"
	print str(counter_unique_no_pass)
	print "Good lines "
	print str(good_lines)
	AF_counter = chrom + "\t" + str(counter_unique_no_pass) + "\t" + str(well_001) + "\t" + str(well_005) + "\t" + str(well_005_plus) + \
				"\t" + str(inova_001) + "\t" + str(inova_005) + "\t" + str(inova_005_plus) + "\n"
	
	print AF_counter
	counterf.write(AF_counter)
	counterf.close()
	#o.close()
	cl.close()
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

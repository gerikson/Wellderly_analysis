"""
Count 36mers

"""

import os, sys, gzip, datetime 

def main(chrom):

	cluster ="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/final_filter_file/filter_wellderly_36mers."+str(chrom)+".txt.gz"

	cl = gzip.open(cluster)


	cluster_dict = dict()


	for line in cl:
		if "#" in line:
			continue
		#keep only the entries with yes	
		elif "YES" in line:
			tp_line = line.split("\t")
			key = tp_line[1] + "_" + tp_line[2] + "_" + tp_line[3]
			cluster_dict[key] = "YES"
			
	cl.close()
	print "cluster dict indexed! Clustered var in "+ chrom + " " +str(len(cluster_dict))

	cf = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/parts/counter_parts_by_chrom.txt")

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
	for l in cf:

		tp_line = l.split("\t")
		if chrom == tp_line[0]:
			number_of_parts = int(tp_line[1]) + 1
			print "chrom found"
			cf.close()
			
			for part in range(1,number_of_parts): 
				part_count += 1
				if part%100 == 0:
					print str(part)

				AF_file = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/parts/AF_parts/"+chrom+".part"+str(part)+".txt.gz")
				for line in AF_file:
					counter += 1
					tp_line = line.strip().split("\t")
					dict_key = tp_line[1]+"_"+tp_line[3]+"_"+tp_line[4]
					try:
						
						if cluster_dict[dict_key]:

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
						else:
							continue
					except:
						good_lines = good_lines + 1
				AF_file.close()
			break
		else:
			continue


	#o.write(file_buffer)
	file_buffer = ""
	#o.flush()
	print "Total lines"
	print str(counter)
	print "New Filtered out lines"
	print str(counter_unique_no_pass)
	print "Good lines "
	print str(good_lines)
	
	counter_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_count_filters/final/36mer_filters_counter.txt"
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

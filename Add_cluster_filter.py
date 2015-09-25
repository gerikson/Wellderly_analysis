"""
Combine Cluster filter wth other filters

"""

import os, sys, gzip, datetime 

def main(chrom):

	cluster = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/clustered_variants/clustered." + chrom + ".vcf.gz"
	filterfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/filter_file/filters." + chrom + ".vcf.gz"
	output = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/filter_file/VQHIGH_white_cluster/VQHIGH_white_cluster." + chrom + ".vcf.gz"
	cl = gzip.open(cluster)
	fil = gzip.open(filterfile)
	o = gzip.open(output, 'w')
	
	o.write("#CHROM\tPOS\tID\tREF\tALT\tVQHIGH\tVQHIGH_IN_WHITE\tClustered\n")
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

	file_buffer = ""
	counter = 0
	counter_unique_no_pass = 0
	good_lines = 0
	for line in fil:
		counter = counter + 1
		if "#" in line:
			continue
		elif "YES" in line:
			file_buffer = file_buffer + line
		else:
			tp_line = line.split("\t")
			key = tp_line[1] + "_" + tp_line[3] + "_" + tp_line[4]
			if key in cluster_dict.keys():
				counter_unique_no_pass = counter_unique_no_pass + 1
				l = line.strip()
				final = l + "\t\t\tYES\n"
				file_buffer = file_buffer + final
			else:
				good_lines = good_lines + 1
				l = line.strip()
				final = l + "\t\t\t\n"
				file_buffer = file_buffer + final

		if counter%100 == 0:

			o.write(file_buffer)
			file_buffer = ""
			o.flush()

		if counter%10000 ==0:
			print "Total lines"
			print str(counter)
			print "Good lines"
			print str(counter_unique_no_pass)
			sys.stdout.flush()

	o.write(file_buffer)
	file_buffer = ""
	o.flush()
	print "Total lines"
	print str(counter)
	print "New Filtered out lines"
	print str(counter_unique_no_pass)
	print "Good lines "
	print str(good_lines)
	sys.stdout.flush()
	o.close()
	cl.close()
	fil.close()




if __name__ == '__main__':

    print "Python Version: " + sys.version
    #args = parse_args()
    #import time
    #start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    #end = time.time()
    #print 'This took {0} seconds'.format(end - start)

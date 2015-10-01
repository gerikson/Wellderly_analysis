"""
Counts repeats by AF

"""
import os, sys, gzip, datetime


def main(chrom):

	input_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/wellderly_inova.VQHIGH.0.95white."+str(chrom)+".vcf.gz"
	counter_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/final_filter_file/missing_genotype_counter.txt"
	filter_file ="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/final_filter_file/missingGeno_filterFile."+str(chrom)+".txt.gz"
	
	f = gzip.open(input_filename)
	c = open(counter_file, "a")
	filt = gzip.open(filter_file, "w")

	counter = 0
	good_lines = 0


	total_wellderly = 0
	total_inova = 0
	filt.write("Chrom\tbegin\tRef\tAlt\tMissingFilter\n")
	filter_block = ""
	for line in f:

		if line[:2] == "##":
			print "header"
		elif line[:1] == "#":
			tp_line = line.split("\t")
			for i in tp_line:
				if i.endswith("DID"):
					total_wellderly += 1
				elif i.endswith("ASM"):
					total_inova += 1
			print "Total Wellderly " + str(total_wellderly)
			print "Total inova " + str(total_inova)
		else:

			counter += 1
			if counter%100 == 0:
				filt.write(filter_block)
				filter_block = ""
				filt.flush()

			if counter%10000 == 0:
				print datetime.datetime.now().time()
				print "total lines " + str(counter)
				print "Good lines " + str(good_lines)
				sys.stdout.flush()


			tp_line = line.strip().split()
			filter_block = filter_block + tp_line[0]+"\t"+tp_line[1]+"\t"+tp_line[3]+"\t"+tp_line[4]
			
			well_missing = 0
			inova_missing = 0
			for index, gen in enumerate(tp_line[9:]):
				index = index + 9

				if gen[0] == '.' or gen[2] == '.':
					if index < 529:
						well_missing += 1
					else:
						inova_missing += 1

			if well_missing > 51 or inova_missing >68:
				#print "WELLDERLY missing " + str(well_missing)
				#print "inova_missing " + str(inova_missing) 
				filter_block = filter_block + "\tYES\n" 
				continue
			else:
				filter_block = filter_block + "\t\n" 
				good_lines += 1
			

	
	print "Total lines " + str(counter)
	print "Total good lines " + str(good_lines)

	AF_counter = chrom + "\t" + str(counter) + "\t" + str(good_lines)
	print AF_counter
	c.write(AF_counter)
	c.close()
	filt.write(filter_block)
	filter_block = ""
	filt.flush()
	
	f.close()
	filt.close()


if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)
"""
Extract all of the variants that have ANY clustered genotypes

"""

import os, sys, gzip, datetime 




DID_start = 0
Inova_start = 0

total_DID = 0
total_ASM = 0



def main(chrom):
	infile = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_withCorrectGenotypes/wellderly_inova." + chrom + ".vcf.gz")
	outfile = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_ALL_clustered/ALL_clustered."+chrom+".txt.gz", "w")


	outfile.write("#CHROM\tPOS\tID\tREF\tALT\n")

	counter = 0
	filter_counter = 0

	buffe = ""
	for line in infile:

		if line[:2] == "##":
			print "header"


		elif line[0] == "#":
			print "header"

		else:
			counter += 1
			line = line.strip()
			if filter_counter%100 == 0:

				outfile.write(buffe)
				buffe = ""
				outfile.flush()

			if counter%100000 ==0:
				print "Total lines"
				print str(counter)
				print "Filtered lines"
				print str(filter_counter)
				print datetime.datetime.now().time()
				sys.stdout.flush()

			#If line ends with "::" no variant at that position was clustered
			if line[-2:] == "::":
				continue
			else:
				filter_counter+=1
				tp_line = line.split("\t")
				filter_line = "\t".join(tp_line[:5])
				buffe = buffe+filter_line+"\n"



	outfile.write(buffe)
	infile.close()
	outfile.close()
	print "Total lines"
	print str(counter)
	print "filtered lines"
	print str(filter_counter)
	print "End"
	counter_file = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_ALL_clustered/ALL_clustered_counter.txt", "a")
	counter_file.write(chrom+"\t"+str(counter)+"\t"+str(filter_counter)+"\n")
	counter_file.close()


if __name__ == '__main__':

    print "Python Version: " + sys.version
    #args = parse_args()
    #import time
    #start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    #end = time.time()
    #print 'This took {0} seconds'.format(end - start)

"""
Extract related individuals

"""

import os, sys, gzip, datetime 

def main(chrom):
	#infile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_allFilters_36kmer_snpsOnly_AF0.01/vcf_nokmer_snps_AF0.01." + chrom + ".vcf.gz"
	#outfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_allFilters_36kmer_snpsOnly_AF0.01/vcf_nokmer_snps_AF0.01.noRelated." + chrom + ".vcf.gz"
	
	#For the riskogram
	infile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/filtered_snps.withHead.txt"
	outfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/filtered_snps.noRelated.withHead.txt"
	relatedfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_allFilters_36kmer_snpsOnly_AF0.01/eliminate_individuals.txt"

	w = open(relatedfile)
	related_id = []

	header = []
	#index the whites
	ln = w.readline()
	ln = ln.strip()
	#Create dictionary here instead
	whites_id = ln.split("\t")


	#i = gzip.open(infile)
	#o = gzip.open(outfile, 'w')

	i = open(infile)
	o = open(outfile, 'w')



	counter_white = 0
	counter = 0

	buffe = ""

	
	for line in i:

		if line[:2] == "##":
			print "header"
			o.write(line)

		#Extract only the white individuals ID
		elif line[0] == "#":
			#o.write(line)
			line = line.strip()
			header = line.split("\t")
			final_header = "\t".join(header[:9])
			for index, gen in enumerate(header[9:]):
				index = index +9
				#Extract whites only
				if gen in whites_id:
					continue
				else:
					final_header = final_header + "\t" + gen
					counter_white += 1
			final_header = final_header + "\n"
			o.write(final_header)
			print "# of whites " + str(counter_white)

		else:
			counter += 1
			line = line.strip()
			tp_line = line.split("\t")

			white_line = "\t".join(tp_line[:9])
			#count_id = 0
			for index, gen in enumerate(tp_line[9:]):
				index = index +9
				indiv_name = header[index]
				if indiv_name in whites_id:
					continue
				else:
					#count_id += 1
					white_line = white_line + "\t" + gen
			buffe = buffe + white_line + "\n"


			if counter%100 == 0:
				o.write(buffe)
				buffe = ""
				o.flush()

			if counter%10000 ==0:
				print "Total lines"
				print str(counter)
				sys.stdout.flush()

	o.write(buffe)
	i.close()
	o.close()
	#clus.close()
	print "Total lines"
	print str(counter)
	print "End"


if __name__ == '__main__':

    print "Python Version: " + sys.version
    #args = parse_args()
    #import time
    #start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    #end = time.time()
    #print 'This took {0} seconds'.format(end - start)

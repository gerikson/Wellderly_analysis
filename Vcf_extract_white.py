"""
Extract lines whites only
Remove variants that were not present in the whites

"""

import os, sys, gzip, datetime 

def main(chrom):
	infile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_withCorrectGenotypes/wellderly_inova." + chrom + ".vcf.gz"
	outfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_whiteOnly/wellderly_inova.0.95white." + chrom + ".vcf.gz"
	whitesfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_withCorrectGenotypes/White_0.95.txt"

	filterfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_whiteOnly/filter.0.95white." + chrom + ".vcf.gz"
	counter_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_whiteOnly/Whites_only_counter.txt"


	w = open(whitesfile)
	whites_id = []

	header = []
	#index the whites
	ln = w.readline()
	ln = ln.strip()
	#Create dictionary here instead
	whites_id = ln.split()


	i = gzip.open(infile)
	o = gzip.open(outfile, 'w')
	#clus = gzip.open(clusteredfile, 'w')
	filt =gzip.open(filterfile, 'w')
	counterf = open(counter_file, 'a')
	#counterf.write("CHROM\tOriginalLines\tFilteredLines\n")
	filt.write("#CHROM\tPOS\tID\tREF\tALT\tNot_white\n")

	counter_white = 0
	counter = 0
	good_counter = 0
	buffe = ""
	buffe_filter = ""
	
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
					final_header = final_header + "\t" + gen
					counter_white += 1
			final_header = final_header + "\n"
			o.write(final_header)
			print "# of whites " + str(counter_white)

		else:
			counter += 1
			line = line.strip()
			tp_line = line.split("\t")

			filter_line = "\t".join(tp_line[:5])

			filter_line = filter_line + "\t"
			white_line = "\t".join(tp_line[:9])
			#count_id = 0
			for index, gen in enumerate(tp_line[9:]):
				index = index +9
				indiv_name = header[index]
				if indiv_name in whites_id:
					#count_id += 1
					white_line = white_line + "\t" + gen
			
			#check if the new line with whites only have any alternate alleles
			#(if there are any alternate alleles they have VQHIGH or VQLOW)
			if 'VQHIGH' in white_line or 'VQLOW' in white_line:	
				good_counter += 1
				buffe = buffe + white_line + "\n"
				#Don't add anything
				filter_line = filter_line + "\t\n"
			else:
				filter_line = filter_line + "\tYES\n"


			buffe_filter = buffe_filter + filter_line

			if good_counter%100 == 0:
				o.write(buffe)
				buffe = ""
				filt.write(buffe_filter)
				buffe_filter = ""
				o.flush()
				filt.flush()

			if counter%10000 ==0:
				print "Total lines"
				print str(counter)
				print "Good lines"
				print str(good_counter)
				sys.stdout.flush()

	o.write(buffe)
	i.close()
	o.close()
	filt.write(buffe_filter)
	#clus.close()
	print "Total lines"
	print str(counter)
	print "Good lines"
	print str(good_counter)
	print "End"
	t = chrom + "\t" + str(counter) + "\t" + str(good_counter) + "\n"
	counterf.write(t)
	counterf.close()
	filt.close()




if __name__ == '__main__':

    print "Python Version: " + sys.version
    #args = parse_args()
    #import time
    #start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    #end = time.time()
    #print 'This took {0} seconds'.format(end - start)

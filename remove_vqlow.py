"""
Extract lines that don't have any variant called VQHIGH
Extract only white individuals, if there are no VQHIGH in the white dataset, remove those entries

"""

import os, sys, gzip, datetime 

def main(chrom):
	infile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_withCorrectGenotypes/wellderly_inova." + chrom + ".vcf.gz"
	outfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/wellderly_inova.VQHIGH.0.95white." + chrom + ".vcf.gz"
	whitesfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_withCorrectGenotypes/White_0.95.txt"
	#clusteredfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/clustered_variant/clustered." + chrom + ".vcf.gz"
	filterfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/filters." + chrom + ".vcf.gz"
	counter_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/VQHIGH_filters_counter.txt"


	w = open(whitesfile)
	wh = []

	header = {}
	#index the whites
	ln = w.readline()
	ln = ln.strip()
	wh= ln.split()
	
	whites_id = {}
	for i in wh:
		whites_id[i] = "Y"


	i = gzip.open(infile)
	o = gzip.open(outfile, 'w')
	#clus = gzip.open(clusteredfile, 'w')
	filt =gzip.open(filterfile, 'w')
	counterf = open(counter_file, 'a')
	#counterf.write("CHROM\tOriginalLines\tFilteredLines\n")
	filt.write("#CHROM\tPOS\tID\tREF\tALT\tVQHIGH\tVQHIGH_IN_WHITE\n")

	counter = 0
	good_counter = 0
	buffe = ""
	buffe_filter = ""
	DID_start = 0
	Inova_start = 0
	for line in i:

		if line[:2] == "##":
			print "header"
			o.write(line)

		#Extract only the white individuals ID
		elif line[0] == "#":
			#o.write(line)
			line = line.strip()
			h = line.split("\t")

			#Put the data into a dictionary first instead of an array
			for index, he in  enumerate(h):
				header[index] = he

			final_header = "\t".join(h[:9])

			for index, gen in enumerate(h[9:]):
				index = index +9
				#Extract whites only
				if gen in whites_id.keys():
					final_header = final_header + "\t" + gen

			final_header = final_header + "\n"
			o.write(final_header)

		else:
			counter += 1
			line = line.strip()
			tp_line = line.split("\t")

			filter_line = "\t".join(tp_line[:5])

			if 'VQHIGH' in line:
				filter_line = filter_line + "\t"
				#check the last 2 elements of the line, if it ends in '::' no variants were clustered
				#NO CLUSTER VARIANTS
				#if line[-2:] == "::":

				white_line = "\t".join(tp_line[:9])
				#count_id = 0
				for index, gen in enumerate(tp_line[9:]):
					index = index +9
					indiv_name = header[index]

					#USE a dictionary instead of an array, for the indexing
					try: 
						if whites_id[indiv_name] == "Y":
							#count_id += 1
							white_line = white_line + "\t" + gen
					except:
						continue
				#check if the new line with whites only have any VQHIGHs
				if 'VQHIGH' in white_line:	
					good_counter += 1
					buffe = buffe + white_line + "\n"
					#Dont add anything
					filter_line = filter_line + "\t\n"
				else:
					filter_line = filter_line + "\tYES\n"

			else:
				filter_line = filter_line + "\tYES\t\n"

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

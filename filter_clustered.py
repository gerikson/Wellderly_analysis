"""
Filter clustered variants present in 10 percent of wellderly/inova samples

"""

import os, sys, gzip, datetime 
counter_alleles_inova = 0
counter_alleles_wellderly = 0

well_001 = 0
well_005 = 0
well_005_plus = 0

inova_001 = 0
inova_005 = 0
inova_005_plus = 0


DID_start = 0
Inova_start = 0

total_DID = 0
total_ASM = 0

def count_AF(line):
	tp_line = line.split("\t")
	alleles = tp_line[4].split(",")

	total_well_alleles = 0
	total_inova_alleles = 0
	total_well_alt = 0
	total_inova_alt = 0
	dict_of_alleles_well={}
	dict_of_alleles_inova={}

	well_AF = 0
	inova_AF = 0
	if len(alleles) > 1:
		#print "len alleles > 1"
		#create a dictionary that would count the number of alleles
		for index, i in enumerate(alleles):
			#add one to the index, alleles start at 1
			ind = index + 1
			dict_of_alleles_well[ind] = 0
			dict_of_alleles_inova[ind] = 0


	for index, gen in enumerate(tp_line[9:]):
		index = index + 9
		if len(alleles) == 1:
			#tp_gen = gen[:3]
			#First Allele
			if gen[0] == "0":
				if index < 529:
					total_well_alleles += 1
				else:
					total_inova_alleles += 1
			elif gen[0] == "1":
				if index < 529:
					total_well_alleles += 1
					total_well_alt +=1
				else:
					total_inova_alleles += 1
					total_inova_alt +=1
			#Second Allele
			if gen[2] == "0":
				if index < 529:
					total_well_alleles += 1
				else:
					total_inova_alleles += 1
			elif gen[2] == "1":
				if index < 529:
					total_well_alleles += 1
					total_well_alt +=1
				else:
					total_inova_alleles += 1
					total_inova_alt +=1
		
		else:
			#print "Multi alleles"
			#First Allele
			if gen[0] != ".":
				if gen[0] == "0":
					if index < 529:
						total_well_alleles += 1
					else:
						total_inova_alleles += 1
				
				else:

					if index < 529:
						total_well_alleles += 1
						#Increment the dictionary value for that allele
						try:
							dict_of_alleles_well[int(gen[0])] +=1
						except:
							print "gen[0] " + str(gen[0])
							print "alleles " + tp_line[4]
					else:
						total_inova_alleles += 1
						try:
							dict_of_alleles_inova[int(gen[0])] +=1
						except:
							print "gen[0] " + str(gen[0])
							print "alleles " + tp_line[4]

			#Second Allele
			if gen[2] != ".":
				if gen[2] == "0":
					if index < 529:
						total_well_alleles += 1
					else:
						total_inova_alleles += 1
				else:
					if index < 529:
						total_well_alleles += 1
						#Increment the dictionary value for that allele
						try:
							dict_of_alleles_well[int(gen[2])] +=1
						except:
							print "gen[0] " + str(gen[2])
							print "alleles " + tp_line[4]
					else:
						total_inova_alleles += 1
						try:
							dict_of_alleles_inova[int(gen[2])] +=1
						except:
							print "gen[0] " + str(gen[2])
							print "alleles " + tp_line[4]

			inova_max_value = max(dict_of_alleles_inova.values())
			well_max_value = max(dict_of_alleles_well.values())

	try:
		if len(alleles) == 1:
			well_AF = float(total_well_alt)/float(total_well_alleles)
		else:
			# For multi allelic values remove the highest AF from 1
			well_AF = 1 - float(well_max_value)/float(total_well_alleles)

		#print str(AF)
	except:
		well_AF = 0.0

	try:
		if len(alleles) == 1:
			inova_AF = float(total_inova_alt)/float(total_inova_alleles)
		else:
			inova_AF = 1 - float(inova_max_value)/float(total_inova_alleles)
		#print str(AF)
	except:
		inova_AF = 0.0

	if well_AF < 0.01:
		well_001 += 1
	elif well_AF < 0.05:
		well_005 += 1
	else:
		well_005_plus += 1

	if inova_AF < 0.01:
		inova_001 += 1
	elif inova_AF < 0.05:
		inova_005 += 1
	else:
		inova_005_plus +=1



def main(chrom):
	infile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/wellderly_inova.VQHIGH.0.95white." + chrom + ".vcf.gz"
	outfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_clustered/clustered." + chrom + ".vcf.gz"
	counter_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_clustered/cluster_filters_counter.txt"
	filter_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_clustered/filters.clustered." + chrom + ".vcf.gz"

	inf = gzip.open(infile)
	o = gzip.open(outfile, 'w')
	filt = gzip.open(filter_file, 'w')
	counterf = open(counter_file, 'a')

	#counterf.write("CHROM\tOriginalLines\tFilteredLines\n")
	filt.write("#CHROM\tPOS\tID\tREF\tALT\tCluster\n")
	counterf.write("\t\t\tWellderly\t\t\tInova\n")
	counterf.write("AF\tTotal Var\t Filtered Var\t<0.01\t0.01-0.05\t>0.05\t<0.01\t0.01-0.05\t>0.05\n")

	counter = 0
	filter_counter = 0
	buffe_filter = ""
	buffe = ""
	for line in inf:

		if line[:2] == "##":
			print "header"
			o.write(line)

		#Extract only the white individuals ID
		elif line[0] == "#":
			print "header"
			o.write(line)
		else:
			counter += 1
			line = line.strip()
			tp_line = line.split("\t")

			filter_line = "\t".join(tp_line[:5])

			if line[-2:] == "::":
				filter_line = filter_line + "\t\n"
			else:
				
				clusters = line.split(":")
				well_clust = []
				inova_clust = []
				try:
					count_elements = len(clusters)
					#print "Elements counted " + str(count_elements)
					well_clust = clusters[count_elements-1].split(",")
					inova_clust = clusters[count_elements-2].split(",")
					
					#print "Length clustered inova: " + str(len(inova_clust))
				except:
					print "Except"

				#If more then 10% in either wellderly or inova clusetered together
				if len(well_clust) > 60 or len(inova_clust) >150:
					#print "found"
					filter_counter = filter_counter + 1
					filter_line = filter_line + "\tYES\n"
					count_AF(line)
				else:
					filter_line = filter_line + "\t\n"
					buffe = buffe + line

			buffe_filter = buffe_filter + filter_line

			if filter_counter%100 == 0:

				filt.write(buffe_filter)
				buffe_filter = ""
				filt.flush()
				o.write(buffe)
				buffe = ""
				o.flush()

			if counter%10000 ==0:
				print "Total lines"
				print str(counter)
				print "Good lines"
				print str(filter_counter)
				sys.stdout.flush()


	filt.write(buffe_filter)
	o.write(buffe)
	#clus.close()
	print "Total lines"
	print str(counter)
	print "Good lines"
	print str(filter_counter)
	print "End"
	AF_counter = chrom + "\t" + str(counter) + "\t" + str(counter_unique_no_pass) + "\t" + str(well_001) + "\t" + str(well_005) + "\t" + str(well_005_plus) + \
					"\t" + str(inova_001) + "\t" + str(inova_005) + "\t" + str(inova_005_plus)
	
	print AF_counter
	counterf.write(AF_counter)
	counterf.close()
	inf.close()
	o.close()
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

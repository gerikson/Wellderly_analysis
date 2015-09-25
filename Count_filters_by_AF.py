"""
Count the filters depending of the AF
"""

import os, sys, gzip, datetime 

def main(chrom):

	infile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/wellderly_inova.VQHIGH.0.95white." + chrom + ".vcf.gz"
	
	counter_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/counter_by_alleleFreq.txt"

	i = gzip.open(infile)
	counterf = open(counter_file, 'a')

	#counterf.write("\tWellderly\tt\t\tInova\n")
	#counterf.write("AF\t<0.01\t0.01-0.05\t>0.05\t<0.01\t0.01-0.05\t>0.05\n")

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

	counter = 0

	for line in i:

		if line[:2] == "##":
			print "header"
			#o.write(line)

		#Extract only the white individuals ID
		elif line[0] == "#":
			#o.write(line)
			line = line.strip()
			header = line.split("\t")
			for index, indiv in enumerate(header):
				print indiv
				if indiv.endswith("DID") and DID_start == 0:
					DID_start = index
				elif indiv.endswith("ASM") and Inova_start == 0:
					Inova_start = index

				if indiv[:-3] == "DID":
					total_DID += 1
				elif indiv[:-3] == "ASM":
					total_ASM += 1

			print "Total Wellderly: " + str(total_DID)
			print "Total Inova: " + str(total_ASM)
			print "Wellderly start at: " + str(DID_start)
			print "Inova starts at: " + str(Inova_start)

		else:
			counter += 1
			line = line.strip()
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
				


			
			if counter%10000 ==0:
				print "Total lines"
				print str(counter)
				sys.stdout.flush()
			

	print "Total lines"
	print str(counter)
	print "End"
	
	AF_counter = chrom + "\t" + str(well_001) + "\t" + str(well_005) + "\t" + str(well_005_plus) + \
					"\t" + str(inova_001) + "\t" + str(inova_005) + "\t" + str(inova_005_plus)
	
	print AF_counter
	counterf.write(AF_counter)
	counterf.close()
	i.close()
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

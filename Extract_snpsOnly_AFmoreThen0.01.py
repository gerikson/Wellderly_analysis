"""
Extracts separately the snps with AF >0.01 and the delins with AF>0.01

"""

import os, sys, gzip, datetime 

def main(chrom):


	infile="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_clustered_repeats_homopoly_etc_missing_cov/v1_wellderly_inova.VQHIGH.0.95white.nocluster.repeats.etc.missing.cov."+str(chrom)+".vcf.gz"
	snp_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_AFmore0.01/vcf_snps_AF0.01."+str(chrom)+".vcf.gz"
	delins_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_delins_AFmore0.01/vcf_delins_AF0.01."+str(chrom)+".vcf.gz"
	counter_snps="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_AFmore0.01/snps_counter_AF0.01.txt"
	counter_delins="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_delins_AFmore0.01/delins_counter_AF0.01.txt"
	counterfile_byAF="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_AFmore0.01/counter_file_by_AF.txt"

	i = gzip.open(infile)
	#snpfile = gzip.open(snp_file, "w")
	#delinsfile = gzip.open(delins_file, "w")
	count_snps = open(counter_snps, "a")
	count_delins = open(counter_delins, "a")
	counterf = open(counterfile_byAF, "a")

	counter = 0
	counter_snps_AF001 = 0
	counter_delins_AF001 = 0
	
	well_001 = 0
	well_005 = 0
	well_005_plus = 0

	inova_001 = 0
	inova_005 = 0
	inova_005_plus = 0

	#snp_buffer = ""
	#delin_buffer = ""

	for line in i:
		counter += 1
		if line[:1] == "#":
			#snpfile.write(line)
			#delinsfile.write(line)
			continue
		else:
			tp_line = line.split("\t")
			key = tp_line[1] + "_" + tp_line[3] + "_" + tp_line[4]
			is_this_snp = True

			alleles = tp_line[4].split(",")

			#verify is this variant is a snp or a delins etc
			if len(tp_line[3]) > 1:
				is_this_snp = False

			for al in alleles:
				if len(al) > 1:
					is_this_snp = False
					break

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

			
			if well_AF > 0.01 or inova_AF > 0.01:
				if is_this_snp:
					#snp_buffer = snp_buffer+line
					counter_snps_AF001 += 1
				else:
					#delin_buffer = delin_buffer+line
					counter_delins_AF001 += 1
			
		'''
		if counter%100 == 0:
			snpfile.write(snp_buffer)
			snp_buffer = ""
			snpfile.flush()

			delinsfile.write(delin_buffer)
			delin_buffer = ""
			delinsfile.flush()
		'''
		if counter%10000 ==0:
			print "Total lines"
			print str(counter)
			sys.stdout.flush()

	'''
	snpfile.write(snp_buffer)
	snp_buffer = ""
	snpfile.flush()

	delinsfile.write(delin_buffer)
	'''
	print "Final Total lines"
	print str(counter)

	AF_counter = chrom + "\t" + str(counter) + "\t" + str(well_001) + "\t" + str(well_005) + "\t" + str(well_005_plus) + \
				"\t" + str(inova_001) + "\t" + str(inova_005) + "\t" + str(inova_005_plus)
	
	print AF_counter
	counterf.write(AF_counter)
	counterf.close()
	count_snps.write(chrom+"\t"+str(counter)+"\t"+str(counter_snps_AF001)+"\n")
	count_snps.close()
	count_delins.write(chrom+"\t"+str(counter)+"\t"+str(counter_delins_AF001)+"\n")
	count_delins.close()
	i.close()
	#snpfile.close()
	#delinsfile.close()





if __name__ == '__main__':

    print "Python Version: " + sys.version
    #args = parse_args()
    #import time
    #start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    #end = time.time()
    #print 'This took {0} seconds'.format(end - start)

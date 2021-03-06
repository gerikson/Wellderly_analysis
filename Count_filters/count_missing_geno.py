"""
Counts missing genotype
Updated to exclude the related individuals
Set VQLOW to missing
And extract the counts again

"""
import os, sys, gzip, datetime


def main(chrom):

	input_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/wellderly_inova.VQHIGH.0.95white."+str(chrom)+".vcf.gz"
	counter_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_count_filters/missing_genotype.txt"
	#filter_file ="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_clustered_repeats_homopoly_etc_missing/filter_missing_genotype.txt.gz"
	filter_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/final_filter_file/filter_wellderly_missingGeno"+str(chrom)+".txt.gz"
	
	relatedfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_allFilters_36kmer_snpsOnly_AF0.01/eliminate_individuals.txt"

	w = open(relatedfile)
	header = []
	#index the whites
	ln = w.readline()
	ln = ln.strip()
	#Create dictionary here instead
	related_id = ln.split("\t")
	w.close()
	header = []
	f = gzip.open(input_filename)
	filt = gzip.open(filter_file, "w")

	counter = 0
	good_lines = 0
	missing_counter = 0

	well_001 = 0
	well_005 = 0
	well_005_plus = 0

	inova_001 = 0
	inova_005 = 0
	inova_005_plus = 0

	total_wellderly = 0
	total_inova = 0
	filt.write("Chrom\tbegin\tRef\tAlt\tMissingFilter\n")
	filter_block = ""
	for line in f:

		if line[:2] == "##":
			print "header"
			continue
		elif line[:1] == "#":

			line = line.strip()
			header = line.split("\t")
			final_header = "\t".join(header[:9])
			counter_white = 0
			for index, gen in enumerate(header[9:]):
				index = index +9
				#Extract whites only
				if gen in related_id:
					continue
				else:
					final_header = final_header + "\t" + gen
					counter_white += 1
			final_header = final_header + "\n"
			print "# of related " + str(counter_white)
			continue
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


			#Remove related individuals
			noRelated_line = tp_line[:9]
			for index, gen in enumerate(tp_line[9:]):
				index = index +9
				indiv_name = header[index]
				if indiv_name in related_id:
					continue
				else:
					noRelated_line.append(gen)

			#print "Removed related" 

			#Replace VQLOW and count missing
			well_missing = 0
			inova_missing = 0
			for index,i in enumerate(noRelated_line):
				if index > 8:
					var = i.split(":")
					geno = var[0]
					geno_split = geno.split("/")
					if "VQLOW" in i:
						#print "VQLOW"
						if var[3]=="VQLOW":
						    geno_split[0]='.'
						if var[4]=="VQLOW":
						    geno_split[1]='.'
						var[0]="/".join(geno_split)
						#print "new geno " + ":".join(var)
						noRelated_line[index] = ":".join(var)

					if geno_split[0] == '.' or geno_split[1] == '.':
						if index < 511:
							well_missing += 1
						else:
							inova_missing += 1





			'''
			#This was the original version
			for index, gen in enumerate(tp_line[9:]):
				index = index + 9
				#geno = gen.split(":")
				#if "." in geno[0]:
				if gen[0] == '.' or gen[2] == '.':
					if index < 529:
						well_missing += 1
					else:
						inova_missing += 1
			#print "wellderly missing " + str(well_missing)
			#print "inova missing " + str(inova_missing)
			'''

			if well_missing > 51 or inova_missing >68:
				filter_block = filter_block + "\tYES\n" 

				missing_counter += 1

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

				continue
			else:
				filter_block = filter_block + "\t\n" 
				good_lines += 1
				#filter_block = filter_block + "\t"
			


	print "Total missing counter " + str(missing_counter)

	
	print "Total lines " + str(counter)
	print "Total good lines " + str(good_lines)
	c = open(counter_file, "a")
	AF_counter = chrom + "\t" + str(counter) + "\t" + str(missing_counter) + "\t" + str(well_001) + "\t" + str(well_005) + "\t" + str(well_005_plus) + \
			"\t" + str(inova_001) + "\t" + str(inova_005) + "\t" + str(inova_005_plus) + "\n"

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
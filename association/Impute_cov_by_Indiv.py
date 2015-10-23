"""
Impute the coverage information 

"""
import os, sys, gzip, datetime


def main():

	ch = "1"
	input_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/smallest_pvalues.correct.sorted.vcf"
	output_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/smallest_pvalues.correct.sorted.inputedCov.vcf"
	
	coverage_inova="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/MedianCov_snpsOfInterest_inova.txt.header"
	coverage_wellderly="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/MedianCov_snpsOfInterest_wellderly.txt.header"


	f = open(input_filename)
	outf = open(output_file, "w")
	

	cov_file = open(coverage_inova)
	covereage_in = {}
	print "Start creating dictionary"
	counter = 0
	head = []
	for i in cov_file:
		counter += 1
		#If this is the first line
		if counter == 1:
			head = i.strip().split(",")
			print "lenght header: " + str(len(head))
		else:
			tp_line = i.split(":")
			var = tp_line[1]
			tp_var = var.strip().split(",")
			
			var_dict={}
			var_id=""
			for index, v in enumerate(tp_var):
				if index == 0:
					t_v = v.strip().split("_")
					var_id = t_v[0]
				else:
					#Put data in a dictionary inside the dictionary
					w = head[index].strip()
					var_dict[w]=v

			covereage_in[var_id]=var_dict
			
	print "Dictionary created inova " + str(len(covereage_in))

	cov_file.close()

	cov_file = open(coverage_wellderly)
	print "Start creating dictionary"
	counter = 0
	head = []
	covereage_well = {}
	for i in cov_file:
		counter += 1
		#If this is the first line
		if counter == 1:
			head = i.strip().split(",")
			for index, iv in enumerate(head):
				new_iv=iv+"-DID"
				head[index]=new_iv
			print "lenght header: " + str(len(head))
		else:
			tp_line = i.split(":")
			var = tp_line[1]
			tp_var = var.strip().split(",")			
			var_dict={}
			var_id=""
			for index, v in enumerate(tp_var):
				if index == 0:
					t_v = v.strip().split("_")
					var_id = t_v[0]
				else:
					#Put data in a dictionary inside the dictionary
					w = head[index].strip()
					var_dict[w]=v

			covereage_well[var_id]=var_dict
			
	print "Dictionary created wellderly "+ str(len(covereage_well))

	cov_file.close()
	
	#Test
	#print covereage_well["128220474_1"]["GS000015806-DID"]

	counter = 0
	found = 0
	total_wellderly = 0
	total_inova = 0
	header_vcf = []

	for line in f:
		cov_well = ""
		cov_in = ""

		if line[:2] == "##":
			outf.write(line)
			print "header"
		elif line[:1] == "#":
			outf.write(line)
			tp_line = line.strip().split("\t")
			#header = tp_line[9:]
			header_vcf = tp_line
			#print header_vcf
			for ind, i in enumerate(tp_line):
				if i.endswith("DID"):
					total_wellderly += 1
				elif i.endswith("ASM"):					
					#print "start index of inova: " + str(ind)					
					total_inova += 1
					#continue

			print "Total Wellderly " + str(total_wellderly) 
			print "Total inova " + str(total_inova)
		else:

			counter += 1

			found_var = False
			if counter%10000 == 0:
				print datetime.datetime.now().time()
				print "total lines " + str(counter)
				sys.stdout.flush()


			tp_line = line.strip().split()

			#final_line = tp_line[:9]

			
			for index, i in enumerate(tp_line):

				if index > 8:
					var = tp_line[index].split(":") 
					indiv = header_vcf[index]

					v_i = str(tp_line[1])
					try:
						c = covereage_well[v_i][indiv]
						found_var = True
						try:
							var[5] = c.strip()
							tp_line[index] = ":".join(var)
							found_var = True
						except:
							tp_line[index] = ":".join(var)+":"+c.strip()
							found_var = True

					except:
						try:
							c = covereage_in[v_i][indiv]
							found_var = True
							try:
								var[5] = c.strip()
								tp_line[index] = ":".join(var)
								found_var = True
							except:
								tp_line[index] = ":".join(var)+":"+c.strip()
								found_var = True
						except:
							continue
						

			final_line = "\t".join(tp_line) + "\n"
			#print final_line
			outf.write(final_line)
			if found_var:
				found += 1
			











			'''
			#Extract VQLOW counts
			wellderly_VQLOW = 0
			inova_VQLOW = 0


			#if this is a different chromosome
			if tp_line[0] != ch:
				ch = tp_line[0]
				print "New Chrom " + ch
				coverage_inova="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/inova_median_coverage/MediansCompiled/medians_chrm_"+str(ch)+".csv"
				coverage_wellderly="/gpfs/group/torkamani/bhuvan/wellderly/coverage/CoverageInfo/MediansCompiled/medians_chrm_"+str(ch)+".csv"
			    

				f = open(input_filename)
				outf = open(output_file, "a")
				
				cov_file = open(coverage_wellderly)
				covereage_well = {}
				print "Start creating dictionary"
				for i in cov_file:
				    ln = i.split(",")
				    covereage_well[ln[0]] = ln[1]
				print "Dictionary created wellderly"
				cov_file.close()

				cov_file = open(coverage_inova)
				covereage_in = {}
				print "Start creating dictionary"
				for i in cov_file:
				    ln = i.split(",")
				    covereage_in[ln[0]] = ln[1]
				print "Dictionary created inova"
				cov_file.close()

			pos = tp_line[1] + "_" + str(len(tp_line[3]))

			#Extract the coverage from the dictionary

			try:
			    cov_well = covereage_well[pos]
			except:
			    print "coverage not found wellderly"
			    print tp_line[:5]

			try:
			    cov_in = covereage_in[pos]
			except:
			    print "coverage not found wellderly"
			    print tp_line[:5]

			alleles = tp_line[4].split(",")

			total_well_alleles = 0
			total_inova_alleles = 0
			total_well_alt = 0
			total_inova_alt = 0
			dict_of_alleles_well={}
			dict_of_alleles_inova={}

			well_AF = 0
			inova_AF = 0
			missing_well=0
			missing_inova=0
			if len(alleles) > 1:
				#print "len alleles > 1"
				#create a dictionary that would count the number of alleles
				for index, i in enumerate(alleles):
					#add one to the index, alleles start at 1
					ind = index + 1
					dict_of_alleles_well[ind] = 0
					dict_of_alleles_inova[ind] = 0


			for index, gen in enumerate(tp_line[9:]):
				#Extract VQLOW counts
				VQ_info = gen.split(":")
				if VQ_info[3]=="VQLOW":
					if index < 520:
						wellderly_VQLOW += 1
					else:
						inova_VQLOW += 1

				if VQ_info[4]=="VQLOW":
					if index < 520:
						wellderly_VQLOW += 1
					else:
						inova_VQLOW += 1

				index = index + 9
				if len(alleles) == 1:
					#tp_gen = gen[:3]
					#First Allele
					if gen[0] == "0":
						if index < 520:
							total_well_alleles += 1
						else:
							total_inova_alleles += 1
					elif gen[0] == "1":
						if index < 520:
							total_well_alleles += 1
							total_well_alt +=1
						else:
							total_inova_alleles += 1
							total_inova_alt +=1
					if gen[0] == ".":
						if index < 520:
							missing_well += 1
						else:
							missing_inova += 1
					#Second Allele
					if gen[2] == "0":
						if index < 520:
							total_well_alleles += 1
						else:
							total_inova_alleles += 1
					elif gen[2] == "1":
						if index < 520:
							total_well_alleles += 1
							total_well_alt +=1
						else:
							total_inova_alleles += 1
							total_inova_alt +=1

					if gen[2] == ".":
						if index < 520:
							missing_well += 1
						else:
							missing_inova += 1
				
				else:
					#print "Multi alleles"
					#First Allele
					if gen[0] == ".":
						if index < 520:
							missing_well += 1
						else:
							missing_inova += 1
					else:
						if gen[0] == "0":
							if index < 520:
								total_well_alleles += 1
							else:
								total_inova_alleles += 1
						
						else:

							if index < 520:
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
					if gen[2] == ".":
						if index < 520:
							missing_well += 1
						else:
							missing_inova += 1
					else:
						if gen[2] == "0":
							if index < 520:
								total_well_alleles += 1
							else:
								total_inova_alleles += 1
						else:
							if index < 520:
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

			final_line = "\t".join(tp_line[:5]) + "\t" +str(well_AF) + "\t" + str(inova_AF) + "\t" + cov_well.strip() + "\t" + cov_in.strip() + "\t" + str(missing_well) + "\t"+ str(missing_inova) + "\t"+ str(wellderly_VQLOW) + "\t" + str(inova_VQLOW)+"\n"
			outf.write(final_line)
			print final_line


			
	
	print "Total lines " + str(counter)
	'''

	print "Total varinats " + str(counter)
	print "Found variants " + str(found)
	outf.close()
	f.close()



if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main()
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)
"""
From the cypher annotations counts:
Total # of SNPs (% novel)  
     # Nonsense (% novel)
     # Splice-Site (% novel)
     # Missense (% novel)
     # Synonymous (% novel)
     # Intronic (% novel)
     # Intergenic (% novel)
Total # of INS (% novel)
Total # of DEL (% novel)
     # Frameshift (% novel)
     In-frame (% novel)
     Splice-Site (% novel)
     Intronic (% novel)
     Intergenic (% novel)


"""
import os, sys, gzip, datetime


def main(chrom):

	#plink_filename_freq_wellderly=open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/Plink_files/"+str(chrom)+".wellderly.frq")
	#plink_filename_freq_inova=open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/Plink_files/"+str(chrom)+".inova.frq")	
	
	#AF_file = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/AF/final_"+chrom+".txt.gz")
	AF_file = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/FINAL_AF/final_"+chrom+".txt.gz")

	annotation_files = []
	annotation_file1 = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/byChrom/" +str(chrom) +".txt.gz"
	annotation_file2 = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/byChrom_file2/" +str(chrom) +".txt.gz"
	annotation_files.append(annotation_file2)
	annotation_files.append(annotation_file1)

	#Index the frequency file based on 0 index
	freq_dict = {}
	snps_count = 0
	ins_count = 0
	del_count = 0
	other_count = 0
	
	snps_count_inova = 0
	ins_count_inova = 0
	del_count_inova = 0
	other_count_inova = 0

	counter = 0
	counter_dictionary = {}
	inova_counter_dictionary = {}
	thousand_dictionary = {}

	counter_dictionary["Total_snps_original_data"] = [0,0,0]
	counter_dictionary["Total_ins_original_data"] = [0,0,0]
	counter_dictionary["Total_del_original_data"] = [0,0,0]
	counter_dictionary["Total_other_original_data"] = [0,0,0]

	inova_counter_dictionary["Total_snps_original_data_inova"] = [0,0,0]
	inova_counter_dictionary["Total_ins_original_data_inova"] = [0,0,0]
	inova_counter_dictionary["Total_del_original_data_inova"] = [0,0,0]
	inova_counter_dictionary["Total_other_original_data_inova"] = [0,0,0]

	for line in AF_file:

		counter += 1	
		tp_line = line.strip().split()
		#Check to see if this is a snp
		if len(tp_line[3]) == 1 and len(tp_line[4]) == 1:
			begin = int(tp_line[1]) - 1
			dict_key = str(begin)+"_"+tp_line[3]
			dict_value_array = "_".join(tp_line[4:])
			freq_dict[dict_key] = dict_value_array
			if float(tp_line[5]) > 0.0:
				snps_count += 1
				counter_dictionary = add_counts_original_data(float(tp_line[5]), counter_dictionary, "Total_snps_original_data")
			if float(tp_line[6]) > 0.0:
				snps_count_inova +=1
				inova_counter_dictionary = add_counts_original_data(float(tp_line[6]), inova_counter_dictionary, "Total_snps_original_data_inova")
		#Check to see if this is snp with multiallele
		elif len(tp_line[3]) == 1 and ',' in tp_line[4]:

			alleles = tp_line[4].split(",")
			this_is_snp = True
			for i in alleles:
				#If the lenght of the allele is bigger the 1 this is a insertions
				if len(i) > 1:
					this_is_snp = False

			#if this is a snp
			if this_is_snp:
				'''
				if float(tp_line[5]) > 0.0:
					snps_count += 1
				if float(tp_line[6]) > 0.0:
					snps_count_inova +=1
				'''
				if float(tp_line[5]) > 0.0:
					snps_count += 1
					counter_dictionary = add_counts_original_data(float(tp_line[5]), counter_dictionary, "Total_snps_original_data")
				if float(tp_line[6]) > 0.0:
					snps_count_inova +=1
					inova_counter_dictionary = add_counts_original_data(float(tp_line[6]), inova_counter_dictionary, "Total_snps_original_data_inova")

				begin = int(tp_line[1]) - 1
				dict_key = str(begin)+"_"+tp_line[3]
				#dict_value_array = []
				#dict_value_array.append(tp_line[4:])
				dict_value_array = "_".join(tp_line[4:])
				freq_dict[dict_key] = dict_value_array

			#Else this is a insertion
			else:
				#or maybe other
				if tp_line[4][0] == tp_line[3]:

					'''
					if float(tp_line[5]) > 0.0:
						ins_count += 1
					if float(tp_line[6]) > 0.0:
						ins_count_inova +=1
					'''
					if float(tp_line[5]) > 0.0:
						ins_count += 1
						counter_dictionary = add_counts_original_data(float(tp_line[5]), counter_dictionary, "Total_ins_original_data")
					if float(tp_line[6]) > 0.0:
						ins_count_inova +=1
						inova_counter_dictionary = add_counts_original_data(float(tp_line[6]), inova_counter_dictionary, "Total_ins_original_data_inova")

					begin = tp_line[1]
					#Don't even have to trim the alts, we will just look to see if we find it
					dict_key=begin+"_-"
					#dict_value_array = []
					#dict_value_array.append(tp_line[4:])
					dict_value_array = "_".join(tp_line[4:])
					freq_dict[dict_key] = dict_value_array
				else:
					print "indel "+line
					'''
					if float(tp_line[5]) > 0.0:
						other_count += 1
					if float(tp_line[6]) > 0.0:
						other_count_inova +=1
					'''
					if float(tp_line[5]) > 0.0:
						other_count += 1
						counter_dictionary = add_counts_original_data(float(tp_line[5]), counter_dictionary, "Total_other_original_data")
					if float(tp_line[6]) > 0.0:
						other_count_inova +=1
						inova_counter_dictionary = add_counts_original_data(float(tp_line[6]), inova_counter_dictionary, "Total_other_original_data_inova")

					begin = int(tp_line[1]) - 1
					#Don't even have to trim the alts, we will just look to see if we find it
					dict_key=str(begin)+"_"+tp_line[3]
					#dict_value_array = []
					#dict_value_array.append(tp_line[4:])
					dict_value_array = "_".join(tp_line[4:])
					freq_dict[dict_key] = dict_value_array

		#this is a insertion (OR delins if the first letter doens't match)
		elif len(tp_line[3]) == 1 and len(tp_line[4]) > 1:
			#Verify that this is a insertion and not a delins, first letter
			#of the alt is same as the ref
			if tp_line[4][0] == tp_line[3]:
				'''
				if float(tp_line[5]) > 0.0:
					ins_count += 1
				if float(tp_line[6]) > 0.0:
					ins_count_inova +=1
				'''
				if float(tp_line[5]) > 0.0:
					ins_count += 1
					counter_dictionary = add_counts_original_data(float(tp_line[5]), counter_dictionary, "Total_ins_original_data")
				if float(tp_line[6]) > 0.0:
					ins_count_inova +=1
					inova_counter_dictionary = add_counts_original_data(float(tp_line[6]), inova_counter_dictionary, "Total_ins_original_data_inova")

				begin = tp_line[1]
				#Don't even have to trim the alts, we will just look to see if we find it
				dict_key=begin+"_-"
				#dict_value_array = []
				#dict_value_array.append(tp_line[4:])
				dict_value_array = "_".join(tp_line[4:])
				freq_dict[dict_key] = dict_value_array
			else:
				#Treat it as delin
				print "indel "+line
				'''
				if float(tp_line[5]) > 0.0:
					other_count += 1
				if float(tp_line[6]) > 0.0:
					other_count_inova +=1
				'''

				if float(tp_line[5]) > 0.0:
					other_count += 1
					counter_dictionary = add_counts_original_data(float(tp_line[5]), counter_dictionary, "Total_other_original_data")
				if float(tp_line[6]) > 0.0:
					other_count_inova +=1
					inova_counter_dictionary = add_counts_original_data(float(tp_line[6]), inova_counter_dictionary, "Total_other_original_data_inova")

				begin = int(tp_line[1]) - 1
				#Don't even have to trim the alts, we will just look to see if we find it
				dict_key=str(begin)+"_"+tp_line[3]
				#dict_value_array = []
				#dict_value_array.append(tp_line[4:])
				dict_value_array = "_".join(tp_line[4:])
				freq_dict[dict_key] = dict_value_array

		#this is a deletion
		elif len(tp_line[3]) > 1 and len(tp_line[4]) == 1:
			

			#This also migth be a delin, verify the first character
			if tp_line[4] == tp_line[3][0]:
				print "del "+line
				'''
				if float(tp_line[5]) > 0.0:
					del_count += 1
				if float(tp_line[6]) > 0.0:
					del_count_inova +=1
				'''
				if float(tp_line[5]) > 0.0:
					del_count += 1
					counter_dictionary = add_counts_original_data(float(tp_line[5]), counter_dictionary, "Total_del_original_data")
				if float(tp_line[6]) > 0.0:
					del_count_inova +=1
					inova_counter_dictionary = add_counts_original_data(float(tp_line[6]), inova_counter_dictionary, "Total_del_original_data_inova")

				begin = tp_line[1]
				ref = tp_line[3][1:]
				dict_key = begin + "_" + ref
				dict_value_array = "-_"+"_".join(tp_line[5:])
				freq_dict[dict_key] = dict_value_array
			else:
				#This is a delin
				print "indel "+line
				'''
				if float(tp_line[5]) > 0.0:
					other_count += 1
				if float(tp_line[6]) > 0.0:
					other_count_inova +=1
				'''

				if float(tp_line[5]) > 0.0:
					other_count += 1
					counter_dictionary = add_counts_original_data(float(tp_line[5]), counter_dictionary, "Total_other_original_data")
				if float(tp_line[6]) > 0.0:
					other_count_inova +=1
					inova_counter_dictionary = add_counts_original_data(float(tp_line[6]), inova_counter_dictionary, "Total_other_original_data_inova")

				begin = int(tp_line[1])-1 
				dict_key = str(begin) + "_" + tp_line[3]
				dict_value_array = tp_line[4] + "_"+"_".join(tp_line[5:])
				freq_dict[dict_key] = dict_value_array
		#Do some trimming for indels
		else:
			print "indel "+line
			'''
			if float(tp_line[5]) > 0.0:
				other_count += 1
			if float(tp_line[6]) > 0.0:
				other_count_inova +=1
			'''
			if float(tp_line[5]) > 0.0:
				other_count += 1
				counter_dictionary = add_counts_original_data(float(tp_line[5]), counter_dictionary, "Total_other_original_data")
			if float(tp_line[6]) > 0.0:
				other_count_inova +=1
				inova_counter_dictionary = add_counts_original_data(float(tp_line[6]), inova_counter_dictionary, "Total_other_original_data_inova")

			begin = int(tp_line[1]) - 1
			ref = tp_line[3]
			var = tp_line[4]
			final_index = 0

			start = 0
			for i in range(0,min( len(var),len(ref)) ):
				if var[i] == ref[i]:
					start += 1
				else:
					break

			ref = ref[start:]

			if start > 0 :
				begin = begin + start

			dict_key = str(begin)+"_"+ref

			#dict_value_array.append(tp_line[4:])
			dict_value_array = "_".join(tp_line[4:])
			freq_dict[dict_key] = dict_value_array			

	print "Total snps " + str(snps_count)
	print "Total ins " + str(ins_count)
	print "Total del " + str(del_count)
	print "Total other " + str(other_count)
	'''
	counter_dictionary["Total_snps_original_data"] = snps_count
	counter_dictionary["Total_ins_original_data"] = ins_count
	counter_dictionary["Total_del_original_data"] = del_count
	counter_dictionary["Total_other_original_data"] = other_count
	'''
	print "Total snps inova " + str(snps_count_inova)
	print "Total ins inova " + str(ins_count_inova)
	print "Total del inova " + str(del_count_inova)
	print "Total other inova " + str(other_count_inova)
	'''
	inova_counter_dictionary["Total_snps_original_data_inova"] = snps_count_inova
	inova_counter_dictionary["Total_ins_original_data_inova"] = ins_count_inova
	inova_counter_dictionary["Total_del_original_data_inova"] = del_count_inova
	inova_counter_dictionary["Total_other_original_data_inova"] = other_count_inova
	'''
	for key in counter_dictionary:
		val = counter_dictionary[key]
		print key
		for i in val:
			print str(i)

	for key in inova_counter_dictionary:
		val = inova_counter_dictionary[key]
		print key
		for i in val:
			print str(i)


	#temp_counter = 0
	line_counter = 0
	for annot_file in annotation_files:
		print annot_file
		afile = gzip.open(annot_file)
		for line in afile:
			line_counter += 1
			if line_counter % 100000 == 0:
				print str(line_counter)

			if counter == 100:
				print "break"
				break
			tp_line = line.strip().split("\t")
			#Verify if this is even found in the wellderly
			#If found return the allele frequency of wellderly, inova and boolean if this is novel
			wel, inova, thous, freq_dict = check_presence_in_filtered_dataset(tp_line, tp_line[8], freq_dict)
			#print "Wellderly AF " + str(wel)
			#print "Inova AF " + str(inova)
			#print "Thousand genomes AF " + str(thous)
			if wel > 0.0  or inova >0.0:

				splice_site = tp_line[10].replace('///','')
				splice_site =splice_site.replace('-','')

				if tp_line[3] == "snp":
					#Count how many snps were found in the annotatio
					if wel > 0.0:
						counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_snps" , wel, thousand_dictionary, thous)
					if inova > 0.0:
						inova_counter_dictionary, thousand_dictionary = add_counts(inova_counter_dictionary, "inova_snps" , inova, thousand_dictionary, thous)

					
					if 'Nonsense' in tp_line[7]:
						if wel > 0.0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_Nonsense" , wel, thousand_dictionary, thous)
						if inova > 0.0:
							inova_counter_dictionary, thousand_dictionary = add_counts(inova_counter_dictionary, "inova_Nonsense" , inova, thousand_dictionary, thous)


					elif not splice_site == '':
						if wel > 0.0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_splice_site" , wel, thousand_dictionary, thous)
						if inova > 0.0:
							inova_counter_dictionary, thousand_dictionary = add_counts(inova_counter_dictionary, "inova_splice_site" , inova, thousand_dictionary, thous)

					elif 'Nonsynonymous' in tp_line[7]:
						if wel > 0.0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_Nonsynonymous" , wel, thousand_dictionary, thous)
						if inova > 0.0:
							inova_counter_dictionary, thousand_dictionary = add_counts(inova_counter_dictionary, "inova_Nonsynonymous" , inova, thousand_dictionary, thous)

					elif 'Synonymous' in tp_line[7]:
						if wel > 0.0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_Synonymous" , wel, thousand_dictionary, thous)
						if inova > 0.0:
							inova_counter_dictionary, thousand_dictionary = add_counts(inova_counter_dictionary, "inova_Synonymous" , inova, thousand_dictionary, thous)

					elif 'Intron' in tp_line[9]:
						if wel > 0.0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_Intron" , wel, thousand_dictionary, thous)
						if inova > 0.0:
							inova_counter_dictionary, thousand_dictionary = add_counts(inova_counter_dictionary, "inova_Intron" , inova, thousand_dictionary, thous)

					else:
						#print "intergenic? "+ tp_line[9]
						if wel > 0.0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_Intergenic" , wel, thousand_dictionary, thous)
						if inova > 0.0:
							inova_counter_dictionary, thousand_dictionary = add_counts(inova_counter_dictionary, "inova_Intergenic" , inova, thousand_dictionary, thous)


				else:
					if tp_line[3] == "ins":
						if wel > 0.0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_ins" , wel, thousand_dictionary, thous)
						if inova > 0.0:
							inova_counter_dictionary, thousand_dictionary = add_counts(inova_counter_dictionary, "inova_ins" , inova, thousand_dictionary, thous)

					elif tp_line[3] == "del":
						if wel > 0.0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_del" , wel, thousand_dictionary, thous)
						if inova > 0.0:
							inova_counter_dictionary, thousand_dictionary = add_counts(inova_counter_dictionary, "inova_del" , inova, thousand_dictionary, thous)

					else:
						if wel > 0.0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_indel" , wel, thousand_dictionary, thous)
						if inova > 0.0:
							inova_counter_dictionary, thousand_dictionary = add_counts(inova_counter_dictionary, "inova_indel" , inova, thousand_dictionary, thous)


					if 'Frameshift' in tp_line[7]:
						if wel > 0.0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_complex_frameshift" , wel, thousand_dictionary, thous)
						if inova > 0.0:
							inova_counter_dictionary, thousand_dictionary = add_counts(inova_counter_dictionary, "inova_complex_frameshift" , inova, thousand_dictionary, thous)


					elif 'In_Frame' in tp_line[7]:
						if wel > 0.0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_complex_inFrame" , wel, thousand_dictionary, thous)
						if inova > 0.0:
							inova_counter_dictionary, thousand_dictionary = add_counts(inova_counter_dictionary, "inova_complex_inFrame" , inova, thousand_dictionary, thous)

					elif not splice_site == '':
						if wel > 0.0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_complex_SpliceSite" , wel, thousand_dictionary, thous)
						if inova > 0.0:
							inova_counter_dictionary, thousand_dictionary = add_counts(inova_counter_dictionary, "inova_complex_SpliceSite" , inova, thousand_dictionary, thous)


					elif 'Intron' in tp_line[9]:
						if wel > 0.0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_complex_Intron" , wel, thousand_dictionary, thous)
						if inova > 0.0:
							inova_counter_dictionary, thousand_dictionary = add_counts(inova_counter_dictionary, "inova_complex_Intron" , inova, thousand_dictionary, thous)
					
					else:
						if wel > 0.0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_complex_Intergenic" , wel, thousand_dictionary, thous)
						if inova > 0.0:
							inova_counter_dictionary, thousand_dictionary = add_counts(inova_counter_dictionary, "inova_complex_Intergenic" , inova, thousand_dictionary, thous)
				
		afile.close()

	counter_line = ""
	for key in counter_dictionary:
		val = counter_dictionary[key]
		counter_line = counter_line+ chrom +"\t"+key
		print key
		for i in val:
			counter_line = counter_line + "\t" + str(i)
			print str(i)
		counter_line = counter_line + "\n"

	for key in inova_counter_dictionary:
		val = inova_counter_dictionary[key]
		counter_line = counter_line+ chrom +"\t"+key
		print key
		for i in val:
			counter_line = counter_line + "\t" + str(i)
			print str(i)
		counter_line = counter_line + "\n"

	counter_file = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/table2.counter.final", 'a')
	counter_file.write(counter_line)
	counter_file.close()
	
	'''
	thousand_line = ""
	for key in thousand_dictionary:
		val = thousand_dictionary[key]
		data_val = counter_dictionary[key]
		thousand_line = thousand_line + chrom +"\t"+ key
		print key
		for index, i in enumerate(val):
			data = data_val[index]
			try:
				perc = float(i)/float(data)
			except:
				perc = 0
			thousand_line = thousand_line + "\t" + str(perc) 
			print str(perc)
		thousand_line = thousand_line + "\n"


	counter_file = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/table2.counter.thousandG", 'a')
	counter_file.write(thousand_line)
	counter_file.close()
	'''
	
	#Check to see if there are any variants that weren't found in annotation
	counter_not_found = 0
	for dict_key in freq_dict.keys():
		counter_not_found += 1
		print "Not found"
		print "Dict_key " + dict_key
		print "Value " + freq_dict[dict_key]

	print "Total variants not found in the annotation " + str(counter_not_found)
	print "Total variants at start " + str(counter)

def check_presence_in_filtered_dataset(tp_line, thousand_af, freq_dict):
	welld_found = 0.0
	inova_found = 0.0
	thousand_novel = 0
	
	dict_key = tp_line[1]+"_"+tp_line[4]
	#print dict_key
	annotation_alt = tp_line[5]
	try:
		found_something = freq_dict[dict_key]
		tp_found = found_something.split("_")
		#print "found something " + str(tp_found[0])

		
		alt = tp_found[0]
		#print "Alt allele " + str(alt)

		if annotation_alt in alt:
			#print "annotation alt " + annotation_alt
			#print "alt " + alt

			inova_found = float(tp_found[2])
			welld_found = float(tp_found[1]) 
			
			#print "Wellderly AF " + str(welld_found)
			#print "Inova AF " + str(inova_found)
			#remove this entry from the dictionary so the multiple alleles
			#won't be annotated multiple times
			del freq_dict[dict_key]
		
		#print "Wellderly AF " + str(welld_found)
		#print "Inova AF " + str(inova_found)

	except:
		welld_found = 0.0
		inova_found = 0.0

	#check 1000 genomes only if found in either wellderly or inova
	novel=False
	if welld_found >0.0 or inova_found>0.0:
		#print "Welld AF " + str(welld_found)
		#print "Inova AF " + str(inova_found)
		
		#if thousand novel column has a valid entry it would be fload
		try:
			thousand_novel = float(thousand_af)
		#except this is a novel variant
		except:
			thousand_novel = 0.0

		
		if thousand_novel == 0.0:
			novel = True

	return welld_found, inova_found, novel, freq_dict

def add_counts_original_data(alle_freq, counts_dict, dict_key):

	count_array = counts_dict[dict_key]

	if alle_freq < 0.01:
		new_freq = count_array[0] + 1
		count_array[0] = new_freq

	elif alle_freq < 0.05:
		new_freq = count_array[1] + 1
		count_array[1] = new_freq

	else:
		new_freq = count_array[2] + 1
		count_array[2] = new_freq

	counts_dict[dict_key] = count_array
	return counts_dict


def add_counts(counts_dict, dict_key, alle_freq, thousandg_dict, novel):
	#Extract the count array if this is not the first round
	try:
		count_array = counts_dict[dict_key]
	except:
		count_array = [0,0,0]

	#Extract the novel array
	try:
		novel_array = thousandg_dict[dict_key]
	except:
		novel_array = [0,0,0]			

	if alle_freq < 0.01:
		new_freq = count_array[0] + 1
		count_array[0] = new_freq

		if novel:
			novel_count = novel_array[0]+1
			novel_array[0] = novel_count

	elif alle_freq < 0.05:
		new_freq = count_array[1] + 1
		count_array[1] = new_freq

		if novel:
			novel_count = novel_array[1]+1
			novel_array[1] = novel_count
	else:
		new_freq = count_array[2] + 1
		count_array[2] = new_freq

		if novel:
			novel_count = novel_array[2]+1
			novel_array[2] = novel_count

	counts_dict[dict_key] = count_array
	thousandg_dict[dict_key] = novel_array

	return counts_dict, thousandg_dict


if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)
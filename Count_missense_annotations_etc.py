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

	plink_filename_freq_wellderly=open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/Plink_files/"+str(chrom)+".wellderly.frq")
	plink_filename_freq_inova=open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/Plink_files/"+str(chrom)+".inova.frq")	


	annotation_files = []
	annotation_file1 = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/byChrom/" +str(chrom) +".txt.gz"
	annotation_file2 = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/byChrom_file2/" +str(chrom) +".txt.gz"
	annotation_files.append(annotation_file2)
	annotation_files.append(annotation_file1)

	#Index the vcf file first 0 based 
	wellderly_freq_dict = {}
	wellderly_snps_count = 0
	wellderly_ins_count = 0
	wellderly_del_count = 0
	wellderly_other_count = 0

	counter = 0

	for line in plink_filename_freq_wellderly:
		counter += 1	
		#Is this the header
		if counter == 1:
			continue
		else:
			tp_line = line.strip().split()
			#Verify if AF==0 don't store it
			if float(tp_line[4]) == 0.0:
				continue
			else:
				#Is this a snp
				if len(tp_line[2]) == 1 and len(tp_line[3]) == 1:
					wellderly_snps_count +=1
					pos = tp_line[1].strip().split("-")
					key = pos[0].strip() +"_"+tp_line[2].strip()+"_"+tp_line[3].strip()
					#print key
					wellderly_freq_dict[key] = tp_line[4]
					#print key
				#Is this a ins or del, begin position is the same
				elif len(tp_line[2]) == 1 and len(tp_line[3]) > 1:
					wellderly_del_count +=1
					pos = tp_line[1].strip().split("-")
					key = pos[0] +"_"+tp_line[3][1:]+"_-"
					#print "Key deletion " + key
					wellderly_freq_dict[key] = tp_line[4]
				elif len(tp_line[2]) > 1 and len(tp_line[3]) == 1:
					wellderly_ins_count +=1
					pos = tp_line[1].strip().split("-")
					key = pos[0] +"_-_"+tp_line[2][1:]
					#print "Key insertion " + key
					wellderly_freq_dict[key] = tp_line[4]
				#This is delins
				elif len(tp_line[2]) > 1 and len(tp_line[3]) > 1:
					wellderly_del_count +=1
					pos = tp_line[1].strip().split("-")
					key = pos[0] +"_"+tp_line[2][1:]+"_"+tp_line[3][1:]
					#print "Key delins " + key
					wellderly_freq_dict[key] = tp_line[4]
				else:
					print "No idea wtf is this"
					print line


	print "total snps: "+str(wellderly_snps_count)
	print "Total ins: "+str(wellderly_ins_count)
	print "Total del "+str(wellderly_del_count)
	print "Total other "+str(wellderly_other_count)
	print "Size of welderly dictionary " + str(len(wellderly_freq_dict))
	#NOTE!!! need to differentiate between del and ins from annotation file not from plink file
	#Plink files alleles can be flipped
	plink_filename_freq_wellderly.close()

	#Index inova frequency
	inova_freq_dict = {}
	inova_snps_count = 0
	inova_ins_count = 0
	inova_del_count = 0
	inova_other_count = 0

	counter = 0
	for line in plink_filename_freq_inova:
		counter += 1	
		#Is this the header
		if counter == 1:
			continue
		else:
			tp_line = line.strip().split()
			#Verify if AF==0 don't store it
			if float(tp_line[4]) == 0.0:
				continue
			else:
				#Is this a snp
				if len(tp_line[2]) == 1 and len(tp_line[3]) == 1:
					inova_snps_count +=1
					pos = tp_line[1].strip().split("-")
					key = pos[0].strip() +"_"+tp_line[2].strip()+"_"+tp_line[3].strip()
					inova_freq_dict[key] = tp_line[4]
				#Is this a ins or del, begin position is the same
				elif len(tp_line[2]) == 1 and len(tp_line[3]) > 1:
					inova_del_count +=1
					pos = tp_line[1].strip().split("-")
					key = pos[0] +"_"+tp_line[3][1:]+"_-"
					#print "Key deletion " + key
					inova_freq_dict[key] = tp_line[4]
				elif len(tp_line[2]) > 1 and len(tp_line[3]) == 1:
					inova_ins_count +=1
					pos = tp_line[1].strip().split("-")
					key = pos[0] +"_-_"+tp_line[2][1:]
					#print "Key insertion " + key
					inova_freq_dict[key] = tp_line[4]
				#This is delins
				elif len(tp_line[2]) > 1 and len(tp_line[3]) > 1:
					inova_del_count +=1
					pos = tp_line[1].strip().split("-")
					key = pos[0] +"_"+tp_line[2][1:]+"_"+tp_line[3][1:]
					#print "Key delins " + key
					inova_freq_dict[key] = tp_line[4]
				else:
					print "No idea wtf is this"
					print line

	print "total inova snps: "+str(inova_snps_count)
	print "Total inova ins: "+str(inova_ins_count)
	print "Total inova del "+str(inova_del_count)
	print "Total inova other "+str(inova_other_count)
	print "Size of inova dictionary " + str(len(inova_freq_dict))
	#NOTE!!! need to differentiate between del and ins from annotation file not from plink file
	#Plink files alleles can be flipped
	plink_filename_freq_inova.close()

	counter_dictionary = {}
	thousand_dictionary = {}
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
			wel, inova, thous = check_presence_in_filtered_dataset(tp_line, tp_line[8], wellderly_freq_dict, inova_freq_dict)
			
			if wel > 0  or inova >0:

				splice_site = tp_line[10].replace('///','')
				splice_site =splice_site.replace('-','')

				if tp_line[3] == "snp":
					#Count how many snps were found in the annotatio
					if wel > 0:
						counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_snps" , wel, thousand_dictionary, thous)
					if inova > 0:
						counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "inova_snps" , inova, thousand_dictionary, thous)

					
					if 'Nonsense' in tp_line[7]:
						if wel > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_Nonsense" , wel, thousand_dictionary, thous)
						if inova > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "inova_Nonsense" , inova, thousand_dictionary, thous)


					elif not splice_site == '':
						if wel > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_splice_site" , wel, thousand_dictionary, thous)
						if inova > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "inova_splice_site" , inova, thousand_dictionary, thous)

					elif 'Nonsynonymous' in tp_line[7]:
						if wel > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_Nonsynonymous" , wel, thousand_dictionary, thous)
						if inova > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "inova_Nonsynonymous" , inova, thousand_dictionary, thous)

					elif 'Synonymous' in tp_line[7]:
						if wel > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_Synonymous" , wel, thousand_dictionary, thous)
						if inova > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "inova_Synonymous" , inova, thousand_dictionary, thous)

					elif 'Intron' in tp_line[9]:
						if wel > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_Intron" , wel, thousand_dictionary, thous)
						if inova > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "inova_Intron" , inova, thousand_dictionary, thous)

					elif 'Exon' in tp_line[9] or 'Intron' in tp_line[9]:
						continue
					#if it was neither Exon or Intron in location, this must be intergenic	
					else:
						#print "intergenic? "+ tp_line[9]
						if wel > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_Intergenic" , wel, thousand_dictionary, thous)
						if inova > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "inova_Intergenic" , inova, thousand_dictionary, thous)


				else:
					if tp_line[3] == "ins":
						if wel > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_ins" , wel, thousand_dictionary, thous)
						if inova > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "inova_ins" , inova, thousand_dictionary, thous)

					elif tp_line[3] == "del":
						if wel > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_del" , wel, thousand_dictionary, thous)
						if inova > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "inova_del" , inova, thousand_dictionary, thous)


					if 'Frameshift' in tp_line[7]:
						if wel > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_complex_frameshift" , wel, thousand_dictionary, thous)
						if inova > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "inova_complex_frameshift" , inova, thousand_dictionary, thous)


					elif 'In_Frame' in tp_line[7]:
						if wel > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_complex_inFrame" , wel, thousand_dictionary, thous)
						if inova > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "inova_complex_inFrame" , inova, thousand_dictionary, thous)

					elif not splice_site == '':
						if wel > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_complex_SpliceSite" , wel, thousand_dictionary, thous)
						if inova > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "inova_complex_SpliceSite" , inova, thousand_dictionary, thous)


					elif 'Intron' in tp_line[9]:
						if wel > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_complex_Intron" , wel, thousand_dictionary, thous)
						if inova > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "inova_complex_Intron" , inova, thousand_dictionary, thous)
					
					elif 'Exon' in tp_line[9] or 'Intron' in tp_line[9]:
						continue

					else:
						if wel > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "well_complex_Intergenic" , wel, thousand_dictionary, thous)
						if inova > 0:
							counter_dictionary, thousand_dictionary = add_counts(counter_dictionary, "inova_complex_Intergenic" , inova, thousand_dictionary, thous)
					

					


		afile.close()

	counter_line = ""
	for key in counter_dictionary:
		val = counter_dictionary[key]
		counter_line = counter_line+ chrom +"\t"+key
		print key
		for i in val:
			counter_line = counter_line + "\t" + str(i)
			print i
		counter_line = counter_line + "\n"

	counter_file = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/table2.counter.final", 'a')
	counter_file.write(counter_line)
	counter_file.close()
	
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
			print i
		thousand_line = thousand_line + "\n"


	counter_file = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/table2.counter.thousandG", 'a')
	counter_file.write(thousand_line)
	counter_file.close()


def check_presence_in_filtered_dataset(tp_line, thousand_af, wellderly_freq_dict, inova_freq_dict):
	welld_found = 0
	inova_found = 0
	thousand_novel = 0
	
	#chrom = tp_line[1]
	#if this is a snp, use the end position:
	if tp_line[3] == "snp":
		key = tp_line[2]+"_"+tp_line[4]+"_"+tp_line[5]
		#print key
		key2 = tp_line[2]+"_"+tp_line[5]+"_"+tp_line[4]
	else:

		#deletions and insertions were already switched, live key1 and key2 as they are
		key = tp_line[1]+"_"+tp_line[4]+"_"+tp_line[5]
		#print key
		key2 = tp_line[1]+"_"+tp_line[4]+"_"+tp_line[5]

	try:
		#Looking for second key first, more probably by plink
		welld_found = float(wellderly_freq_dict[key2])
		#print "welld found! " + val
		#welld_found = 1
		#Remove key, in case of duplicates this key won't be found again
		del wellderly_freq_dict[key2]


	except:
		try:
			welld_found = float(wellderly_freq_dict[key])
			#print "welld found! " + val
			#welld_found = 1
			del wellderly_freq_dict[key]
		except:
			welld_found = 0.0


	try:
		inova_found = float(inova_freq_dict[key2])
		#print "Inova found! " + val
		#inova_found = 1
		del inova_freq_dict[key2]
	except:
		try:
			inova_found = float(inova_freq_dict[key])
			#print "Inova found! " + val
			#inova_found = 1
			del inova_freq_dict[key]
		except:
			inova_found = 0.0

	#check 1000 genomes only if found in either wellderly or inova
	novel=False
	if welld_found >0.0 or inova_found>0.0:

		#if thousand novel column has a valid entry it would be fload
		try:
			thousand_novel = float(thousand_af)
		#except this is a novel variant
		except:
			thousand_novel = 0.0

		
		if thousand_novel == 0.0:
			novel = True

	return welld_found, inova_found, novel

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
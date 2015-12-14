"""
Extract cat1 annotations

"""

import os, sys, gzip, datetime
import cPickle
import resource
import gc
import time


def main(chrom):
	ch = chrom[3:]
	print "chrom is: " + ch

	AF_file = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/FINAL_AF/final_"+chrom+".txt.gz")
	#AF_file = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/FINAL_wellderly_all_filters_withHWEandAD.noZeroAF"+chrom+".vcf.gz")

	annotation_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/cat1_final/cat1_part"
	#outputfile = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/cat1_final/results_annot_cat1_plus_vcf.txt", "a")
	outputfile = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/cat1_final/AF_results_annot_cat1_plus_vcf.txt", "a")

	freq_dict = {}
	freq_dict_AF = {}
	snps_count=0
	ins_count=0
	del_count=0
	other_count=0
	header = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/cat1_final/vcf_header", "w")
	
	'''
	for line in AF_file:

		if line[0] == "#":
			header.write(line)
			continue
	
		tp_line = line.strip().split()
		dict_value_array = []
		#Check to see if this is a snp
		if len(tp_line[3]) == 1 and len(tp_line[4]) == 1:
			begin = int(tp_line[1]) - 1
			dict_key = str(begin)+"_"+tp_line[3]
			#print dict_key
			
			dict_value_array.append(tp_line[4])
			dict_value_array.append(line)
			freq_dict[dict_key] = dict_value_array
			snps_count += 1
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
				snps_count += 1
				begin = int(tp_line[1]) - 1
				dict_key = str(begin)+"_"+tp_line[3]
				dict_value_array.append(tp_line[4])
				dict_value_array.append(line)
				freq_dict[dict_key] = dict_value_array

			#Else this is a insertion
			else:
				#or maybe other
				if tp_line[4][0] == tp_line[3]:
					ins_count += 1
					begin = tp_line[1]
					#Don't even have to trim the alts, we will just look to see if we find it
					dict_key=begin+"_-"
					dict_value_array.append(tp_line[4])
					dict_value_array.append(line)
					freq_dict[dict_key] = dict_value_array
				else:
					#Treat it as delin
					other_count += 1
					begin = int(tp_line[1]) - 1
					#Don't even have to trim the alts, we will just look to see if we find it
					dict_key=str(begin)+"_"+tp_line[3]
					dict_value_array.append(tp_line[4])
					dict_value_array.append(line)
					freq_dict[dict_key] = dict_value_array

		#this is a insertion (OR delins if the first letter doens't match)
		elif len(tp_line[3]) == 1 and len(tp_line[4]) > 1:
			#Verify that this is a insertion and not a delins, first letter
			#of the alt is same as the ref
			if tp_line[4][0] == tp_line[3]:
				ins_count += 1
				begin = tp_line[1]
				#Don't even have to trim the alts, we will just look to see if we find it
				dict_key=begin+"_-"
				dict_value_array.append(tp_line[4])
				dict_value_array.append(line)
				freq_dict[dict_key] = dict_value_array
			else:
				#Treat it as delin
				other_count += 1
				begin = int(tp_line[1]) - 1
				#Don't even have to trim the alts, we will just look to see if we find it
				dict_key=str(begin)+"_"+tp_line[3]
				dict_value_array.append(tp_line[4])
				dict_value_array.append(line)
				freq_dict[dict_key] = dict_value_array

		#this is a deletion
		elif len(tp_line[3]) > 1 and len(tp_line[4]) == 1:
			

			#This also migth be a delin, verify the first character
			if tp_line[4] == tp_line[3][0]:
				#This is insertion
				del_count +=1
				begin = tp_line[1]
				ref = tp_line[3][1:]
				dict_key = begin + "_" + ref
				dict_value_array.append("-")
				dict_value_array.append(line)
				freq_dict[dict_key] = dict_value_array
			else:
				#This is a delin
				other_count += 1
				begin = int(tp_line[1])-1
				dict_key = str(begin) + "_" + tp_line[3]
				dict_value_array.append(tp_line[4])
				dict_value_array.append(line)
				freq_dict[dict_key] = dict_value_array
		#Do some trimming for indels
		else:
			other_count += 1
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

			dict_value_array.append(tp_line[4])
			dict_value_array.append(line)
			freq_dict[dict_key] = dict_value_array			
	'''
	counter = 0
	
	for line in AF_file:

		counter += 1	
		tp_line = line.strip().split()
		#Check to see if this is a snp
		if len(tp_line[3]) == 1 and len(tp_line[4]) == 1:
			begin = int(tp_line[1]) - 1
			dict_key = str(begin)+"_"+tp_line[3]
			dict_value_array = "_".join(tp_line[4:])
			freq_dict[dict_key] = dict_value_array
			snps_count += 1
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
				snps_count += 1
				begin = int(tp_line[1]) - 1
				dict_key = str(begin)+"_"+tp_line[3]
				dict_value_array = "_".join(tp_line[4:])
				freq_dict[dict_key] = dict_value_array

			#Else this is a insertion
			else:
				#or maybe other
				if tp_line[4][0] == tp_line[3]:
					ins_count += 1
					begin = tp_line[1]
					#Don't even have to trim the alts, we will just look to see if we find it
					dict_key=begin+"_-"
					dict_value_array = "_".join(tp_line[4:])
					freq_dict[dict_key] = dict_value_array
				else:
					#Treat it as delin
					other_count += 1
					begin = int(tp_line[1]) - 1
					#Don't even have to trim the alts, we will just look to see if we find it
					dict_key=str(begin)+"_"+tp_line[3]
					dict_value_array = "_".join(tp_line[4:])
					freq_dict[dict_key] = dict_value_array

		#this is a insertion (OR delins if the first letter doens't match)
		elif len(tp_line[3]) == 1 and len(tp_line[4]) > 1:
			#Verify that this is a insertion and not a delins, first letter
			#of the alt is same as the ref
			if tp_line[4][0] == tp_line[3]:
				ins_count += 1
				begin = tp_line[1]
				#Don't even have to trim the alts, we will just look to see if we find it
				dict_key=begin+"_-"
				dict_value_array = "_".join(tp_line[4:])
				freq_dict[dict_key] = dict_value_array
			else:
				#Treat it as delin
				other_count += 1
				begin = int(tp_line[1]) - 1
				#Don't even have to trim the alts, we will just look to see if we find it
				dict_key=str(begin)+"_"+tp_line[3]
				dict_value_array = "_".join(tp_line[4:])
				freq_dict[dict_key] = dict_value_array

		#this is a deletion
		elif len(tp_line[3]) > 1 and len(tp_line[4]) == 1:
			

			#This also migth be a delin, verify the first character
			if tp_line[4] == tp_line[3][0]:
				del_count +=1
				begin = tp_line[1]
				ref = tp_line[3][1:]
				dict_key = begin + "_" + ref
				dict_value_array = "-_"+"_".join(tp_line[5:])
				freq_dict[dict_key] = dict_value_array
			else:
				#This is a delin
				other_count += 1
				begin = int(tp_line[1])-1 
				dict_key = str(begin) + "_" + tp_line[3]
				dict_value_array = tp_line[4] + "_"+"_".join(tp_line[5:])
				freq_dict[dict_key] = dict_value_array
		#Do some trimming for indels
		else:
			other_count += 1
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
	

	header.close()
	AF_file.close()
	print "Total snps " + str(snps_count)
	print "Total ins " + str(ins_count)
	print "Total del " + str(del_count)
	print "Total other " + str(other_count)

	total_cat1 = 0
	cat1_found = 0
	multialleles = 0
	for i in range(1,4):

		annot = annotation_file +str(i)+".txt"
		annotFile = open(annot)
		for line in annotFile:
			line = line.strip()
			total_cat1 +=1
			tp_line = line.split("\t")
			#check if it's the same chromosome
			#file 1 and 2 have different indexes then 3
			if tp_line[1] == chrom or tp_line[0] == chrom:
				print "same chrom"
				if i == 3:
					dict_key = tp_line[1]+"_"+tp_line[4]
				else:
					dict_key = tp_line[2]+"_"+tp_line[5]
				#print dict_key
				if i == 3:
					annotation_alt = tp_line[5]
				else:
					annotation_alt = tp_line[6]
				
				try:
					found_something = freq_dict[dict_key]
					print "found something " + found_something
					
					#print "found something " + str(tp_found[0])
					#Version 1, vcf file
					'''
					alt = found_something[0]

					if annotation_alt in alt:
						cat1_found += 1
						if "," in alt:
							multialleles += 1
						if i == 3:
							#print tp_line[20:]
							final_line=found_something[1].strip()+"\t"+"\t".join(tp_line[:6])+"\t"+"\t".join(tp_line[20:]) + "\n"
						else:
							print tp_line[1:7]
							final_line=found_something[1].strip()+"\t"+"\t".join(tp_line[1:7])+"\t"+"\t".join(tp_line[8:53]) + "\t"+"\t".join(tp_line[53:108])+"\n"

						vcf_data = found_something[1].strip().split("\t")
						print "lenght of vcf file " + str(len(vcf_data))
						f_line = final_line.split("\t")
						print "lenght of final_line " + str(len(f_line)) 
						outputfile.write(final_line)
					'''
					#version 2, AF file
					tp_found = found_something.split("_")
					alt = tp_found[0]
					print "alt is "+ alt
					print "annotation alt is "+str(annotation_alt)
					if annotation_alt in alt:
						cat1_found += 1
						if "," in alt:
							multialleles += 1
						if i == 3:
							#print tp_line[20:]
							final_line=dict_key + "\t" + "\t".join(tp_found)+"\t"+"\t".join(tp_line[:6])+"\t"+"\t".join(tp_line[20:]) + "\n"
						else:
							print tp_line[1:7]
							final_line=dict_key + "\t" +"\t".join(tp_found)+"\t"+"\t".join(tp_line[1:7])+"\t"+"\t".join(tp_line[8:53]) + "\t"+"\t".join(tp_line[53:108])+"\n"

						vcf_data = found_something[1].strip().split("\t")
						print "lenght of vcf file " + str(len(vcf_data))
						f_line = final_line.split("\t")
						print "lenght of final_line " + str(len(f_line)) 
						outputfile.write(final_line)
					
				except:
					continue
		
		annotFile.close()

	outputfile.close()
	print "Total cat1 " + str(total_cat1)
	print "Total found "+ str(cat1_found)
	#counterfile = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/cat1_final/counter_file.txt", "a")
	counterfile = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/cat1_final/AF_counter_file.txt", "a")

	counter_f = chrom + "\t"+str(total_cat1)+"\t"+str(cat1_found)+"\t"+str(multialleles)+"\n"
	counterfile.write(counter_f)
	counterfile.close()

if __name__ == '__main__':

    print "Python Version: " + sys.version

    print "start"
    main(sys.argv[1])
    print "done"

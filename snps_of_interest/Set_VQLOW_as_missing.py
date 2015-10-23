"""
Remove VQLOW set as missing

"""

import os, sys, gzip, datetime 

def main():
	'''
	snp_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/desease_snps-corected.txt"
	vcf_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/filtered_snps.noRelated.withHead.withCov.vcf"
	#vcf_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/unfiltered_snps.txt"
	out_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/final_snp_file.txt"
	#not_found_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/snps_notFound_in_filtered_dataset.txt"
	#not_found_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/snps_notFound_in_unfiltered_dataset.txt"
	excluded_snps="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/other_excluded_snps.txt"
	'''
	snp_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/disease_snps_alzeimers_corrected.txt"
	vcf_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/filtered_snps_alzeimers.txt"
	out_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/final_snp_file.alzeimers.txt"
	excluded_snps="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/other_excluded_snps.alzeimerstxt"

	snps = open(snp_file)
	inp = open(vcf_file)
	out = open(out_file, "w")
	nf = open(excluded_snps, "w")

	dict_of_snps = {}
	for line in snps:
		tp_line = line.split("\t")
		dict_of_snps[tp_line[2]] = line

	counter = 0

	dict_of_AF = {}
	easy_count = 0
	bigger_length_ref = 0
	not_found = 0
	for line in inp:
		counter += 1

		if line.startswith("#"):
			out.write(line)
			continue

		line = line.replace("VQLOW",".")
		tp_line = line.split("\t")


		#Extract the snp from the vcf file, compare if the REF and Alt match 
		try:
			begin = tp_line[1]
			variant_line = dict_of_snps[begin]
			
			origin_line = variant_line.split("\t")
			if origin_line[3] == tp_line[3] and origin_line[4] == tp_line[4]: #and origin_line[4] == tp_line[4]:
				easy_count += 1
				out.write(line)
				continue
			#Effect allele is the reference, calculate AF of the reference in the wellderly and inova
			elif origin_line[4] == tp_line[3] and origin_line[3] == tp_line[4]:
				easy_count +=1
				out.write(line)
				continue
			elif "," in tp_line[4]:
				new_alt = tp_line[4].split(",")
				if origin_line[3] == tp_line[3]: #and origin_line[4] == tp_line[4]:
					for i in new_alt:
						if i == origin_line[4]:
							#print "multy alleles"
							#print origin_line
							#print tp_line[:9]
							easy_count += 1
							out.write(line)
							continue

				#Effect allele is the reference, calculate AF of the reference in the wellderly and inova
				elif origin_line[4] == tp_line[3]:
					for i in new_alt:
						if i == origin_line[3]:
							
							#print "multy alleles"
							#print origin_line
							#print tp_line[:9]
							easy_count += 1
							out.write(line)
							continue

			else:

				#out.write(line)
				nf.write(line)
				print "Not found"
				print origin_line
				not_found += 1
				continue
			
		#Variants not found in the wellderly dataset after filtering 
		except:
			print "Not found"
			print origin_line
			not_found += 1
			nf.write(line)



	print "Lines found " + str(counter)
	print "Easy lines " + str(easy_count)
	print "Bigger lenth ref " + str(bigger_length_ref)
	print "Not found " + str(not_found)
	snps.close()
	inp.close()
	out.close()
	#nf.close()



if __name__ == '__main__':

    print "Python Version: " + sys.version
    #args = parse_args()
    #import time
    #start = time.time()
    main()
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    #end = time.time()
    #print 'This took {0} seconds'.format(end - start)

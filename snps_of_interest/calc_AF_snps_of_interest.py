"""
Calculate AF snps of interest

"""

import os, sys, gzip, datetime 

def main():

	snp_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/desease_snps-corected.txt"
	vcf_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/filtered_snps.txt"
	#vcf_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/unfiltered_snps.txt"
	out_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/final_snp_file.txt"
	#not_found_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/snps_notFound_in_filtered_dataset.txt"
	#not_found_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/snps_notFound_in_unfiltered_dataset.txt"
	excluded_snps="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/other_excluded_snps.txt"
	
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
	for line in inp:
		counter += 1

		tp_line = line.split("\t")

		#Extract the snp from the vcf file, compare if the REF and Alt match 
		try:
			begin = tp_line[1]
			variant_line = dict_of_snps[begin]
			
			origin_line = variant_line.split("\t")
			if origin_line[3] == tp_line[3] and origin_line[4] == tp_line[4]: #and origin_line[4] == tp_line[4]:
				easy_count += 1
				out.write(line)

			#Effect allele is the reference, calculate AF of the reference in the wellderly and inova
			elif origin_line[4] == tp_line[3] and origin_line[3] == tp_line[4]:
				easy_count +=1
				out.write(line)
			else:
				#print tp_line[:7]
				#bigger_length_ref += 1
				nf.write(line)
				
			
		#Variants not found in the wellderly dataset after filtering 
		except:
			nf.write(line)



	print "Lines found " + str(counter)
	print "Easy lines " + str(easy_count)
	print "Bigger lenth ref " + str(bigger_length_ref)
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

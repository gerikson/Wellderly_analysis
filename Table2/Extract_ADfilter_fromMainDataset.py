"""
Extract the variants that were filtered out by the allele depth 
from the main dataset
"""
import os, sys, gzip, datetime


def main(chrom):



	filtered_out_rare_variants = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_rare_variants/vcf_rare_variant_filtered/"+str(chrom)+".vcf.gz")
	almost_final_dataset = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/sanity_check_wellderly_all_filters_withHWE.noZeroAF"+chrom+".vcf.gz")
	final_dataset = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/FINAL_wellderly_all_filters_withHWEandAD.noZeroAF"+chrom+".vcf.gz", "w")
	af_file = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/AF_noAD_filter/final_"+chrom+".txt.gz")

	header = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/header.txt")
	
	for line in header:
		final_dataset.write(line)

	header.close()


	postfilter_dictionary = {}

	postfilter_counter = 0
	#Index the filtered out rare variants
	for line in filtered_out_rare_variants:
		if line[0] == "#":
			continue
		else:
			postfilter_counter += 1
			tp_line = line.strip().split("\t")
			dict_key = tp_line[1]+"_"+tp_line[3]+"_"+tp_line[4]
			postfilter_dictionary[dict_key] = "Y"
	
	filtered_out_rare_variants.close()
	print "Variants to be removed counter " + str(postfilter_counter)
	print "Size of the variants to be removed dictionary " + str(len(postfilter_dictionary))

	#Index the AFs
	AF_dictionary = {}
	for line in af_file:
		af_array = []
		tp_line = line.strip().split("\t")
		dict_key = tp_line[1]+"_"+tp_line[3]+"_"+tp_line[4]
		af_array.append(float(tp_line[5]))
		af_array.append(float(tp_line[6]))
		AF_dictionary[dict_key] = af_array

	af_file.close()
	print "Size of the AF dictionary " + str(len(AF_dictionary))

	#copy to file the variants that were not found in the postfilter dictionary
	prefilter_counter = 0
	removed_variants_counter = 0
	found_variants = 0
	write_buffer = ""
	bdictionary_of_variants_to_be_removed = {}

	well_001 = 0
	well_005 = 0
	well_005_plus = 0

	inova_001 = 0
	inova_005 = 0
	inova_005_plus = 0


	for line in almost_final_dataset:

		prefilter_counter += 1			
		if found_variants%100 == 0:
			final_dataset.write(write_buffer)
			write_buffer = ""
			final_dataset.flush()

		if prefilter_counter%10000 == 0:
			print datetime.datetime.now().time()
			print "total lines " + str(prefilter_counter)
			print "found variants " + str(found_variants)
			print "variants to be filtered out " + str(removed_variants_counter)
			sys.stdout.flush()


		tp_line = line.strip().split("\t")
		dict_key = tp_line[1]+"_"+tp_line[3]+"_"+tp_line[4]
		#check to see if this is one of the variants that we need to filter out
		try:
			val = postfilter_dictionary[dict_key]
			removed_variants_counter += 1
			continue
		except:
			found_variants += 1
			write_buffer = write_buffer + line
			#This needs to have a key, if it doesn't we are fucked
			af_array = AF_dictionary[dict_key]
			well_AF = float(af_array[0])
			inova_AF = float(af_array[1])

			if well_AF > 0.0:
				if well_AF < 0.01:
				    well_001 += 1
				elif well_AF < 0.05:
				    well_005 += 1
				else:
				    well_005_plus += 1

			if inova_AF > 0.0:
				if inova_AF < 0.01:
				    inova_001 += 1
				elif inova_AF < 0.05:
				    inova_005 += 1
				else:
				    inova_005_plus +=1

	final_dataset.write(write_buffer)
	final_dataset.close()
	print "Good variants " + str(found_variants)
	print "Removed variants " + str(removed_variants_counter)

	'''
	counter_file = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_count_filters/final/Counter_by_AF_after_ALL_filters.txt", "a")

	AF_counter = chrom + "\t" + str(found_variants) + "\t" + str(removed_variants_counter) + "\t" + str(well_001) + "\t" + str(well_005) + "\t" + str(well_005_plus) + "\t" + str(inova_001) + "\t" + str(inova_005) + "\t" + str(inova_005_plus) + "\n"
	
	#add AF counter
	counter_file.write(AF_counter)
	counter_file.close()
	'''
	almost_final_dataset.close()


if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)
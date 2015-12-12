"""
Extract the variants that were filtered out byt the allele depth filter

I will need to extract them from the entire main dataset
"""
import os, sys, gzip, datetime


def main(chrom):


	rare_variants_postfilter =gzip.open("/gpfs/group/stsi/data/mrueda/VCF/"+str(chrom)+".0_01.dp.vcf.gz")
	rare_variants_prefilter = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_rare_variants/vcf_rareVariant."+str(chrom)+".vcf.gz")
	filtered_out_rare_variants = gzip.open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_rare_variants/vcf_rare_variant_filtered/"+str(chrom)+".vcf.gz","w")


	postfilter_dictionary = {}

	postfilter_counter = 0
	for line in rare_variants_postfilter:
		if line[0] == "#":
			continue
		else:
			postfilter_counter += 1
			tp_line = line.strip().split("\t")
			dict_key = tp_line[1]+"_"+tp_line[3]+"_"+tp_line[4]
			postfilter_dictionary[dict_key] = "Y"
	
	rare_variants_postfilter.close()
	print "Postfilter counter " + str(postfilter_counter)
	print "Size of the postfilder dictionary " + str(len(postfilter_dictionary))

	#copy to file the variants that were not found in the postfilter dictionary
	prefilter_counter = 0
	removed_variants_counter = 0
	found_variants = 0
	write_buffer = ""
	dictionary_of_variants_to_be_removed = {}
	for line in rare_variants_prefilter:

		if line[0] == "#":
			continue
		else:
			prefilter_counter += 1			
			if removed_variants_counter%100 == 0:
				filtered_out_rare_variants.write(write_buffer)
				write_buffer = ""
				filtered_out_rare_variants.flush()

			if prefilter_counter%10000 == 0:
				print datetime.datetime.now().time()
				print "total lines " + str(prefilter_counter)
				print "found variants " + str(found_variants)
				print "variants to be filtered out " + str(removed_variants_counter)
				sys.stdout.flush()


			tp_line = line.strip().split("\t")
			dict_key = tp_line[1]+"_"+tp_line[3]+"_"+tp_line[4]
			try:
				val = postfilter_dictionary[dict_key]
				found_variants += 1
				continue
			except:
				removed_variants_counter += 1
				write_buffer = write_buffer + line
				dictionary_of_variants_to_be_removed[dict_key] = "Y"

	#copy last lines to file
	filtered_out_rare_variants.write(write_buffer)
	filtered_out_rare_variants.close()
	write_buffer = ""
	rare_variants_prefilter.close()

	print "total lines " + str(prefilter_counter)
	print "found variants " + str(found_variants)
	print "size of dictionary of postfilter variants " + str(len(postfilter_dictionary))

	print "variants to be filtered out " + str(removed_variants_counter)
	print "size of dictinary to be filtered out " + str(len(dictionary_of_variants_to_be_removed))	

	#Get rid of post filter dictionary, we need only dictionary_of_variants_to_be_removed
	postfilter_dictionary = {}



if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)
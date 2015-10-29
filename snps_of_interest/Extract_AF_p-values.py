"""
Extract snps of interest from final association

"""
import os, sys, gzip, datetime


def main():


	in_file="/gpfs/group/stsi/data/nwineing/wellderly/grs/all_weights.csv"
	association="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/final_assoc/test-results.assoc.logistic"
	outfile="/gpfs/group/stsi/data/nwineing/wellderly/grs/all_weights_AF_p-value_v2.csv"

	inf = open(in_file)
	assf = open(association)
	outf = open(outfile, "w")


	snps_dict = {}
	counter = 0
	for i in inf:
		counter += 1
		tp_line = i.split(",")
		dict_key=tp_line[1]+"_"+tp_line[2]
		print "original dict key " + dict_key
		snps_dict[dict_key]=i

	inf.close()

	print "total vars of interest: " + str(counter)
	print "dict dimension " + str(len(snps_dict))

	counter = 0
	found_lines = 0
	for line in assf:
		counter += 1
		if counter == 1:
			continue
		tp_line = line.strip().split("\t")
		#dict_key=tp_line[1].strip()+"_"+tp_line[2].strip()
		dict_key=tp_line[1].strip()+"_"+tp_line[2].strip()
		#print "after dict key" + dict_key

		try:
			var = snps_dict[dict_key]

			#extract the affect allele
			coord = tp_line[0].split("-")

			#final_line = var.strip()+","+coord[3]+","+tp_line[4] +","+tp_line[5]+","+tp_line[3]+"\n"
			final_line = var.strip()+","+"-"+","+"-" +","+tp_line[5]+","+tp_line[3]+"\n"
			print "Variant found! "+ final_line
			found_lines += 1
			outf.write(final_line)
		except:
			continue


	print "Total lines " + str(counter)
	print "var found " + str(found_lines) 

	assf.close()
	outf.close()


if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main()
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)
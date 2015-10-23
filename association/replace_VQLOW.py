"""
Calculate AF snps of interest

"""

import os, sys, gzip, datetime 

def main():


	vcf_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_allFilters_36kmer_snpsOnly_AF0.01/final_vcf_nokmer_snps_AF0.01.noRelated.vcf.gz"
	out_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/vcf_noVQHIGH.vcf.gz"

	inp = gzip.open(vcf_file)
	out = gzip.open(out_file, "w")


	buffe = ""
	counter = 0
	print "Start"
	for line in inp:
		counter += 1
		if counter%100 == 0:

			out.write(buffe)
			buffe = ""
			out.flush()

		if counter%10000 ==0:
			print "Total lines"
			print str(counter)
			sys.stdout.flush()

		line = line.replace("VQLOW",".")
		buffe = buffe + line

		'''
		if line[:1] == "#"
			print "header"
			out.write(line)
		else:
		'''
	out.write(buffe)
	buffe = ""
	out.flush()





	print "Lines found " + str(counter)

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

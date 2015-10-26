"""
Fix vcf file, remove all of the #

"""

import os, sys, gzip, datetime 

def main():


	vcf_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/final_assoc/final_wellderly_inova_AF0.05.vcf.gz"
	out_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/final_assoc/final_wellderly_inova_AF0.05.fixed.vcf.gz"
	inp = gzip.open(vcf_file)
	out = gzip.open(out_file, "w")


	buffe = ""
	counter = 0
	header_line=0
	print "Start"
	for line in inp:
		counter += 1
		if line[0] == "#":
			header_line += 1
			print str(counter)
			continue
		if counter%100 == 0:

			out.write(buffe)
			buffe = ""
			out.flush()

		if counter%10000 ==0:
			print "Total lines"
			print str(counter)
			sys.stdout.flush()

		buffe = buffe + line

	
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

"""
Replace VQLOW with missing

"""

import os, sys, gzip, datetime 

def main(sample):


	vcf_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_allFilters_36kmer_snpsOnly_AF0.01/final_vcf_nokmer_snps_AF0.01.noRelated."+str(sample)+".vcf.gz"
	out_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_final_allFilters_noVQLOW_snpsOnly_AF0.01/final_vcf_noVQLOW_snps_AF0.01.noRelated."+str(sample)+".vcf.gz"

	inp = gzip.open(vcf_file)
	out = gzip.open(out_file, "w")


	buffe = ""
	counter = 0
	counter_with_VQLOW = 0
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
		
		if line[:1] == "#":
			print "header"
			out.write(line)
		

		elif "VQLOW" in line:
			counter_with_VQLOW += 1
			tp_line = line.strip().split("\t")
			for index,i in enumerate(tp_line):
				if "VQLOW" in i:
					#print "VQLOW"
					var = i.split(":")
					geno = var[0]
					geno_split = geno.split("/")
					if var[3]=="VQLOW":
						geno_split[0]='.'
					if var[4]=="VQLOW":
						geno_split[1]='.'
					var[0]="/".join(geno_split)
					#print "new geno " + ":".join(var)
					tp_line[index] = ":".join(var)
			
			line = "\t".join(tp_line) + "\n" 

		buffe = buffe + line

	out.write(buffe)
	buffe = ""
	out.flush()





	print "Lines found " + str(counter)
	print "Counter with VQLOW " + str(counter_with_VQLOW)
	inp.close()
	out.close()
	#nf.close()



if __name__ == '__main__':

    print "Python Version: " + sys.version
    #args = parse_args()
    #import time
    #start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    #end = time.time()
    #print 'This took {0} seconds'.format(end - start)

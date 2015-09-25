"""
Extract one individual from merged data
"""

import os, sys, gzip, datetime 

def main(sample):

	CG_file = "/gpfs/group/stsi/data/projects/wellderly/annai/200Genomes.recalibrated.vcf.gz"
	chr21_cg = "/gpfs/group/stsi/data/projects/wellderly/annai/chr21.ILLUMINA."+sample+".gz"


	o = gzip.open(chr21_cg, "wb")

	cg_data = gzip.open(CG_file)

	buffer_of_lines = ""
	count_column = 0
	total_variants = 0
	total_homo_ref = 0
	total_missing = 0
	valid_variants = 0

	chrom = "chrM"
	for line in cg_data:
		if line[:2] == "##":
			o.write(line)
			print "header"
		elif line[:1] == "#":
			#o.write(line)
			print line
			line = line.strip()
			tp_line = line.split("\t")
			
			for col in tp_line:
				if col == "LP6005831-DNA_B01":
				#if col == "GS000026835-DID":
					print "found"
					#print str(count_column)
					fl = "\t".join(tp_line[:9]) + "\t" + col + "\n"
					o.write(fl)
					break
				count_column = count_column + 1
			
			print "Going though entries, column to extract: " + str(count_column)

		else:
			#print line[:4]
			if line[:5] == "chr21":
				total_variants = total_variants + 1
				if total_variants%100:
					o.write(buffer_of_lines)
					buffer_of_lines = ""
					#o.flush()

				if total_variants%10000 == 0:
					print str(total_variants)
					#sys.stdout.flush()

				line = line.strip()
				tp_line = line.split("\t")

				if 'PASS' in tp_line[6]:
					#total_variants = total_variants + 1
					indiv = tp_line[count_column]

					gen = indiv.split(":")
					if gen[0] == "./.":
						total_missing = total_missing + 1
						continue
					elif gen[0] == "0/0":
						total_homo_ref = total_homo_ref + 1
						continue
					else:
						valid_variants = valid_variants + 1
						buffer_of_lines = buffer_of_lines + "\t".join(tp_line[:9]) + "\t" + indiv + "\n"
						#buffer_of_lines = buffer_of_lines + "\t".join(tp_line[:8]) + "\t" + gen[0] + "\n"
						


	o.write(buffer_of_lines)
	print "total variants " + str(total_variants)
	#print "total homo ref " + str(total_homo_ref)
	#print "total missing " + str(total_missing)
	print "valid varinats " + str(valid_variants)

	cg_data.close()
	o.close()


if __name__ == '__main__':

    print "Python Version: " + sys.version
    #args = parse_args()
    #import time
    #start = time.time()
    sample = 'LP6005831-DNA_B01'
    main(sample)
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    #end = time.time()
    #print 'This took {0} seconds'.format(end - start)

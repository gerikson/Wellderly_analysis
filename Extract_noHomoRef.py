"""
QC combined data from annai
"""

import os, sys, gzip, datetime 

def main():
	combinedfile = "/gpfs/group/stsi/data/projects/wellderly/annai/merged-recode.vcf.gz"
	#CG_file = "/gpfs/group/stsi/data/projects/wellderly/annai/wellderly_601CG_chr21.vcf.gz"
	#illumina_file = "/gpfs/group/stsi/data/projects/wellderly/annai/200Genomes.recalibrated.vcf.gz"
	#chr21_cg = "/gpfs/group/stsi/data/projects/wellderly/annai/chr21.cg.gz"
	CG_file = "/gpfs/group/stsi/data/projects/wellderly/annai/chr21.illumina.gz"
	out = "/gpfs/group/stsi/data/projects/wellderly/annai/chr21.nohomoRef.ILLUMINA.noLow.cg.gz"

	combined = gzip.open(combinedfile)
	#chr21 = gzip.open(chr21_cg, "wb")
	o = gzip.open(out, "wb")

	combined_dict = dict()
	
	print "Creating dictionary combined data..."
	total_merged=0
	total_real_geno=0
	half_allele_count=0
	homo_ref = 0
	hetero_ref = 0
	homo_alt = 0
	other = 0
	for line in combined:
		if line.startswith("#"):
			print "header"
		else:
			total_merged = total_merged + 1
			line = line.strip()
			l = line.split("\t")
			dict_key = l[1] + "_" + l[3] + "_" + l[4]
			combined_dict[dict_key] = l[9]
			if l[9] == ".":
				continue
			else:
				total_real_geno = total_real_geno + 1
				if l[9] == "./1":
					half_allele_count = half_allele_count + 1
				elif l[9] == "0/0":
					homo_ref = homo_ref + 1
				elif l[9] == "0/1" or l[9] == "1/0":
					hetero_ref = hetero_ref + 1
				elif l[9] == "1/1":
					homo_alt = homo_alt + 1
				else:
					other = other + 1

	print "Total merged " + str(total_merged)
	print "Total real geno " + str(total_real_geno)
	print "Half allele count " + str(half_allele_count)
	print "Homo ref " + str(homo_ref)
	print "Hetero ref " + str(hetero_ref)
	print "Homo alt " + str(homo_alt)
	print "Other " + str(other)

	cg_data = gzip.open(CG_file)
	total_variants = 0
	total_homo_ref = 0
	total_missing = 0

	for line in cg_data:
		if line[:2] == "##":
			#chr21.write(line)
			o.write(line)
			print "header"
		elif line[:1] == "#":
			o.write(line)
			print line
			line = line.strip()
			tp_line = line.split("\t")
			'''
			for col in tp_line:
				if col == "GS000026835-DID":
					print "found"
					#print str(count_column)
					fl = "\t".join(tp_line[:8]) + "\t" + col + "\n"
					o.write(fl)
					break
				#count_column = count_column + 1
			'''
			print "Going though entries"

		else:
			total_variants = total_variants + 1
			line = line.strip()
			tp_line = line.split("\t")

			if tp_line[6] == 'PASS':
				#ID of column of interest
				#gen = tp_line[count_column].split(":")
				gen = tp_line[8].split(":")
				if gen[0] == "./.":
					total_missing = total_missing + 1
					continue
				elif gen[0] == "0/0":
					total_homo_ref = total_homo_ref + 1
				else:
					o.write(line + "\n")

				if total_variants%10000 == 0:
					print str(total_variants)


	print "total variants " + str(total_variants)
	print "total homo ref " + str(total_homo_ref)
	print "total missing " + str(total_missing)


	combined.close()
	cg_data.close()
	o.close()


if __name__ == '__main__':

    print "Python Version: " + sys.version
    #args = parse_args()
    #import time
    #start = time.time()
    main()
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    #end = time.time()
    #print 'This took {0} seconds'.format(end - start)

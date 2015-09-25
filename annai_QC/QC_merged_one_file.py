"""
Extract one individual from merged data
"""

import os, sys, gzip, datetime 

def main(sample):

	merged_file = "/gpfs/group/stsi/data/projects/wellderly/annai/new_data/merged-chr22-explicit-half-call-recode.vcf.gz"
	outfile = "/gpfs/group/stsi/data/projects/wellderly/annai/new_data/concordant_merged_illumina_cg_LP6005831-DNA_B01_GS000026835-DID.txt"
	i = gzip.open(merged_file)
	mat = open(outfile, "w")

	count_column = 0
	illumina_col = 0
	cg_col = 0

	total_variants = 0

	#mat.write("Complete Genomics GS000026835-DID \n ")
	mat.write("ILLUMINA LP6005831-DNA_B01 \n ")
	mat.write("\t0/1\t1/1\t./1\tother\t.\n")

	print "Extracting data from CG file..."
	count_column=0
	not_found = 0
	found =0
	not_equal = 0
	total_variants = 0
	total_real_variants_cg = 0
	diff_geno = 0
	weird = 0

	#counters
	homo_homo = 0
	homo_hetero = 0
	homo_homoalt = 0
	homo_half = 0
	homo_other = 0
	homo_missing = 0

	het_homo = 0
	het_hetero = 0
	het_homoalt = 0
	het_half = 0
	het_other = 0
	het_missing = 0

	homohet_homo = 0
	homohet_hetero = 0
	homohet_homoalt = 0
	homohet_half = 0
	homohet_other = 0
	homohet_missing = 0

	half_homo = 0
	half_hetero = 0
	half_homoalt = 0
	half_half = 0
	half_other = 0
	half_missing = 0

	other_homo = 0
	other_hetero = 0
	other_homoalt = 0
	other_half = 0
	other_other = 0
	other_missing = 0

	missing_homo = 0
	missing_hetero = 0
	missing_homoalt = 0
	missing_half = 0
	missing_other = 0
	missing_missing = 0
	same_geno = 0

	missing = 0
	for line in i:
		if line[:2] == "##":
			#o.write(line)
			print "header"
		elif line[:1] == "#":
			#o.write(line)
			print line
			line = line.strip()
			tp_line = line.split("\t")
			
			for col in tp_line:
				if col == "GS000026835-DID":
				#if col == "GS000020592-DID":
					print "found"
					#cg_col = count_column
					illumina_col = count_column
				#elif col == "LP6005830-DNA_A01":
				elif col == "LP6005831-DNA_B01":
					cg_col = count_column
					#illumina_col = count_column

				count_column = count_column + 1
			
			print "Illumina column " + str(illumina_col)
			print "CG column " + str(cg_col)
			print "Going though entries, column to extract: " + str(count_column)

		else:


			if total_variants%10000 == 0:
				print str(total_variants)
				#sys.stdout.flush()


			total_variants = total_variants + 1
			line = line.strip()
			tp_line = line.split("\t")


			if tp_line[illumina_col] == '.' and tp_line[cg_col] == '.':
				missing = missing + 1
			else:
				total_real_variants_cg = total_real_variants_cg + 1

				if tp_line[illumina_col] == tp_line[cg_col]:

					if tp_line[cg_col] == '0/0':
						homo_homo = homo_homo + 1
					elif tp_line[cg_col] == '0/1' or tp_line[cg_col] == '1/0':
						het_hetero = het_hetero + 1
					elif tp_line[cg_col] == '1/1':
						homohet_homoalt = homohet_homoalt + 1
					elif tp_line[cg_col] == './1' or tp_line[cg_col] == '1/.':
						half_half = half_half + 1
					else:
						other_other = other_other + 1
					same_geno = same_geno + 1
				else:
					diff_geno = diff_geno + 1
					if tp_line[cg_col] == "0/0":
						if tp_line[illumina_col] == '0/1' or tp_line[illumina_col] == '1/0':
							homo_hetero = homo_hetero + 1
						elif tp_line[illumina_col] == '1/1':
							homo_homoalt = homo_homoalt + 1
						elif tp_line[illumina_col] == './1' or tp_line[illumina_col] == '1/.':
							homo_half = homo_half + 1
						else:
							homo_other = homo_other + 1
					elif tp_line[cg_col] == "0/1" or tp_line[cg_col] == "1/0" :
						if tp_line[illumina_col] == '0/0':
							het_homo = het_homo + 1
						elif tp_line[illumina_col] == '0/1' or tp_line[illumina_col] == '1/0':
							het_hetero = het_hetero + 1
						elif tp_line[illumina_col] == '1/1':
							het_homoalt = het_homoalt + 1
						elif tp_line[illumina_col] == './1' or tp_line[illumina_col] == '1/.':
							het_half = het_half + 1
						elif tp_line[illumina_col] == ".":
							het_missing = het_missing + 1
						else:
							het_other = het_other + 1
					elif tp_line[cg_col] == "1/1":
						if tp_line[illumina_col] == '0/0':
							homohet_homo = homohet_homo + 1
						if tp_line[illumina_col] == '0/1' or tp_line[illumina_col] == '1/0':
							homohet_hetero = homohet_hetero + 1
						elif tp_line[illumina_col] == '1/1':
							homohet_homoalt = homohet_homoalt + 1
						elif tp_line[illumina_col] == './1' or tp_line[illumina_col] == '1/.':
							homohet_half = homohet_half + 1
						elif tp_line[illumina_col] == ".":
							homo_missing = homo_missing + 1
						else:
							homohet_other = homohet_other + 1
					elif tp_line[cg_col] == "./1" or tp_line[cg_col] == "1/.":
						if tp_line[illumina_col] == '0/0':
							half_homo = homohet_homo + 1
						elif tp_line[illumina_col] == '0/1' or tp_line[illumina_col] == '1/0':
							half_hetero = half_hetero + 1
						elif tp_line[illumina_col] == '1/1':
							half_homoalt = half_homoalt + 1
						elif tp_line[illumina_col] == './1' or tp_line[illumina_col] == '1/.' :
							half_half = half_half + 1
						elif tp_line[illumina_col] == ".":
							half_missing = half_missing + 1
						else:
							half_other = half_other + 1
					elif tp_line[cg_col] == ".":
						if tp_line[illumina_col] == '0/0':
							missing_homo = homohet_homo + 1
						elif tp_line[illumina_col] == '0/1' or tp_line[illumina_col] == '1/0':
							missing_hetero = half_hetero + 1
						elif tp_line[illumina_col] == '1/1':
							missing_homoalt = half_homoalt + 1
						elif tp_line[illumina_col] == './1' or tp_line[illumina_col] == '1/.' :
							missing_half = half_half + 1
						elif tp_line[illumina_col] == ".":
							missing_missing = half_missing + 1
						else:
							missing_other = half_other + 1

			if total_variants%10000 == 0:
				print str(total_variants)


	print "total variants " + str(total_variants)
	print "total_real_variants_cg " + str(total_real_variants_cg)
	print "total found:" + str(found)
	print "total found same geno " + str(same_geno)
	print "total diff_geno " + str(diff_geno)
	print "total_not_found " + str(not_found) 
	print "weird merged = ./1 and CG_data = 0/0 count " + str(weird)
	
	print "homo_homo " + str(homo_homo)
	print "homo hetero " + str(homo_hetero)
	print "homo homoalt " + str(homo_homoalt)
	print "homo half " + str(homo_half)
	print "homo other "+ str(homo_other)
	print "homo missing " + str(homo_missing)

	print "het het " + str(het_homo)
	print "het hetero " + str(het_hetero)
	print "het homoalt " + str(het_homoalt)
	print "het half " + str(het_half)
	print "het other "+ str(het_other)
	print "het missing " + str(het_missing)

	print "homohet homo " + str(homohet_homo)
	print "homohet hetero " + str(homohet_hetero)
	print "homohet homoalt " + str(homohet_homoalt)
	print "homohet half " + str(homohet_half)
	print "homohet other " + str(homohet_other)
	print "homohet missing " + str(homohet_missing)

	print "half homo " + str(half_homo)
	print "half hetero " + str(half_hetero)
	print "half homoalt " + str(half_homoalt)
	print "half half " + str(half_half)
	print "half other " + str(half_other)
	print "half missing " + str(half_missing)

	print "other homo " + str(other_homo)
	print "other hetero " + str(other_hetero)
	print "other homoalt " + str(other_homoalt)
	print "other half " + str(other_half)
	print "other other " + str(other_other)
	print "other missing " + str(other_missing)
	
	concordant = float(same_geno)/float(total_real_variants_cg)
	print "concordant rate " + str(concordant)
	mat.write("0/1\t"+str(het_hetero) + "\t"+str(het_homoalt)+"\t"+str(het_half)+"\t" +str(het_other)+"\t"+str(het_missing)+"\n")
	mat.write("1/1\t"+str(homohet_hetero) + "\t"+str(homohet_homoalt)+"\t"+str(homohet_half)+"\t" +str(homohet_other)+"\t"+str(homohet_missing)+"\n")
	mat.write("1/.\t"+str(half_hetero) + "\t"+str(half_homoalt)+"\t"+str(half_half)+"\t" +str(half_other)+"\t"+str(half_missing)+"\n")
	mat.write("Other\t"+str(other_hetero) + "\t"+str(other_homoalt)+"\t"+str(other_half)+"\t" +str(other_other)+"\t"+str(other_missing)+"\n")
	mat.write(".\t"+str(missing_hetero) + "\t"+str(missing_homoalt)+"\t"+str(missing_half)+"\t" +str(missing_other)+"\t"+str(missing_missing)+"\n")
	mat.write("Concordant rate " + str(concordant))


	print "total variants " + str(total_variants)
	print "total missing " + str(missing)

	
	i.close()

	mat.close()

if __name__ == '__main__':

    print "Python Version: " + sys.version
    #args = parse_args()
    #import time
    #start = time.time()
    sample = 'GS000026835-DID'
    main(sample)
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    #end = time.time()
    #print 'This took {0} seconds'.format(end - start)

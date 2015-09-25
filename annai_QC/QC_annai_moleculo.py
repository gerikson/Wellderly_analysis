"""
QC combined data from annai
"""

import os, sys, gzip, datetime 

def main():
	combinedfile = "/gpfs/group/stsi/data/projects/wellderly/annai/merged-recode.vcf.gz"
	#CG_file = "/gpfs/group/stsi/data/projects/wellderly/annai/wellderly_601CG_chr21.vcf.gz"
	#illumina_file = "/gpfs/group/stsi/data/projects/wellderly/annai/200Genomes.recalibrated.vcf.gz"
	#chr21_cg = "/gpfs/group/stsi/data/projects/wellderly/annai/chr21.cg.gz"
	CG_file = "/gpfs/group/stsi/data/projects/wellderly/annai/chr21.nohomoRef.ILLUMINA.noLow.cg.gz"

	combined = gzip.open(combinedfile)
	#chr21 = gzip.open(chr21_cg, "wb")

	combined_dict = dict()
	
	print "Creating dictionary combined data..."
	total_merged=0
	total_real_geno=0
	half_allele_count=0
	homo_ref = 0
	hetero_ref = 0
	homo_alt = 0
	for line in combined:
		if line.startswith("#"):
			print "header"
		else:
			total_merged = total_merged + 1
			line = line.strip()
			l = line.split("\t")
			#dict_key = l[1] + "_" + l[3] + "_" + l[4]

			if l[9] == ".":
				continue
			else:
				dict_key = l[1] + "_" + l[3]
				combined_dict[dict_key] = l[9]
				total_real_geno = total_real_geno + 1
				if l[9] == "./1":
					half_allele_count = half_allele_count + 1
				elif l[9] == "0/0":
					homo_ref = homo_ref + 1
				elif l[9] == "0/1" or l[9] == "1/0":
					hetero_ref = hetero_ref + 1
				elif l[9] == "1/1":
					homo_alt = homo_alt + 1

	print "Total merged " + str(total_merged)
	print "Total real geno " + str(total_real_geno)
	print "Half allele count " + str(half_allele_count)
	print "Homo ref " + str(homo_ref)
	print "Hetero ref " + str(hetero_ref)
	print "Homo alt " + str(homo_alt)

	cg_data = gzip.open(CG_file)


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
	same_geno = 0

	low_vqslod = 0
	for line in cg_data:
		if line[:2] == "##":
			#chr21.write(line)
			print "header"
		elif line[:1] == "#":
			#chr21.write(line)
			print line
			line = line.strip()
			tp_line = line.split("\t")
			for col in tp_line:
				if col == "GS000026835-DID":
					print "found"
					print str(count_column)
					#fl = "\t".join(tp_line[:5]) + "\t" + col + "\n"
					#chr21.write(fl)
					break
				count_column = count_column + 1
			print "Going though entries"

		else:
			total_variants = total_variants + 1
			line = line.strip()
			tp_line = line.split("\t")

			#ID of column of interest
			#gen = tp_line[count_column].split(":")
			gen = tp_line[8].split(":")
			if gen[0] == "./.":
				continue
			else:
				total_real_variants_cg = total_real_variants_cg + 1

				#write only the info of interest
				#fl = "\t".join(tp_line[:8]) + "\t" + tp_line[count_column] + "\n"
				#chr21.write(fl)
				
				#dict_key = tp_line[1] + "_" + tp_line[3] + "_" + tp_line[4]
				dict_key = tp_line[1] + "_" + tp_line[3]
				if dict_key in combined_dict.keys():
					found = found + 1
					if combined_dict[dict_key] == gen[0]:
						if gen[0] == '0/0':
							homo_homo = homo_homo + 1
						elif gen[0] == '0/1' or gen[0] == '1/0':
							het_hetero = het_hetero + 1
						elif gen[0] == '1/1':
							homohet_homoalt = homohet_homoalt + 1
						elif gen[0] == './1' or gen[0] == '1/.':
							half_half = half_half + 1
						else:
							other_other = other_other + 1
						same_geno = same_geno + 1
					else:
						diff_geno = diff_geno + 1
						if gen[0] == "0/0":
							if combined_dict[dict_key] == '0/1' or combined_dict[dict_key] == '1/0':
								homo_hetero = homo_hetero + 1
							elif combined_dict[dict_key] == '1/1':
								homo_homoalt = homo_homoalt + 1
							elif combined_dict[dict_key] == './1' or combined_dict[dict_key] == '1/.':
								homo_half = homo_half + 1
							else:
								homo_other = homo_other + 1
						elif gen[0] == "0/1" or gen[0] == "1/0" :
							if combined_dict[dict_key] == '0/0':
								het_homo = het_homo + 1
							elif combined_dict[dict_key] == '0/1' or combined_dict[dict_key] == '1/0':
								het_hetero = het_hetero + 1
							elif combined_dict[dict_key] == '1/1':
								het_homoalt = het_homoalt + 1
							elif combined_dict[dict_key] == './1' or combined_dict[dict_key] == '1/.':
								het_half = het_half + 1
							else:
								het_other = het_other + 1
						elif gen[0] == "1/1":
							if combined_dict[dict_key] == '0/0':
								homohet_homo = homohet_homo + 1
							if combined_dict[dict_key] == '0/1' or combined_dict[dict_key] == '1/0':
								homohet_hetero = homohet_hetero + 1
							elif combined_dict[dict_key] == '1/1':
								homohet_homoalt = homohet_homoalt + 1
							elif combined_dict[dict_key] == './1' or combined_dict[dict_key] == '1/.':
								homohet_half = homohet_half + 1
							else:
								homohet_other = homohet_other + 1
						elif gen[0] == "./1" or gen[0] == "1/.":
							if combined_dict[dict_key] == '0/0':
								half_homo = homohet_homo + 1
							elif combined_dict[dict_key] == '0/1' or combined_dict[dict_key] == '1/0':
								half_hetero = half_hetero + 1
							elif combined_dict[dict_key] == '1/1':
								half_homoalt = half_homoalt + 1
							elif combined_dict[dict_key] == './1' or combined_dict[dict_key] == '1/.' :
								half_half = half_half + 1
							else:
								half_other = half_other + 1
				else:
					not_found = not_found + 1
					if gen[0] == "0/0":
						homo_missing = homo_missing + 1
					if gen[0] == '0/1' or gen[0] == "1/0":
						het_missing = het_missing + 1
					elif gen[0] == '1/1':
						homohet_missing = homohet_missing + 1
					elif gen[0] == './1' or gen[0] == "1/.":
						half_missing = half_missing + 1
					else:
						other_missing = other_missing + 1

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


	print "low_vqslod " + str(low_vqslod)

	combined.close()
	cg_data.close()


if __name__ == '__main__':

    print "Python Version: " + sys.version
    #args = parse_args()
    #import time
    #start = time.time()
    main()
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    #end = time.time()
    #print 'This took {0} seconds'.format(end - start)

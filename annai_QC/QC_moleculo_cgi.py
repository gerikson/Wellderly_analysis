"""
QC combined data from annai
"""

import os, sys, gzip, datetime 

def main():

	CG_file = "/gpfs/group/stsi/data/projects/wellderly/annai/chr21.GS000026835-DID.cg.gz"
	illumina_file = "/gpfs/group/stsi/data/projects/wellderly/annai/chr21.ILLUMINA.LP6005831-DNA_B01.gz"
	
	cg_illumina_matrix = "/gpfs/group/stsi/data/projects/wellderly/annai/cg_illumina_matrix_cgFirst_NOvqhigh.txt"

	mat = open(cg_illumina_matrix, "w")
	mat.write("Complete Genomics GS000026835-DID \n ")
	#mat.write("ILLUMINA LP6005831-DNA_B01 \n ")
	mat.write("\t0/1\t1/1\t./1\tother\tmissing\n")
	cg = gzip.open(CG_file)

	illumina = gzip.open(illumina_file)
	combined_dict = dict()
	
	print "Creating dictionary combined data..."
	total_merged=0
	total_real_geno=0
	half_allele_count=0
	homo_ref = 0
	hetero_ref = 0
	homo_alt = 0
	other = 0
	for line in cg:
		if line.startswith("#"):
			print "header"
		else:
			total_merged = total_merged + 1
			line = line.strip()
			l = line.split("\t")
			
			geno = l[9].split(":")
			g = geno[0]
			if g == "./." :
				continue
			elif 'VQHIGH' in geno:
				alt = l[4].split(',')
				a = alt[0]
				#dict_key = l[1] + "_" + l[3] + "_" + a
				dict_key = l[1] + "_" + l[3]
				combined_dict[dict_key] = g
				total_real_geno = total_real_geno + 1
				if g == "./1" or g == "1/.":
					half_allele_count = half_allele_count + 1
				elif g == "0/0":
					homo_ref = homo_ref + 1
				elif g == "0/1" or g == "1/0":
					hetero_ref = hetero_ref + 1
				elif g == "1/1":
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

	#cg_data = gzip.open(CG_file)
	cg_data = gzip.open(illumina_file)

	print "Extracting data from CG file..."
	count_column=0
	not_found = 0
	found =0
	not_equal = 0
	total_variants = 0
	total_real_variants_cg = 0
	diff_geno = 0
	weird = 0
	chrom = "chrM"

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
	for line in illumina:
		if line[:2] == "##":
			#chr21.write(line)
			print "header"
		elif line[:1] == "#":
			print line
			'''
			line = line.strip()
			tp_line = line.split("\t")
			for col in tp_line:
				if col == "LP6005831-DNA_B01":
					print "found"
					print str(count_column)
					fl = "\t".join(tp_line[:5]) + "\t" + col + "\n"
					break
				count_column = count_column + 1
			print "Going though entries"
			'''
		else:
			total_variants = total_variants + 1
			line = line.strip()
			tp_line = line.split("\t")

			gen = tp_line[9].split(":")
			if gen[0] == "./.":
				continue
			else:
				total_real_variants_cg = total_real_variants_cg + 1

				dict_key = tp_line[1] + "_" + tp_line[3]
				if dict_key in combined_dict.keys():
					#print "CGI: " + combined_dict[dict_key]
					#print "Illumina " + gen[0]
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
					#print dict_key
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
	print "total found   " + str(found)
	print "total found same geno  " + str(same_geno)
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

	cg.close()
	cg_data.close()
	mat.write("0/1\t"+str(het_hetero) + "\t"+str(het_homoalt)+"\t"+str(het_half)+"\t" +str(het_other)+"\t"+str(het_missing)+"\n")
	mat.write("1/1\t"+str(homohet_hetero) + "\t"+str(homohet_homoalt)+"\t"+str(homohet_half)+"\t" +str(homohet_other)+"\t"+str(homohet_missing)+"\n")
	mat.write("1/.\t"+str(half_hetero) + "\t"+str(half_homoalt)+"\t"+str(half_half)+"\t" +str(half_other)+"\t"+str(half_missing)+"\n")
	mat.write("Other\t"+str(other_hetero) + "\t"+str(other_homoalt)+"\t"+str(other_half)+"\t" +str(other_other)+"\t"+str(other_missing)+"\n")

	mat.close()
	#chr21.close()

if __name__ == '__main__':

    print "Python Version: " + sys.version
    #args = parse_args()
    #import time
    #start = time.time()
    main()
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    #end = time.time()
    #print 'This took {0} seconds'.format(end - start)

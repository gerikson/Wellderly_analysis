"""
Extract variants that have either cohort a p-value <0.00001

"""
import os, sys, gzip, datetime


def main():

	hwe_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/final_assoc/plink.hwe"
	assoc_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/final_assoc/test-results.assoc.logistic.short"
	
	out_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/final_assoc/final_association.txt"
	wellderly_freq ="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/final_assoc/Wellderly_AF.frq"
	inova_freq = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/final_assoc/Inova_AF.frq"
	
	hwe = open(hwe_file)
	outf = open(out_file, "w")
	snps = open(assoc_file)

	well_freq = open(wellderly_freq)
	in_freq = open(inova_freq)

	#Create a dictionary of the hwe bellow threshold and another
	#dict for the ones over threshold 
	
	remove_dict = {}
	all_dict = {}
	counter =0
	removed_var=0
	for i in hwe:
		counter += 1
		if counter == 1:
			continue
		else:
			tp_l = i.strip().split()
			#print "Lenght line: "+ str(len(tp_l))
			p_val = float(tp_l[8])

			#Keep only the variants that are bellow 1e-4
			if p_val <0.0001:
				d_k = tp_l[1].strip()
				#print d_k
				remove_dict[d_k] = p_val
				removed_var += 1
			else:
				d_k = tp_l[1].strip()
				#add which one is this, annafected (UNAFF) or affected (AFF) or ALL
				fdk = tp_l[1].strip()+"_" +tp_l[2].strip()
				#print d_k
				all_dict[fdk] = p_val

	hwe.close()

	print "total hwe: " + str(counter)
	print "total removed: " + str(removed_var)
	print "dict dimension " + str(len(all_dict))
	
	
	#Dict of the wellderly freq
	#welld_dict = {}
	#all_welld_dict = {}
	counter =0
	removed_var=0
	for i in well_freq:
		counter += 1
		if counter == 1:
			continue
		else:
			tp_l = i.strip().split()
			#print "Lenght line: "+ str(len(tp_l))
			try:
				freq = float(tp_l[4])

				#Get rid of the variants that are bellow 0.05 or >0.95
				if freq <0.05 or freq >0.95 :
					d_k = tp_l[1].strip()
					#print d_k
					#welld_dict[d_k] = freq 
					removed_var += 1
				else:
					d_k = tp_l[1].strip()
					fdk = d_k+"_well"
					all_dict[fdk] = freq
			except:
				print i
				#continue

	well_freq.close()

	print "total well: " + str(counter)
	print "total removed: " + str(removed_var)
	print "dict dimension " + str(len(all_dict))

	#Dict of the inova freq
	#inova_dict = {}
	#all_inova_dict = {}
	counter =0
	removed_var=0
	for il in in_freq:
		counter += 1
		if counter == 1:
			continue
		else:
			tp_l = il.strip().split()
			#print "Lenght line: "+ str(len(tp_l))
			try:
				freq = float(tp_l[4])

				#Get rid of the variants that are bellow 0.05 or >0.95
				if freq <0.05 or freq >0.95 :
					d_k = tp_l[1].strip()
					#print d_k
					#inova_dict[d_k] = freq 
					removed_var += 1
				else:
					d_k = tp_l[1].strip()
					fdk = d_k+"_in"
					all_dict[fdk] = freq
			except:
				print i
				#continue

	in_freq.close()
	
	print "total inova: " + str(counter)
	print "total removed: " + str(removed_var)
	print "dict dimension " + str(len(all_dict))

	var_emited = 0
	total_lines = 0
	for line in snps:
		total_lines += 1
		l = line.strip().split()
		var = l[0].strip()
		#print var
		if total_lines == 1:
			final_line = "SNP\tCHR\tBP\tP\twellAF\tinovaAF\tP_hwe_inova\tP_hwe_well\tP_hwe_all\n"
			outf.write(final_line)
		else:
			try:
				#If the data is in the remove dictionary will just 
				#coninue to the next variant
				new_line = remove_dict[var]
				var_emited += 1
				continue
			except:
				#final_line=""
				#Extract the AF wellderly and the inova if there is no such freq 
				#continue to the next iteration 
				tl = l[0].split("-")
				final_line=l[0] +"\t"+tl[0]+"\t"+tl[1] + "\t" + l[1]
				try:
					fdk = var+"_well"
					freq=all_dict[fdk]

					final_line=final_line + "\t" +str(freq)
				except:
					#print fdk
					continue
				try:
					fdk = var+"_in"
					freq=all_dict[fdk]
					final_line=final_line + "\t" +str(freq)
				except:
					#print fdk
					continue

				
				#Extract the p_values by all separatelly UNAFF, AFF and ALL
				try:
					p_inova=var+"_UNAFF"
					p_val = all_dict[p_inova]

					final_line = final_line+"\t"+str(p_val)
					#outf.write(final_line)
				except:

					final_line = final_line + "\t-"
					#outf.write(final_line)
				try:
					p_inova=var+"_AFF"
					p_val = all_dict[p_inova]
					final_line = final_line + "\t" + str(p_val)
					#outf.write(final_line)
				except:

					final_line = final_line + "\t-"
					#outf.write(final_line)

				try:
					p_inova=var+"_ALL"
					p_val = all_dict[p_inova]
					final_line = final_line + "\t" + str(p_val)+"\n"
					#outf.write(final_line)
				except:

					final_line = final_line + "\t-\n"
					#outf.write(final_line)
				outf.write(final_line)


	print "Total lines " + str(total_lines)
	print "var emited " + str(var_emited) 
	outf.close()
	snps.close()
	hwe.close()



if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main()
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)
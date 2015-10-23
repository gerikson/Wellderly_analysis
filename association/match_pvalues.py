"""
Match p_values

"""
import os, sys, gzip, datetime


def main():

	snp_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/smallest_pvalues.correct.txt"
	output_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/smallest_pvalues.combined.txt"
	input_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/vcf_file_smallest_pvalues_AF.vcf"
	no_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/smallest_pvalues.not_found.txt"
	
	f = open(input_file)
	outf = open(output_file, "w")
	snps = open(snp_file)

	nf = open(no_file, "w")

	data_dict = {}

	for i in f:
		tp_l = i.strip().split("\t")
		l_key =tp_l[0]+"_"+tp_l[1]
		data_dict[l_key] = i

	f.close()

	for line in snps:
		l = line.strip().split("\t")
		tp_l = l[0].split("_")
		d_k =  tp_l[0]+ "_" + tp_l[1]
		final_line=""
		try:
			new_line = data_dict[d_k]
			new_split_line = new_line.split("\t")

			#Verify if the allelese match even though they might be inversed
			
			if (tp_l[2] == new_split_line[3] or tp_l[2] == new_split_line[4]) or (tp_l[3] == new_split_line[3] or tp_l[3] == new_split_line[4]):
				final_line = line.strip() + "\t" + new_line.strip() + "\n"
				outf.write(final_line)
			else:
				print "Alleles don't match"
				nf.write(line)

		except:
			print "Var not found"
			nf.write(line)


	outf.close()
	snps.close()
	nf.close()



if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main()
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)
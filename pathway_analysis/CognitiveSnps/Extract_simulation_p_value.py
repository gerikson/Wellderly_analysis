'''
Extract smallest p-value by gene, this is the version 2 with each snp assigned to one single gene
'''
#import numpy as np

import os, sys, gzip, datetime
print 'start'
print datetime.datetime.now().time()


def main(sample):

	plink_results="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/simulation/plink_results/"
	sim_pvalue="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/cognitive_snps/simulation/"

	array_indexes = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/cognitive_snps/snp_position_array")

	out=sim_pvalue+str(sample)+".sim_snps_p_value.txt"

	outputf = open(out,"w")

	extract_pValue_column="awk '{print $9}' "+plink_results+"results"+str(sample)+".assoc.logistic >" + sim_pvalue+str(sample)+".p-values.txt"
	os.system(extract_pValue_column)

	input_file = sim_pvalue+str(sample)+".p-values.txt"

	p_value_array = []
	with open(input_file) as f:

	    for line in f:
	    	line = line.strip()
	    	p_value_array.append(line)



	print "Size of the p_value array " + str(len(p_value_array))


	for line in array_indexes:
		index = int(line.strip())
		pval = p_value_array[index]
		outputf.write(pval+"\n")


	os.system("rm "+sim_pvalue+str(sample)+".p-values.txt")
	outputf.close()
	array_indexes.close()

if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main(sys.argv[1])
    end = time.time()





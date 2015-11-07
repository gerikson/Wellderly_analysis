'''
Extract smallest p-value by gene, this is the version 2 with each snp assigned to one single gene
'''
#import numpy as np

import os, sys, gzip, datetime
print 'start'
print datetime.datetime.now().time()


def main(sample):

	#plink_results="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/simulation/plink_results/"
	#min_pvalue="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/simulation/min_per_gene/"

	real_plink="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/results.assoc.logistic"
	min_pvalue="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/REAL_DATA_smallest_p_value_per_gene.txt"


	gene_coord_file=open("/gpfs/home/gerikson/wellderly/resources/gene_names_coordinates_plink_final.sorted.txt")
	out=min_pvalue

	outputf = open(out,"w")

	extract_pValue_column="awk '{print $9}' "+real_plink+ " >" + min_pvalue+".temp"
	os.system(extract_pValue_column)

	input_file = min_pvalue+".temp"

	p_value_array = []
	with open(input_file) as f:

	    for line in f:
	    	line = line.strip()
	    	if line != "P":
	    		p_value_array.append(float(line))



	print "Size of the p_value array " + str(len(p_value_array))


	for line in gene_coord_file:
		tp_line = line.strip().split("\t")
		#If only one snp in this gene, extract p_value of that gene
		if int(tp_line[2]) == int(tp_line[3]):
			smallest_p_value_per_gene = p_value_array[int(tp_line[2])]
			final_line = tp_line[0] + "\t" + tp_line[1] + "\t" + str(smallest_p_value_per_gene) + "\n"
		else:
			smallest_p_value_per_gene = min(p_value_array[int(tp_line[2]):int(tp_line[3])])
			final_line = tp_line[0] + "\t" + tp_line[1] + "\t" + str(smallest_p_value_per_gene) + "\n"
		#print final_line
		outputf.write(final_line)

	#os.system("rm "+min_pvalue+str(sample)+".p-values.txt")
	outputf.close()
	gene_coord_file.close()

if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main(sys.argv[1])
    end = time.time()





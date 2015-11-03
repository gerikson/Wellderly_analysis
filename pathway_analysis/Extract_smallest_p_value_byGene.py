'''
Extract smallest p-value by gene
'''
#import numpy as np

import os, sys, gzip, datetime
print 'start'
print datetime.datetime.now().time()


def main(sample):

	plink_results="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/simulation/plink_results/"
	min_pvalue="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/simulation/min_per_gene/"


	#working_dir="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/backup_plink_files/"
	#out=working_dir+"smallest_p_value_per_gene.txt"
	gene_coord_file=open("/gpfs/home/gerikson/wellderly/resources/gene_names_coordinates_plink.sorted.txt")
	out=min_pvalue+str(sample)+".min_p-values_per_gene.txt"

	outputf = open(out,"w")

	extract_pValue_column="awk '{print $9}' "+plink_results+"results"+str(sample)+".assoc.logistic >" + min_pvalue+str(sample)+".p-values.txt"
	#FOR REAL DATA
	#extract_pValue_column="awk '{print $9}' "+working_dir+"plink.results.assoc.logistic >" + working_dir+".p-values.txt"
	os.system(extract_pValue_column)

	'''
	#remove the first line, no header necessary
	remove_head = "tail -n +2"+working_dir+"p_values.txt >"+working_dir+"p_values.noHead.txt"
	os.system(remove_head)
	#Replace new line with tab, easier to load in array
	replace = "tr '\n' '\t' < "+working_dir+"p_values.noHead.txt >"+working_dir+"p_values.line.txt" 
	os.system(replace)
	'''

	#input_file = working_dir+"p_values.txt"
	input_file = min_pvalue+str(sample)+".p-values.txt"

	p_value_array = []
	with open(input_file) as f:
	    #p_value_array = [map(float, line.split("\t")) for line in f]
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

	os.system("rm "+min_pvalue+str(sample)+".p-values.txt")
	outputf.close()
	gene_coord_file.close()

if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)




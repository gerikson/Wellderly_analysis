'''
Remove duplicate values 

'''
#import numpy as np

import os, sys, gzip, datetime, math



def main():


	results = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/final_results/RESULTS_negLog10_noLD.noChrom")
	#final_results = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/final_results/FINAL_FINAL_RESULTS_p_values", "w")
	log_results = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/final_results/FINAL_FINAL_RESULTS_negLog10_exon.noChrom.noLD.rnk", "w")


	real_pVal_dict = {}


	for line in results:
		tp_line = line.strip().split("\t")
		#gene = tp_line[0].split("_")
		#dict_key = gene[1]
		dict_key = tp_line[0]
		#Verify if this gene was already present
		try:
			p_value = real_pVal_dict[dict_key]
			if float(tp_line[1]) < float(p_value):
				print "samller value"
				real_pVal_dict[dict_key] = tp_line[1]
		except:
			real_pVal_dict[dict_key] = tp_line[1]

	results.close()

	print "Size of real P-value dict " + str(len(real_pVal_dict))

	print 'start'
	print datetime.datetime.now().time()




	#Store data from dictionary to file
	for dict_key in real_pVal_dict:
		p_value = real_pVal_dict[dict_key]

		#print p_value
		#final_results.write(dict_key+"\t"+str(p_value)+"\n")
		logp = (-1)*(math.log10(float(p_value)))
		log_results.write(dict_key+"\t"+str(logp)+"\n")

	#final_results.close()
	log_results.close()



	print 'End'
	print datetime.datetime.now().time()








if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main()
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)

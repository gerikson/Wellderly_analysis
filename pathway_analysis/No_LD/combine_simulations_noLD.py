'''
Go trough each simultaion file, per each gene if the p_value of the simulation
is lower (more significant) then the real p-value, the gene gets +1 points

final_p-value = (points + 1)/(simulations +1)

For 10,0000 simulations final p-value is

final_p-value = (points + 1)/(10,001)

'''
#import numpy as np

import os, sys, gzip, datetime, math



def main():
	#real_p_values = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/REAL_DATA_smallest_p_value_per_gene.txt")
	real_p_values = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/REAL_DATA_smallest_p_value_per_gene_exon_noLD.txt")
	path_to_input_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/simulation/min_per_gene/"

	results = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/final_results/RESULTS_p_values_noLD","w")
	log_results = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/final_results/RESULTS_negLog10_noLD","w")


	#Create dictionary of real p_values
	#And a dictionary with genes as the key but epty for now that will
	#be simulation results
	real_pVal_dict = {}
	sim_results_dict = {}

	for line in real_p_values:
		tp_line = line.strip().split("\t")
		dict_key = tp_line[0]+"_"+tp_line[1]
		real_pVal_dict[dict_key] = tp_line[2]
		sim_results_dict[dict_key] = 0
	real_p_values.close()

	print "Size of real P-value dict " + str(len(real_pVal_dict))
	print "Size of simulation dict " + str(len(sim_results_dict))
	
	print 'start'
	print datetime.datetime.now().time()
	for sample in range(1,10001):
		if sample % 100 == 0:
			print str(sample)

		sim_file = open(path_to_input_file+str(sample)+".min_p-values_per_gene.txt")
		for line in sim_file:
			tp_line = line.strip().split("\t")
			try:
				dict_key = tp_line[0]+"_"+tp_line[1]
				#I don't use 'try' to get dict entry 
				#Every single line should have a existing gene in the real p-value file
				#If it doesn't we got problems
				#Verify if the simulation p-value is smaller the the real p_value
				if float(tp_line[2]) <float(real_pVal_dict[dict_key]):
					sim_count = sim_results_dict[dict_key]
					sim_results_dict[dict_key] = sim_count+1
			except:
				print "line " + line
				print "sim_file " + str(sample)+".min_p-values_per_gene.txt"
		sim_file.close()

	#Store data from dictionary to file
	for dict_key in sim_results_dict:
		count_sim = sim_results_dict[dict_key]
		#print count_sim
		p_value = float(count_sim+1)/float(10001)
		#print p_value
		results.write(dict_key+"\t"+str(p_value)+"\n")
		logp = (-1)*(math.log10(p_value))
		log_results.write(dict_key+"\t"+str(logp)+"\n")

	results.close()
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
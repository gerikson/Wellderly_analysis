'''
Start calculating the minimum of each gene for all 10k jobs of the simulation
'''

import os, sys, gzip, datetime


print 'start'
print datetime.datetime.now().time()

counter = 0
QSUB = "qsub -q workq -M gerikson@scripps.edu -l mem=8G -l cput=9600:00:00 -l walltime=500:00:00 "
jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/extract_min_pvalue/"


for sample in range(1,10001):
	if sample%100 == 0:
		print str(sample)


	command = "python /gpfs/home/gerikson/scripts/Wellderly_scripts/GitHub/pathway_analysis/Extract_smallest_p_value_byGene.py " + str(sample)
	jobfile = jobs_folder + str(sample) + ".job"         
	outjob = open(jobfile, 'w')
	outjob.write("#!/bin/csh\n")                    
	outjob.write("#PBS -S /bin/bash\n")
	outjob.write("#PBS -l nodes=1:ppn=1\n")
	outjob.write("#PBS -l mem=8gb\n")
	outjob.write("#PBS -l walltime=500:00:00\n")
	outjob.write("#PBS -l cput=9600:00:00\n")
	outjob.write("#PBS -m n\n")
	outjob.write(command + "\n")
	outjob.close()  
	execute = QSUB + ' -e '+ jobs_folder + str(sample) + '.job.err -o ' + jobs_folder + str(sample)  + '.job.out ' + jobfile
	#print execute 

	sys.stdout.flush()
	clustnum = os.popen(execute, 'r')
	jobnum = clustnum.readline().strip()
	clustnum.close()





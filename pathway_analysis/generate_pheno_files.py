'''
Generate *.pheno files and start plink jobs
'''

import os, sys, gzip, datetime
id_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/id.pheno"
pheno_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/val.pheno"
print 'start'
print datetime.datetime.now().time()

counter = 0
QSUB = "qsub -q stsi -M gerikson@scripps.edu -l mem=8G -l cput=9600:00:00 -l walltime=500:00:00 "
jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/simulation."


for sample in range(1001,10000):
	if sample%100 == 0:
		print str(sample)
	file_to_shuffle="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/val.pheno"
	pheno_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/simulation/pheno_files/temp"+str(sample)+".pheno"
	pheno_shuffled = pheno_file+".shuff"
	shuffle_command = "shuf "+file_to_shuffle+" >"+pheno_shuffled
	#print shuffle_command
	os.system(shuffle_command)
	#combine the id file with the *.pheno file
	combined_shuffled_temp = pheno_file+".shuff.temp"
	combine_command ="paste "+id_file+" "+pheno_shuffled+" | column -s $' ' -t >"+combined_shuffled_temp
	#print combine_command
	os.system(combine_command)
	final_shuffled="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/simulation/pheno_files/plink"+str(sample)+".pheno"
	sed_command="sed 's/\t/ /g' "+combined_shuffled_temp+" >"+final_shuffled
	#print sed_command
	os.system(sed_command)

	os.system("rm "+pheno_shuffled)
	os.system("rm "+combined_shuffled_temp)

	plink_results="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/simulation/plink_results/results"+str(sample)
	working_folder="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/"
	plink_input=working_folder+"plink"
	pca_file=working_folder+"plink-PCA.eigenvec"

	plink_command = "/gpfs/home/nwineing/plink --bfile "+plink_input+" --logistic hide-covar --pheno "+final_shuffled+" --covar "+pca_file+" --allow-no-sex --out "+plink_results
	jobfile = jobs_folder + str(sample) + ".job"         
	outjob = open(jobfile, 'w')
	outjob.write("#!/bin/csh\n")                    
	outjob.write("#PBS -S /bin/bash\n")
	outjob.write("#PBS -l nodes=1:ppn=1\n")
	outjob.write("#PBS -l mem=8gb\n")
	outjob.write("#PBS -l walltime=500:00:00\n")
	outjob.write("#PBS -l cput=9600:00:00\n")
	outjob.write("#PBS -m n\n")
	outjob.write("cd "+working_folder+"\n")
	outjob.write(plink_command + "\n")
	outjob.close()  
	execute = QSUB + ' -e '+ jobs_folder + str(sample) + '.job.err -o ' + jobs_folder + str(sample)  + '.job.out ' + jobfile
	#print execute 

	sys.stdout.flush()
	clustnum = os.popen(execute, 'r')
	jobnum = clustnum.readline().strip()
	clustnum.close()





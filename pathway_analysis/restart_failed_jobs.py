'''
Find the jobs that did not succesfully complete and restart the simulation for those jobs
'''
import os, sys, gzip, datetime
id_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/id.pheno"
pheno_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/val.pheno"
print 'start'
print datetime.datetime.now().time()

counter = 0
QSUB = "qsub -q workq -M gerikson@scripps.edu -l mem=8G -l cput=9600:00:00 -l walltime=500:00:00 "
jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/simulation."

def start_job(sample):
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


#find all files that have expected size:
#find -size 348M >jobs.succesfully.completed

#remove the begin and end of the file for easier sorting
#sed 's/.\/results//g' jobs.succesfully.completed >jobs.succesfully.completed.fixed
#sed 's/.assoc.logistic//g' jobs.succesfully.completed.fixed >jobs.succesfully.completed.fixed.v2
#sort -k 1n,1 jobs.succesfully.completed.fixed.v2 >jobs.succesfully.completed.sorted


completed_jobs=open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/simulation/plink_results/jobs.succesfully.completed.sorted")

counter = 0
line = completed_jobs.readline()

'''
while line != "":
	counter = counter + 1
	#if this job was in the completed job file continue
	if str(counter) != line.strip():
		#restart the job
		print str(counter)
		sample = str(counter)
		start_job(sample)

	else:
		#Go to next line
		line = completed_jobs.readline()	
'''
start_job("10000")	






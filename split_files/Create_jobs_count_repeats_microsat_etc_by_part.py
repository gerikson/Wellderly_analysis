"""
Count variants excluded by the median coverage

"""

import os, sys, gzip, datetime 

def create_job_file(sample):

	cf = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/parts/counter_parts_by_chrom.txt")

	print 'start'
	print datetime.datetime.now().time()

	counter = 0
	QSUB = "qsub -q workq -M gerikson@scripps.edu -l mem=16G -l cput=9600:00:00 -l walltime=500:00:00 "
	
	#Repeat
	#jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/count_filters/final/split_files/repeat."
	#Homopolymer
	jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/count_filters/final/split_files/homop."	
	#Micrositelite
	#jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/count_filters/final/split_files/microsite."
	#SegDup
	#jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/count_filters/final/split_files/SegDup."
	#RepeatMask
	#jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/count_filters/final/split_files/RepeatMask."


	for l in cf:

		tp_line = l.split("\t")
		if sample == tp_line[0]:
			number_of_parts = int(tp_line[1]) + 1
			print "chrom found"	
			for part in range(1,number_of_parts): 

				if part%100 == 0:
					print str(part)

				chrom = sample
				#command = "python /gpfs/home/gerikson/scripts/Wellderly_scripts/GitHub/split_files/Count_repeats_microsat_etc_by_part.py " + str(sample) + " "+ str(part)
				#command = "python /gpfs/home/gerikson/scripts/Wellderly_scripts/GitHub/split_files/Count_segDup_byParts.py " + str(sample) + " "+ str(part)
				#command = "python /gpfs/home/gerikson/scripts/Wellderly_scripts/GitHub/split_files/Count_microsite_byPart.py " + str(sample) + " "+ str(part)
				#command = "python /gpfs/home/gerikson/scripts/Wellderly_scripts/GitHub/split_files/count_repeatMask_byPart.py " + str(sample) + " "+ str(part)
				command = "python /gpfs/home/gerikson/scripts/Wellderly_scripts/GitHub/split_files/count_homopoly_byPart.py " + str(sample) + " "+ str(part)

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
				execute = QSUB + ' -e '+ jobs_folder + str(part)+ "."+str(chrom) + '.job.err -o ' + jobs_folder + str(part) + "."+str(chrom) + '.job.out ' + jobfile
            

				sys.stdout.flush()
				clustnum = os.popen(execute, 'r')
				jobnum = clustnum.readline().strip()
				clustnum.close()
	cf.close()



if __name__ == '__main__':

	print "Python Version: " + sys.version
	
	create_job_file("chr1")
	create_job_file("chr2")
	create_job_file("chr7")
	'''
	#For homopolymer run all of the data, I think it might have been a mistake in the original dataset
	#RepeatMask gets all chroms by part too, no ideea what I did that last time, ma
	#chr4 for homopolymers completed

	for sample in range(1,23):   
		chrom = "chr"+str(sample)
		create_job_file(chrom)
	'''
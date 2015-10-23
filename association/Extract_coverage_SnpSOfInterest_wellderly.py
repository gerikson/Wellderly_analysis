"""
Extract coverage for snp of interest from the coverage files wellderly or inova

"""
import os, sys, gzip, datetime
import subprocess as sp
#from subprocess import Popen, PIPE
import shlex


def main():

	QSUB = "qsub -q workq -M gerikson@scripps.edu -l mem=4G -l cput=9600:00:00 -l walltime=500:00:00 "
	jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/extract_snp/cov_inova."

	#snp_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/desease_snps-corected.txt"
	snp_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/smallest_pvalues.correct.txt"
	coverage_output="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/Cov_snpsOfInterest_welld.txt"
			
	snps = open(snp_file)
	#out_filt = open(output_filename, 'w')
	#out_unfilt = open(unfiltered_output, 'w')
	counter = 0
	for line in snps:
		counter += 1
		print str(counter)
		tp_line = line.strip().split("\t")
		ch = tp_line[0].split("_")
		chrom = ch[0]
		start_position = str(int(ch[1]) - 1)

		var = start_position +"_1"
		file_pattern = "coverages_for_chrm_"+str(chrom)+"_*"
		#command = "zcat " + filtered_filename + " | awk '{if ($2 == "+start_position+") print $0}' >>"+output_filename

		command = "find . -name " + '"'+file_pattern+'" | xargs grep -E '+"'"+var+"' >>" +coverage_output
		print command
		jobfile = jobs_folder + str(var) + ".job"         
		outjob = open(jobfile, 'w')
		outjob.write("#!/bin/bash\n")                    
		outjob.write("#PBS -S /bin/bash\n")
		outjob.write("#PBS -l nodes=1:ppn=1\n")
		outjob.write("#PBS -l mem=8gb\n")
		outjob.write("#PBS -l walltime=500:00:00\n")
		outjob.write("#PBS -l cput=9600:00:00\n")
		outjob.write("#PBS -m n\n")
		outjob.write("cd /gpfs/group/torkamani/bhuvan/wellderly/coverage/CoverageInfo/Whites\n")
		#outjob.write("cd /gpfs/group/stsi/data/projects/wellderly/GenomeComb/inova_coverage_byIndiv/Whites\n")
		outjob.write(command + "\n")
		outjob.close()  
		execute = QSUB + ' -e '+ jobs_folder + str(var) + '.job.err -o ' + jobs_folder + str(var)  + '.job.out ' + jobfile
		sys.stdout.flush()
		clustnum = os.popen(execute, 'r')
		jobnum = clustnum.readline().strip()
		print jobnum
		clustnum.close()


	snps.close()


if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main()
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)
"""
Extract snps of interest from vcf file

"""
import os, sys, gzip, datetime
import subprocess as sp
#from subprocess import Popen, PIPE
import shlex


def run(cmd):
  """Runs the given command locally and returns the output, err and exit_code."""
  if "|" in cmd:    
    cmd_parts = cmd.split('|')
  else:
    cmd_parts = []
    cmd_parts.append(cmd)
  i = 0
  p = {}
  for cmd_part in cmd_parts:
    cmd_part = cmd_part.strip()
    if i == 0:
      p[i]=Popen(shlex.split(cmd_part),stdin=None, stdout=PIPE, stderr=PIPE)
    else:
      p[i]=Popen(shlex.split(cmd_part),stdin=p[i-1].stdout, stdout=PIPE, stderr=PIPE)
    i = i +1
  (output, err) = p[i-1].communicate()
  exit_code = p[0].wait()

  return str(output), str(err), exit_code


def main():

	QSUB = "qsub -q stsi -M gerikson@scripps.edu -l mem=8G -l cput=9600:00:00 -l walltime=500:00:00 "
	jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/extract_snp/alzeimers."

	snp_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/final_ALL_snps-corrected.txt"
	#snp_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/disease_snps_alzeimers_corrected.txt"
	#output_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/filtered_snps_alzeimers.txt"
	#unfiltered_output="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/unfiltered_alzeimers.txt"
	output_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/final_ALL_snps-corrected.vcf"

	snps = open(snp_file)
	#out_filt = open(output_filename, 'w')
	#out_unfilt = open(unfiltered_output, 'w')
	counter = 0
	for line in snps:
		counter += 1
		print str(counter)
		tp_line = line.strip().split("\t")
		chrom = tp_line[1]
		start_position = tp_line[2]
		filtered_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/final_assoc/final_wellderly_inova_AF0.05.fixed.withHead.vcf.gz"
		#filtered_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_allFilters_36kmer_snpsOnly_AF0.01/final_vcf_nokmer_snps_AF0.01.noRelated.chr"+str(chrom)+".vcf.gz"
		#unfiltered_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/wellderly_inova.VQHIGH.0.95white.chr"+str(chrom)+".vcf.gz"
		
		
		command = "zcat " + filtered_filename + " | awk '{if ($2 == "+start_position+") print $0}' >>"+output_filename
		jobfile = jobs_folder + str(start_position) + "filtered.job"         
		outjob = open(jobfile, 'w')
		outjob.write("#!/bin/bash\n")                    
		outjob.write("#PBS -S /bin/bash\n")
		outjob.write("#PBS -l nodes=1:ppn=1\n")
		outjob.write("#PBS -l mem=8gb\n")
		outjob.write("#PBS -l walltime=500:00:00\n")
		outjob.write("#PBS -l cput=9600:00:00\n")
		outjob.write("#PBS -m n\n")
		outjob.write(command + "\n")
		outjob.close()  
		execute = QSUB + ' -e '+ jobs_folder + str(start_position) + 'filtered.job.err -o ' + jobs_folder + str(start_position)  + 'filtered.job.out ' + jobfile
		sys.stdout.flush()
		clustnum = os.popen(execute, 'r')
		jobnum = clustnum.readline().strip()
		print jobnum
		clustnum.close()

		'''
		command2 = "zcat " + unfiltered_filename + " | awk '{if ($2 == "+start_position+") print $0}' >>"+unfiltered_output
		jobfile2 = jobs_folder + str(start_position) + "unfiltered.job"         
		outjob = open(jobfile2, 'w')
		outjob.write("#!/bin/bash\n")                    
		outjob.write("#PBS -S /bin/bash\n")
		outjob.write("#PBS -l nodes=1:ppn=1\n")
		outjob.write("#PBS -l mem=4gb\n")
		outjob.write("#PBS -l walltime=500:00:00\n")
		outjob.write("#PBS -l cput=9600:00:00\n")
		outjob.write("#PBS -m n\n")
		outjob.write(command2 + "\n")
		outjob.close()  
		execute = QSUB + ' -e '+ jobs_folder + str(start_position) + 'unfiltered.job.err -o ' + jobs_folder + str(start_position)  + 'unfiltered.job.out ' + jobfile2
		sys.stdout.flush()
		clustnum = os.popen(execute, 'r')
		jobnum = clustnum.readline().strip()
		print jobnum
		clustnum.close()
		'''
		
		'''
		#check the filtered file first
		command = "zcat "+filtered_filename
		p1 = sp.Popen(shlex.split(command), stdout=sp.PIPE)
		command2 = " awk '{if ($2 == "+start_position+") print $0}'"
		p2 = sp.Popen(shlex.split(command2) , stdin=p1.stdout, stdout=sp.PIPE)
		output, error = p2.communicate()
		out_filt.write(output)

		output = ""
		#check the unfiltered file
		command = "zcat "+unfiltered_filename
		p1 = sp.Popen(shlex.split(command), stdout=sp.PIPE)
		command2 = " awk '{if ($2 == "+start_position+") print $0}'"
		p2 = sp.Popen(shlex.split(command2) , stdin=p1.stdout, stdout=sp.PIPE)
		output, error = p2.communicate()
		out_unfilt.write(output)
		'''
	snps.close()
	#out_filt.close()
	#out_unfilt.close()

if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main()
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)
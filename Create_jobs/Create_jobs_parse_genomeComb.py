import os, sys, datetime
"""
This needs module load python/2.6.5 in .bash_profile
"""
def create_job_file(sample):

    #command = "cg compar2vcf /gpfs/group/stsi/data/projects/wellderly/GenomeComb/inova_results/mcompar."+str(sample)+ ".txt.gz /gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_results/genome_comb_results_chr"+str(sample)+".vcf"
    command = "python /gpfs/home/gerikson/scripts/Wellderly_scripts/GitHub/parse_genome_comb.py chr" + str(sample)
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
    print execute 

    sys.stdout.flush()
    clustnum = os.popen(execute, 'r')
    jobnum = clustnum.readline().strip()
    clustnum.close()


print 'start'
print datetime.datetime.now().time()

counter = 0
QSUB = "qsub -q stsi -M gerikson@scripps.edu -l mem=8G -l cput=9600:00:00 -l walltime=500:00:00 "
jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/parse_all/all_data_to_vcf."

#create_job_file("22")

for sample in range(3,23):    
    create_job_file(sample)

create_job_file("X")
create_job_file("Y")
create_job_file("M")
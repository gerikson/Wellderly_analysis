import os, sys, datetime

def create_job_file(sample):

    chrom = "chr"+str(sample)
    command = "python /gpfs/home/gerikson/scripts/Wellderly_scripts/GitHub/split_files/split_files.py "+chrom
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
    print "done"


print 'start'
print datetime.datetime.now().time()

counter = 0
QSUB = "qsub -q workq -M gerikson@scripps.edu -l mem=16G -l cput=9600:00:00 -l walltime=500:00:00 "
jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/count_filters/final/splitf."



for sample in range(1,23):       
    create_job_file(sample)

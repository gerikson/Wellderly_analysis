import os, sys, datetime

def create_job_file(sample):

    #command = "python /gpfs/home/gerikson/scripts/Wellderly_scripts/GitHub/split_annotation_by_chrom.py " + str(sample)
    command = "python /gpfs/home/gerikson/scripts/Wellderly_scripts/GitHub/split_by_chrom_v1.py " + str(sample)
    #jobfile = jobs_folder +"file3.job"  
    jobfile = jobs_folder +"file2.job"        
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
    #execute = QSUB + ' -e '+ jobs_folder + 'file3.job.err -o ' + jobs_folder + 'file3.job.out ' + jobfile
    execute = QSUB + ' -e '+ jobs_folder + 'file2.job.err -o ' + jobs_folder + 'file2.job.out ' + jobfile
    print execute 

    sys.stdout.flush()
    clustnum = os.popen(execute, 'r')
    jobnum = clustnum.readline().strip()
    clustnum.close()


print 'start'
print datetime.datetime.now().time()

counter = 0
QSUB = "qsub -q stsi -M gerikson@scripps.edu -l mem=8G -l cput=9600:00:00 -l walltime=500:00:00 "
jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/annotation"

#sample = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/batch_WHspocvQya_smallvars_vypdIYlUQa_FINAL.txt.gz"
sample="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/batch_WHspocvQya_smallvars_z2LKFEMXob_FINAL.txt.gz"
#sample="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/annotations/wellderly_inova_annotations_part3.txt.gz"
create_job_file(sample)
import os, sys, datetime

def create_job_file(sample):

    #command = "python /gpfs/home/gerikson/scripts/Wellderly_scripts/GitHub/split_files/count_clustered.py chr" + str(sample)
    #command = "python /gpfs/home/gerikson/scripts/Wellderly_scripts/GitHub/split_files/Count_missing.py chr" + str(sample)
    #command = "python /gpfs/home/gerikson/scripts/Wellderly_scripts/GitHub/split_files/Initial_high_quality_variants_count.py chr" + str(sample)
    #command = "python /gpfs/home/gerikson/scripts/Wellderly_scripts/GitHub/split_files/count_coverage.py chr" + str(sample)
    #command = "python /gpfs/home/gerikson/scripts/Wellderly_scripts/GitHub/split_files/Count_repeats_microsat_etc.py chr" + str(sample)
    #command = "python /gpfs/home/gerikson/scripts/Wellderly_scripts/GitHub/split_files/Count_hardyWeinberg.py chr" + str(sample)
    #command = "python /gpfs/home/gerikson/scripts/Wellderly_scripts/GitHub/split_files/count_36mers.py chr" + str(sample)
    command = "python /gpfs/home/gerikson/scripts/Wellderly_scripts/GitHub/split_files/Count_final_after_allFilters.py chr" + str(sample)
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
QSUB = "qsub -q stsi -M gerikson@scripps.edu -l mem=32G -l cput=9600:00:00 -l walltime=500:00:00 "
#jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/count_filters/final/split_files/cluster."
#jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/count_filters/final/split_files/cat1.AF"
#jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/count_filters/final/split_files/high_qual.count"
#jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/count_filters/final/split_files/missing.count"
#jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/count_filters/final/split_files/Coverage.count"
#jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/count_filters/final/split_files/repeat.count"
#jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/count_filters/final/split_files/homopolymer.count"
#jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/count_filters/final/split_files/microsat.count"
#jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/count_filters/final/split_files/SegDup.count"
#jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/count_filters/final/split_files/RepeatMask.count"
#jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/count_filters/final/split_files/HWE.count"
#jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/count_filters/final/split_files/36mer.count"
jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/count_filters/final/split_files/ALLfilters_noAD.count"


#for sample in range(1,23):   
for sample in range(1,23):   
    create_job_file(sample)
'''
create_job_file(4)
create_job_file(5)
create_job_file(6)
create_job_file(6)
create_job_file(7)
'''
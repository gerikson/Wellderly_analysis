import os, sys, datetime

def create_job_file():
    #Counter file
    cf = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/parts/counter_parts_by_chrom.txt")
    chrom_dict = {}

    for line in cf:

        tp_line = line.split("\t")
        print "chrom " + tp_line[0]
        chrom = tp_line[0]
        number_of_parts = int(tp_line[1]) + 1

        for sample in range(1,number_of_parts): 

            if sample%100 == 0:
                print str(sample)
                
            command = "python /gpfs/home/gerikson/scripts/Wellderly_scripts/GitHub/split_files/Extract_AF_on_the_prefilter_dataset_parts.py "+chrom + " " + str(sample)
            jobfile = jobs_folder + str(sample) + "."+str(chrom)+".job"         
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
            execute = QSUB + ' -e '+ jobs_folder + str(sample)+ "."+str(chrom) + '.job.err -o ' + jobs_folder + str(sample) + "."+str(chrom) + '.job.out ' + jobfile
            
            sys.stdout.flush()
            clustnum = os.popen(execute, 'r')
            jobnum = clustnum.readline().strip()
            clustnum.close()
            


print 'start'
print datetime.datetime.now().time()

counter = 0
QSUB = "qsub -q workq -M gerikson@scripps.edu -l mem=16G -l cput=9600:00:00 -l walltime=500:00:00 "
jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/count_filters/final/AF_by_parts/AF_part"

     
create_job_file()

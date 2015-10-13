import os, sys, datetime
"""
vcf to plink
"""
def create_job_file(sample):

    infile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_AFmore0.01/vcf_snps_AF0.01.chr"+str(sample)+".vcf.gz"
    outfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/chr."+str(sample)
    command = "/gpfs/home/nwineing/plink --vcf "+infile+" --double-id --vcf-half-call m --make-bed --out "+outfile
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
jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/association/vcf_to_plink."

#create_job_file("22")

for sample in range(1,23):    
    create_job_file(sample)

create_job_file("X")
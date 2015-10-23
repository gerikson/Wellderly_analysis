"""
Extract snps of interest from vcf file

"""
import os, sys, gzip, datetime
import subprocess as sp
#from subprocess import Popen, PIPE
import shlex


def main():

    QSUB = "qsub -q workq -M gerikson@scripps.edu -l mem=4G -l cput=9600:00:00 -l walltime=500:00:00 "
    jobs_folder = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/jobfolder/extract_snp/assoc_snps."

    snp_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/smallest_pvalues.plus.txt"
    vcf_file="final_vcf_nokmer_snps_AF0.01.noRelated.chr"

    out_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/smallest_pvalues.correct.vcf"
            
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
        start_position = ch[1]
        #start_position = tp_line[1]

        file_pattern = vcf_file + chrom+".vcf.gz"
        #command = "zcat " + filtered_filename + " | awk '{if ($2 == "+start_position+") print $0}' >>"+output_filename
        #command = "find . -name " + '"'+file_pattern+'" | xargs grep -E '+"'"+var+"' >>" +coverage_output
        #command = "zcat "+file_pattern+" | head -50 | awk '{ if ($2 == "+start_position+") print $0 }' >>"+out_file
        command = "zcat "+ file_pattern + " | grep '" + start_position + "' >>"+out_file
        #command = "find . -name "+file_pattern+' | xargs grep -E '+"'"+start_position+"' >>" +out_file
        print command
        jobfile = jobs_folder + str(start_position) + ".job"         
        outjob = open(jobfile, 'w')
        outjob.write("#!/bin/bash\n")                    
        outjob.write("#PBS -S /bin/bash\n")
        outjob.write("#PBS -l nodes=1:ppn=1\n")
        outjob.write("#PBS -l mem=8gb\n")
        outjob.write("#PBS -l walltime=500:00:00\n")
        outjob.write("#PBS -l cput=9600:00:00\n")
        outjob.write("#PBS -m n\n")
        outjob.write("cd /gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_allFilters_36kmer_snpsOnly_AF0.01/\n")
        outjob.write(command + "\n")
        outjob.close()  
        execute = QSUB + ' -e '+ jobs_folder + str(start_position) + '.job.err -o ' + jobs_folder + str(start_position)  + '.job.out ' + jobfile
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
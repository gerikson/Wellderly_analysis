"""
Extract snps of interest from vcf file

"""
import os, sys, gzip, datetime
import subprocess as sp
#from subprocess import Popen, PIPE
import shlex

'''
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
'''


def main():

	
	snp_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/desease_snps.txt"
	output_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/filtered_snps.txt"
	unfiltered_output="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_snps_of_interest/unfiltered_snps.txt"
			
	snps = open(snp_file)
	out_filt = open(output_filename, 'w')
	out_unfilt = open(unfiltered_output, 'w')
	counter = 0
	for line in snps:
		counter += 1
		print str(counter)
		tp_line = line.strip().split("\t")
		chrom = tp_line[1]
		start_position = tp_line[2]
		filtered_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_clustered_repeats_homopoly_etc_missing_cov/v1_wellderly_inova.VQHIGH.0.95white.nocluster.repeats.etc.missing.cov.chr"+str(chrom)+".vcf.gz"
		unfiltered_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/wellderly_inova.VQHIGH.0.95white.chr"+str(chrom)+".vcf.gz"

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

	snps.close()
	out_filt.close()
	out_unfilt.close()

if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main()
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)
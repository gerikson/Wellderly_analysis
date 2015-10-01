"""
Extract whole genome coverage and exome coverage from CG data
"""

import os
import cPickle

Root = '/gpfs/group/stsi/data/projects/wellderly'

output_wg = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/MedianCovereage_wellderly/total_cov_wg_and_exome.txt"
white_wg = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/MedianCovereage_wellderly/whiteOnly_wg_and_exome.txt"
whitesfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_withCorrectGenotypes/White_0.95.txt"

outp = open(output_wg, 'w')
outp_wh = open(white_wg, 'w')
wh = open(whitesfile)
whites_id = {}

for i in wh:
	wht = i.split()
	for j in wht:
		print "White ID " + j.strip()
		whites_id[j.strip()] = j.strip()


DirecDict = {}
total_white = 0
for batchNum in ['01','02','03','04','05','06','07','08','09','10','11']:
	batchName   = 'batch'+batchNum
	batchFolder = os.path.join(Root,batchName)
	
	print batchName

	for dirs, subdirs, files in os.walk(batchFolder):
		#print files
		for f in files:
			if 'summary' in f:
				#print f
				t_dir = dirs.split("/")
				#print t_dir[7] + "\t" + f
				final_line = t_dir[8] + "\t" + f
				final_file = dirs+"/"+f
				print final_file
				#opfile = open(f)
				opfile = open(final_file)

				for line in opfile:
					if 'Gross mapping yield' in line:
						#print line
						tp = line.split('\t')
						whole_genome_cov = float(tp[2].strip())/3.2
						#print str(whole_genome_cov)
						final_line = final_line + "\t" + str(whole_genome_cov)
					if 'Exome coverage' in line and 'weightSumSequenceCoverage >= 10x' in line:
						print line
						exome_cov = line.split("\t")
						print exome_cov[2]
						final_line = final_line + "\t" + exome_cov[2] + "\n"
				opfile.close()
				try:
					id_white = t_dir[8].strip()
					print id_white
					if whites_id[id_white] == id_white:
						total_white += total_white
						outp.write(final_line)
						outp_wh.write(final_line)
				except:
					outp.write(final_line)


outp.close()
outp_wh.close()
print "Total white " + str(total_white)
print "End!"
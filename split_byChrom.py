import os, sys, gzip, datetime


infile=open("/gpfs/home/gerikson/wellderly/filter_data/RepeatMasker/RepeatMasker_alldata.txt")

outfilename="/gpfs/home/gerikson/wellderly/filter_data/RepeatMasker/byChrom_ALL/RepeatMasker."

chrom="chr1"
name=outfilename+chrom+".txt"
outfile=open(name,"w")

line=infile.readline()
while line.strip() != "":
	
	if "#bin" in line:
		line=infile.readline()
		continue

	templine=line.split("\t")
	if chrom == templine[5]:
		'''
		if templine[11] == 'Simple_repeat' or templine[11] == 'Satellite' or templine[11] == 'Low_complexity':
			outfile.write(line)
		'''
		outfile.write(line)
	else:
		break
		'''
		outfile.close()
		chrom=templine[5]
		name=outfilename+chrom+".txt"
		outfile=open(name,"w")
		if templine[11] == 'Simple_repeat' or templine[11] == 'Satellite' or templine[11] == 'Low_complexity':
			outfile.write(line)
		'''
	line=infile.readline()

infile.close()
outfile.close()		


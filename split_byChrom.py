import os, sys, gzip, datetime


#infile=open("/gpfs/home/gerikson/wellderly/filter_data/RepeatMasker/RepeatMasker_alldata.txt")
infile=open("/gpfs/home/gerikson/wellderly/filter_data/36kmer/unique_36mer.txt")
#outfilename="/gpfs/home/gerikson/wellderly/filter_data/RepeatMasker/byChrom_ALL/RepeatMasker."
outfilename="/gpfs/home/gerikson/wellderly/filter_data/36kmer/byChrom/unique_36mer"
chrom="chr1"
name=outfilename+chrom+".txt"
outfile=open(name,"w")
counter =0
print "herro?"
line=infile.readline()
while line.strip() != "":
	counter += 1
	#print str(counter)
	'''
	if "#bin" in line:
		line=infile.readline()
		continue
	'''
	if counter%100000 == 0:
		print str(counter)

	templine=line.split("\t")
	if chrom == templine[0]:
		'''
		if templine[11] == 'Simple_repeat' or templine[11] == 'Satellite' or templine[11] == 'Low_complexity':
			outfile.write(line)
		'''
		outfile.write(line)
		line=infile.readline()

	else:
		print "New Chrom"
		outfile.close()
		chrom=templine[0]
		name=outfilename+chrom+".txt"
		outfile=open(name,"w")
		outfile.write(line)
		line=infile.readline()
		'''
		outfile.close()
		chrom=templine[5]
		name=outfilename+chrom+".txt"
		outfile=open(name,"w")
		if templine[11] == 'Simple_repeat' or templine[11] == 'Satellite' or templine[11] == 'Low_complexity':
			outfile.write(line)
		'''

infile.close()
outfile.close()		


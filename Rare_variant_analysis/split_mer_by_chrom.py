import os, sys, gzip, datetime



infile=open("/gpfs/home/gerikson/wellderly/filter_data/24mer/unique_24mer.txt")
outfilename="/gpfs/home/gerikson/wellderly/filter_data/34mer/byChrom/unique_24mer"
chrom="chr1"
name=outfilename+chrom+".txt"
outfile=open(name,"w")
counter =0
print "herro?"
line=infile.readline()
while line.strip() != "":
	counter += 1
	if counter%100000 == 0:
		print str(counter)

	templine=line.split("\t")
	if chrom == templine[0]:

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


infile.close()
outfile.close()		


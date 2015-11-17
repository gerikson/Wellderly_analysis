import os, sys, gzip, datetime



infile=open("/gpfs/home/gerikson/wellderly/filter_data/24mer/24mer.unique.txt")
outfilename="/gpfs/home/gerikson/wellderly/filter_data/24mer/byChrom/unique_24mer"
chrom="chr1"
name=outfilename+chrom+".txt"
outfile=open(name,"w")
counter =0
print "herro?"
line=infile.readline()
buffe = ""
while line.strip() != "":
	counter += 1
	if counter%100000 == 0:
		print str(counter)
	if counter%100 == 0:
		outfile.write(buffe)
		buffe = ""
		outfile.flush()

	templine=line.strip().split("\t")
	try:
		if chrom == templine[0]:
			buffe = buffe+line
			line=infile.readline()

		else:
			print "New Chrom"
			outfile.write(buffe)
			buffe = ""
			outfile.flush()
			outfile.close()
			chrom=templine[0]
			name=outfilename+chrom+".txt"
			outfile=open(name,"w")
			outfile.write(line)
			line=infile.readline()
	except:
		continue

outfile.write(buffe)
buffe = ""
outfile.flush()
infile.close()
outfile.close()		


import os, sys, gzip, datetime


infile=open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/exclude_All_repeats/final_vcf_allChrom_snps_AF0.01-0.05_MAF.bim.short")

outfilename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/exclude_All_repeats/byChrom/chr"

chrom="1"
name=outfilename+chrom+".txt"
outfile=open(name,"w")

line=infile.readline()
while line.strip() != "":
	line = line.strip()
	templine=line.split("_")
	if chrom==templine[0]:
		new_l = "chr" + "\t".join(templine) + "\n"
		outfile.write(new_l)
	else:
		outfile.close()
		templine=line.split("_")
		chrom=templine[0]
		name=outfilename+chrom+".txt"
		outfile=open(name,"w")
		new_l = "chr" + "\t".join(templine) + "\n"
		outfile.write(new_l)

	line=infile.readline()

infile.close()
outfile.close()		


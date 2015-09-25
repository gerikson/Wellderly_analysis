"""
Extract Variants in/adjacent to homopolmyer repeat minimum length of 6bp.
date: June 25th, 2015
author: gerikson
"""

import os, sys, gzip, datetime

def open_chrom(j):
	print "New chrom " + str(j)
	infilename="/gpfs/group/stsi/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/chr"+str(j)+".fa"
	infile=open(infilename, 'r')
	outfilename="/gpfs/home/gerikson/wellderly/filter_data/Homopolymer/byChrom/homopolymer.chr"+str(j)+".txt"
	outfile=open(outfilename, 'w')

	counter=0
	start_counter=0

	amino_acid=""
	repeat_lenght=0

	for line in infile:
		if line[0]!='>':
			for i in line.strip():
				#print i
				counter=counter+1
				if i.lower()==amino_acid:
					#If this is part of a bigger repeat
					if repeat_lenght>0:
						repeat_lenght=repeat_lenght+1
						#print str(repeat_lenght)
					#if this is a newrepeat
					else:
						repeat_lenght=2
						#That means we already had 2 aa minus the adjusent 
						start_counter=counter-3
				else:
					amino_acid=i.lower()
					if repeat_lenght>6:
						homopolmyer="chr"+str(j)+"\t"+str(start_counter)+"\t"+str(counter)+"\n"
						outfile.write(homopolmyer)
						repeat_lenght=0
						#print "Polymer found! " + homopolmyer
					else:
						repeat_lenght=0


	infile.close()
	outfile.close()


for j in range(1,22):
	open_chrom(j)

open_chrom("X")
open_chrom("Y")




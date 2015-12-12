import os, sys, gzip

def main(sample):
	path_to_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/wellderly_inova.VQHIGH.0.95white."+str(sample)+".vcf.gz"
	infile = gzip.open(path_to_file)
	chunksize = 10000
	fid = 1
	os.system("mkdir /gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/parts/"+sample)
	outfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/parts/"+sample+"/"+sample+".part"
	ft = outfile+str(fid) +".txt.gz"
	print ft
	f = gzip.open(ft, "w")
	counter_file_by_chrom = open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/parts/counter_parts_by_chrom.txt", "a")
	block = ""
	counter = 0
	for line in infile:
		if line[0] == "#":
			continue
		block = block+line
		counter += 1
		if counter%100 == 0:
			f.write(block)
			block = ""

		if counter%10000 == 0:
			f.close()
			fid = fid + 1
			print str(fid)
			outfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/parts/"+sample+"/"+sample+".part"
			ft = outfile+str(fid) +".txt.gz"
			f = gzip.open(ft, "w")
		
	f.write(block)
	block = ""
	f.close()

	infile.close()
	counter_file_by_chrom.write(sample + "\t" + str(fid) + "\t"+ str(counter) + "\n")
	counter_file_by_chrom.close()
	count_line = sample + "\t" + str(fid) + "\t"+ str(counter) + "\n"
	print count_line
	print "DONE!"

main(sys.argv[1])
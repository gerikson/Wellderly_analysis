"""
Remove the variants that have more then 10% missing genotypes

"""
import os, sys, gzip, datetime


def main(chrom):

	input_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_clustered_repeats_homopoly_etc/wellderly_inova.VQHIGH.0.95white.nocluster.repeats.etc."+str(chrom)+".vcf.gz"
	output_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_clustered_repeats_homopoly_etc_missing/wellderly_inova.VQHIGH.0.95white.nocluster.repeats.etc.missing."+str(chrom)+".vcf.gz"
	header_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_clustered_repeats_homopoly_etc_missing/head.txt"
	counter_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_clustered_repeats_homopoly_etc_missing/counter_wellderly_inova.VQHIGH.0.95white.nocluster.repeats.etc.missing.txt"

	o = gzip.open(output_filename, 'w')
	h = open(header_filename)
	for line in h:
		l = line.strip()
		o.write(l + "\n")
	h.close()

	
	c = open(counter_file, "a")
	print 'Calculating...'

	counter = 0
	good_lines = 0
	block = ""
	f = gzip.open(input_filename)

	for line in f:

		if line[:1] == "#":
			print "header"
			o.write(line)

		else:

			counter += 1

			if counter%100 == 0:
				o.write(block)
				block = ""
				o.flush()

			if counter%10000 == 0:
				print datetime.datetime.now().time()
				print "total lines " + str(counter)
				print "Good lines " + str(good_lines)
				sys.stdout.flush()

			tp_line = line.strip().split()
			
			well_missing = 0
			inova_missing = 0

			for index, gen in enumerate(tp_line[9:]):
				index = index + 9
				#geno = gen.split(":")
				#if "." in geno[0]:
				if gen[0] == '.' or gen[2] == '.':
					if index < 529:
						well_missing += 1
					else:
						inova_missing += 1
			#print "wellderly missing " + str(well_missing)
			#print "inova missing " + str(inova_missing)

			if well_missing > 51 or inova_missing >68:
				continue
			else:
				block = block + line
				good_lines += 1
			

	o.write(block)

	print "Total lines " + str(counter)
	print "Total good lines " + str(good_lines)

	c.write(chrom + "\t" + str(counter) + "\t" + str(good_lines)+ "\n")
	c.close()
	f.close()
	o.close()


if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)
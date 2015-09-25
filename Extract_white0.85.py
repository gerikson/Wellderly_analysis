"""
Extract whites only from the already parsed out VQHIGH variants

Extract variants that have missing/unknow genotypes in more then 10 percent wellderly or inova
Extract variants that cave median covereage <10 percent or >100 
There are variants that don't have any covereage at all, keep those for now

TO DO: deal with the clustered variants, not many really 

"""
import os, sys, gzip, datetime 

def main(chrom):
	infile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_withCorrectGenotypes/"+chrom+".vqhigh.vcf.gz"
	outfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_withCorrectGenotypes/"+chrom+".filtered0.85.vcf.gz"
	whites = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_withCorrectGenotypes/white_0.85.txt"
	
	i = gzip.open(infile)
	o = gzip.open(outfile, 'w')
	w = open(whites)
	whites_id = []

	header = []
	#index the whites
	ln = w.readline()
	ln = ln.strip()
	whites_id = ln.split()
	#whites_id.append(ln)
	

	print "Total whites " + str(len(whites_id))
	print whites_id
	counter = 0
	good_counter = 0
	buffe = ""
	
	for line in i:
		if line[:2] == "##":
			print "header"
			o.write(line)
		elif line[0] == "#":
			o.write(line)
			line = line.strip()
			header = line.split("\t")
			#print header
		else:
			counter += 1
			if counter%10000 ==0:
				print "Total lines"
				print str(counter)
				print "Good lines"
				print str(good_counter)
				sys.stdout.flush()

			if good_counter%100 == 0:
				o.write(buffe)
				buffe=""
				o.flush()

			line = line.strip()
			line = line.split("\t")
			#line_for_cluster = tp_line.split(":")
			#print "cluster data"
			#print line_for_cluster[-2]

			missing_well=0
			missing_inova=0

			median_covereage=0
			total_cov_count = 0

			white_line = []
			for index, gen in enumerate(line[9:]):
				index = index +9
				#Extract whites only
				indiv_name = header[index]
				#print indiv_name
				if indiv_name in whites_id:
					white_line.append(gen)

					tp_gen = gen.split(":")
					if "." in tp_gen[0]:
						if indiv_name[-3] == 'DID':
							missing_well += 1
						else:
							missing_inova += 1

					try:
						if '?' not in tp_gen[5]:
							try:
								median_covereage += int(tp_gen[5])
								total_cov_count += 1
							except:
								print "Bad covereage " + tp_gen[5]
					except:
						continue

			#If we have more then 10% missing:

			if missing_well > 57 or missing_inova >74:
				#print "too many missing"
				continue
			
			#Verify covereage
			elif total_cov_count > 0:
				if (median_covereage/total_cov_count) <10 or (median_covereage/total_cov_count) > 100:
					#print "Weird covereage " + str(median_covereage/total_cov_count)
					continue
			else:
				#print line[:9]
				good_counter += 1
				buffe = buffe + "\t".join(line[:9]) + "\t" + "\t".join(white_line) + "\n"

					

				#print header[index]
				#print header
				#break
			'''
			if good_counter%100 == 0:
				o.write(buffe)
				buffe = ""
			'''


	o.write(buffe)
	w.close()
	i.close()
	o.close()
	print "Total lines"
	print str(counter)
	print "Good lines"
	print str(good_counter)
	print "End"




if __name__ == '__main__':

    print "Python Version: " + sys.version
    #args = parse_args()
    import time
    #start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    #end = time.time()
    #print 'This took {0} seconds'.format(end - start)
    #main()
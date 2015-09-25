def main():
    infile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/plink_data_analysis/wellderly.pheno"
    exclude_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_withCorrectGenotypes/white_0.85.txt"
    outfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/plink_data_analysis/wellderly.pheno.excluded0.85.txt"

    exclude = []
    with open(exclude_file) as e:
    	ln = e.readline()
    	ln = ln.strip()
    	exclude = ln.split()
    	print "length of exclude array" + str(len(exclude))

    not_found_counter = 0
    with open(infile) as i:
    	with open(outfile, 'w') as o:
	    	for line in i:
	    		line = line.strip()
	    		l = line.split()
	    		if l[0] in exclude:
	    			o.write("\t".join(l)+"\n")
	    			#not_found_counter += 0
	    			#o.write(l[0]+"\t"+l[1]+"\n")
	    		else:
	    			not_found_counter += 0

if __name__ == '__main__':

    #print "Python Version: " + sys.version
    main()
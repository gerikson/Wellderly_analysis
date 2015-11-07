""" 
Extract position of snp based on start position
"""

from myvariant.variant import Variant

#infile="/Users/gerikson/Desktop/Sublime/Wellderly_scripts/GitHub/snps_of_interest/disease_snps_alzeimers.txt"
#outfile="/Users/gerikson/Desktop/Sublime/Wellderly_scripts/GitHub/snps_of_interest/disease_snps_alzeimers_corrected.txt"
infile="/Users/gerikson/Desktop/CognitiveSNPs/snpID.txt"
outfile="/Users/gerikson/Desktop/CognitiveSNPs/snpID_positiona.txt"
inf = open(infile)
outf = open(outfile, 'w')
diff_coordinates = 0
counter = 0
#desease = ""
for line in inf:
	#print line
	#p_line = line.split("\t")
	#rsID = tp_line[0]

	rsID = line.strip()
	#print rsID
	try:
		counter += 1
		results = Variant.find_by(q=rsID)
		for r in results:
			#print r._id 
			begin = r.dbsnp['hg19']['start']
			temp_id = r._id
			temp_id1 = temp_id.split(":")
			chrom = temp_id1[0]
			outf.write(rsID + "\t" + chrom+"\t" +str(begin)+"\n")

	except:
		print "Not Found " + rsID
		outf.write(rsID)
		continue

#print desease
print "Total variants found " + str(counter)
print "Total diff coordinates: " + str(diff_coordinates)
inf.close()
outf.close()
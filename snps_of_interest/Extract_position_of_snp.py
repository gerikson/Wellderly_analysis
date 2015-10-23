""" 
Extract position of snp based on start position
"""

from myvariant.variant import Variant

infile="/Users/gerikson/Desktop/Sublime/Wellderly_scripts/GitHub/snps_of_interest/disease_snps_alzeimers.txt"
outfile="/Users/gerikson/Desktop/Sublime/Wellderly_scripts/GitHub/snps_of_interest/disease_snps_alzeimers_corrected.txt"
inf = open(infile)
outf = open(outfile, 'w')
diff_coordinates = 0
counter = 0
#desease = ""
for line in inf:
	#print line
	tp_line = line.split("\t")
	rsID = tp_line[0]
	#print rsID
	try:
		counter += 1
		results = Variant.find_by(q=rsID)
		for r in results:
			#print r._id 
			begin = r.dbsnp['hg19']['start']
			if str(begin) != tp_line[2]:
				tp_line[2] = str(begin)
				corrected_line = "\t".join(tp_line)
				outf.write(corrected_line)
				diff_coordinates +=1
			else:
				corrected_line = "\t".join(tp_line)
				outf.write(corrected_line)
	
	except:
		print "Not Found " + rsID
		corrected_line = "\t".join(tp_line)
		outf.write(corrected_line)
		continue
	
#print desease
print "Total variants found " + str(counter)
print "Total diff coordinates: " + str(diff_coordinates)
inf.close()
outf.close()
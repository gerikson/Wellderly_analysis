""" 
Extract position of snp based on start position
"""

from myvariant.variant import Variant

infile="/Users/gerikson/Desktop/CognitiveSNPs/RegionalPlots/RegionalPlotData_v2.txt"
outfile="/Users/gerikson/Desktop/CognitiveSNPs/RegionalPlots/RegionalPlotData_correctSnps.txt"
inf = open(infile)
outf = open(outfile, 'w')
diff_coordinates = 0
counter = 0

dict_pos = {}
for line in inf:

	var = line.strip().split("\t")
	dict_pos[var[2]] = "Y"

inf.close()

counter = 0
results = Variant.find_by(q="chr6:26491708-26530420")
for r in results:
	counter += 1

	try:
		begin = r.dbsnp['hg19']['start']
		found_snp = dict_pos[str(begin)]
		rsID = r.dbsnp["rsid"]
		print rsid
	except:
		print "shit"
		continue

print "Counter " + str(counter)
inf.close()
outf.close()
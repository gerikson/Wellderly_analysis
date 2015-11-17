'''
Generate an array of the snp position
Total snps in snp dictionary 134
Total snps found 84
'''

import os, sys, gzip, datetime


#association_snps=open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/cognitive_snps/real_data.plink.association")
#outfile=open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/cognitive_snps/snp_position_array")
found_snps="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/cognitive_snps/found_snpID_position"
snp_position_file=open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/cognitive_snps/snpID_position.sorted.txt")
missing_snps=open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/cognitive_snps/missing_snpID", "w")
# Prep
print 'Prepping snp position file...'

snp_dict={}

with open(found_snps) as f:
    for line in f:
        items = line.strip().split('\t')


        snp_dict[items[0]] = "Y"

counter = 0
found = 0
for line in snp_position_file:


    counter +=1
    tp_line = line.strip().split()

    #check if the begin position is one the cognitive snps
    try:
        array_of_data = snp_dict[tp_line[0]]

    except:
        missing_snps.write(tp_line[0]+"\n")
        continue

print "Total snps in snp dictionary " + str(len(snp_dict))
print "Total snps found " + str(found)
snp_position_file.close()
missing_snps.close()



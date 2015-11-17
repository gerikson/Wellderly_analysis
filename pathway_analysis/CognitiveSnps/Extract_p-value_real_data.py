'''
Generate an array of the snp position
Total snps in snp dictionary 134
Total snps found 84
'''

import os, sys, gzip, datetime


association_snps=open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/cognitive_snps/real_data.plink.association")
#outfile=open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/cognitive_snps/snp_position_array", "w")
found_snps=open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/cognitive_snps/found_snpID_position_p-value", "w")
snp_position_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/cognitive_snps/snpID_position.sorted.txt"

# Prep
print 'Prepping snp position file...'

snp_dict={}

with open(snp_position_file) as f:
    for line in f:
        items = line.strip().split('\t')
        array_of_data = []
        array_of_data.append(items[0])
        array_of_data.append(items[1])

        snp_dict[items[2]] = array_of_data

counter = 0
found = 0
for line in association_snps:
    if counter == 0:
        counter += 1
        print "header"
        continue
    else:
        #Add 1 to the counter, will deduct from array index if we actually 
        #find the snp in the cognitive snp list
        counter +=1
        tp_line = line.strip().split()
        chrom=str(tp_line[0].strip())
        begin_poz = tp_line[2]
        #check if the begin position is one the cognitive snps
        try:
            array_of_data = snp_dict[begin_poz]
            rsID = array_of_data[0]
            ch = array_of_data[1]
            #Make sure this position is on the same chromosome
            if chrom == ch:
                found +=1
                #Extracting 1 becouse arrays are 0 based easier for next step
                #outfile.write(str(counter-1)+"\n")
                found_snps.write(rsID+"\t"+ch+"\t"+begin_poz+"\t"+tp_line[8]+"\n")
                
        except:
            continue

print "Total snps in snp dictionary " + str(len(snp_dict))
print "Total snps found " + str(found)
association_snps.close()
#outfile.close()
found_snps.close()



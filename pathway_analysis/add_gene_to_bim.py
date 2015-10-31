'''
Add genes to plink file
'''

import os, sys, gzip, datetime


chrom="chr1"

infile=open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/plink.bim")
outfile=open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/plink.withGene.bim", "w")

gene_file="/gpfs/home/gerikson/wellderly/resources/byChrom/Genes_"+chrom+".sorted.txt"


# Prep
print 'Prepping gene file...'
gene_begin=[]
gene_end=[]
gene_dict={}

with open(gene_file) as f:
    for line in f:
        items = line.split('\t')
        gene_begin.append(int(items[1]))
        gene_end.append(int(items[2]))
        gene_dict[int(items[1])] = items[0]

ch="1"
for line in infile:
    #print ch
    tp_line = line.strip().split("\t")
    new_ch=str(tp_line[0].strip())

    if new_ch == ch:
        begin =int(tp_line[3])
        final_string=tp_line[0]+"\t"+tp_line[1]
        gene_count=0
        for index, i in enumerate(gene_begin):
            if begin >= i:
                if begin <= gene_end[index]:
                    gene_count+=1
                    if gene_count == 1:
                        #print "Gene "+ gene_dict[i]
                        final_string=final_string+"_"+gene_dict[i]
                    else:
                        final_string=final_string+"///"+gene_dict[i]

    else:
        print "new chrom " + tp_line[0]
        ch=tp_line[0]
        gene_file="/gpfs/home/gerikson/wellderly/resources/byChrom/Genes_chr"+ch+".sorted.txt"
        gene_begin=[]
        gene_end=[]
        gene_dict={}
        gene_count=0
        with open(gene_file) as f:
            for line in f:
                items = line.split('\t')
                gene_begin.append(int(items[1]))
                gene_end.append(int(items[2]))
                gene_dict[int(items[1])] = items[0]

        begin =int(tp_line[3])
        final_string=tp_line[0]+"\t"+tp_line[1]
        for index, i in enumerate(gene_begin):
            if begin >= i:
                if begin <= gene_end[index]:
                    gene_count+=1
                    if gene_count == 1:
                        #print "Gene "+ gene_dict[i]
                        final_string=final_string+"_"+gene_dict[i]
                    else:
                        final_string=final_string+"///"+gene_dict[i]


    final_string=final_string+"\t"+tp_line[2]+"\t"+tp_line[3]+"\t"+tp_line[4]+"\t"+tp_line[5]+"\n"
    outfile.write(final_string)


infile.close()
outfile.close()


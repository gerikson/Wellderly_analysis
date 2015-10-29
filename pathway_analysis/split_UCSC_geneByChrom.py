'''
Extracted the kgXref genes and mapped them to the UCSC known genes using the mysql database:
mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -N -A -e 'select kgXref.kgID, kgXref.geneSymbol,knownGene.name,knownGene.chrom,knownGene.txStart,knownGene.txEnd from kgXref, knownGene where knownGene.name=kgXref.kgID' >genes_positions.txt

Split the gene file by Chromosomes, for the genes with multiple transcripts keep the 
minimum begin and maximum end, deduct 100,000 from begin and add 100,000 to end

Results:
Counter 82960
Good chrom counter 78827
Total genes 29256
'''

import os, sys, gzip, datetime

good_chroms =["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                "chr21","chr22","chrX","chrY","chrM"]
chrom="chr1"
infile=open("/gpfs/home/gerikson/wellderly/resources/genes_positions.txt")

outfilename="/gpfs/home/gerikson/wellderly/resources/byChrom/Genes_"+chrom+".txt"

outfile=open(outfilename,"a")
counter =0
good_chrom_counter = 0
gene_dict={}
total_genes = 0
for line in infile:
    counter += 1
    line = line.strip()
    tp_line = line.split("\t")
    if tp_line[3] in good_chroms:
        good_chrom_counter+=1
    else:
        continue
    #If this is a different chromosome store the previous chromosome dictionary
    #Extract 100KB from begin and 100KB to end
    #and open a new file for the new chrom
    if tp_line[3] != chrom:
        print chrom + "\t" + str(len(gene_dict))
        total_genes = total_genes + len(gene_dict)
        for value in gene_dict:
            coord = gene_dict[value]
            tp_coord = coord.split("_")
            begin = int(tp_coord[0])-100000
            end = int(tp_coord[1])+100000
            final_line=value+"\t" + str(begin)+"\t"+str(end)+"\t"+tp_coord[2]+"\t"+tp_coord[3]+"\n"
            outfile.write(final_line)
            
        #Clear the doctionary        
        gene_dict = {}
        outfile.close()
        chrom = tp_line[3]
        outfilename="/gpfs/home/gerikson/wellderly/resources/byChrom/Genes_"+chrom+".txt"
        outfile=open(outfilename,"a")

    #Extract gene name
    gene = tp_line[1]

    try:
        pos = gene_dict[gene]
        #verify if the new coordinates are different then the existing ones
        #If begin is smaller then previous or end is bigger then previous set the
        #new coordinates and update the dictionary
        tp_pos = pos.split("_")
        beg = tp_pos[0]
        end = tp_pos[1]
        if int(tp_line[4]) < int(beg):
            #print "new begin! " +tp_line[4]
            #print "old begin " + beg
            beg = tp_line[4]

        if int(tp_line[5]) > int(end):
            #print "new end! " +tp_line[5]
            #print "old end " + end
            end = tp_line[5]

        pos = beg+"_"+end+"_"+tp_pos[2]+"_"+tp_pos[3]+"///"+tp_line[0]
        gene_dict[gene]=pos
        #print "Gene succesfully updated"
    except:
        pos=tp_line[4]+"_"+tp_line[5] + "_"+tp_line[3]+"_"+tp_line[0]
        gene_dict[gene]=pos

#Store last chromosome
print chrom + "\t" +str(len(gene_dict))
total_genes = total_genes + len(gene_dict)
for value in gene_dict:
    coord = gene_dict[value]
    tp_coord = coord.split("_")
    begin = int(tp_coord[0])-100000
    end = int(tp_coord[1])+100000
    final_line=value+"\t" + str(begin)+"\t"+str(end)+"\t"+tp_coord[2]+"\t"+tp_coord[3]+"\n"
    outfile.write(final_line)

print "Counter "+str(counter)
print "Good chrom counter "+str(good_chrom_counter)
print "Total genes " + str(total_genes)
infile.close()
outfile.close() 

for sample in range(1,23):  
    infile= "/gpfs/home/gerikson/wellderly/resources/byChrom/Genes_chr"+str(sample)+".txt" 
    sortedfile= "/gpfs/home/gerikson/wellderly/resources/byChrom/Genes_chr"+str(sample)+".sorted.txt" 
    os.system("sort -k2,2n "+infile+" > "+ sortedfile)

    


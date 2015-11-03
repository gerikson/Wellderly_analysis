'''
Extract gene name and the first and last index of the bim file snp that has it
3799647 plink.withGene.bim
881747 no_assigned_gene.txt
'''

import os, sys, gzip, datetime


chrom="1"
infile=open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/plink.withGene.bim")
outfilename="/gpfs/home/gerikson/wellderly/resources/gene_names_coordinates_plink.txt"



outfile=open(outfilename,"w")

counter =0

gene_dict={}

for line in infile:
    line = line.strip()
    tp_line = line.split("\t")

    #Verify if same chromosome, if different save gene_dict fo file
    if tp_line[0] != chrom:
        for g in gene_dict:
            coord = gene_dict[g]
            outfile.write(chrom + "\t" + g + "\t" + str(coord[0])+"\t" + str(coord[1]) + "\n")
        chrom = tp_line[0]
        gene_dict={}


    genes = tp_line[1].split("_")

    try:
        gen = genes[1].strip().split("///")


        for g in gen:
            #Verify if this gene is already present
            #If present update the end coordinate, start coordinate doesn't need updating
            try:
                coords = gene_dict[g]
                coords[1] = counter
                gene_dict[g] = coords
            except:
                coords = [counter,counter]
                gene_dict[g] = coords
    except:
        print line

    counter += 1

for g in gene_dict:
    coord = gene_dict[g]
    outfile.write(chrom + "\t" + g + "\t" + str(coord[0])+"\t" + str(coord[1]) + "\n")
chrom = tp_line[0]
gene_dict={}
infile.close()
outfile.close()

'''
Extract gene name and the first and last index of the bim file snp that has it
3799647 plink.withGene.bim
881747 no_assigned_gene.txt

Reassign genes, only 1 gene per snp
'''

import os, sys, gzip, datetime


chrom="1"
infile=open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/plink.withGene.bim")
single_gene_bim=open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/pathway_analysis/plink.withSingleGene.bim", "w")
ch = 'chr'+chrom
gene_file="/gpfs/home/gerikson/wellderly/resources/byChrom/Genes_"+ch+".sorted.txt"
overlapping_gene_file=open("/gpfs/home/gerikson/wellderly/resources/overlapping_genes-exons.txt","a")

outfilename="/gpfs/home/gerikson/wellderly/resources/gene_names_coordinates_plink_final_exons.txt"

gene_coord_dict = {}
with open(gene_file) as f:
    for line in f:
        items = line.strip().split('\t')
        gene_pos=[]
        gene_pos.append(int(items[1]))
        gene_pos.append(int(items[2]))
        gene_coord_dict[items[0]] = gene_pos

outfile=open(outfilename,"w")

counter = -1

gene_dict={}
wtf_gene_dict={}
wtf_gene = 0
total_overlapping_genes = 0

for line in infile:
    
    counter += 1
    final_gene = ""
    final_snp_coord = ""
    line = line.strip()
    tp_line = line.split("\t")

    #Verify if same chromosome, if different save gene_dict fo file
    if tp_line[0] != chrom:
        for g in gene_dict:
            coord = gene_dict[g]
            outfile.write(chrom + "\t" + g + "\t" + str(coord[0])+"\t" + str(coord[1]) + "\n")
        
        #Extract the overlapping genes to file
        for gene in overlapping_genes:
            coordinates = gene_coord_dict[gene]
            begin = int(coordinates[0])+100000
            end = int(coordinates[1])-100000
            overlapping_gene_file.write(chrom + "\t" + gene+"\t"+str(begin)+"\t"+str(end)+"\n")
        #Clear dictionary 
        overlapping_genes = {}

        chrom = tp_line[0]
        #Load the next chromosome coordinates
        ch = 'chr'+chrom
        gene_file="/gpfs/home/gerikson/wellderly/resources/byChrom/Genes_"+ch+".sorted.txt"
        gene_coord_dict = {}
        with open(gene_file) as f:
            for line in f:
                items = line.strip().split('\t')
                gene_pos=[]
                gene_pos.append(int(items[1]))
                gene_pos.append(int(items[2]))
                gene_coord_dict[items[0]] = gene_pos


        gene_dict={}


    genes = tp_line[1].split("_")
    final_snp_coord = genes[0]
    #Was this snp associated with any genes previously
    try:
        gen = genes[1].strip().split("///")
    except:
        single_gene_bim.write(line+"\n")
        #print line
        continue
    #if this snp has only one gene, we don't have to do anything,
    #just update the coordinates
    #Verify if this gene is already present
    #If present update the end coordinate, start coordinate doesn't need updating
    if len(gen) == 1:  
        #print "Only one gene assigned to this snps"
        final_gene = gen[0]  
        try:
            coords = gene_dict[final_gene]
            coords[1] = counter
            gene_dict[final_gene] = coords
        except:
            coords = [counter,counter]
            gene_dict[final_gene] = coords        
    #if there multiple genes that overlap with this
    else:

        final_gene = ""
        min_distance = []
        gene_with_coords = []
        #Extract the snp position
        pos = tp_line[1].split("-")
        snppos = int(pos[1])
        #go through all of the genes and extract the minimum distace to either begin or end


        snp_inside_gene = 0
        count_overlaping_genes = 0
        temp_overlap = ""
                
        overlapping_genes = {}
        for g in gen:
            try:
                coords = gene_coord_dict[g]
            except:
                wtf_gene+=1
                #I previouslu removed some genes that had same name on 2 different areast of genome
                #Those genes were already in the generated plink file
                wtf_gene_dict[g] = wtf_gene
                continue

            begin = int(coords[0]) + 100000
            end = int(coords[1]) - 100000
            #If this snps is actually inside the gene and not adjacent, just pic the gene
            if snppos >= begin and snppos <= end:
                #print "This is inside single gene"

                #Uncoment this if you want to keep al overlaping genes
                #Or exclude them, whatever
                
                final_gene = g
                if temp_overlap != "":
                    overlapping_genes[temp_overlap]="Y"
                    overlapping_genes[g]="Y"
                temp_overlap = g
                
                snp_inside_gene += 1
                count_overlaping_genes += 1



        if snp_inside_gene == 0:
            continue
        #if the snp is inside gene and there is only one gene overlapping it
        elif snp_inside_gene == 1 :
            try:
                coords = gene_dict[final_gene]
                coords[1] = counter
                gene_dict[final_gene] = coords
            except:
                coords = [counter,counter]
                gene_dict[final_gene] = coords
        else:
            #if multiple genes overlap the same snp, pick the one where the snp
            #is closest to the begin
            dist_to_begin_array = []
            gene_array = []
            for gene in overlapping_genes:
                coords = gene_coord_dict[gene]
                begin = int(coords[0]) + 100000
                dist_to_begin = abs(int(snppos) - begin)
                dist_to_begin_array.append(dist_to_begin)
                gene_array.append(gene)

            min_dist_to_begin = min(dist_to_begin_array)

            #Extract the index of the dist to begin and it's gene
            final_gene = ""
            for index, val in enumerate(dist_to_begin_array):
                if val == min_dist_to_begin:
                    final_gene = gene_array[index]

            #Edit the coordinates of the final gene
            try:
                coords = gene_dict[final_gene]
                coords[1] = counter
                gene_dict[final_gene] = coords
            except:
                coords = [counter,counter]
                gene_dict[final_gene] = coords

            if count_overlaping_genes > 1:
                #print "overlapping genes " + str(count_overlaping_genes) 
                total_overlapping_genes += 1
                #print overlapping_genes

    single_gene_bim.write(tp_line[0]+"\t"+final_snp_coord+"_"+final_gene+"\t"+tp_line[2]+"\t"+tp_line[3]+"\t"+tp_line[4]+"\t"+tp_line[5]+"\n")





for g in gene_dict:
    coord = gene_dict[g]
    outfile.write(chrom + "\t" + g + "\t" + str(coord[0])+"\t" + str(coord[1]) + "\n")
for gene in overlapping_genes:
    coordinates = gene_coord_dict[gene]
    begin = int(coordinates[0])+100000
    end = int(coordinates[1])-100000
    overlapping_gene_file.write(chrom + "\t" + gene+"\t"+str(begin)+"\t"+str(end)+"\n")

#xtract the overlapping genes and coordinates to file
print "Wtf genes : "+ str(wtf_gene)
print "Lenght of wtf gene: " + str(len(wtf_gene_dict))
print "Total genes: " + str(counter)
print "Total overlapping genes: " + str(total_overlapping_genes)

chrom = tp_line[0]
gene_dict={}
infile.close()
outfile.close()
single_gene_bim.close()

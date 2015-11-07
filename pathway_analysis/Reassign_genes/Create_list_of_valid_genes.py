'''
Extracted the kgXref genes and mapped them to the UCSC known genes using the mysql database:
mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -N -A -e 'select kgXref.kgID, kgXref.geneSymbol,knownGene.name,knownGene.chrom,knownGene.txStart,knownGene.txEnd from kgXref, knownGene where knownGene.name=kgXref.kgID' >genes_positions.txt

Split the gene file by Chromosomes, for the genes with multiple transcripts keep the 
minimum begin and maximum end, deduct 100,000 from begin and add 100,000 to end

Results:
Counter 82960
Good chrom counter 78827
Total genes 29256

# after emiting the genes with UCSC transcripts that don't even overlap
Counter 82960
Good chrom counter 78827
Total genes 28316

#After emiting removing all of the genes that are not present in the 
5th column of the kgXref column

Unique genes in kgXref file:
zcat kgXref.txt.gz | awk -F "\t" '{print $5}' | sort -u | wc -l
28517

Our results (after eliminating the non overlaping - same gene name)
Counter 82960
Good chrom counter 78827
Total genes 28251

'''

import os, sys, gzip, datetime

def is_this_a_real_kgXref_gene(gene):
    try:
        g = kgXref_real_genes[gene]
        return True
    except:
        return False

#Is this one of the gene areas that don't even overlap
def is_this_not_a_problem_gene(gene):
    try:
        prob_gen = problem_gene[gene]
        return False
    except:
        return True



good_chroms =["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                "chr21","chr22","chrX","chrY","chrM"]
chrom="chr1"
infile=open("/gpfs/home/gerikson/wellderly/resources/genes_positions.txt")
outfilename="/gpfs/home/gerikson/wellderly/resources/byChrom/Genes_"+chrom+".txt"

kgXref_file = gzip.open("/gpfs/home/gerikson/wellderly/resources/kgXref.txt.gz")

kgXref_real_genes = {}
#Extract only the genes that have values in the 5th column
for line in kgXref_file:
    line = line.strip().split("\t")
    real_gene = line[4]
    if real_gene.strip() == "":
        continue
    else:
        kgXref_real_genes[real_gene] = "Y"

kgXref_file.close()
print "Real genes found " + str(len(kgXref_real_genes))


#Insert the problem genes in a dictionary
#These genes should not be kept since the overlap multiple regions
p_gene=open("/gpfs/home/gerikson/wellderly/resources/problem_genes_v1.txt")
problem_gene={}
for line in p_gene:
    line = line.strip()
    problem_gene[line] = "Y"
print "Problem genes found " + str(len(problem_gene))


outfile=open(outfilename,"a")

counter =0
good_chrom_counter = 0
gene_dict={}
total_genes = 0
unreal_genes = 0
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
        #print chrom + "\t" + str(len(gene_dict))
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
    gene = gene.replace(" ","_")
    #If this gene present in the 5th column of kgXref 
    #and is not one of the non overlapping genes?
    if is_this_a_real_kgXref_gene(gene):

        if is_this_not_a_problem_gene(gene):
            #print str(is_this_a_real_kgXref_gene(gene))
            #print str(is_this_not_a_problem_gene(gene))

            try:
                pos = gene_dict[gene]
                #verify if the new coordinates are different then the existing ones
                #If begin is smaller then previous or end is bigger then previous set the
                #new coordinates and update the dictionary
                tp_pos = pos.split("_")
                beg = tp_pos[0]
                end = tp_pos[1]

                '''
                #Verify if the new begin is bigger then the end +100,000
                #That means the 2 transcripts are not overlapping should treat it differently:
                #This looks like 2 separate genes, need to me kept twice

                if int(tp_line[4]) >(int(end)+100000):
                    
                    print gene
                    print pos
                    problem_gene.write(gene+"\n")
                    
                    pos=tp_line[4]+"_"+tp_line[5] + "_"+tp_line[3]+"_"+tp_line[0]
                    gene_dict[gene]=pos
                '''

                #if this overlap ends  before the begin of existing overlap then break
                if int(tp_line[5]) < (int(beg)-100000):
                    #remove that gene, it's non overlaping, FAKE - same names, different locations
                    del gene_dict[gene]
                    #print "too small"
                    continue
                elif int(tp_line[4]) > (int(end)+100000):
                    #remove that gene, it's non overlaping, FAKE - same names, different locations
                    del gene_dict[gene]
                    #print "too big"
                    continue


                #else:   
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

    else:
        unreal_genes += 1
        #print str(is_this_a_real_kgXref_gene(gene))
        #print str(is_this_not_a_problem_gene(gene))

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
print "Genes that are not real " + str(unreal_genes)
infile.close()
outfile.close() 
p_gene.close()

for sample in range(1,23):  
    infile= "/gpfs/home/gerikson/wellderly/resources/byChrom/Genes_chr"+str(sample)+".txt" 
    sortedfile= "/gpfs/home/gerikson/wellderly/resources/byChrom/Genes_chr"+str(sample)+".sorted.txt" 
    os.system("sort -k2,2n "+infile+" > "+ sortedfile)

    


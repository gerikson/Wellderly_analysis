"""
Filter by HWE

"""
import os, sys, gzip, datetime


def extract_AF(tp_line):

    alleles = tp_line[4].split(",")
    
    total_well_alleles = 0
    total_inova_alleles = 0
    total_well_alt = 0
    total_inova_alt = 0
    dict_of_alleles_well={}
    dict_of_alleles_inova={}
    
    well_AF = 0
    inova_AF = 0
    if len(alleles) > 1:
        #print "len alleles > 1"
        #create a dictionary that would count the number of alleles
        for index, i in enumerate(alleles):
            #add one to the index, alleles start at 1
            ind = index + 1
            dict_of_alleles_well[ind] = 0
            dict_of_alleles_inova[ind] = 0


    for index, gen in enumerate(tp_line[9:]):
        index = index + 9
        if len(alleles) == 1:
            #tp_gen = gen[:3]
            #First Allele
            if gen[0] == "0":
                if index < 521:
                    total_well_alleles += 1
                else:
                    total_inova_alleles += 1
            elif gen[0] == "1":
                if index < 521:
                    total_well_alleles += 1
                    total_well_alt +=1
                else:
                    total_inova_alleles += 1
                    total_inova_alt +=1
            #Second Allele
            if gen[2] == "0":
                if index < 521:
                    total_well_alleles += 1
                else:
                    total_inova_alleles += 1
            elif gen[2] == "1":
                if index < 521:
                    total_well_alleles += 1
                    total_well_alt +=1
                else:
                    total_inova_alleles += 1
                    total_inova_alt +=1
        
        else:
            #print "Multi alleles"
            #First Allele
            if gen[0] != ".":
                if gen[0] == "0":
                    if index < 521:
                        total_well_alleles += 1
                    else:
                        total_inova_alleles += 1
                
                else:

                    if index < 521:
                        total_well_alleles += 1
                        #Increment the dictionary value for that allele
                        try:
                            dict_of_alleles_well[int(gen[0])] +=1
                        except:
                            print "gen[0] " + str(gen[0])
                            print "alleles " + tp_line[4]
                    else:
                        total_inova_alleles += 1
                        try:
                            dict_of_alleles_inova[int(gen[0])] +=1
                        except:
                            print "gen[0] " + str(gen[0])
                            print "alleles " + tp_line[4]

            #Second Allele
            if gen[2] != ".":
                if gen[2] == "0":
                    if index < 521:
                        total_well_alleles += 1
                    else:
                        total_inova_alleles += 1
                else:
                    if index < 521:
                        total_well_alleles += 1
                        #Increment the dictionary value for that allele
                        try:
                            dict_of_alleles_well[int(gen[2])] +=1
                        except:
                            print "gen[0] " + str(gen[2])
                            print "alleles " + tp_line[4]
                    else:
                        total_inova_alleles += 1
                        try:
                            dict_of_alleles_inova[int(gen[2])] +=1
                        except:
                            print "gen[0] " + str(gen[2])
                            print "alleles " + tp_line[4]

            inova_max_value = max(dict_of_alleles_inova.values())
            well_max_value = max(dict_of_alleles_well.values())

    try:
        if len(alleles) == 1:
            well_AF = float(total_well_alt)/float(total_well_alleles)
        else:
            # For multi allelic values remove the highest AF from 1
            well_AF = 1 - float(well_max_value)/float(total_well_alleles)

        #print str(AF)
    except:
        well_AF = 0.0

    try:
        if len(alleles) == 1:
            inova_AF = float(total_inova_alt)/float(total_inova_alleles)
        else:
            inova_AF = 1 - float(inova_max_value)/float(total_inova_alleles)
        #print str(AF)
    except:
        inova_AF = 0.0

    return well_AF, inova_AF

def main(chrom):
    ch = chrom[3:]
    print "chrom is: " + ch
    


    workingdir="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/"+str(chrom) + "/"
    #os.system("mkdir "+workingdir)
    input_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/"+str(chrom)+"/wellderly_all_filters_"+str(chrom)+".vcf.gz"
    output_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/wellderly_all_filters_withHWE"+str(chrom)+".vcf.gz"
    counter_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/counter.txt"
    pheno_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/final_assoc/final_vcf_allChrom_snps_AF0.01.pheno"


    c = open(counter_file, "a")

    '''    
    os.system("cd "+workingdir)
    plink_command="/gpfs/home/nwineing/plink --vcf "+input_filename+" --double-id --vcf-half-call m --make-bed --out "+workingdir+chrom
    os.system(plink_command)
    
    make_snpID="awk '{print $1, $4"+'"-"'+"$5, $3, $4, $5, $6}' < "+workingdir+chrom+".bim > "+workingdir+chrom+".temp"
    os.system(make_snpID)
    
    rm_snp = "rm "+workingdir+chrom+".bim"
    os.system(rm_snp)
    os.system("mv "+workingdir+chrom+".temp "+workingdir+chrom+".bim")
    
    run_hwe = "/gpfs/home/nwineing/plink --bfile "+workingdir+chrom+" -hardy --allow-no-sex --pheno "+pheno_file +" --out "+workingdir+"plink"
    os.system(run_hwe)
    print "End hwe extraction"
    '''
    

    hwe_file=open(workingdir+"plink.hwe")
    remove_dict = {}
    counter = 0
    removed_var = 0
    for i in hwe_file:
        counter += 1
        if counter == 1:
            continue
        else:
            tp_l = i.strip().split()
            #print "Lenght line: "+ str(len(tp_l))
            p_val = float(tp_l[8])

            #Keep only the variants that are bellow 1e-4
            if p_val <0.0001:
                d_k = tp_l[1].strip()
                #print d_k
                remove_dict[d_k] = p_val
                removed_var += 1

    hwe_file.close()
    print "Total lines in the hwe file "+ str(counter)
    print "Total lines with hwe < 0.0001 "+ str(removed_var)
    counter = 0
    good_lines = 0
    block = ""
    f = gzip.open(input_filename)
    o = gzip.open(output_filename, 'w')
    counter = 0
    good_lines = 0

    well_001 = 0
    well_005 = 0
    well_005_plus = 0

    inova_001 = 0
    inova_005 = 0
    inova_005_plus = 0

    total_wellderly = 0
    total_inova = 0

    for line in f:

        if line[0] == "#":
            print "header"
            o.write(line)
            continue

        counter += 1

        if counter%100 == 0:
        	o.write(block)
        	#filt.write(filter_block)
        	block = ""
        	#filter_block = ""
        	o.flush()

        if counter%10000 == 0:
        	print datetime.datetime.now().time()
        	print "total lines " + str(counter)
        	print "Good lines " + str(good_lines)
        	sys.stdout.flush()


        tp_line = line.strip().split("\t")
        dict_key = tp_line[1]+"-"+tp_line[3]
        #If this variant is in the hwe datatset
        try:
            val = remove_dict[dict_key]
            print "removed variant"
            continue
        except:
            good_lines += 1
            well_AF, inova_AF = extract_AF(tp_line)
            if well_AF < 0.01:
                well_001 += 1
            elif well_AF < 0.05:
                well_005 += 1
            else:
                well_005_plus += 1

            if inova_AF < 0.01:
                inova_001 += 1
            elif inova_AF < 0.05:
                inova_005 += 1
            else:
                inova_005_plus +=1
            block = block+line

    
    f.close()
    o.write(block)
    block = ""
    o.flush()
    o.close()
    AF_counter = chrom + "\t" + str(counter) + "\t" + str(good_lines) + "\t" + str(well_001) + "\t" + str(well_005) + "\t" + str(well_005_plus) + \
        "\t" + str(inova_001) + "\t" + str(inova_005) + "\t" + str(inova_005_plus) + "\n"

    print AF_counter
    c.write(AF_counter)
    c.close()
    print "End!"
    
if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)
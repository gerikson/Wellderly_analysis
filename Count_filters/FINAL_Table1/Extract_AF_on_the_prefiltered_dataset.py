"""
Extract AF on the prefiltered dataset, 
after removing the related individuals
and after setting VQLOW as missing

"""
import os, sys, gzip, datetime
import cPickle
import resource
import gc
import time

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
        final_index = 0
        for index, i in enumerate(alleles):
            final_index = index
            dict_of_alleles_well[index] = 0
            dict_of_alleles_inova[index] = 0

        #add the last index
        dict_of_alleles_well[final_index+1] = 0
        dict_of_alleles_inova[final_index+1] = 0

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

    if total_well_alleles > 0:
        if len(alleles) == 1:
            well_AF = float(total_well_alt)/float(total_well_alleles)
        else:
            # For multi allelic values remove the highest AF from 1
            well_AF = 1 - (float(well_max_value)/float(total_well_alleles))

        #print str(AF)
    else:
        well_AF = 0.0

    if total_inova_alleles > 0:
        if len(alleles) == 1:
            inova_AF = float(total_inova_alt)/float(total_inova_alleles)
        else:
            inova_AF = 1 - (float(inova_max_value)/float(total_inova_alleles))
        #print str(AF)
    else:
        inova_AF = 0.0

    return well_AF, inova_AF

def main(chrom):
    ch = chrom[3:]
    print "chrom is: " + ch

    input_filename = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/wellderly_inova.VQHIGH.0.95white."+chrom+".vcf.gz"
    output_filename = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_noRelated/wellderly_inova.VQHIGH.0.95white.noRelated"+chrom+".vcf.gz"
    AF_file = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_noRelated_AF/AF_"+chrom+".gz"

    relatedfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_allFilters_36kmer_snpsOnly_AF0.01/eliminate_individuals.txt"
    
    w = open(relatedfile)
    header = []
    #index the whites
    ln = w.readline()
    ln = ln.strip()
    #Create dictionary here instead
    related_id = ln.split("\t")
    w.close()


    block = ""
    AF_block = ""
    noRelated_block = ""
    f = gzip.open(input_filename)
    o = gzip.open(output_filename, 'w')
    af_file = gzip.open(AF_file, 'w')

    counter = 0
    good_lines = 0

    well_001 = 0
    well_005 = 0
    well_005_plus = 0

    inova_001 = 0
    inova_005 = 0
    inova_005_plus = 0

    header = []
    final_header = ""
    for line in f:

        if line[:2] == "##":
            print "header"
            o.write(line)
            continue

        #Extract only the white individuals ID
        elif line[0] == "#":
            line = line.strip()
            header = line.split("\t")
            final_header = "\t".join(header[:9])
            counter_white = 0
            for index, gen in enumerate(header[9:]):
                index = index +9
                #Extract whites only
                if gen in related_id:
                    continue
                else:
                    final_header = final_header + "\t" + gen
                    counter_white += 1
            final_header = final_header + "\n"
            o.write(final_header)
            print "# of norelated " + str(counter_white)
            continue


        counter += 1

        if counter%1000 == 0:
            mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            print "memory usage before:", mem
            o.write(noRelated_block)
            noRelated_block = ""
            o.flush()

            af_file.write(AF_block)
            noRelated_block = ""
            af_file.flush()
            print "memory usage after:", mem

        if counter%10000 == 0:
            print datetime.datetime.now().time()
            print "total lines " + str(counter)
            print "Good lines " + str(good_lines)
            sys.stdout.flush()


        tp_line = line.strip().split("\t")
        #Remove related individuals
        noRelated_line = tp_line[:9]
        for index, gen in enumerate(tp_line[9:]):
            index = index +9
            indiv_name = header[index]
            if indiv_name in related_id:
                continue
            else:
                noRelated_line.append(gen)
        
        #Set VQLOW as missing
        for index,i in enumerate(noRelated_line):
            if index > 8:
                var = i.split(":")
                geno = var[0]
                geno_split = geno.split("/")
                if "VQLOW" in i:
                    #print "VQLOW"
                    if var[3]=="VQLOW":
                        geno_split[0]='.'
                    if var[4]=="VQLOW":
                        geno_split[1]='.'
                    var[0]="/".join(geno_split)
                    #print "new geno " + ":".join(var)
                    noRelated_line[index] = ":".join(var)

        #noRelated_block = noRelated_block + "\t".join(noRelated_line) + "\n"
        well_AF, inova_AF = extract_AF(noRelated_line)

        #print "Wellderly_af " + str(well_AF)
        #print "Inova AF" +str(inova_AF)


        #If there are no alleles, don't keep it
        if well_AF == 0.0 and inova_AF == 0.0:

            continue
        else:

            noRelated_block = noRelated_block + "\t".join(noRelated_line) + "\n"

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

            AF_block = AF_block + "\t".join(tp_line[:5]) + "\t" +str(well_AF) + "\t" +str(inova_AF) + "\n"

            good_lines += 1
            continue



    o.write(noRelated_block)
    noRelated_block = ""
    o.flush()

    af_file.write(AF_block)
    noRelated_block = ""
    af_file.flush()

    c=open("/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_noRelated_AF/counter.txt", "a")
    AF_counter = chrom + "\t" + str(counter) + "\t" + str(good_lines) + "\t" + str(well_001) + "\t" + str(well_005) + "\t" + str(well_005_plus) + "\t" + str(inova_001) + "\t" + str(inova_005) + "\t" + str(inova_005_plus) + "\n"

    print AF_counter
    c.write(AF_counter)
    c.close()

    print "Total lines " + str(counter)
    print "Total good lines " + str(good_lines)

    f.close()
    o.close()
    af_file.close()
    print "DONE! "

if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)
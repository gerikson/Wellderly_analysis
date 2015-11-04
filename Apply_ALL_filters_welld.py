"""
Apply rest of filters to the wellderly dataset
* inova coverage
* 36kmer
* remove related individuals
* replace VQLOW with missing
* remove vars with >10%missing Geno (VQLOW might have added some more variants to missing)
* If there are no alleles due to eliminating VQLOW and related don't keep the variant


"""
import os, sys, gzip, datetime

class Overlap(object):
 
        def __init__(self,arr1Start,arr2Start,arr1End,arr2End):
 
        ## Two sets of start and end sites. Each start and end define the
        ## various islands
 
                self.a1s = arr1Start
                self.a2s = arr2Start
                self.a1e = arr1End
                self.a2e = arr2End


        def main(self):
 
                a1s = self.a1s
                a2s = self.a2s
                a1e = self.a1e
                a2e = self.a2e
 
                overlap_array = self.overlapCounts(a1s,a2s,a1e,a2e)
                return overlap_array
 

       	
       	#Overlaping one begin and end with and array of overlaps (a2s and a2e)
        def overlapCounts(self,a1s,a2s,a1e,a2e):
       
 
                #AllOverlaps = []
                #i = 0
                #while i < len(a1s):
                       
                totalOverlap = 0
                #beg1 = a1s[i]
                #end1 = a1e[i]
                beg1 = int(a1s)
                end1 = int(a1e)


                indx = self.insert_pos(a2s,a2s,beg1)
                if indx == "ERROR":
                    return False

                try:
                        beg2 = a2s[indx]
                        end2 = a2e[indx]
                except IndexError:
                       
                        indx = indx-1
                        beg2 = a2s[indx]
                        end2 = a2e[indx]
               
                overlap = self.overlapNumber(beg1,end1,beg2,end2)
                totalOverlap += overlap

                flagF, indxF = True, indx
                flagB, indxB = True, indx

                while flagF:

                        overlap = False
                        indxF += 1
                        try:
                                beg2 = a2s[indxF]
                                end2 = a2e[indxF]
                        except IndexError:
                                flagF = False
                                break

                        overlap = self.overlapNumber(beg1,end1,beg2,end2)
                        if overlap:
                                return True
                                #totalOverlap += overlap
                        else:
                                flagF = False

                while flagB:

                        overlap = False
                        indxB -= 1
                       
                        if indxB < 0:
                                break
                        try:
                                beg2 = a2s[indxB]
                                end2 = a2e[indxB]
                        except IndexError:
                                flagB = False
                                break
                       
                        overlap = self.overlapNumber(beg1,end1,beg2,end2)
                        if overlap:
                                return True
                                #totalOverlap += overlap
                        else:
                                flagB = False
                       
                       
 
                        #i += 1
 
                #print (AllOverlaps)
                return False   
 
        def insert_pos(self,arrI,arr,value):
 
        ## arrI - to get the final index
        ## arr  - the array to work on, same as arrI
       
            if len(arr) == 1:
            	if value > arr[0]:
            		return arrI.index(arr[0])+1
            	else:
            		return arrI.index(arr[0])
            else:
            	arr_mid = len(arr)/2
            	arrLow  = arr[:arr_mid]
            	arrHigh = arr[arr_mid:]
                #print arrLow
            	#try:
            	if value>arrLow[-1]:
            		out = self.insert_pos(arrI,arrHigh,value)
            	else:
            		out = self.insert_pos(arrI,arrLow,value)

                return out


               
        def overlapPercent(self,s1,e1,s2,e2):
 
        ## The function which calculates the overlap between two isalnds
        ##
        ## seq1  s1---------------------e1
        ## seq2             s2---------------------e2
        ##
        ## seq1  s1----------------------------e1
        ## aeq2            s2--------------e2
        ##
        ## seq1           s1--------------------e1
        ## aeq2    s2----------------e2
        ##
        ## seq1            s1-----------e1
        ## aeq2    s2-----------------------------e2
        ##
                
				if (s1<=s2) and (e1<e2) and (e1>=s2):
					return True
				elif (s1<=s2) and (e1>=e2):
					return True
				elif (s1>s2) and (s1<=e2) and (e1>=e2):
					return True
				elif (s1>s2) and (e1<e2):
					return True
				else:
					return False

 				'''
                seqL1 = float(e1-s1+1)
                seqL2 = float(e2-s2+1)
 
                if (s1<=s2) and (e1<e2) and (e1>=s2):
                        overlap = float(e1-s2+1)
                        ovrPerc = min((overlap/seqL1),(overlap/seqL2))*100.0
                elif (s1<=s2) and (e1>=e2):
                        overlap = seqL2
                        ovrPerc = (overlap/seqL1)*100.0
                elif (s1>s2) and (s1<=e2) and (e1>=e2):
                        overlap = float(e2-s1+1)
                        ovrPerc = min((overlap/seqL1),(overlap/seqL2))*100.0           
                elif (s1>s2) and (e1<e2):
                        overlap = seqL1
                        ovrPerc = (overlap/seqL2)*100.0
                else:
                        ovrPerc = 0.0
 
                return ovrPerc 
                '''
 



        def overlapNumber(self,s1,e1,s2,e2):
 
			seqL1 = float(e1-s1+1)
			seqL2 = float(e2-s2+1)
				

			if (s1<=s2) and (e1<e2) and (e1>=s2):
			        overlap = float(e1-s2+1)
			elif (s1<=s2) and (e1>=e2):
			        overlap = seqL2
			elif (s1>s2) and (s1<=e2) and (e1>=e2):
			        overlap = float(e2-s1+1)
			elif (s1>s2) and (e1<e2):
			        overlap = seqL1
			else:
			        overlap = 0.0
			return overlap   
 
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
    
    kmer_file="/gpfs/home/gerikson/wellderly/filter_data/36kmer/byChrom/unique_36mer" + str(chrom) + ".txt"
    coverage_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/inova_median_coverage/MediansCompiled/medians_chrm_"+str(ch)+".csv"  
    

    input_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly_clustered_repeats_homopoly_etc_missing_cov/v1_wellderly_inova.VQHIGH.0.95white.nocluster.repeats.etc.missing.cov."+str(chrom)+".vcf.gz"
    workingdir="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/"+str(chrom)
    os.system("mkdir "+workingdir)
    output_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/"+str(chrom)+"/wellderly_all_filters_"+str(chrom)+".vcf.gz"
    #output_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_rare_variants/vcf_rareVariant."+str(chrom)+".vcf.gz"
    
    counter_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/counter.txt"
    relatedfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_allFilters_36kmer_snpsOnly_AF0.01/eliminate_individuals.txt"

    w = open(relatedfile)
    header = []
    #index the whites
    ln = w.readline()
    ln = ln.strip()
    #Create dictionary here instead
    related_id = ln.split("\t")
    w.close()

    c = open(counter_file, "a")
    kmer_begin = []
    kmer_end = []



    print 'Prepping kmer...'
    with open(kmer_file) as f:
        for line in f:
            items = line.split('\t')
            kmer_begin.append(int(items[1]))
            kmer_end.append(int(items[2]))
    

    cov_file = open(coverage_file)
    covereage = {}
    print "Start creating dictionary"
    for i in cov_file:
        ln = i.split(",")
        covereage[ln[0]] = ln[1]
    print "Dictionary created"
    cov_file.close()


    print 'Calculating...'

    counter = 0
    good_lines = 0


    block = ""
    f = gzip.open(input_filename)
    o = gzip.open(output_filename, 'w')

    for line in f:

        if line[:2] == "##":
            print "header"
            o.write(line)
            continue

        #Extract only the white individuals ID
        elif line[0] == "#":
            #o.write(line)
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
            print "# of related " + str(counter_white)
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
        
        #Check inova coverage
        pos = tp_line[1] + "_" + str(len(tp_line[3]))

        #Extract the coverage from the dictionary
        cov = "0.0"
        try:
            cov = covereage[pos]
        except:
            print "cov. not found"


        if float(cov) < 10.0 or float(cov) > 110.0:
            continue
        else:
            #print "Inova coverage passed"
            begin = int(tp_line[1])
            end = int(tp_line[1]) + 1
            #Overlapping with the kmers that have value of 1
            #If overlap found, good! Unique kmer
            repOverlap = Overlap(begin, kmer_begin, end, kmer_end)
            overlapRep = repOverlap.main()
            if overlapRep:
                #print "kmer passed"
                #Remove related individuals
                noRelated_line = tp_line[:9]
                for index, gen in enumerate(tp_line[9:]):
                    index = index +9
                    indiv_name = header[index]
                    if indiv_name in related_id:
                        continue
                    else:
                        #count_id += 1
                        noRelated_line.append(gen)

                #print "Removed related" 

                #Replace VQLOW and count missing
                well_missing = 0
                inova_missing = 0
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
                   
                        if geno_split[0] == '.' or geno_split[1] == '.':
                            if index < 511:
                                well_missing += 1
                            else:
                                inova_missing += 1
                #print "Removed VQLOW and counter missing"


                if well_missing > 51 or inova_missing >68:
                    continue
                else:
                    block = block + "\t".join(noRelated_line) + "\n"
                    good_lines += 1
                    continue
                    '''
                    well_AF, inova_AF = extract_AF(noRelated_line)
                    #print "Well AF " + str(well_AF)
                    #print "Inova AF " + str(inova_AF)

                    #If there are no alleles, don't keep it
                    if well_AF == 0.0 and inova_AF == 0.0:
                        #print noRelated_line
                        continue
                    elif well_AF <= 0.05 and inova_AF <= 0.05:
                        #print "Good line"
                        # Extract only AF <0.05 in either cohort
                        block = block + "\t".join(noRelated_line) + "\n"
                        good_lines += 1
                        continue
                    '''
    		

    #filt.write(filter_block)
    o.write(block)

    #print "Total g counter " + str(repeatMask_counter)

    print "Total lines " + str(counter)
    print "Total good lines " + str(good_lines)
    c.write(chrom+"\t"+str(counter)+"\t"+str(good_lines)+"\n")

    c.close()
    f.close()
    o.close()


if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)
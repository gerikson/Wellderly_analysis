"""
Counts variants removed by 36mers
If variant overlaps with a 36mers of value 1, continue, else filter out

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
 

def main(chrom):

	ch = chrom[3:]
	print "chrom is: " + ch
	input_filename="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_filtered_VQHIGH_whiteOnly/wellderly_inova.VQHIGH.0.95white."+str(chrom)+".vcf.gz"
	kmer_file="/gpfs/home/gerikson/wellderly/filter_data/36kmer/byChrom/unique_36mer" + str(chrom) + ".txt"

	counter_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_count_filters/filter_36mers_wellderly.txt"
	filter_file="/gpfs/group/stsi/data/projects/wellderly/GenomeComb/final_filter_file/filter_wellderly_36mers."+str(chrom)+".txt.gz"
	
	inputf = gzip.open(input_filename)
	c = open(counter_file, "a")
	filt = gzip.open(filter_file, "w")
	repeatMask_begin = []
	repeatMask_end = []

	count_36mer = 0
	print 'Prepping 36mers...'
	with open(kmer_file) as f:
		for line in f:
			count_36mer += 1
			items = line.split('\t')
			repeatMask_begin.append(int(items[1]))
			repeatMask_end.append(int(items[2]))
    
	print "Total repeats " + str(count_36mer)
	print 'Calculating...'

	counter = 0
	good_lines = 0
	missing_counter = 0

	well_001 = 0
	well_005 = 0
	well_005_plus = 0

	inova_001 = 0
	inova_005 = 0
	inova_005_plus = 0

	total_wellderly = 0
	total_inova = 0
	filt.write("Chrom\tbegin\tRef\tAlt\t36mer\n")
	filter_block = ""


	for line in inputf:

		if line[:2] == "##":
			print "header"
		elif line[:1] == "#":
			tp_line = line.split("\t")
			for i in tp_line:
				if i.endswith("DID"):
					total_wellderly += 1
				elif i.endswith("ASM"):
					total_inova += 1
			print "Total Wellderly " + str(total_wellderly)
			print "Total inova " + str(total_inova)
		else:

			counter += 1
			if counter%100 == 0:
				filt.write(filter_block)
				filter_block = ""
				filt.flush()

			if counter%10000 == 0:
				print datetime.datetime.now().time()
				print "total lines " + str(counter)
				print "Good lines " + str(good_lines)
				sys.stdout.flush()


			tp_line = line.strip().split("\t")
			filter_block = filter_block + tp_line[0]+"\t"+tp_line[1]+"\t"+tp_line[3]+"\t"+tp_line[4]
			begin = int(tp_line[1])
			end = int(tp_line[1]) + 1

			#If variant overlaps with a 36mers of value 1, continue, else filter out
			repOverlap = Overlap(begin, repeatMask_begin, end, repeatMask_end)
			overlapRep = repOverlap.main()
			if overlapRep:
				filter_block = filter_block + "\t\n" 
				good_lines += 1
				continue
			else:
				filter_block = filter_block + "\tYES\n" 

				missing_counter += 1

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
							if index < 529:
								total_well_alleles += 1
							else:
								total_inova_alleles += 1
						elif gen[0] == "1":
							if index < 529:
								total_well_alleles += 1
								total_well_alt +=1
							else:
								total_inova_alleles += 1
								total_inova_alt +=1
						#Second Allele
						if gen[2] == "0":
							if index < 529:
								total_well_alleles += 1
							else:
								total_inova_alleles += 1
						elif gen[2] == "1":
							if index < 529:
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
								if index < 529:
									total_well_alleles += 1
								else:
									total_inova_alleles += 1
							
							else:

								if index < 529:
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
								if index < 529:
									total_well_alleles += 1
								else:
									total_inova_alleles += 1
							else:
								if index < 529:
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

				continue

			


	print "Total missing counter " + str(missing_counter)

	
	print "Total lines " + str(counter)
	print "Total good lines " + str(good_lines)

	AF_counter = chrom + "\t" + str(counter) + "\t" + str(missing_counter) + "\t" + str(well_001) + "\t" + str(well_005) + "\t" + str(well_005_plus) + \
			"\t" + str(inova_001) + "\t" + str(inova_005) + "\t" + str(inova_005_plus) + "\n"

	print AF_counter
	c.write(AF_counter)
	c.close()
	filt.write(filter_block)
	filter_block = ""
	filt.flush()
	filt.close()


if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main(sys.argv[1])
    #main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)
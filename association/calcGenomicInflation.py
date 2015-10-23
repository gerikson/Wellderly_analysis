import os, sys, gzip, datetime
'''
636259 final_vcf_allChrom_snps_AF0.01-MAF_LD_pruned.bim
9106745 final_vcf_allChrom_snps_AF0.01-results.assoc.logistic
'''

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
 


def median(mylist):
    sorts = sorted(mylist)
    length = len(sorts)
    if not length % 2:
        return (sorts[length / 2] + sorts[length / 2 - 1]) / 2.0
    return sorts[length / 2]

def main():
    ch = "1"

    #RepeatMasker_filename="/gpfs/home/gerikson/wellderly/filter_data/RepeatMasker/byChrom_ALL/RepeatMasker.chr"+ch+".txt"

    input_filename = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/final_vcf_allChrom_snps_AF0.01-results.assoc.logistic"
    LD_pruned = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_association/final_vcf_allChrom_snps_AF0.01-MAF_LD_pruned.bim"

    LD_dict = {}
    with open(LD_pruned) as f:
        for line in f:
            tp_line = line.strip().split("\t")
            LD_dict[tp_line[1]] = "Y"

    print "Variants loaded"

    print 'Prepping RepeatMask...'
    count_rep = 0
    
    '''
    repeatMask_begin = []
    repeatMask_end = []

    with open(RepeatMasker_filename) as f:
        for line in f:
            items = line.split('\t')
            repeatMask_begin.append(int(items[6]))
            repeatMask_end.append(int(items[7]))
            count_rep += 1
    '''

    input_f = open(input_filename)
    p_value_array = []
    counter = 0
    repeatMask_counter = 0
    for line in input_f:
        counter += 1
        if counter %10000 == 0:
            print str(counter)
        tp_line = line.strip().split()
        coord = tp_line[1]
        coord_data = coord.split("_")
        #print coord_data
        #print coord
        try:
            if LD_dict[tp_line[1]]:
                #print tp_line[8]
                if tp_line[8] == "NA":
                    continue
                else:
                    '''
                    new_ch = coord_data[0]
                    if new_ch != ch:
                        print "New Chrom! "
                        ch = new_ch
                        RepeatMasker_filename="/gpfs/home/gerikson/wellderly/filter_data/RepeatMasker/byChrom_ALL/RepeatMasker.chr"+ch+".txt"

                        with open(RepeatMasker_filename) as rep:
                            for lineRep in rep:
                   
                                items = lineRep.split('\t')
                                repeatMask_begin.append(int(items[6]))
                                repeatMask_end.append(int(items[7]))
                                count_rep += 1
                    #Add RepeatMask
                    begin = int(coord_data[1])
                    #print str(begin)
                    end = int(begin) + 1
                    #print str(end)
                    repOverlap = Overlap(begin, repeatMask_begin, end, repeatMask_end)
                    overlapRep = repOverlap.main()

                    if overlapRep:
                        repeatMask_counter += 1
                        continue
                    else:
                        #print "Yay"
                    '''
                    p_value_array.append(float(tp_line[8]))
                    #except:
                    #    print tp_line[8]
        except:
            continue


    print "Lenght of p value array: " + str(len(p_value_array))
    print "Lenght of LD dictionary: " + str(len(LD_dict))
    #print "Total eliminated due to repeats: " + str(repeatMask_counter)
    median_p_value = median(p_value_array)

    #> qchisq(.4937,df=1, lower.tail=FALSE)/qchisq(0.5,df=1,lower.tail=FALSE)
    #[1] 1.02971
    print "Median p-value: " + str(median_p_value)

    input_f.close()


if __name__ == '__main__':

    print "Python Version: " + sys.version

    import time
    start = time.time()
    main()
    end = time.time()
    #print 'This took {0} seconds'.format(end - start)
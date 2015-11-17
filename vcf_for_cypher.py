import gzip, sys, datetime

def main(chrom):
    #infile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_inova_withCorrectGenotypes/wellderly_inova."+chrom+".vcf.gz"
    #outfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/data_for_cypher/wellderly_inova_ noGenotypes_"+chrom+".vcf.gz"
    infile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_wellderly_ALL_filters/"+chrom+"/wellderly_all_filters_"+chrom+".vcf"
    outfile = "/gpfs/group/stsi/data/projects/wellderly/GenomeComb/vcf_annotation/"+chrom+"/wellderly_all_filters_"+chrom+".vcf"
    counter = 0
    buffered = ""
    i = open(infile)
    o = open(outfile, 'w')
    #with gzip.open(infile, 'rb') as i:
    #	with gzip.open(outfile, 'wb') as o:
    for line in i:
        if line[:2] == "##":
            print "header"
            o.write(line)
        elif line[0] == "#":
            tp_line = line.strip().split("\t")
            final_line = "\t".join(tp_line[:10]) + "\n"
            print final_line
            o.write(final_line)

        else:
            counter = counter + 1
            tp_line = line.strip().split("\t")
            
            final_line = "\t".join(tp_line[:10]) + "\n"
            buffered = buffered + final_line
            #print final_line

            if counter%100 == 0:
                o.write(buffered)
                buffered = ""
                o.flush()

            if counter%10000 == 0:
                print datetime.datetime.now().time()
                print str(counter)
                sys.stdout.flush()


    o.write(buffered)
    o.flush()
    i.close()
    o.close()            


if __name__ == '__main__':

    print "Python Version: " + sys.version
    main(sys.argv[1])

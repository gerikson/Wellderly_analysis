#PBS -S /bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l mem=8gb
#PBS -l walltime=240:00:00
#PBS -l cput=9600:00:00
#PBS -m n

set -e
set -o pipefail

# set resources and parameters
source /gpfs/group/stsi/methods/annotation/sg-adviser/variant_calling/bin.exomeMultiSample/parameters.txt
exome=$EXOM/SureSelect.hg19.chr1.bed

# -v args 
echo "args: -v proj_path=$proj_path,sid=$sid"

echo "proj_path=$proj_path"
echo "sid=$sid"
echo "exome=$exome"
echo

cd "$proj_path/bam"

# merge sorted bam files (all lanes/rg) for chr1 only: $sid.chr1.sorted.bam

for f in $sid.*.chr1.sorted.bam; do
    bams="$bams I=$f"
done

$JAVA -Xmx4g -Djava.io.tmpdir=$PBSREMOTEDIR -jar $PIC/MergeSamFiles.jar \
        $bams O=${sid}.chr1.sorted.bam SO=coordinate AS=true \
        CREATE_INDEX=false VALIDATION_STRINGENCY=SILENT
$SAM index ${sid}.chr1.sorted.bam

# generate simple coverage stats: alignable region = "chr1: SureSelect exome"

# calculate non-duplicate percentage
read exon_span <<<$( awk '{span+=($3-$2)}END{print span}' $exome )
read dup_sum <<<$( $SAM depth -b $exome ${sid}.chr1.sorted.bam | awk '{sum+=$3}END{print sum}' )

# calculate percentage of reads in exome and out of exome


total_reads=$( $SAM view -c ${sid}.chr1.sorted.bam)
echo "total_reads=$total_reads"

echo "$SAM view ${sid}.chr1.sorted.bam" > read_count.txt
head -c -1 -q read_count.txt $EXOME_COORD > exome.sh
exome_reads=$( bash exome.sh | wc -l) 

let out_exome=$total_reads-$exome_reads

echo "exome reads=$exome_reads"
echo "out_exome=$out_exome"

# calculate mean insert size of library
# read ins_size <<<$( $SAM view ${sid}.chr1.dedup.bam | awk '$9>0&&$9<600{sum+=$9;cnt++}END{print sum/cnt}' )

$JAVA -Xmx4g -Djava.io.tmpdir=$PBSREMOTEDIR -jar $PIC/CollectInsertSizeMetrics.jar \
        I=${sid}.chr1.dedup.bam O=${sid}.chr1.insert_size.txt \
        H=${sid}.chr1.insert_hs.pdf VALIDATION_STRINGENCY=SILENT
read ins_size <<<$( grep '^MEDIAN_INSERT_SIZE' -A1 ${sid}.chr1.insert_size.txt | tail -n+2 |awk '{print $5}' )

# oupt stats for coverage, nonduplicates, and mean_insert_size


( $SAM depth -b $exome ${sid}.chr1.dedup.bam | \
awk -v sid=$sid -v exon_span=$exon_span -v dup_sum=$dup_sum -v ins_size=$ins_size -v exome_reads=$exome_reads \
-v total_reads=$total_reads -v out_exome=$out_exome '
{depth_sum+=$3} 
$3>=10{c10++} 
END {
  depth=sprintf("%4.1f",depth_sum/exon_span)
  r10=sprintf("%4.1f",100*c10/exon_span)
  nondup=sprintf("%4.1f",100*depth_sum/dup_sum)
  read_exome=sprintf("%4.1f",100*exome_reads/total_reads)
  read_out=sprintf("%4.1f",100*out_exome/total_reads)
  print "sampleID\tmean_coverage\tten_reads\%\tnonduplicate\%\tmean_insert_size\treads_in_exome\%\treads_out_of_exome\%";
  print sid "\t\t" depth "\t\t" r10 "\t\t" nondup "\t\t" ins_size "\t\t" read_exome "\t\t" read_out
}'
) > ../results/${sid}.coverage.txt

exit


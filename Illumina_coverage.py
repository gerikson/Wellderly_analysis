
#count the sum of depth
samtools depth '/projects/wellderly/illumina_bams_forEdico/LP6005833-DNA_D01.bam' | awk '{sum+=$3}END{print sum}' 

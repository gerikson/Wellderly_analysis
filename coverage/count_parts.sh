cd '/gpfs/group/torkamani/bhuvan/wellderly/coverage/SnpFiles'

for i in `seq 1 22`;
do
    t='chrm_'$i'_part'
    temp=`ls | grep $t | wc -l`
    echo chrm $i $temp 		
done

declare -a arr=("X" "Y" "M")
for i in "${arr[@]}";
do
    t='chrm_'$i'_part'
    temp=`ls | grep $t | wc -l`
    echo chrm $i $temp 		
done


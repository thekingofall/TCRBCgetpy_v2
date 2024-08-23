
cat $1 |while read id;
do echo $id
arr=($id)
fq1=${arr[0]}
fq2=${arr[1]}

echo $fq1
echo $fq2
rename $fq1 $fq2 *$2
done 
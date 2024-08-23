cat $1|while read id;
do
arr=($id);
fq1=${arr[0]};
fq2=${arr[1]}；
echo $fq1;
mkdir -p $fq1."mixdata"
/home/maolp/mao/Project/20211110_xcc_TCR/2.cleandata/mixcr-3.0.13/mixcr  analyze amplicon -s  hsa --starting-material RNA --5-end no-v-primers --3-end c-primers --adapters adapters-present $fq1".cdrdata/"$fq1".FQ2.t1.cf.fastq" $fq1".mixdata/"$fq1
done
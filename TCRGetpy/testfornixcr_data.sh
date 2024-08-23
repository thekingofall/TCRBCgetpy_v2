ls *_1.fq.gz > 1
ls *_2.fq.gz > 2
ls *_2.fq.gz |cut -d"/" -f 8|cut -d"." -f 1  > 0
paste 0 1 2  > config.clean 
cat config.clean |while read id;
do
arr=($id);
arr=($id)
fq2=${arr[2]}
fq1=${arr[1]}
sample=${arr[0]}
mkdir -p $sample."mixdata"
/home/maolp/mao/Project/20211110_xcc_TCR/2.cleandata/mixcr-3.0.13/mixcr  analyze amplicon -s  hsa --starting-material RNA --5-end no-v-primers --3-end c-primers --adapters adapters-present $fq1 $fq2 $sample."mixdata"/$sample
done
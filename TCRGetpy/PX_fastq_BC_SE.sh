rm -rf *split
seqkit split2 -1 *1.clean.fq -2 *2.clean.fq -l 1G -j 10
cp ../PXTCR01_main.py *split
cd *split

rm -rf testrun.red
rm -rf testrun.red.completed
ls *A_1* >1
ls *A_2* >2
paste 1 2 >test.clean
cat test.clean|while read id;
do 
echo $id
arr=($id)
fq1=${arr[0]}
fq2=${arr[1]}
echo $fq1
echo $fq2
echo "python PXTCR01_main.py --FQ1 $fq1 --FQ2  $fq2 --Module f " >> testrun.red

done

ParaFly -c testrun.red -CPU 10
echo "now merge file"
cat *FQ1*  > merge.f1
cat *FQ2*  > merge.f2
echo "now merge fastq end "

mv merge* ../
cd ..
name=$(ls *_1.clean.fq)
echo ${name%_*}
n2=${name%_*}
mv merge.f1 $n2".FQ1.fastq"
mv merge.f2 $n2".FQ2.fastq"
rm -rf *split
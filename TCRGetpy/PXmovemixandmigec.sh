echo $1
echo $2
echo $num
fold="migec".$2.$1".checkfold"
rm -rf $fold
echo $fold
mkdir -p $fold
for i in *cdrdata;
do 
echo $i
echo $i/*$1.$2.clean.txt
(cp -r $i/*$1.$2.clean.txt $fold)


done


fold2="mixcr".$2".checkfold"
echo $fold2
#rm -rf fold2
mkdir -p $fold2
for i in *mixcrdata;
do 
echo $i
echo $i/*$2.txt
(cp -r $i/*$2.txt $fold2)
done

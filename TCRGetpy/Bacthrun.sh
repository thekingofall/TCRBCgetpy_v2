if [ $# -lt 2 ];
then
echo "the argument is less than 3:<< bash thisrun.sh <get|plot|statplot|movefile|> <BCfile> <foldname><Module b |bc>"
else


for fq1 in *_1.clean.fq;
do 

fq2=${fq1%%_1.clean.fq}"_2.clean.fq";

fq=${fq1%%_1.clean.fq}

if [ $1 == "get" ]; 
then
python PXTCR01_main.py --FQ1 $fq1 --FQ2 $fq2 --Module $4 --BC $2

python  PXTCR01_main.py --FQ1 $fq1 --FQ2 $fq2 --Statwhat a --BC $2
python PXTCR01_main.py --Statindex a1 --BC $2 --Num 3
echo "hello $2"; 
elif [ $1 == "plot" ]
then
TRAfile="TRAdir."$3
TRBfile="TRBdir."$3
TRAB_raw_file="TRAB_raw_dir."$3
rm -rf $TRAfile
rm -rf $TRBfile
mkdir -p $TRAfile
mkdir -p $TRBfile
mkdir -p $TRAB_raw_file


for d in *cdrdata;do echo $d ;(cp -r $d/*TRA.clean.txt $TRAfile/);done
python PXTCR01_main.py --Statindex plot --DF  $TRAfile/ --SV TRA_clean_plot

for j in *cdrdata;do echo $j ;(cp -r $j/*TRB.clean.txt  $TRBfile/);done
python PXTCR01_main.py --Statindex plot --DF  $TRBfile/ --SV TRB_clean_plot
for h in *cdrdata;do echo $h ;(cp -r $h/*TRA,TRB.cdrblast.txt  $TRAB_raw_file/);done
python PXTCR01_main.py --Statindex plot --DF  $TRAB_raw_file/ --SV TRAB_raw_plot


elif [ $1 == "statplot" ]
then

python  PXTCR01_main.py --FQ1 $fq1 --FQ2 $fq2 --Statwhat a --BC $2
for i in {1,2,3};do echo $i;python PXTCR01_main.py --Statindex a1 --BC $2 --Num $i ;done

time=$(date "+%Y-%m-%d-%H")
echo $time

TRAfile="TRAdir."$3.$time
TRBfile="TRBdir."$3.$time
TRAB_raw_file="TRAB_raw_dir."$3.$time
rm -rf $TRAfile
rm -rf $TRBfile
mkdir -p $TRAfile
mkdir -p $TRBfile
mkdir -p $TRAB_raw_file


for d in *cdrdata;do echo $d ;(cp -r $d/*TRA.clean.txt $TRAfile/);done
python PXTCR01_main.py --Statindex plot --DF  $TRAfile/ --SV TRA_clean_plot

for j in *cdrdata;do echo $j ;(cp -r $j/*TRB.clean.txt  $TRBfile/);done
python PXTCR01_main.py --Statindex plot --DF  $TRBfile/ --SV TRB_clean_plot
for h in *cdrdata;do echo $h ;(cp -r $h/*TRA,TRB.cdrblast.txt  $TRAB_raw_file/);done
python PXTCR01_main.py --Statindex plot --DF  $TRAB_raw_file/ --SV TRAB_raw_plot


elif [ $1 == "movefile" ]
then
TRAfile="TRAdir."$3
TRBfile="TRBdir."$3
TRAB_raw_file="TRAB_raw_dir."$3
rm -rf $TRAfile
rm -rf $TRBfile
rm -rf $TRAB_raw_file
mkdir -p $TRAfile/dir.t{1..3}
mkdir -p $TRBfile/dir.t{1..3}
mkdir -p $TRAB_raw_file/dir.t{1..3}


for d in *cdrdata;do echo $d ;(cp -r $d/*t1.TRA.clean.txt $TRAfile/dir.t1);(cp -r $d/*t2.TRA.clean.txt $TRAfile/dir.t2);(cp -r $d/*t3.TRA.clean.txt $TRAfile/dir.t3);done
#python PXTCR01_main.py --Statindex plot --DF  $TRAfile/ --SV TRA_clean_plot


for d in *cdrdata;do echo $d ;(cp -r $d/*t1.TRB.clean.txt $TRBfile/dir.t1);(cp -r $d/*t2.TRB.clean.txt $TRBfile/dir.t2);(cp -r $d/*t3.TRB.clean.txt $TRBfile/dir.t3);done
#python PXTCR01_main.py --Statindex plot --DF  $TRAfile/ --SV TRA_clean_plot

#python PXTCR01_main.py --Statindex plot --DF  $TRBfile/ --SV TRB_clean_plot
for h in *cdrdata;do echo $h ;(cp -r $h/*t1TRA,TRB.cdrblast.txt  $TRAB_raw_file/dir.t1);(cp -r $h/*t2TRA,TRB.cdrblast.txt  $TRAB_raw_file/dir.t2);(cp -r $h/*t3TRA,TRB.cdrblast.txt  $TRAB_raw_file/dir.t3);done
else 
echo  "wanring:nothing happend!pelase choose get or plot"


fi; 
done
fi
# cd $TRAfile
# rename P7M1S $1  P7M1S*
# cd ..
# cd $TRBfile
# rename P7M1S $1  P7M1S*

#mkdir TRBtest;for i in XCC*;do echo $i; cp -r $i/*TRBdir/*  TRBtest/;done

#mkdir PXDIR
#mv P7* PXDIR
# rename PXM1S XCC02
#rename P7M1S XCC02 *txt

import os,sys
BC=sys.argv[1]
foldname=sys.argv[2]
runfile=" "
arr=open(BC,"r")
for line in arr: 
    fqname2=line.rstrip().split("\t")[0]
    fqfile=fqname2+"_mixcrdata"
    os.system("mkdir -p "+fqfile)
    fqtest1=fqname2+".FQ1.fastq"
    fqtest2=fqname2+".FQ2.fastq"
    fqseq=line.rstrip().split("\t")[1]
    print(fqname2)
    print(fqseq)
    os.system("java -jar /home/maolp/mao/Project/20211110_xcc_TCR/2.cleandata/mixcr-3.0.13/mixcr.jar  analyze amplicon -s hsa --starting-material DNA --5-end no-v-primers --3-end c-primers --adapters adapters-present "+fqtest1+" "+fqtest2 +" "+fqfile+"/"+foldname+"_"+fqname2)
    
    


import sys,os
arr=open(sys.argv[1],"r")
fq1=sys.argv[2]
fq2=sys.argv[3]
os.system("rm -rf run_data2.txt")

for line in arr: 
    fqname=line.rstrip().split("\t")[0]+".name"
    fqseq=line.rstrip().split("\t")[1]
    print(fqname)
    print(fqseq)
    os.system("cat "+fq1+" |grep -B1 "+fqseq+"|grep @|awk '{print $1}'> "+fqname)

    os.system("echo \"python PXTCR03_BCsplit.py "+str(fqname)+" " +str(fq1)+"\">>run_data2.txt")
    os.system("echo \"python PXTCR03_BCsplit.py "+str(fqname)+" " +str(fq2)+"\">>run_data2.txt")

os.system("ParaFly -c run_data2.txt -CPU 8")
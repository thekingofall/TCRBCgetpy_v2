
import sys, os,glob,re

data=sys.argv[1]
dataw=open("stat_file_name.txt","w")
for i in glob.glob(data+"/*"):
    if re.findall("FQ1.fastq",i):
        dataw.write(i+"\n")

        print(i)
        name=i.split("/")[-1].split(".")[0]
        stat_file_row=os.popen("cat "+i+"| wc -l ").read() 
        Fq_data_count=round(int(stat_file_row)/4,0)
        print(name+"\t"+str(Fq_data_count))
        dataw.write(name+"\t"+str(Fq_data_count)+"\n")
    elif re.findall("1.clean.fq.gz",i) or re.findall("R1_001.fastq.gz",i):
        dataw.write(i+"\n")
        print(i)
        name=i.split("/")[-1].split(".")[0]+"_clean"
        stat_file_row=os.popen("zcat "+i+"| wc -l ").read()
        Fq_data_count=round(int(stat_file_row)/4,0)
        print(name+"\t"+str(Fq_data_count))
        dataw.write(name+"\t"+str(Fq_data_count)+"\n")

stat_file_name=glob.glob("*FQ1.fastq")



# name=stat_file_name.split(".")[0]
for line in stat_file_name:
    name=line.split(".")[0]
    stat_file_row=os.popen("cat "+line+"| wc -l ").read()
    Fq_data_count=round(int(stat_file_row)/4,0)

    print(name+"\t"+str(Fq_data_count))
    dataw.write(name+"\t"+str(Fq_data_count)+"\n")
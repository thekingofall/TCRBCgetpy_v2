import os,sys

import pandas as pd
import re
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
inputnamefile=sys.argv[1]
R1=sys.argv[1]
R2=sys.argv[2]
samseq=sys.argv[3]
name1=str(R1).split(".")[0].split("_1")[0]
rawname=name1+"R1.name"
def tran_fast(rawname=rawname,dat1=R1,dat2=R2,samseq=samseq):

    os.system("cat "+dat1+" |grep -B1 "+samseq+"|grep @|awk '{print $1}'>> "+ rawname )
    #os.system("cat "+dat2+" |grep -B1 "+samseq+"|grep @|awk '{print $1}'>> "+ rawname )

tran_fast()


inputfile=sys.argv[2]


def BCsplit_BC(inputnamefile,inputfile,FQ):
  print(inputnamefile)

  df_id=pd.read_csv(inputnamefile,names=["name"])

  #print(name)
  name_list=set(df_id["name"].values)
  outputname=str(inputnamefile).split(".")[0]
  out_handle = open(outputname+"."+FQ+".fastq", "w+")
  
  index=0
  for title,seq,qual in FastqGeneralIterator(open(inputfile)):
    # if index==10:
    #   break
    title1="@"+title.split(" ")[0]
    
        #print(title1)
    if title1 in name_list:
      seq1=seq.split(str(samseq))
      if len(seq1)==2:
        umi=seq1[1][0:16]
        if len(umi)>1 and str(umi[0])=="T":
          UMI=umi
          seq_record="\n".join(["@"+title+"   BC   UMI:"+UMI, seq, "+", qual])
          out_handle.write(seq_record)
          out_handle.write("\n")
    index +=1
  print(index)
  out_handle.close()


  
def BCsplit_SE(inputnamefile,inputfile,FQ):
  print(inputnamefile)

  df_id=pd.read_csv(inputnamefile,names=["name"])

  #print(name)
  name_list=set(df_id["name"].values)
  outputname=str(inputnamefile).split(".")[0]
  out_handle = open(outputname+"."+FQ+".fastq", "w+")
  index=0
  for title,seq,qual in FastqGeneralIterator(open(inputfile)):
    # if index==10:
    #   break
    title1="@"+title.split(" ")[0]
    #print(title1)
    if title1 in name_list:
      seq_record="\n".join(["@"+title, seq, "+", qual])
      out_handle.write(seq_record)
      out_handle.write("\n")
  

    index +=1
  print(index)
  out_handle.close()

BCsplit_BC(inputnamefile=rawname,inputfile=R1,FQ="FQ1_1BCtem")
os.system("cat *FQ1_1BCtem* |grep @ |awk '{print $1}' > FQ1_1BCtem.name")
rawname2="FQ1_1BCtem.name"
BCsplit_SE(inputnamefile=rawname2,inputfile=R2,FQ="FQ2_2SQtem")

# BCsplit(inputnamefile=rawname,inputfile=R1,FQ="FQ1_tem")

# BCsplit(inputnamefile=rawname,inputfile=R1,FQ="FQ1_tem")
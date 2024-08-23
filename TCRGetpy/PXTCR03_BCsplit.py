
import os
import pandas as pd
import sys
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
inputnamefile=sys.argv[1]
print(inputnamefile)
inputfile=sys.argv[2]
FQ=inputfile.split(".")[1]
print(FQ)
def BCsplit(inputnamefile,inputfile,FQ):

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

BCsplit(inputnamefile=inputnamefile,inputfile=inputfile,FQ=FQ)

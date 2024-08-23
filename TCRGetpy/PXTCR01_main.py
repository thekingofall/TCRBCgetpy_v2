
import re,os,sys
import datetime
import argparse
import pandas as pd
import gzip
import configparser

# from Bio import SeqIO
# from Bio.SeqIO.QualityIO import FastqGeneralIterator
# from  Atestforarg import  testfordata



parser = argparse.ArgumentParser(description="usage:<FQ1><FQ2><TCR>")
Filename=parser.add_argument_group('Filname',"Fastq name")
TCRargument=parser.add_argument_group('TCRargument',"Some argument about Migec")
Modulearg=parser.add_argument_group('Module argument',"Some argument about Module")
#Filename.add_argument("--FQ1", type=str, help="Fastq R1",required=True)
Filename.add_argument("--FQ1", type=str, help="Fastq R1")
Filename.add_argument("--FQ2", type=str, help="Fastq R2")
Filename.add_argument("--BC", type=str, help="file with BC name and BC seq")
Filename.add_argument("--DF", type=str, help="fold with want to get their stat")
Filename.add_argument("--SV", type=str, help="fold with want to Save DF stat")
Filename.add_argument("--ini", type=str, help="cogfig.ini")
Filename.add_argument("--FN", type=str, help="filename you want to define",default="Filename")
TCRargument.add_argument("--TCR", type=str, help="TCR types:<TRA|TRB|TRA,TRB|...>,defalut:TRA,TRB",default="TRA,TRB")
TCRargument.add_argument("--Num", type=str, help="TCR num,default:3",default="3")

Modulearg.add_argument("--Module",type=str, help="Some Module you want to run:"+"\n "+"a: run all;f: only for transform ")
Modulearg.add_argument("--Statwhat",type=str, help="Some data you want to havce a stat:"+"\n "+"a:run all;a1:run some ;p:only for primary data",default="p")
Modulearg.add_argument("--BothSeq",type=str, help="Some same Seq in your data:"+"\n "+"default:CAGTGGTATCAACGCAGAG;Others seq you need",default="CAGTGGTATCAACGCAGAG")
Modulearg.add_argument("--Statindex",type=str, help=" plot:what folder you want to have a plot and stat ,but if want to use this module ,you are better run the Statwhat module before")
args= parser.parse_args()

def transformer(FRdata,FRdata2,WRdata1,WRdata2,same_seq="CAGTGGTATCAACGCAGAG"):
    while True:
        line1=FRdata.readline().rstrip().decode('utf-8')
        if not line1:
        #if index==10:
            break
        line2=FRdata.readline().rstrip().decode('utf-8')
        lineUMI_before=line2.split(same_seq) #这是共有barcode序列
        line3=FRdata.readline().rstrip().decode('utf-8')
        line4=FRdata.readline().rstrip().decode('utf-8')
        line2_1=FRdata2.readline().rstrip().decode('utf-8') 
        line2_2=FRdata2.readline().rstrip().decode('utf-8')
        line2_3=FRdata2.readline().rstrip().decode('utf-8')
        line2_4=FRdata2.readline().rstrip().decode('utf-8')  
        if len(lineUMI_before)==2 and len(lineUMI_before[1][0:23])>1:
            umi=lineUMI_before[1][0:23]
            # print(umi)
            # print(len(umi))
            if len(umi)==23 and str(umi[0])=="T" and str(umi[5])=="T": 
                lineUMI=umi
                line4_start=len(lineUMI_before[0])
                line4_end=line4_start+len(lineUMI)
                score=line4[line4_start:line4_end]
                data1=line1+"\t"+"BC"+" "+"UMI:"+lineUMI+":"+score+"\n"+line2+"\n"+line3+"\n"+line4
                #print(data1)
                WRdata1.write(data1+"\n")
                data2=line2_1+"\t"+"SE"+" "+"UMI:"+lineUMI+":"+score+"\n"+line2_2+"\n"+line2_3+"\n"+line2_4
                #print(data2)
                WRdata2.write(data2+"\n")
        else:
            continue


    WRdata1.close()
    WRdata2.close()

def migecrun(name,num=str(3),TCRtypes="TRA,TRB"):
        # migecpath=os.getcwd()
        os.system("java -jar   "+ migecpath + " Assemble -m "+num+" -q 15 --mask 0:1 --filter-collisions  "+ name+".FQ1.fastq  " + name+".FQ2.fastq  "+ name+".cdrdata")   
        os.system("java -jar   "+ migecpath + " CdrBlast -a -R "+TCRtypes+" "+name+".cdrdata/"+name+".FQ2.t"+num+".cf.fastq  "+name+".cdrdata/"+name+"_t"+num+"TRA,TRB.cdrblast.txt")

# def migecrun(name,num=str(3),TCRtypes="TRA,TRB",SPECIES="HomoSapiens"):
#         # migecpath=os.getcwd()
#         os.system("java -jar   "+ migecpath + " Assemble -m  "+num+" -q 15 --mask 0:1 --filter-collisions  "+ name+".FQ1.fastq  " + name+".FQ2.fastq  "+ name+".cdrdata")   
#         os.system("java -jar   "+ migecpath + " CdrBlast -S "+SPECIES+"  -a -R "+TCRtypes+" "+name+".cdrdata/"+name+".FQ2.t"+num+".cf.fastq  "+name+".cdrdata/"+name+"_t"+num+"TRA,TRB.cdrblast.txt")

def migecrun1(name,num=str(3),TCRtypes="TRA,TRB"):
        os.system("java -jar   "+ migecpath + " Assemble -m "+num+" -q 15 --mask 0:1 --filter-collisions  "+ name+".FQ1.fastq  " + name+".FQ2.fastq  "+ name+".cdrdata")
def migecrun2(name,num=str(3),TCRtypes="TRA,TRB"):        
        os.system("java -jar   "+ migecpath + " CdrBlast -a -R "+TCRtypes+" "+name+".cdrdata/"+name+".FQ2.t"+num+".cf.fastq  "+name+".cdrdata/"+name+"_t"+num+"TRA,TRB.cdrblast.txt")
        # os.system("java -jar   "+ migecpath + " CdrBlast -S "+SPECIES+"  -a -R "+TCRtypes+" "+name+".cdrdata/"+name+".FQ2.t"+num+".cf.fastq  "+name+".cdrdata/"+name+"_t"+num+"TRA,TRB.cdrblast.txt")
def migecrun3(name,num=str(3),TCRtypes="TRA,TRB",SPECIES="HomoSapiens"):        
        os.system("java -jar   "+ migecpath + " CdrBlast -a -R "+TCRtypes+" "+name+".cdrdata/"+name+".FQ2.t"+num+".cf.fastq  "+name+".cdrdata/"+SPECIES+"_"+name+"_t"+num+"TRA,TRB.cdrblast.txt")
        # os.system("java -jar   "+ migecpath + " CdrBlast -S "+SPECIES+"  -a -R "+TCRtypes+" "+name+".cdrdata/"+name+".FQ2.t"+num+".cf.fastq  "+name+".cdrdata/"+name+"_t"+num+"TRA,TRB.cdrblast.txt")

def migecrunfast(name,num=str(3),TCRtypes="TRA,TRB"):
    migeccpu=open("run_migec_cpu.txt","a+")
    command="java -jar   "+ migecpath + " Assemble -m "+num+" -q 15 --mask 0:1 --filter-collisions  "+ name+".FQ1.fastq  " + name+".FQ2.fastq  "+ name+".cdrdata"
    migeccpu.write(command+"\n")


        

def migecrunfast(name,num=str(3),TCRtypes="TRA,TRB"):
    #os.system("rm  -rf run_migec_cpu.txt")
    migeccpu=open("run_migec_cpu.txt","a+")
    command="java -jar   "+ migecpath + " Assemble -m "+num+" -q 15 --mask 0:1 --filter-collisions  "+ name+".FQ1.fastq  " + name+".FQ2.fastq  "+ name+".cdrdata"
    migeccpu.write(command+"\n")

def migecfast(BC2,n=4,TCRtypes="TRA,TRB"):
    print("Now fast migec")
    arr=open(BC2,"r")
    os.system("rm  -rf run_migec_cpu.txt")

    for line in arr: 
        
        _fqname3=line.rstrip().split("\t")[0]
        #_fqseq3=line.rstrip().split("\t")[1]
        
        #print(fqname2)
        #print(fqseq)

        for i in range(1,n):
            num2=str(i)
            print(i)
            migecrunfast(name=_fqname3,num=num2,TCRtypes="TRA,TRB") 
    os.system("ParaFly -c run_migec_cpu.txt -CPU 24")
    return(_fqname3)   




def trans_data(name):
    WR1=open(name+".FQ1.fastq","w+")
    WR2=open(name+".FQ2.fastq","w+")
    transformer(FRdata=gzip.open(args.FQ1,"r"),FRdata2=gzip.open(args.FQ2,"r"),WRdata1=WR1,WRdata2=WR2,same_seq=Sameseq)
    WR3=open(name+".FQ1.fastq","a")
    WR4=open(name+".FQ2.fastq","a")
    transformer(FRdata=gzip.open(args.FQ2,"r"),FRdata2=gzip.open(args.FQ1,"r"),WRdata1=WR3,WRdata2=WR4,same_seq=Sameseq)


def draw_picture(name,num=str(3),TCRtypes="TRA,TRB"):
    run="Rscript   "+ Drawpiepath + "  "+ name+".cdrdata/"+name+"_t"+num+"TRA,TRB.cdrblast.txt     "+name+".cdrdata/"+name+"_"+TCRtypes+"_t"+num
    print(run)
    os.system(run)

def draw_picture2(name,num=str(3),TCRtypes="TRA,TRB",SPECIES="HomoSapiens"):
    run="Rscript   "+ Drawpiepath + "  "+ name+".cdrdata/"+SPECIES+"_"+name+"_t"+num+"TRA,TRB.cdrblast.txt     "+name+".cdrdata/"+SPECIES+"_"+name+"_"+TCRtypes+"_t"+num
    print(run)
    os.system(run)

def testfordata(R,endnum=10,Tseq="T....T....T....TCTTGGG"):
    FR=open(R,"r")

    WR1=open("FQ1.fastq","a")
    index=0
    while True:
        line1=FR.readline().rstrip().decode('utf-8')
        if index==endnum:
            break
        line2=FR.readline().rstrip().decode('utf-8')
        line3=FR.readline().rstrip().decode('utf-8')
        line4=FR.readline().rstrip().decode('utf-8')
        lineUMI_before=re.findall(Tseq,line2)
        if len(lineUMI_before)!=0:
            print(lineUMI_before[0])
            line2_dist=line2.split(lineUMI_before[0])
            line4_start=len(line2_dist[0])
            line4_end=line4_start+len(lineUMI_before[0])
            score=line4[line4_start:line4_end]
            #print(line2_dist)
            print(line4_start)
            print(line4_end)
            print(score)

        index+=1
    return index   


def BC_split(BC,fq1,fq2):
    os.system("rm -rf run_data2.txt")
    os.system("rm -rf run_data2.txt.completed")
    arr=open(BC,"r")

    for line in arr: 
        print(line)
        fqname=line.rstrip().split("\t")[0]+".name"
        fqseq=line.rstrip().split("\t")[1]
        print(fqname)
        print(fqseq)
        os.system("cat "+fq1+" |grep -B1 "+fqseq+"|grep @|awk '{print $1}'> "+fqname)

        os.system("echo \"python   "+  PXTCR03_BCsplitpath + " "+str(fqname)+" " +str(fq1)+"\">>run_data2.txt")
        os.system("echo \"python   "+  PXTCR03_BCsplitpath + " "+str(fqname)+" " +str(fq2)+"\">>run_data2.txt")

    os.system("ParaFly -c run_data2.txt -CPU 24")
    return(fqname)
    # run_command1="python PXTCR02_OS.py "+ BC +" "+fq1+" " +fq2+""
    # os.system(run_command1)
def count_line(stat_file_name):
    _stat_file_row=os.popen("cat "+stat_file_name+"| wc -l ").read()
    _Fq_data_count=round(int(_stat_file_row)/4)
    _lsdata=list()
    _data1=print('the reads of {} is {}'.format(stat_file_name,_Fq_data_count))
    _lsdata.append(_Fq_data_count)
    _lsdata.append(_data1)
    return _lsdata
def stat_module(FQone=args.FQ1,stat_file="stat_file",stat_what=args.Statwhat):
    _all_stat_num=list()
    # os.system("echo 'Reads in "+FQone+"'")

    # val1=os.popen("cat "+FQone+"| wc -l ").read()

    fq_count=count_line(FQone)[0]
    # fastq_count_val1=os.popen("cat "+fq1+"| wc -l ").read()
    # fastq_count=round(int(fastq_count_val1)/4)
    fastq_count=count_line(fq1)[0]
    _all_stat_num.append(fq_count)
    _all_stat_num.append(fastq_count)
    fastq_pop=str(fastq_count/fq_count*100)+"%"
    stat_file=open(stat_file,"a+")
    stat_file.write("File"+","+"Reads_count"+","+"All_Percentage"+"\n")
    stat_file.write(FQone +","+str(fq_count)+","+"100%"+"\n")
    stat_file.write(fq1 +","+str(fastq_count)+","+str(fastq_pop)+"\n")
    stat_file.close()
    return _all_stat_num


def tran_fast(fname,dat1,dat2,samseq,R=args.FQ1,R2=args.FQ2):
    rawname=fname+".name"

    os.system("cat "+dat1+" |grep -B1 "+samseq+"|grep @|awk '{print $1}'>> "+ rawname )
    os.system("cat "+dat2+" |grep -B1 "+samseq+"|grep @|awk '{print $1}'>> "+ rawname )
    os.system("python   "+  PXTCR03_BCsplitpath + " "+rawname+" "+ R)
    os.system("python   "+  PXTCR03_BCsplitpath + " "+rawname+" "+ R2)



if __name__=='__main__':   
    starttime = datetime.datetime.now()
    if not os.path.exists(str(args.ini)) :
        print("config.ini file not exist")
    else:
        print("ini file exist")
        con = configparser.ConfigParser()
        con.read(args.ini,encoding='utf-8')
        
        software=dict(con.items("software"))
        # print(software)
        # print(software["migec"].rstrip('"').split('"')[-1])
        migecpath=os.path.join(str(os.getcwd()),software["migec"].rstrip('"').split('"')[-1])
        Drawpiepath=os.path.join(str(os.getcwd()),software["drawpie"].rstrip('"').split('"')[-1])
        PXTCR03_BCsplitpath=os.path.join(str(os.getcwd()),software["pxtcr03_bcsplit"].rstrip('"').split('"')[-1])
        PXTCR02_OSpath=os.path.join(str(os.getcwd()),software["pxtcr02_os"].rstrip('"').split('"')[-1])
        Draw_immu_index=os.path.join(str(os.getcwd()),software["draw_immu_index"].rstrip('"').split('"')[-1])

        # print(Drawpiepath)
        # print(migecpath)
        
        R=args.FQ1
        R2=args.FQ2
        TCRtypes=args.TCR
        num=str(args.Num)
        
        name=str(R).split(".")[0].split("_1")[0]
        name=args.FN
        
        Stattext=str(name)+"_stat.csv"
        print(name)
        print(Stattext)
        verbosity=args.Module
        fq1=name+".FQ1.fastq"
        fq2=name+".FQ2.fastq"
        BC=args.BC
        Sameseq=args.BothSeq
        Statornot=args.Statwhat
        #print(name)
        print("Now ,Start transform")


        if verbosity == "a":
            print("this module whill help you get the BC file,then Split the BC as you need ,then run a migec in three times ")
            trans_data(name=name)
            migecrun1(name=name,num=num,TCRtypes=TCRtypes)
            migecrun2(name=name,num=num,TCRtypes=TCRtypes)
            draw_picture(name,num=num,TCRtypes=TCRtypes)
            arr=open(BC,"r")
            for line in arr: 
                fqname2=line.rstrip().split("\t")[0]
                fqseq=line.rstrip().split("\t")[1]
                print(fqname2)
                print(fqseq)
                for i in range(1,4):
                    num2=str(i)
                    print(i)
                    migecrun(name=fqname2,num=num2,TCRtypes=TCRtypes)
                    
                    draw_picture(fqname2,num=num2,TCRtypes=TCRtypes)
        elif verbosity == "f":
            print("this module get the Reads of the FQ file")
            print(name)
            trans_data(name=name)
        # elif verbosity == "sb":
        #     print("this module get the Reads of the FQ file")
        #     print(name)
            # trans_data(name=name)
        elif verbosity=="fs":
            print("this module run the migec")
            print(name)
            trans_data(name=name)
            migecrun1(name=name,num=num,TCRtypes=TCRtypes)
            migecrun2(name=name,num=num,TCRtypes=TCRtypes)
            draw_picture(name,num=num,TCRtypes=TCRtypes)
        
        elif verbosity=="fs2":
            print(name)
            trans_data(name=name)
            migecrun(name=name,num=num,TCRtypes=TCRtypes)
            # migecrun2(name=name,num=num,TCRtypes=TCRtypes)
            draw_picture(name,num=num,TCRtypes=TCRtypes)
        elif verbosity=="fsm":
            print("this module run the migec")
            print(name)
            # trans_data(name=name)
            # migecrun1(name=name,num=num,TCRtypes=TCRtypes)
            migecrun3(name=name,num=num,TCRtypes="TRB",SPECIES="MacacaMulatta")
            draw_picture2(name,num=num,TCRtypes="TRB",SPECIES="MacacaMulatta")
        elif verbosity=="fsma":
            print("this module run the migec")
            print(name)
            trans_data(name=name)
            migecrun1(name=name,num=num,TCRtypes=TCRtypes)
            migecrun3(name=name,num=num,TCRtypes="TRB",SPECIES="MacacaMulatta")
            draw_picture2(name,num=num,TCRtypes="TRB",SPECIES="MacacaMulatta")

        elif verbosity=="s":
            print("this module run the migec")
            migecrun1(name=name,num=num,TCRtypes=TCRtypes)
            migecrun2(name=name,num=num,TCRtypes=TCRtypes)
            draw_picture(name,num=num,TCRtypes=TCRtypes)
        elif verbosity=="t":
            draw_picture(name,num=num,TCRtypes=TCRtypes)
        elif verbosity=="o":
            testfordata(R)
        elif verbosity=="b":
            print("this module whill help you get the BC file,then Split the BC as you need ,then run a migec")
            trans_data(name=name)
            BC_split(BC=BC,fq1=fq1,fq2=fq2)
            
            arr=open(BC,"r")
            for line in arr: 
                fqname2=line.rstrip().split("\t")[0]
                fqseq=line.rstrip().split("\t")[1]
                print(fqname2)
                print(fqseq)
                migecrun1(name=fqname2,num=num,TCRtypes=TCRtypes)
                migecrun2(name=fqname2,num=num,TCRtypes=TCRtypes)
                draw_picture(fqname2,num=num,TCRtypes=TCRtypes)
        elif verbosity=="bc":
            BC_split(BC=BC,fq1=fq1,fq2=fq2)
            arr=open(BC,"r")
            for line in arr: 
                fqname2=line.rstrip().split("\t")[0]
                fqseq=line.rstrip().split("\t")[1]
                print(fqname2)
                print(fqseq)
                fq_in_count=count_line(fqname2+".FQ1.fastq")
                
                migecrun(name=fqname2,num=num,TCRtypes=TCRtypes)
                draw_picture(fqname2,num=num,TCRtypes=TCRtypes)
        elif verbosity=="bm":
            # BC_split(BC=BC,fq1=fq1,fq2=fq2)
            arr=open(BC,"r")
            for line in arr: 
                fqname2=line.rstrip().split("\t")[0]
                fqseq=line.rstrip().split("\t")[1]
                print(fqname2)
                print(fqseq)
                # fq_in_count=count_line(fqname2+".FQ1.fastq")
                print(migecpath)
                
                migecrun(name=fqname2,num=num,TCRtypes=TCRtypes)
                draw_picture(fqname2,num=num,TCRtypes=TCRtypes)
        
        elif verbosity=="b3":
            arr=open(BC,"r")
            for line in arr: 
                fqname2=line.rstrip().split("\t")[0]
                fqseq=line.rstrip().split("\t")[1]
                print(fqname2)
                print(fqseq)
                for i in range(1,4):
                    num2=str(i)
                    print(i)
                    migecrun(name=fqname2,num=num2,TCRtypes=TCRtypes)
                    
                    draw_picture(fqname2,num=num2,TCRtypes=TCRtypes)
                # os.system("mkdir "+fqname2)
                # os.system()
        elif args.Statindex=="a":
            print("this module will help stat the reads of fastq file")
                    
            os.system("rm -rf "+Stattext)
            os.system("touch "+Stattext)
            print(name)
            #print(Stattext)
            data_tem=stat_module(FQone=args.FQ1,stat_file=Stattext)
            print('all is{} ,BC is {},Percenge is {}%'.format(str(data_tem[0]),str(data_tem[1]),str(data_tem[1]/data_tem[0]*100)))
            if args.Statwhat=="a" and args.BC :
                arr=open(BC,"r")
                tem_data=open(Stattext,"a+")
                for line in arr: 
                    fqname2=line.rstrip().split("\t")[0]
                    fqseq=line.rstrip().split("\t")[1]
                    fq_tem=count_line(fqname2+".FQ1.fastq")
                    tem_data.write(fqname2 +","+str(fq_tem[0])+","+str()+str(fq_tem[0]/data_tem[1]*100)+"%"+"\n")
                arr.close()
                tem_data.close()
        elif args.Statindex=="a1" and args.BC:
            print("this module will help you stat the migec file ")
            #import pandas as pd
            num3=str(3)
            num3=num
            arr=open(BC,"r")
            # dat2=pd.DataFrame(columns=["fqname_types","raw_types","clean_atypes","TRA_types","TRB_types"])
            
            dat2=pd.DataFrame(columns=["fqname_types","raw_types","raw_sum","clean_types","clean_sum","clean_TRA_types","clean_TRA_sum","clean_TRB_types","clean_TRB_sum","raw_max_cdr3_seq","clean_max_cdr3_seq",
            'raw_V_types',"clean_V_types","clean_max_Vsegments","clean_max_V_num","clean_J_types","clean_max_J_num","clean_max_Jsegments"])
            for line in arr: 
                fqname2=line.rstrip().split("\t")[0]
                fqseq=line.rstrip().split("\t")[1]
                cdr_tem=pd.read_csv(fqname2+".cdrdata/"+fqname2+"_t"+str(num3)+"TRA,TRB.cdrblast.txt",header=0,sep="\t")
                print('the types of {} raw is {}'.format(fqname2,cdr_tem.shape[0]))
            
                cdr_clean=cdr_tem[cdr_tem["CDR3 amino acid sequence"].str.contains("\?")==False]
                print('the types of clean of {} is {}'.format(fqname2,cdr_clean.shape[0]))
                cdr_clean.to_csv(fqname2+".cdrdata/"+fqname2+"_t"+str(num3)+"TRA,TRB.cdrblast.clean.csv",index=False)
                TRA=cdr_clean[cdr_clean["V segments"].str.contains("TRA")]
                print('the types of TRA in {} is {} '.format(fqname2,TRA.shape[0]))
                TRA.to_csv(fqname2+".cdrdata/"+fqname2+"_t"+str(num3)+".TRA.clean.txt",index=False,sep="\t")
                TRB=cdr_clean[cdr_clean["V segments"].str.contains("TRB")]
                print('the types of TRB in {} is {} '.format(fqname2,TRB.shape[0]))
                TRB.to_csv(fqname2+".cdrdata/"+fqname2+"_t"+str(num3)+".TRB.clean.txt",index=False,sep="\t")
                # dat={'fqname_types':[fqname2],'raw_types':[str(cdr_tem.shape[0])],'clean_types':[cdr_clean.shape[0]],'TRA_types':[TRA.shape[0]],'TRB_types':[TRB.shape[0]]}
                # dat=pd.DataFrame(dat)   
                dat={
                'fqname_types':[fqname2],
                'raw_types':[str(cdr_tem.shape[0])],
                'raw_sum':[list(cdr_tem.iloc[:,[0]].sum())[0]],
                'clean_types':[cdr_clean.shape[0]],
                'clean_sum':[list(cdr_clean.iloc[:,[0]].sum())[0]],
                'clean_TRA_types':[TRA.shape[0]],
                'clean_TRA_sum':[list(TRA.iloc[:,[0]].sum())[0]],
                'clean_TRB_types':[TRB.shape[0]],
                'clean_TRB_sum':[list(TRB.iloc[:,[0]].sum())[0]],
                'raw_max_cdr3_seq':cdr_tem.iloc[:,[3]].iloc[0][0],
                'clean_max_cdr3_seq':cdr_clean.iloc[:,[3]].iloc[0][0],
                'raw_V_types':[len(cdr_tem["V segments"].unique())],
                "clean_V_types":[len(cdr_clean["V segments"].unique())],
                "clean_max_Vsegments":[cdr_tem["V segments"].value_counts().index[0]],
                "clean_max_V_num":[cdr_tem["V segments"].value_counts()[0]],
                "clean_J_types":[len(cdr_clean["J segments"].unique())],
                "clean_max_J_num":[cdr_tem["J segments"].value_counts()[0]],
                "clean_max_Jsegments":[cdr_tem["J segments"].value_counts().index[0]]
                }
                dat=pd.DataFrame(dat)
                dat2.loc[line,:]=dat.iloc[0,:]
            dat2.to_csv(name+num3+".BC_stat.csv")
        elif args.Statindex=="plot":
            folder_this=args.DF
            folder_this_name=args.SV

            os.system("Rscript " +Draw_immu_index  +" "+folder_this+" " + folder_this_name)

        elif verbosity=="rf":
            migecfast(BC2=BC,n=4,TCRtypes="TRA,TRB")

            arr2=open(BC,"r")
            
            for line in arr2: 
                fqname2=line.rstrip().split("\t")[0]
                fqseq=line.rstrip().split("\t")[1]
                print(fqname2)
                print(fqseq)
                for i in range(1,4):
                    num2=str(i)
                    print(i)
                    migecrun2(name=fqname2,num=num2,TCRtypes=TCRtypes)
                    
                    draw_picture(fqname2,num=num2,TCRtypes=TCRtypes)
        endtime = datetime.datetime.now()
        print (endtime - starttime)
        print("END")





            #python PXTCR01_main.py --Statindex plot --DF TRAtest/ --SV TRAtest2

                # print(cdr_tem.info())
                # print(cdr_tem.describe())




    #python PXTCR01_main.py --FQ1 XCC002_BDDP210001798-1A_1.clean.fq --FQ2 XCC002_BDDP210001798-1A_2.clean.fq --BC BCname.txt --Module b
    #just for stat       
    #python PXTCR01_main.py --FQ1 XCC001_BDDP210001797-1A_1.clean.fq --BC BCname.txt --Statwhat a --Module stat
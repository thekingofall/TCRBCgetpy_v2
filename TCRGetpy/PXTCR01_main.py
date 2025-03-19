import re,os,sys
import datetime
import argparse
import pandas as pd
import gzip
import configparser
import glob

# from Bio import SeqIO
# from Bio.SeqIO.QualityIO import FastqGeneralIterator
# from  Atestforarg import  testfordata


### conda activate  TCRRNA_java8
parser = argparse.ArgumentParser(description="usage:<FQ1><FQ2><TCR>")
Filename=parser.add_argument_group('Filname',"Fastq name")
TCRargument=parser.add_argument_group('TCRargument',"Some argument about Migec")
Modulearg=parser.add_argument_group('Module argument',"Some argument about Module")
#Filename.add_argument("--FQ1", type=str, help="Fastq R1",required=True)
Filename.add_argument("--FQ1", type=str, help="Fastq R1")
Filename.add_argument("--FQ2", type=str, help="Fastq R2")
Filename.add_argument("--FQDir", type=str, help="Directory containing paired FASTQ files")
Filename.add_argument("--BC", type=str, help="file with BC name and BC seq")
Filename.add_argument("--DF", type=str, help="fold with want to get their stat")
Filename.add_argument("--SV", type=str, help="fold with want to Save DF stat")
Filename.add_argument("--ini", type=str, help="cogfig.ini")
Filename.add_argument("--FN", type=str, help="filename you want to define",default="Filename")
Filename.add_argument("--OutDir", type=str, help="Output directory for all results", default=".")
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

def migecrun1(name, output_dir=".", num=str(3), TCRtypes="TRA,TRB"):
    """Run MIGEC Assemble step with specified output directory"""
    # Ensure output directory exists
    ensure_dir(output_dir)
    
    # Create file paths
    fq1_path = get_full_path(output_dir, name+".FQ1.fastq")
    fq2_path = get_full_path(output_dir, name+".FQ2.fastq")
    cdr_dir = get_full_path(output_dir, name+".cdrdata")
    
    # Ensure the cdrdata directory exists
    ensure_dir(cdr_dir)
    
    cmd = f"java -jar {migecpath} Assemble -m {num} -q 15 --mask 0:1 --filter-collisions {fq1_path} {fq2_path} {cdr_dir}"
    print(f"Running command: {cmd}")
    os.system(cmd)

def migecrun2(name, output_dir=".", num=str(3), TCRtypes="TRA,TRB"):
    """Run MIGEC CdrBlast step with specified output directory"""
    # Ensure output directory exists
    ensure_dir(output_dir)
    
    # Create file paths
    cdr_dir = get_full_path(output_dir, name+".cdrdata")
    input_file = get_full_path(cdr_dir, name+f".FQ2.t{num}.cf.fastq")
    output_file = get_full_path(cdr_dir, name+f"_t{num}TRA,TRB.cdrblast.txt")
    
    cmd = f"java -jar {migecpath} CdrBlast -a -R {TCRtypes} {input_file} {output_file}"
    print(f"Running command: {cmd}")
    os.system(cmd)

def migecrun3(name, output_dir=".", num=str(3), TCRtypes="TRA,TRB", SPECIES="HomoSapiens"):
    """Run MIGEC CdrBlast step with species specification and output directory"""
    # Ensure output directory exists
    ensure_dir(output_dir)
    
    # Create file paths
    cdr_dir = get_full_path(output_dir, name+".cdrdata")
    input_file = get_full_path(cdr_dir, name+f".FQ2.t{num}.cf.fastq")
    output_file = get_full_path(cdr_dir, f"{SPECIES}_{name}_t{num}TRA,TRB.cdrblast.txt")
    
    cmd = f"java -jar {migecpath} CdrBlast -a -R {TCRtypes} {input_file} {output_file}"
    print(f"Running command: {cmd}")
    os.system(cmd)

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

def ensure_dir(directory):
    """确保目录存在，如果不存在则创建"""
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory

def get_full_path(base_dir, filename):
    """获取完整文件路径"""
    return os.path.join(base_dir, filename)

def trans_data(name, output_dir=".", r1_file=None, r2_file=None):
    """
    Transform data from paired FASTQ files to the format required for downstream analysis
    
    Parameters:
    - name: sample name
    - output_dir: directory where output files will be saved
    - r1_file: R1 FASTQ file path (if None, use args.FQ1)
    - r2_file: R2 FASTQ file path (if None, use args.FQ2)
    """
    # Ensure output directory exists
    ensure_dir(output_dir)
    
    # Create output file paths
    fq1_path = get_full_path(output_dir, name+".FQ1.fastq")
    fq2_path = get_full_path(output_dir, name+".FQ2.fastq")
    
    WR1=open(fq1_path, "w+")
    WR2=open(fq2_path, "w+")
    
    # Use provided file paths if given, otherwise fall back to args
    r1 = r1_file if r1_file is not None else args.FQ1
    r2 = r2_file if r2_file is not None else args.FQ2
    
    # Check if files are gzipped based on extension
    r1_is_gzip = r1.endswith('.gz')
    r2_is_gzip = r2.endswith('.gz')
    
    # Open files with appropriate method based on compression
    if r1_is_gzip:
        fr1 = gzip.open(r1, "r")
    else:
        fr1 = open(r1, "r")
        
    if r2_is_gzip:
        fr2 = gzip.open(r2, "r")
    else:
        fr2 = open(r2, "r")
    
    transformer(FRdata=fr1, FRdata2=fr2, WRdata1=WR1, WRdata2=WR2, same_seq=Sameseq)
    
    # Close the input files
    fr1.close()
    fr2.close()
    
    # Now do the reverse transformation
    WR3=open(fq1_path, "a")
    WR4=open(fq2_path, "a")
    
    # Reopen files for the reverse transformation
    if r2_is_gzip:
        fr2 = gzip.open(r2, "r")
    else:
        fr2 = open(r2, "r")
        
    if r1_is_gzip:
        fr1 = gzip.open(r1, "r")
    else:
        fr1 = open(r1, "r")
    
    transformer(FRdata=fr2, FRdata2=fr1, WRdata1=WR3, WRdata2=WR4, same_seq=Sameseq)
    
    # Close files again
    fr1.close()
    fr2.close()
    
    return fq1_path, fq2_path


def draw_picture(name, output_dir=".", num=str(3), TCRtypes="TRA,TRB"):
    """Run R script for visualization with specified output directory"""
    # Ensure output directory exists
    ensure_dir(output_dir)
    
    # Create file paths
    cdr_dir = get_full_path(output_dir, name+".cdrdata")
    input_file = get_full_path(cdr_dir, name+f"_t{num}TRA,TRB.cdrblast.txt")
    output_prefix = get_full_path(cdr_dir, name+f"_{TCRtypes}_t{num}")
    
    cmd = f"Rscript {Drawpiepath} {input_file} {output_prefix}"
    print(f"Running command: {cmd}")
    os.system(cmd)

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

def find_paired_fastq_files(directory):
    """
    Search for paired FASTQ files in the given directory.
    Returns a list of tuples (R1_file, R2_file, sample_name)
    """
    # Get all fastq/fq files (including .gz compressed files)
    fastq_files = []
    fastq_files.extend(glob.glob(os.path.join(directory, "*.fastq")))
    fastq_files.extend(glob.glob(os.path.join(directory, "*.fq")))
    fastq_files.extend(glob.glob(os.path.join(directory, "*.fastq.gz")))
    fastq_files.extend(glob.glob(os.path.join(directory, "*.fq.gz")))
    
    # Dictionary to store paired files
    paired_files = []
    
    # Common patterns for R1/R2 or 1/2 designations in filenames
    r1_patterns = ["_1.", "_R1.", "_r1."]
    r2_patterns = ["_2.", "_R2.", "_r2."]
    
    # Find R1 files and their corresponding R2 files
    for f in fastq_files:
        base_name = os.path.basename(f)
        
        for r1_pat, r2_pat in zip(r1_patterns, r2_patterns):
            if r1_pat in base_name:
                r2_file = f.replace(r1_pat, r2_pat)
                
                # Check if R2 file exists
                if r2_file in fastq_files:
                    # Get sample name (everything before _1. or _R1.)
                    sample_name = base_name.split(r1_pat)[0]
                    paired_files.append((f, r2_file, sample_name))
                    break
    
    return paired_files

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
        
        # Handle directory input for auto-pairing
        if args.FQDir:
            paired_files = find_paired_fastq_files(args.FQDir)
            print(f"Found {len(paired_files)} paired FASTQ files in {args.FQDir}")
            
            # Use the specified output directory or create a default one
            output_dir = args.OutDir if args.OutDir != "." else os.path.join(os.getcwd(), "TCR_results")
            ensure_dir(output_dir)
            
            for r1_file, r2_file, sample_name in paired_files:
                print(f"Processing paired files: {os.path.basename(r1_file)} and {os.path.basename(r2_file)}")
                
                # Create sample-specific output directory
                sample_output_dir = os.path.join(output_dir, sample_name)
                ensure_dir(sample_output_dir)
                
                R = r1_file
                R2 = r2_file
                TCRtypes = args.TCR
                num = str(args.Num)
                
                name = sample_name
                if args.FN != "Filename":  # If user specified a custom filename
                    name = args.FN
                
                Stattext = os.path.join(sample_output_dir, f"{name}_stat.csv")
                print(f"Sample name: {name}")
                print(f"Output statistics will be saved to: {Stattext}")
                
                verbosity = args.Module
                BC = args.BC
                Sameseq = args.BothSeq
                Statornot = args.Statwhat
                
                print("Now, Start transform")
                
                if verbosity == "a":
                    print("this module whill help you get the BC file,then Split the BC as you need ,then run a migec in three times ")
                    fq1_path, fq2_path = trans_data(name=name, output_dir=sample_output_dir, r1_file=r1_file, r2_file=r2_file)
                    migecrun1(name=name, output_dir=sample_output_dir, num=num, TCRtypes=TCRtypes)
                    migecrun2(name=name, output_dir=sample_output_dir, num=num, TCRtypes=TCRtypes)
                    draw_picture(name, output_dir=sample_output_dir, num=num, TCRtypes=TCRtypes)
                    if BC:
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
                    fq1_path, fq2_path = trans_data(name=name, output_dir=sample_output_dir, r1_file=r1_file, r2_file=r2_file)
                elif verbosity=="fs":
                    print("this module run the migec")
                    print(name)
                    fq1_path, fq2_path = trans_data(name=name, output_dir=sample_output_dir, r1_file=r1_file, r2_file=r2_file)
                    migecrun1(name=name, output_dir=sample_output_dir, num=num, TCRtypes=TCRtypes)
                    migecrun2(name=name, output_dir=sample_output_dir, num=num, TCRtypes=TCRtypes)
                    draw_picture(name, output_dir=sample_output_dir, num=num, TCRtypes=TCRtypes)
                
                elif verbosity=="fs2":
                    print(name)
                    fq1_path, fq2_path = trans_data(name=name, output_dir=sample_output_dir, r1_file=r1_file, r2_file=r2_file)
                    migecrun(name=name,num=num,TCRtypes=TCRtypes)
                    draw_picture(name, output_dir=sample_output_dir, num=num, TCRtypes=TCRtypes)
                elif verbosity=="fsm":
                    print("this module run the migec")
                    print(name)
                    migecrun3(name=name, output_dir=sample_output_dir, num=num, TCRtypes="TRB", SPECIES="MacacaMulatta")
                    draw_picture2(name, num=num, TCRtypes="TRB", SPECIES="MacacaMulatta")
                elif verbosity=="fsma":
                    print("this module run the migec")
                    print(name)
                    fq1_path, fq2_path = trans_data(name=name, output_dir=sample_output_dir, r1_file=r1_file, r2_file=r2_file)
                    migecrun1(name=name, output_dir=sample_output_dir, num=num, TCRtypes=TCRtypes)
                    migecrun3(name=name, output_dir=sample_output_dir, num=num, TCRtypes="TRB", SPECIES="MacacaMulatta")
                    draw_picture2(name, num=num, TCRtypes="TRB", SPECIES="MacacaMulatta")

                elif verbosity=="s":
                    print("this module run the migec")
                    migecrun1(name=name, output_dir=sample_output_dir, num=num, TCRtypes=TCRtypes)
                    migecrun2(name=name, output_dir=sample_output_dir, num=num, TCRtypes=TCRtypes)
                    draw_picture(name, output_dir=sample_output_dir, num=num, TCRtypes=TCRtypes)
                elif verbosity=="t":
                    draw_picture(name, output_dir=sample_output_dir, num=num, TCRtypes=TCRtypes)
                elif verbosity=="o":
                    testfordata(R)
                elif verbosity=="b":
                    print("this module whill help you get the BC file,then Split the BC as you need ,then run a migec")
                    fq1_path, fq2_path = trans_data(name=name, output_dir=sample_output_dir, r1_file=r1_file, r2_file=r2_file)
                    BC_split(BC=BC,fq1=fq1_path,fq2=fq2_path)
                    
                    arr=open(BC,"r")
                    for line in arr: 
                        fqname2=line.rstrip().split("\t")[0]
                        fqseq=line.rstrip().split("\t")[1]
                        print(fqname2)
                        print(fqseq)
                        migecrun1(name=fqname2, output_dir=sample_output_dir, num=num, TCRtypes=TCRtypes)
                        migecrun2(name=fqname2, output_dir=sample_output_dir, num=num, TCRtypes=TCRtypes)
                        draw_picture(fqname2, output_dir=sample_output_dir, num=num, TCRtypes=TCRtypes)
                elif verbosity=="bc":
                    BC_split(BC=BC,fq1=fq1_path,fq2=fq2_path)
                    arr=open(BC,"r")
                    for line in arr: 
                        fqname2=line.rstrip().split("\t")[0]
                        fqseq=line.rstrip().split("\t")[1]
                        print(fqname2)
                        print(fqseq)
                        fq_in_count=count_line(fqname2+".FQ1.fastq")
                        
                        migecrun(name=fqname2,num=num,TCRtypes=TCRtypes)
                        draw_picture(fqname2, output_dir=sample_output_dir, num=num, TCRtypes=TCRtypes)
                elif verbosity=="bm":
                    arr=open(BC,"r")
                    for line in arr: 
                        fqname2=line.rstrip().split("\t")[0]
                        fqseq=line.rstrip().split("\t")[1]
                        print(fqname2)
                        print(fqseq)
                        print(migecpath)
                        
                        migecrun(name=fqname2,num=num,TCRtypes=TCRtypes)
                        draw_picture(fqname2, output_dir=sample_output_dir, num=num, TCRtypes=TCRtypes)
                
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
                            
                            draw_picture(fqname2, output_dir=sample_output_dir, num=num2, TCRtypes=TCRtypes)
                elif args.Statindex=="a":
                    print("this module will help stat the reads of fastq file")
                    
                    os.system("rm -rf "+Stattext)
                    os.system("touch "+Stattext)
                    print(name)
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
                    num3=str(3)
                    num3=num
                    arr=open(BC,"r")
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
            fq1_path, fq2_path = trans_data(name=name, output_dir=".", r1_file=R, r2_file=R2)
            migecrun1(name=name, output_dir=".", num=num, TCRtypes=TCRtypes)
            migecrun2(name=name, output_dir=".", num=num, TCRtypes=TCRtypes)
            draw_picture(name, output_dir=".", num=num, TCRtypes=TCRtypes)
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
                    
                    draw_picture(fqname2, output_dir=".", num=num2, TCRtypes=TCRtypes)
        elif verbosity == "f":
            print("this module get the Reads of the FQ file")
            print(name)
            fq1_path, fq2_path = trans_data(name=name, output_dir=".", r1_file=R, r2_file=R2)
        elif verbosity=="fs":
            print("this module run the migec")
            print(name)
            fq1_path, fq2_path = trans_data(name=name, output_dir=".", r1_file=R, r2_file=R2)
            migecrun1(name=name, output_dir=".", num=num, TCRtypes=TCRtypes)
            migecrun2(name=name, output_dir=".", num=num, TCRtypes=TCRtypes)
            draw_picture(name, output_dir=".", num=num, TCRtypes=TCRtypes)
        
        elif verbosity=="fs2":
            print(name)
            fq1_path, fq2_path = trans_data(name=name, output_dir=".", r1_file=R, r2_file=R2)
            migecrun(name=name,num=num,TCRtypes=TCRtypes)
            draw_picture(name, output_dir=".", num=num, TCRtypes=TCRtypes)
        elif verbosity=="fsm":
            print("this module run the migec")
            print(name)
            migecrun3(name=name, output_dir=".", num=num, TCRtypes="TRB", SPECIES="MacacaMulatta")
            draw_picture2(name, num=num, TCRtypes="TRB", SPECIES="MacacaMulatta")
        elif verbosity=="fsma":
            print("this module run the migec")
            print(name)
            fq1_path, fq2_path = trans_data(name=name, output_dir=".", r1_file=R, r2_file=R2)
            migecrun1(name=name, output_dir=".", num=num, TCRtypes=TCRtypes)
            migecrun3(name=name, output_dir=".", num=num, TCRtypes="TRB", SPECIES="MacacaMulatta")
            draw_picture2(name, num=num, TCRtypes="TRB", SPECIES="MacacaMulatta")

        elif verbosity=="s":
            print("this module run the migec")
            migecrun1(name=name, output_dir=".", num=num, TCRtypes=TCRtypes)
            migecrun2(name=name, output_dir=".", num=num, TCRtypes=TCRtypes)
            draw_picture(name, output_dir=".", num=num, TCRtypes=TCRtypes)
        elif verbosity=="t":
            draw_picture(name, output_dir=".", num=num, TCRtypes=TCRtypes)
        elif verbosity=="o":
            testfordata(R)
        elif verbosity=="b":
            print("this module whill help you get the BC file,then Split the BC as you need ,then run a migec")
            fq1_path, fq2_path = trans_data(name=name, output_dir=".", r1_file=R, r2_file=R2)
            BC_split(BC=BC,fq1=fq1_path,fq2=fq2_path)
            
            arr=open(BC,"r")
            for line in arr: 
                fqname2=line.rstrip().split("\t")[0]
                fqseq=line.rstrip().split("\t")[1]
                print(fqname2)
                print(fqseq)
                migecrun1(name=fqname2, output_dir=".", num=num, TCRtypes=TCRtypes)
                migecrun2(name=fqname2, output_dir=".", num=num, TCRtypes=TCRtypes)
                draw_picture(fqname2, output_dir=".", num=num, TCRtypes=TCRtypes)
        elif verbosity=="bc":
            BC_split(BC=BC,fq1=fq1_path,fq2=fq2_path)
            arr=open(BC,"r")
            for line in arr: 
                fqname2=line.rstrip().split("\t")[0]
                fqseq=line.rstrip().split("\t")[1]
                print(fqname2)
                print(fqseq)
                fq_in_count=count_line(fqname2+".FQ1.fastq")
                
                migecrun(name=fqname2,num=num,TCRtypes=TCRtypes)
                draw_picture(fqname2, output_dir=".", num=num, TCRtypes=TCRtypes)
        elif verbosity=="bm":
            arr=open(BC,"r")
            for line in arr: 
                fqname2=line.rstrip().split("\t")[0]
                fqseq=line.rstrip().split("\t")[1]
                print(fqname2)
                print(fqseq)
                print(migecpath)
                
                migecrun(name=fqname2,num=num,TCRtypes=TCRtypes)
                draw_picture(fqname2, output_dir=".", num=num, TCRtypes=TCRtypes)
        
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
                    
                    draw_picture(fqname2, output_dir=".", num=num2, TCRtypes=TCRtypes)
        elif args.Statindex=="a":
            print("this module will help stat the reads of fastq file")
                    
            os.system("rm -rf "+Stattext)
            os.system("touch "+Stattext)
            print(name)
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
            num3=str(3)
            num3=num
            arr=open(BC,"r")
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
                    migecrun2(name=fqname2, output_dir=".", num=num2, TCRtypes=TCRtypes)
                    
                    draw_picture(fqname2, output_dir=".", num=num2, TCRtypes=TCRtypes)
        endtime = datetime.datetime.now()
        print (endtime - starttime)
        print("END")





            #python PXTCR01_main.py --Statindex plot --DF TRAtest/ --SV TRAtest2

                # print(cdr_tem.info())
                # print(cdr_tem.describe())




    #python PXTCR01_main.py --FQ1 XCC002_BDDP210001798-1A_1.clean.fq --FQ2 XCC002_BDDP210001798-1A_2.clean.fq --BC BCname.txt --Module b
    #just for stat       
    #python PXTCR01_main.py --FQ1 XCC001_BDDP210001797-1A_1.clean.fq --BC BCname.txt --Statwhat a --Module stat

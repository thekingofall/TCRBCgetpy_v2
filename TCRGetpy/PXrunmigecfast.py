import sys,os
BC2=sys.argv[1]
# arr=open(BC,"r")
        
# os.system("rm  -rf run_migec_cpu.txt")
# def migecrunfast(name,num=str(3),TCRtypes="TRA,TRB"):
#     #os.system("rm  -rf run_migec_cpu.txt")
#     migeccpu=open("run_migec_cpu.txt","a+")
#     command="java -jar migec-1.2.9.jar Assemble -m "+num+" -q 15 --mask 0:1 --filter-collisions  "+ name+".FQ1.fastq  " + name+".FQ2.fastq  "+ name+".cdrdata"
#     migeccpu.write(command+"\n")

# for line in arr: 
    
#     fqname2=line.rstrip().split("\t")[0]
#     fqseq=line.rstrip().split("\t")[1]
#     #print(fqname2)
#     #print(fqseq)

#     for i in range(1,4):
#         num2=str(i)
#         print(i)
#         migecrunfast(name=fqname2,num=num2,TCRtypes="TRA,TRB")
def migecrunfast(name,num=str(3),TCRtypes="TRA,TRB"):
    #os.system("rm  -rf run_migec_cpu.txt")
    migeccpu=open("run_migec_cpu.txt","a+")
    command="java -jar migec-1.2.9.jar Assemble -m "+num+" -q 15 --mask 0:1 --filter-collisions  "+ name+".FQ1.fastq  " + name+".FQ2.fastq  "+ name+".cdrdata"
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
migecfast(BC2=BC2)
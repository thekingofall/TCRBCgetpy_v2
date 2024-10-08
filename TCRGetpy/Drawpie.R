library(RColorBrewer)
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2){
    print("Usage:<filename><Give the png a name>")
}else{
    filename=args[1]
    #print(filename)
    pngname=args[2]
    TCR_condition2<-read.table(filename,sep = "\t",header = T)
    TRA=TCR_condition2[which(substring(TCR_condition2$V.segments,1,3)=="TRA"),]
    TRB=TCR_condition2[which(substring(TCR_condition2$V.segments,1,3)=="TRB"),]
    
    
    Runpie<-function(TCR,name){
    count=as.data.frame(table(TCR$Count))
    count$Var1<-as.numeric(as.character(count$Var1))
    count
    count$sum<-count$Var1*count$Freq
    count<-count[order(count$sum,decreasing = T),]
    count$percent<-signif((count$sum/sum(count$sum))*100,2)
    count$TCRname<-paste0("Clone_",count$Var1)
    count$label=paste(count$TCRname, "(", count$percent, "%)")
    write.csv(count,paste0(pngname,"_pure_",name,".csv"))
    png(paste0(pngname,"_",name,".png"),width=600*4,height=3*600,res=72*3)
    pie(count$sum,labels = count$Var1[1:5],col=colorRampPalette(rev(brewer.pal(10,'Spectral')))(34),border="white",angle = 10,main=paste0(name,"\n","  CDR3 total ",sum(count$sum),"\n","  CDR3 types ", nrow(TCR)))
    legend("left",legend =count$label[1:5] , cex=1.0,fill = colorRampPalette(rev(brewer.pal(10,'Spectral')))(34)[1:4],   border="white",bty = "n")
    dev.off()
    
    }

    if (all(substring(TCR_condition2$V.segments,1,3)=="TRB")){

    Runpie(TRB,name="TCRβ")

    }else if(all(substring(TCR_condition2$V.segments,1,3)=="TRA")){
    
    Runpie(TRA,name="TCRα")
    }else{
    Runpie(TRA,name="TCRα")
    Runpie(TRB,name="TCRβ")
    }

}




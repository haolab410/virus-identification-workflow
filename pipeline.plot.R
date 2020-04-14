library(gplots)
library(geneplotter)
library(optparse)

args <- commandArgs(TRUE)
folder = (args[1])
sample = (args[2])
folder = (paste(folder,'/',sep=''))
print(folder)
print('here')
    
cols<-c("orange","green","blue","red3","gray","purple","black")
meta<-read.table(paste(folder,"/cov.summary.txt",sep=""),check.names=F,comment.char = "",header=T)
index<-sort(meta$"Avg_fold",decreasing=T,index.return=T)
meta<-meta[index$ix,]
cov.th<-50
Accs<-subset(meta,Covered_percent>cov.th,select=c("#ID"))[,1]
#paste(folder,sample,"/",sep="")
  
temp<-read.delim(paste(folder,"/FastViromeExplorer-final-sorted-abundance.tsv",sep=""),check.names=F,comment.char = "",header=T)
index<-match(Accs,temp[,1])
Vnames<-temp[index,2]
  
    #######quality piechart#########
abund<-read.delim(paste(folder,"/abundance.tsv",sep=""),check.names=F,comment.char = "",header=T)
reads.virus<-sum(abund$est_counts)
reads.qc<-read.table(paste(folder,"/reads.qc.summary.txt",sep=""),header=T,sep="\t")
reads.hg<-reads.qc$Reads_mapped
reads.LQ<-reads.qc$Read_QC_Removed
reads.Others<-reads.qc$Reads_Total-reads.hg-reads.LQ-reads.virus
read.distri<-c(reads.virus,reads.hg,reads.Others,reads.LQ)
percent<- paste("(",round(100*read.distri/reads.qc$Reads_Total,digit=2),"%)",sep="")
names(read.distri)<-paste(c("Virus","Host","Others","Low Quality"),percent,sep="")  
pdf(file=paste(folder,"/piechart-reads-qc.pdf",sep=""))
par(mar = c(2,5,2,5))
pie(read.distri, labels = names(read.distri), edges = 200, radius = 1,clockwise = FALSE, density = NULL, angle = 45,col = cols, cex=0.7,border = NULL, lty = NULL, main = NULL)
    
#savepdf(paste(folder,"piechart-reads-qc",sep=""),asp=0.7)
dev.off()
  
###barplotting############
Accs.ditri<-meta[,1]
cov.perc<-meta$Covered_percent
index<-match(Accs.ditri,abund$target_id)
mat.virus<-cbind(abund[index,],cov.perc)
index<-sort(mat.virus$"est_counts",decreasing=T,index.return=T)
mat.virus<-mat.virus[index$ix,]
mat.virus$est_counts<-round(log10(mat.virus$est_counts),digit=2)
temp<-mat.virus$est_counts
temp<-temp/max(temp)
mat.virus$est.counts.norm<-temp
  
temp<-mat.virus$cov.perc
temp<-temp/max(temp)
mat.virus$cov.perc.norm<-temp
  
LeftAxisLabs <- pretty(seq(0, max(mat.virus$est_counts), length.out = 10))
RightAxisLabs <- pretty(seq(0, max(mat.virus$cov.perc), length.out = 10))
LeftAxisAt <- LeftAxisLabs/max(mat.virus$est_counts)
RightAxisAt <- RightAxisLabs/max(mat.virus$cov.perc)
LeftAxisLabs <-paste("1e+0",LeftAxisLabs,sep="")
RightAxisLabs<-paste(RightAxisLabs,"%",sep="")
pdf(file=paste(folder,"/barplot-virus.pdf",sep=""))
par(mar = c(5,5,2,5))
with(mat.virus, plot(target_id,est.counts.norm, type='1',ann=T,axes=FALSE,xlab='',bty="1",ylab="Reads of Coverage",ylim=c(0,max( mat.virus$est.counts.norm)+4.5)))
par(new = T)
if(!(is.na(max(mat.virus$cov.perc.norm)))){
    with(mat.virus, plot(target_id,cov.perc.norm, type='n',axes=FALSE,xlab='',ylab=''))
}
axis(side=4, tick=F,labels=F)
axis(side=1,line=T,tick=F,labels=F)
mtext(side = 4, line = 3, 'Percentage of Coverage',las=0,)
cols.bar<-c("blue", "red")
legend("topright", legend=c("Reads of Coverage", "Percentage of Coverage"),density=c(NA,40),angle=45, fill=cols.bar,border=cols.bar,bty="n")
    
par(new = T)  
#print('here')
if(!(is.na(max(mat.virus$cov.perc.norm)))){
    x<-barplot2(t(as.matrix(mat.virus[,  c( "est.counts.norm","cov.perc.norm")])),beside = T, yaxt = "n", names.arg = mat.virus$target_id, density=c(NA,10),angle=45,ylim=c(0, max(c(LeftAxisAt, RightAxisAt))),xaxt="n",col=cols.bar)
    text(cex=1, x=x[1,]-0.2, y=-0.1, labels=mat.virus$target_id, xpd=TRUE, srt=45)
    axis(4, at = RightAxisAt, labels = RightAxisLabs,las=1)
}
print('here2')
axis(2, at = LeftAxisAt,labels=LeftAxisLabs,las=1,line=-0.5) 
      
    
#savepdf(paste(folder,"barplot-virus",sep=""),width=8)
dev.off()
  
    ####coverage plotting####
num.plot<-length(Accs)
print(num.plot)
if(num.plot > 0){
    pdf(file=(paste(folder,"/coverage-virus.pdf",sep="")))
    par(mfrow=c(num.plot,1))
    for(i in 1:length(Accs))
    {
        data<-read.table(paste(folder,Accs[i],".basecov.txt",sep=""),header=F,sep="\t")
        end<-max(data[,2])
        index<-match(data[,2],0:end)
        x<-0:end
        y<-rep(NA,end)
        y[index]<-data[,3]
        par(mar = c(2,5,2,2))
        label<-paste("Cordinates of ",Vnames[i],"(",Accs[i],")",sep="")
        plot(x=x,y=y,xlab="",ylab="Counts",col=cols[i],type="l",lwd=1.5)
        text(cex=0.5, x=median(x), y=max(y)*0.95, labels=label )
    }
   
    #savepdf(paste(folder,"coverage-virus",sep=""))
    dev.off()
} else {
    pdf(file=(paste(folder,"/coverage-virus.pdf",sep="")))
    #par(pty='s')
    plot.new()
    title('No Virus Found with Coverage Larger than 50%')

    dev.off()
}
  
#plotting.cov.genome(sample, folder)

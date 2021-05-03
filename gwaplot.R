     Show Dotfiles Show Owner/Mode
/projectnb/bs859/spr21/users/tpillars/final_project/
##R --vanilla --args filename plottitle plotfile < gwaplot.R > xx.log
##Assumes PLINK output file with CHR BP and P columns; plots all p-values (need to filter prior to running)

plotscan<-function(pvalue,chrom,physloc,plotname,filename,type="bmp",plotchar=20){
  chrom1<-ifelse(chrom%in%"X","23",ifelse(chrom%in%"Y","24",chrom))
  chrom<-as.numeric(chrom1)
  y<-data.frame(chrom,physloc,pvalue)
  y<-y[order(y$chrom,y$physloc),]
  maxlist<-NULL
  nchrom<-length(unique(chrom))
  for(i in 1:nchrom){
    maxlist<-c(maxlist,max(subset(y,chrom==i)$physloc))}
  chromlist<-data.frame(chr=1:nchrom,max=maxlist)
  chromlist$max<-chromlist$max+500000
  chromlist$cumsum<-NULL
  chromlist$cumsum<-cumsum(as.numeric(chromlist$max))
  chromlist$cumsum2[2:nchrom]<-chromlist$cumsum[1:(nchrom-1)]
  chromlist$cumsum2[1]<-0
  chromlist$addon<-chromlist$cumsum2
  chromlist$chrend<-chromlist$max+chromlist$addon
  chromlist$chrbeg<-chromlist$chrend-chromlist$max
  chromlist$midpt<-chromlist$chrbeg+(chromlist$chrend-chromlist$chrbeg)/2
  y$newloc<-NA
  for (i in 1:dim(chromlist)[1]){
    chrbeg1<-chromlist$chrbeg[i]
    y$newloc<-ifelse(y$chrom==chromlist$chr[i],y$physloc+chrbeg1,y$newloc)}

    y<-y[order(y$newloc),]

  palette(c("darkgreen","blue"))

    if (type == "bmp")
{
	bitmap(paste(filename,".bmp",sep=""),height=6.0, width=15)
}
if(type == "png")
{
	bitmap(paste(filename,".png",sep=""),height=6.0, width=15,type="png16m")
}
if(type == "jpeg")
{
	bitmap(paste(filename,".jpeg",sep=""),height=6.0, width=15,type="jpeg")
} 
if(type == "postscript")
{
	bitmap(paste(filename,".ps",sep=""),height=6.0, width=15,type="psrgb")
} 
if(type == "pdf")
{
	bitmap(paste(filename,".pdf",sep=""),height=6.0, width=15,type="pdfwrite")
} 
if(type == "tiff")
{
	bitmap(paste(filename,".tiff",sep=""),height=6.0, width=15,type="tiffg32d",res=600)
} 

    
  par(mgp=c(2.5,.9,0))

  plot(subset(y$newloc,pvalue<.05),-log10(subset(y$pvalue,y$pvalue<.05)),pch=plotchar,xaxt="n",
  xlab="Chromosome",ylab="-log10(p)",las=1,
  cex.lab=2.5,cex.axis=2,cex.main=2,main=plotname,col=as.numeric(subset(y$chrom,y$pvalue<.05)))
  abline(v=chromlist$chrend,lty=2)
  abline(h=-log10(5e-8),lty=3)
  axis(1,at=chromlist$midpt,labels=chromlist$chr,cex.axis=1.5)
  dev.off()
}

##R --vanilla --args filename plottitle plotfile < gwaplot.R > xx.log
##Assumes PLINK output file with CHR BP and P columns; plots all p-values (need to filter prior to running)

args<-commandArgs(trailingOnly=TRUE)
print(args)
infile<-args[1]
plottitle<-args[2]
plotfile<-args[3]
yy<-read.table(infile,header=T,as.is=T)
yy<-subset(yy,CHR>0)
plotscan(yy$P,yy$CHR,yy$BP,plottitle,plotfile,type="png",plotchar=20)


     Show Dotfiles Show Owner/Mode
/projectnb/bs859/spr21/users/tpillars/final_project/
##run using:
##R --vanilla --args filename plottitle test < qqplot.R > xxx.log
## Assumes PLINK output format, with a TEST and P column


args <- commandArgs(trailingOnly = TRUE)
print(args)
myfile <- args[1]
title<-args[2]
test<-args[3]


qplotpval<-function(xx,title=NULL){
  yy<-sort(subset(xx,!is.na(xx)))
  xlambda<-round(median(qchisq(yy,df=1,lower.tail=FALSE),na.rm=TRUE)/0.455,3)
  xlab1<-paste("Expected -log10(p-value)",title,sep="")
  qq <- (-log10(ppoints(length(yy))))
  yygc<-pchisq(qchisq(yy,1,lower.tail=F)/xlambda,1,lower.tail=F)
  plot(qq, -log10(yy),
       ylab="Observed -log10(p-value)",xlab=xlab1, pch=1,cex.lab=1.75,
	cex.axis=1.75,las=1,col="red")
  points(qq,-log10(yygc),col="green",pch=1)
  abline(0,1)
  text(-log10(min(qq))-1,-log10(0.02),"No GC", col="red",cex=1.5)
  text(-log10(min(qq))-1,-log10(0.01),"GC", col="green",cex=1.5)
title(main=substitute(lambda==xlambda),cex.main=2)

}



mydat<-read.table(myfile,as.is=T,header=T,comment.char="")
if(length(grep("TEST",colnames(mydat))>0))mydat<-subset(mydat,TEST==test)
mydat<-subset(mydat,!is.na(P))

bitmap(paste("qq.",title,".jpeg",sep=""),type="jpeg")
qplotpval(mydat$P,title=title)
dev.off()

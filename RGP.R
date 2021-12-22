#UTF-8 encoding
setwd("D:\\乳腺癌课题 - 已有预后基因\\代码上传\\RGP")####please change this path to the "RGP" file!
load(".\\RGPmodel.RData")
library(limma)
library(survival)
library(survivalROC)
library(survminer)
namedata=c("20685","21653","17705","11121","7390","20711","1456"
           ,"31448","4922","METABRIC")# 10 test sets
pvaluesurvival=c() ###the p values in K-M curves
for (ID in namedata) {
  cox_gene_list=read.table("cox_gene_list.txt", header=T, sep="\t", check.names=F)
  cox_pair_list=read.table("cox_pair_list.txt", header=T, sep="\t", check.names=F)
  expdata=read.table(paste0(ID,"data65.txt"), header=T, sep="\t", check.names=F,row.names = 1)
  out=data.frame()
  for(i in 1:nrow(cox_gene_list))
  {
    num=which(row.names(expdata)[]==cox_gene_list[i,1])
    out=rbind(out,expdata[num,])
    
  }
  #expression data is in "out", then we are going to pair these genes
  tcgaPair=data.frame()
  rt=out
  sampleNum=ncol(rt)
  for(i in 1:(nrow(rt)-1)){
    for(j in (i+1):nrow(rt)){
      pair=ifelse(rt[i,]>rt[j,], 1, 0)
      rownames(pair)=paste0(rownames(rt)[i],"|",rownames(rt)[j])
      tcgaPair=rbind(tcgaPair, pair)
    }
  }
  tcgaOut=rbind(ID=colnames(tcgaPair), tcgaPair)
  tcgaOut=tcgaOut[which(rownames(tcgaOut)[]%in%cox_pair_list[,1]),]
  tcgaOut=t(tcgaOut)
  cli=read.table(paste0(ID,"clinical.txt"), header=T, sep="\t", check.names=F) ####input survival data of these sets
  fustat=vector(length = nrow(tcgaOut))
  futime=vector(length = nrow(tcgaOut))
  for(i in 1:nrow(tcgaOut))
  {
    
    
    num=which(cli[,1]==rownames(tcgaOut)[i])
    fustat[i]=cli[num,3]
    futime[i]=cli[num,2]
    
  }
  tcgaOut=cbind(tcgaOut,fustat,futime)
  tcgaOut[,"futime"]=as.numeric(tcgaOut[,"futime"])/365 ####we used "year" in survival time
  tcgaOut[,"fustat"]=as.numeric(tcgaOut[,"fustat"])
  tcgaOut=as.data.frame(tcgaOut)
  for(i in 1:ncol(tcgaOut))
  {
    tcgaOut[,i]=as.numeric(tcgaOut[,i])
  }
  riskScoreTest=predict(multiCox,type="risk",newdata=tcgaOut)      #we obtained the RGP risk sore of this test set
  medianTrainRisk=median(riskScoreTest)
  riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
  riskfile=cbind(id=rownames(cbind(tcgaOut,riskScoreTest,riskTest)),cbind(tcgaOut,riskScore=riskScoreTest,risk=riskTest))
  write.table(riskfile,
              file=paste0(".\\RGP score\\",ID,"riskTest.txt"),
              sep="\t",
              quote=F,
              row.names=F)#save the RGP risk score of this test set in "RGP score" file.
  
  

  #to draw K-M curves using below codes
  rt=riskfile
  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     title=paste0("GSE",ID),
                     pval=pValue,
                     pval.size=6,
                     conf.int=F,
                     legend.title="risk group",
                     legend.labs=c("high","low"),
                     font.legend=12,
                     fontsize=4,
                     xlab="Time(years)",
                     ylab="Overall survival",
                     break.time.by = 2,
                     palette=c("red","blue"),
                     risk.table=TRUE,
                     risk.table.title="",
                     risk.table.height=.25)
  pdf(file=paste0(".\\survival\\",ID,"survivalTest.pdf"),width=5.5,height=5)#we can get K-M curve of this test set in "survival" file.
  print(surPlot)
  dev.off()
  pvaluesurvival=c(pvaluesurvival,pValue)
  
  
  #to draw ROC curves using below codes
  rocCol=rainbow(3)
  aucText=c()
  
  pdf(file=paste0(".\\ROC\\",ID,"multiROC.pdf"),width=6,height=6)
  
  par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
  
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =3, method="KM")
  plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
       xlab="False positive rate", ylab="True positive rate",
       main=paste0("RGP score in", ID),
       lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
  aucText=c(aucText,paste0("3rd year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
  
  
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =5, method="KM")
  lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],
        lwd = 2)
  aucText=c(aucText,paste0("5th year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
  
  
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =7, method="KM")
  lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],
        lwd = 2)
  aucText=c(aucText,paste0("7th year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
  
  
  
  abline(0,1)
  
  
  legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
  dev.off()
  print(ID)

}
outdata=cbind(namedata,pvaluesurvival)
fix(outdata)





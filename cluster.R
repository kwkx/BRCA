#UTF-8 encoding
workDir="C:\\Users\\zhujie\\Downloads\\BRCA-main\\BRCA-main\\cluster"#please change this path to the "cluster" file!
setwd(workDir)
library(ConsensusClusterPlus)  
library(survival)
library(survminer)
expFile="TCGAdata65.txt"            #input TCGA-BRCA 65 gene expression data 

#input expression data
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)

#to do clustering
maxK=9
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             title=workDir,
                             clusterAlg="km",
                             distance="euclidean",
                             seed=123456,
                             plot="png")


#output clustering labels
clusterNum=3        #according to the CDF curves, k=3 is the optimal cluster number.
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("cluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$cluster))
cluster$cluster=letter[match(cluster$cluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="Cluster.txt", sep="\t", quote=F, col.names=F) #output the cluster labels for every patients in TCGA-BRCA cohort.



clusterFile="Cluster.txt"     #cluster labels.
cliFile="TCGAclinical.txt"               #input survival data.
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365
#merge the cluster labels and survival data.
sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])
#calculate the p value in K-M curve
length=length(levels(factor(rt$cluster)))
diff=survdiff(Surv(futime, fustat) ~ cluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ cluster, data = rt)

#draw the survival curve
bioCol=c("#FF0000","#0066FF","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="cluster",
                   legend.labs=levels(factor(rt[,"cluster"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 1,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)
pdf(file="survival.pdf",onefile = FALSE,width=7,height=5.5)#output the survival curve.
print(surPlot)
dev.off()

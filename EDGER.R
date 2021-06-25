library(edgeR)
library(gplots)
library(plyr)
library(ggplot2)
mtcounts = read.csv('GenecountsMT.csv', row.names = 1)
meta_data = read.csv('metamt.csv')
head(mtcounts)

tol = c("R", "NR", "NR", "R", "R", "NR")
tolcol = c("blue", "red", "red", "blue", "blue", "red")

y<-DGEList(counts=mtcounts, group=tol)
design<- model.matrix(~tol)
keep=NULL

keep<-filterByExpr(y,design)


#keep<-rowSums(cpm(y)>1)>=2
table(keep)
y<- y[keep,]
head(y)
y<-calcNormFactors(y)
head(y)

plotMDS(y,col=tolcol)

y<-estimateDisp(y, design)
results <- exactTest(y)
summary(de<-decideTestsDGE(results))
head(mtcounts)

head(results)
dim(results)
resmt<-topTags(results, n=96)

tail(resmt)
write.table(resmt$table, file = "ExactTest_005_MT.txt", sep="\t")



X1<-rownames(y)
write.table(X1, file="ExactTest_Background_MT.txt", sep="\t")

cpmy<- cpm(y, log=TRUE)
DEGnames<-rownames(resmt)

hcpmy <-subset(cpmy, rownames(cpmy) %in% DEGnames)
hcpmy<-as.matrix(hcpmy)
heatmap.2(hcpmy, scale ="none", col= bluered(100), trace="none",
          density.info = "none", margins=c(9,7), keysize = 1, cexRow=0.2)

head(hcpmy)
deglist <- as.data.frame(hcpmy)



#getting cyp91-3 expression
cpmy<- cpm(y, log=TRUE)

DEGnames<-rownames(y)
hi <-subset(cpmy, rownames(cpmy) %in% DEGnames)
hi<-as.matrix(hi)
head(hi)

write.table(hi, file ='fullgenelistwithnormalizedcounts.txt')

cyp1<-filter(hi, X1 =='rna-XM_026445397.1') 
cyp2<-filter(hi, X1 =='rna-XM_392043.7') 
cypp3<-filter(hi, X1 =='rna-XM_006565484.3')

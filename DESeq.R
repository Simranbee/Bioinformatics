#Data analysis 


#load libraries
library(edgeR)
library(DESeq2)
library(limma)
library(Biobase)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)
library(clusterProfiler)

#input data
just.raw.counts = read.csv("GenecountsB.csv", nrows = 28396, row.names = 1)
rawcounts <- as.matrix(just.raw.counts) 

head(just.raw.counts)
dim(just.raw.counts)
tail(just.raw.counts, 10)

#Input metadata
meta.data = read.csv(file="metab.csv", row.names = 1)
head(meta.data)

#colnames(meta.data)[0] = 'Sample'
meta.data$Tolerance = levels = c('R', 'NR')
#create dataset object

count.data.set = DESeqDataSetFromMatrix(countData=just.raw.counts, 
                                        colData=meta.data, design=~Tolerance)


count.data.set

count.data.set.object <- DESeq(count.data.set)
count.data.set.object

# 1) DATA NORMALIZATION
counts(count.data.set)[1:6,1:6]
count.data.set = estimateSizeFactors(count.data.set)
sizeFactors(count.data.set)       
counts(count.data.set, normalized=TRUE)[1:6,1:6]

boxplot(counts(count.data.set, normalize=TRUE))

vst = varianceStabilizingTransformation(count.data.set)
boxplot(assay(vst))
plotPCA(vst,'Patriline')
plotPCA(vst,'Survival') 
#clustering
d = dist(t(assay(vst)))
h =hclust(d)
plot(h)
k = kmeans(d, centers = 2)
k

#data is clustering, but not by survival, p11, p21 have strange expression
#check to see if this is consistent across organs - it is, but only for P11

# 2) ESTIMATE DISPERSIONS
count.data.set = estimateDispersions(count.data.set)
plotDispEsts(count.data.set)

# 3) APPLY STATISTICS
count.data.set = nbinomWaldTest(count.data.set)
result_table_brain = results(count.data.set)
summary(result_table_brain)
dim(result_df_brain) 
#SHORTCUT for 1-estimate sizefactors, estimatedispersions, nbinomwaldtest:
#shortcut <- DESeq(count.data.set)

View(as.data.frame(result_table_brain))
result_df_brain = as.data.frame(result_table_brain)

#plotted cyp9q1-3
plotCounts(count.data.set, gene='rna-NM_001011613.1', intgroup = 'Tolerance')

#5)FILTERING
#removeNA
sum(complete.cases(result_df_brain))

filter_dfbrain = result_df_brain[complete.cases(result_df_brain),]
dim(filter_dfbrain)


#Filter results
#padj < 0.05 
# log2Foldchange > 1 < -1

sum(filter_dfbrain$padj< 0.05)
filter_dfbrain1 = filter_dfbrain[filter_dfbrain$padj< 0.05,]

abs(filter_dfbrain1$log2FoldChange) > 1
filter3 = filter_dfbrain1[abs(filter_dfbrain1$log2FoldChange) > 1,]
dim(filter3)

#Data visualization
View(filter3)

plotMA(result_table_brain)

# 6) volcano plot
ggplot(filter_dfbrain1, aes(x=log2FoldChange, y=-log10(padj)))+ geom_point()+
  geom_vline(xintercept = 1) + geom_vline(xintercept = -1) +
  geom_hline(yintercept =-log10(0.05))

filter_dfbrain$test = filter_dfbrain$padj< 0.05 & abs(filter_dfbrain$log2FoldChange) > 1


ggplot(filter_dfbrain, aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point(aes(colour =test, alpha = 0.3))+
  scale_colour_manual(values = c('black', 'red')) +
  geom_vline(xintercept = 1, colour = 'green', linetype = 2) + 
  geom_vline(xintercept = -1, colour='green', linetype = 2) +
  geom_hline(yintercept =-log10(0.05), colour ='blue', linetype = 2)+
  xlim(-5, 10) +
  theme(legend.position = 'none')
ggplotly(g)


# 5)HEATMAP->join normalized counts to padj and filter

#ADD COLUMN HEADER TO RESULTFILE
vst


write.csv(result_df_brain, 'resultbrain.txt')

res = read_csv('resultbrain.txt')
head(res)
vst_mat =assay(vst)
write.csv(vst_mat, 'normalizedcountsB.txt')

normcounts = read_csv('normalizedcountsB.txt')
head(normcounts)
annotated_df = left_join(normcounts, res)
View(annotated_df)




#FILTERING


Filter1 = annotated_df[complete.cases(annotated_df),]
dim(filter_dfbrain)
filter2 = Filter1[Filter1$padj< 0.05,]
filter_3 = filter2[abs(filter2$log2FoldChange) > 1,]
dim(filter_3)
head(filter_3)
heat_map <- filter_3[c('X1', 'P11B', 'P1B', 'P21B', 'P7B', 'P8B')]
heat_map
heatmap(heat_map)

#HEATMAP

#GO ANNOTATION
#A- EXPORT GENES
View(filter3)
write.csv(filter3, 'brainFC')
#went to brainFC file and added in Gene as column name


FC = read.csv("brainFC")
FC
DEGene = FC['X']
DEGene
write.csv(DEGene, 'DEGenesBrain.txt', row.names  =FALSE)

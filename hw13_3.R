counts <- read.csv('https://raw.githubusercontent.com/s-a-nersisyan/HSE_bioinformatics_2021/master/seminar13/colon_cancer_tumor_vs_normal_unpaired_counts.tsv', sep = '\t', head = TRUE, )
rownames(counts) <- counts$X
counts$X <- NULL
head(counts)

coldata<- data.frame(x = colnames(counts), tissue = c('tumor','tumor','tumor','tumor','tumor','normal','normal','normal','normal','normal'))
coldata$tissue <- factor(coldata$tissue)

library('DESeq2')
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ tissue)
dds <- DESeq(dds)

resultsNames(dds)
res <- results(dds, name = 'tissue_tumor_vs_normal')
res <- lfcShrink(dds, coef= 'tissue_tumor_vs_normal', type="apeglm")

resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), 
          file="condition_treated_results.csv")
library(BiocManager)
dataDF = read.table('250205_Young_Old_TPM.txt', sep = '\t', header = TRUE, row.names = 1, check.names = F)
metadata = read.table('metadata.txt', sep = '\t', quote = '', stringsAsFactors = F, header = T)
dataDF = dataDF[, metadata$Sample.Name]
dataDF = dataDF[metadata$Sample.Name]

dataDF = data.matrix(dataDF)
dataDF = dataDF[rowSums(dataDF >= 5)>= 9, ]

library(limma)
library(TCC)
tcc = new('TCC', dataDF, group = as.numeric(factor(metadata$Group)))
class(tcc)
show(tcc)
tcc = calcNormFactors(tcc, norm.method = 'tmm', test.method = 'edger', FDR = 0.05)
normDF = getNormalizedData(tcc)
head(normDF)

library(reshape2)
dataDF_melt = melt(dataDF)

colnames(dataDF_melt) = c('gene','sample','count')
dataDF_melt = merge(dataDF_melt, metadata[c('Sample.Name', 'Group')], by.x = 'sample', by.y = 'Sample.Name')

normDF_melt = melt(normDF)
colnames(normDF_melt) = c('gene','sample','count')
normDF_melt = merge(normDF_melt, metadata[c('Sample.Name', 'Group')], by.x = 'sample', by.y = 'Sample.Name')

library(ggplot2)
library(colorRamps)
#Un-normalized 
ggplot(dataDF_melt)+ 
  geom_boxplot(aes(x=sample, y=log2(count+1), fill=Group))+ 
  labs(title = 'counts before normalization') + theme_bw()+ 
  theme(axis.text.x = element_text(angle=90))+
  scale_fill_manual(values = c("O1D" = "#ffd7d6",
                               "O4D" = "#FFA9BA",
                               "O7D" = "#ee7ead",
                               "Y1D" = "#d5edf8",
                               "Y4D" = "#b5bef5",
                               "Y7D" = "#7d8be0"),)
ggsave('Fig1.Boxplot-1.png', width = 7, height = 6)
#Normalized
ggplot(normDF_melt)+ 
  geom_boxplot(aes(x=sample, y=log2(count+1), fill=Group))+ 
  labs(title = 'counts after normalization') + theme_bw()+ 
  theme(axis.text.x = element_text(angle=90))+
  scale_fill_manual(values = c("O1D" = "#ffd7d6",
                               "O4D" = "#FFA9BA",
                               "O7D" = "#ee7ead",
                               "Y1D" = "#d5edf8",
                               "Y4D" = "#b5bef5",
                               "Y7D" = "#7d8be0"),)
ggsave('Fig2.Boxplot-2.png', width = 7, height = 6)

library(ggrepel)
log2DF = log2(normDF + 1)
PCA = prcomp(t(log2DF))
View(PCA)
names(PCA)

#방법1
variances <- PCA$sdev^2
total_variance <- sum(variances)
proportion_variance <- variances / total_variance
cumulative_variance <- cumsum(proportion_variance)
pc1 = round(proportion_variance[1]*100, 2)
pc2 = round(proportion_variance[2]*100, 2)
head(PCA$rotation)
head(sort(abs(PCA$rotation[, "PC1"]), decreasing = TRUE))
pcaDF = merge(PCA$x, metadata[c('Sample.Name', 'Group')], by.x=0, by.y='Sample.Name')
head(pcaDF)
ggplot(pcaDF) + 
  geom_point(aes(x = PC1, y = PC2, color = Group), size = 3)+
  geom_text_repel(aes(x= PC1, y = PC2, label = Row.names), size = 3)+
  labs(x = paste('PC1 (', 
                 round(PCA$sdev[1]/sum(PCA$sdev)*100, 2),
                 '%)', sep = ''),
       y = paste('PC2 (', 
                 round(PCA$sdev[2]/sum(PCA$sdev)*100, 2),
                 '%)', sep = ''))+
  theme_bw()

#방법2
pcaDF = merge(PCA$x, metadata[c('Sample.Name', 'Group')],
              by.x = 0, by.y = 'Sample.Name')
ggplot(pcaDF) +
  geom_point(aes(x = PC1, y = PC2, color = Group), size = 3) +
  geom_hline(yintercept = 0, color = 'red',
             linetype = 'dashed')+
  geom_vline(xintercept = 0, color = 'red',
             linetype = 'dashed')+
  geom_text_repel(aes(x = PC1, y = PC2, label = Row.names), size = 3)+
  labs(x = paste('PC1 (',
                 round(PCA$sdev[1]/sum(PCA$sdev)*100, 2),
                 '%)', sep = ''),
       y = paste('PC2 (',
                 round(PCA$sdev[2]/sum(PCA$sdev)*100, 2),
                 '%)', sep = ''))+
  theme_bw()
ggsave('Fig3.PCAplot.png', width = 6, height = 4)

ggplot(subset(pcaDF, Group %in% c("O1D", "O4D", "O7D")))+
  geom_point(aes(x = PC1, y = PC2, color = Group), size = 3)+
  geom_hline(yintercept = -10, color = 'red',
             linetype = 'dashed')+
  geom_vline(xintercept = 0, color = 'red',
             linetype = 'dashed')+
  geom_text_repel(aes(x = PC1, y = PC2, label = Row.names), size = 3)+
  labs(x = paste('PC1 (',
                 round(PCA$sdev[1]/sum(PCA$sdev)*100, 2),
                 '%)', sep = ''),
       y = paste('PC2 (',
                 round(PCA$sdev[2]/sum(PCA$sdev)*100, 2),
                 '%)', sep = ''))+
  theme_bw()
ggsave('Fig3.PCAplot.old.png', width = 6, height = 4)

ggplot(subset(pcaDF, Group %in% c("Y1D", "Y4D", "Y7D")))+
  geom_point(aes(x = PC1, y = PC2, color = Group), size = 3)+
  geom_hline(yintercept = 15, color = 'red',
             linetype = 'dashed')+
  geom_vline(xintercept = 0, color = 'red',
             linetype = 'dashed')+
  geom_text_repel(aes(x = PC1, y = PC2, label = Row.names), size = 3)+
  labs(x = paste('PC1 (',
                 round(PCA$sdev[1]/sum(PCA$sdev)*100, 2),
                 '%)', sep = ''),
       y = paste('PC2 (',
                 round(PCA$sdev[2]/sum(PCA$sdev)*100, 2),
                 '%)', sep = ''))+
  theme_bw()
ggsave('Fig3.PCAplot.young.png', width = 6, height = 4)

library(ggdendro)
hc = hclust(dist(t(log2DF)), method = 'average')
ggdendrogram(hc, rotate = TRUE)
ggsave('Fig4.dendrogram.png', width = 4, height = 4)

library(factoextra)
fviz_dend(hc, k = 4, rect = TRUE, show_labels = TRUE)

library(pheatmap)
library(colorRamps)
pheatmap(cor(log2DF))
pheatmap(cor(log2DF), 
         main = 'expression correlation', 
         color = 
           colorRampPalette(c('blue', 'white', 'red'))(100),
         treeheight_row = 5, 
         treeheight_col = 5,
         #cutree_rows = 3, 
         #cutree_cols = 4,
         filename = 'Fig5.correlation.png',
         width = 5,
         height = 4)

conTCC = new('TCC', dataDF, 
              group = as.numeric(factor(metadata$Condition)))
conTCC = calcNormFactors(conTCC, norm.method = 'tmm', 
                          test.method = 'deseq2', FDR = 0.05)
conTCC = estimateDE(conTCC, test.method = 'deseq2', FDR = 0.05)
conTCC = getResult(conTCC, sort = T )
conDEG = conTCC[which(abs(conTCC$m.value) >= 1 &
                        conTCC$q.value < 0.05),]
log2DF_conDEG = log2DF[rownames(log2DF) %in% conDEG$gene_id,]
log2DF_conDEG_scale = t(scale(t(log2DF_conDEG)))
pheatmap(log2DF_conDEG_scale, 
         main = 'DEGs by condition',
         color = 
           colorRampPalette(c('blue', 'white', 'red'))(100),
         show_rownames = FALSE,
         cluster_cols = FALSE, 
         clustering_method = 'ward.D', 
         treeheight_row = 10,
         cutree_rows = 4,
         cutree_cols = 6,
         gaps_col = c(3,6,9),
         filename = 'Fig6.DEG_contype.png',
         width = 5,
         height = 5)

tcc <- calcNormFactors(conTCC, norm.method = "tmm", test.method = "edgeR")  # 정규화 수행
tcc <- estimateDE(tcc, test.method = "edgeR", FDR = 0.05)  # DEG 분석 수행
deg_result <- getResult(tcc, sort = TRUE)  # P-value 기준 정렬된 결과
deg_result$max_group_mean <- apply(tcc$count, 1, max)  # 각 유전자의 최대 평균 발현값 계산
deg_result$fold_change <- 2^deg_result$m.value # Fold change 계산 (2^Log2FC)
final_result <- deg_result[, c("gene_id", "max_group_mean", "m.value", "fold_change", "p.value", "q.value")] # 필요한 컬럼만 선택하여 저장
head(final_result)
write.xlsx(final_result, file = "DEG_results.xlsx", rowNames = FALSE)

library(clusterProfiler)
CD4_up = as.numeric(cellTCC$gene_id[cellTCC$m.value <= -2 & cellTCC$q.value <= 0.05])
CD4_enrich = enrichGO(CD4_up, 'org.Mm.eg.db', ont = 'BP', keyType = 'ENTREZID', minGSSize = 5, maxGSSize = 200)
CD4_enrich = CD4_enrich@result
CD4_enrich = CD4_enrich[CD4_enrich$qvalue <= 0.05,]
CD4_enrich = CD4_enrich[order(CD4_enrich$Count, decreasing = c(T)),]
CD4_enrich$ID = factor(CD4_enrich$ID, levels = CD4_enrich$ID)
head(CD4_enrich)
ggplot(CD4_enrich[1:20,])+
  geom_point(aes(x = GeneRatio, y = Description, size = Count,
                 color = qvalue))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_x_discrete(labels = function(x) str_extract(x, "^[0-9]+"))
ggsave('Fig8.CD4-GOenrich.png', width = 8, height = 8)
write.table(CD4_GSEA, 'Table1.CD4-GOenrich.txt', sep = '\t',
            quote = FALSE, col.names = TRUE)

##Extract gene list from GO 
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
go_terms <- c("GO:0002218", "GO:0090716", "GO:0090717", "GO:0002460", "GO:0046631", "GO:0140367", "GO:0001788", "GO:0061760", "GO:0019882", 
              "GO:0050851", "GO:0002404", "GO:0019730", "GO:0140374", "GO:0042113", "GO:0031296", "GO:0030183", "GO:1990117", "GO:0002339",
              "GO:0002263", "GO:0045123", "GO:0006957", "GO:0006958", "GO:0006956", "GO:0001867", "GO:0038178", "GO:0002430", "GO:0050904",
              "GO:0038096", "GO:0002431", "GO:0046629", "GO:0002467", "GO:0002432", "GO:0042386", "GO:0035172", "GO:0002455", "GO:0002434", 
              "GO:0002387", "GO:1902615", "GO:0002433", "GO:0002418", "GO:0002520", "GO:0090713", "GO:0009682", "GO:0002437", "GO:0002220", 
              "GO:0002758", "GO:0002227", "GO:0002366", "GO:0002269", "GO:0050902", "GO:0030595", "GO:0043299", "GO:0001776", "GO:0001909", 
              "GO:0002443", "GO:0002522", "GO:0002523", "GO:0002285", "GO:0031294", "GO:0046651", "GO:0002319", "GO:0035709", "GO:0043379", 
              "GO:0071674", "GO:0002262", "GO:0002274", "GO:0097529", "GO:0030101", "GO:0001779", "GO:0042267", "GO:0002423", "GO:0002423",
              "GO:0002228", "GO:0070236", "GO:0002578", "GO:0050858", "GO:0039532", "GO:0002635", "GO:0045611", "GO:1903707", "GO:0002698",
              "GO:0050777", "GO:0045829", "GO:0002695", "GO:0002686", "GO:0033026", "GO:0060266", "GO:2000524", "GO:0002644", "GO:0034136", 
              "GO:0034144", "GO:0034148", "GO:0034122", "GO:0008228", "GO:0002251", "GO:0002451", "GO:0002458", "GO:0002465", "GO:0002270", 
              "GO:0090720", "GO:0002440", "GO:0002679", "GO:0002536", "GO:0002200", "GO:0002223", "GO:0002286", "GO:0031295", "GO:0030217", 
              "GO:0001913", "GO:0002424", "GO:0042098", "GO:0045058", "GO:0002461", "GO:0002507", "GO:0042092")
combined_genes <- c()
for (go_term in go_terms) {
  genes <- genes_list[[go_term]]
  combined_genes <- c(combined_genes, genes)
}
unique_genes <- unique(combined_genes)
genes_df <- data.frame(Gene = unique_genes)
write.csv(genes_df, "genes_list_combined_unique.csv", row.names = FALSE)

deg_data <- read_excel("DEG_results.xlsx")
genes_list <- read_csv("genes_list_combined_unique.csv")
deg_genes <- deg_data$gene_id
genes_list_unique <- genes_list$Gene
common_deg_data <- deg_data %>% filter(gene_id %in% genes_list_unique)
head(common_deg_data)
write.xlsx(common_deg_data, "common_genes_DEG_results.xlsx")

getwd()
save.image('250205.RData')
load('250205.RData')

# 参考： https://www.jianshu.com/p/8aa995149744
# 参考： https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# 参考：https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

# setwd("F:/projects_on_940xa/2022062200_WN_hESC_LM_and_ChIP_Seq/02_WN_H9_BPK25_LM_Seq/4.1_DEGs_use_featureCounts_and_DESeq2_exon_wtihO/")

# library(DESeq2)
# library(ggplot2)
# library(ggrepel)
# library(pheatmap)

# exp_data <- read.table("0_WN_hESC_H9_BPK25_RNA_Seq_samples.featureCounts.s1.cut.forDESeq2.txt", sep='\t', row.names="Geneid", header=TRUE, check.names = FALSE)
# head(exp_data)

# # H9_0day_BPK25_rep1  H9_0day_BPK25_rep2  H9_0day_BPK25_rep3  H9_0day_DMSO_rep1 H9_0day_DMSO_rep2 H9_0day_DMSO_rep3 H9_4day_BPK25_rep1  H9_4day_BPK25_rep2  H9_4day_BPK25_rep3  H9_4day_DMSO_rep1 H9_4day_DMSO_rep2 H9_4day_DMSO_rep3 H9_9day_BPK25_rep1  H9_9day_BPK25_rep2  H9_9day_BPK25_rep3  H9_9day_DMSO_rep1 H9_9day_DMSO_rep2
# count_matrix <- as.matrix(exp_data[,c("H9_0day_BPK25_rep1","H9_0day_BPK25_rep2","H9_0day_BPK25_rep3","H9_0day_DMSO_rep1","H9_0day_DMSO_rep2","H9_0day_DMSO_rep3","H9_4day_BPK25_rep1","H9_4day_BPK25_rep2","H9_4day_BPK25_rep3","H9_4day_DMSO_rep1","H9_4day_DMSO_rep2","H9_4day_DMSO_rep3","H9_9day_BPK25_rep1","H9_9day_BPK25_rep2","H9_9day_BPK25_rep3","H9_9day_DMSO_rep1","H9_9day_DMSO_rep2")])
# head(count_matrix)

# phenodata_for_all <- read.table("1_WN_hESC_H9_BPK25_RNA_Seq_samples.phenodata.txt", sep='\t', header=TRUE, check.names = FALSE)
# phenodata_for_all$groups = factor(phenodata_for_all$groups, levels = c("d0_BPK25","d0_DMSO","d4_BPK25","d4_DMSO","d9_BPK25","d9_DMSO"))
# phenodata_for_all

# dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = phenodata_for_all, design = ~ groups)
# # 查看过滤前基因数量
# nrow(dds)
# # 57820


# ########################################################################################################
# # 使用不能全部样本中为 0 的条件过滤基因
# dds_filter1 <- dds[ rowSums(counts(dds)) > 1, ]
# # 查看过滤后基因数量
# nrow(dds_filter1)
# # 29788

# dds_filter1 <- dds[ rowSums(counts(dds)) > 3, ]
# # 查看过滤后基因数量
# nrow(dds_filter1)
# # 27142

# # 要求全部样本中的reads总数大于18条（至少要求：一组比较的两类样本的三个重复中每个重复至少3条reads，假定其他样本中不表达也得2*3*3）
# dds_filter1 <- dds[ rowSums(counts(dds)) > 18, ]
# # 查看过滤后基因数量
# nrow(dds_filter1)
# # 22096

# # 
# # 三合一，reads count normalization and find DEGs
# dds_filter1_deseq <- DESeq(dds_filter1)
# # estimating size factors
# # estimating dispersions
# # gene-wise dispersion estimates
# # mean-dispersion relationship
# # final dispersion estimates
# # fitting model and testing

# ## R 提示建议使用 字母、数字、下划线和点作为列名和变量名，这样是安全的。

# # 输出 normalized count 矩阵
# dds_filter1_normalized_counts <- as.data.frame(counts(dds_filter1_deseq, normalized=TRUE))
# head(dds_filter1_normalized_counts)

# # 将 normalized count 合并追加到原始count结果后（不含在全部样本中reads count为0的那些基因）
# # 参照：https://blog.csdn.net/weixin_30751947/article/details/95108476, 
# # merge(x,y, by=0,all.y=TRUE)指将x和y按照行名合并，由by指定按照哪列合并，0代表行名，以y中的数据为主，即y要全部输出，如果x中没有对应行名的数据则用NA代替，如果是all.x=TRUE，则以x为主
# normalized_count_df <- merge(exp_data, dds_filter1_normalized_counts, by=0, all.y=TRUE)
# nrow(normalized_count_df)
# # 22096
# head(normalized_count_df)

# #对列名进行修改
# colnames(normalized_count_df) <- c("gene_id","gene_name","gene_type","H9_0day_BPK25_rep1","H9_0day_BPK25_rep2","H9_0day_BPK25_rep3","H9_0day_DMSO_rep1","H9_0day_DMSO_rep2","H9_0day_DMSO_rep3","H9_4day_BPK25_rep1","H9_4day_BPK25_rep2","H9_4day_BPK25_rep3","H9_4day_DMSO_rep1","H9_4day_DMSO_rep2","H9_4day_DMSO_rep3","H9_9day_BPK25_rep1","H9_9day_BPK25_rep2","H9_9day_BPK25_rep3","H9_9day_DMSO_rep1","H9_9day_DMSO_rep2","H9_0day_BPK25_rep1.normalized","H9_0day_BPK25_rep2.normalized","H9_0day_BPK25_rep3.normalized","H9_0day_DMSO_rep1.normalized","H9_0day_DMSO_rep2.normalized","H9_0day_DMSO_rep3.normalized","H9_4day_BPK25_rep1.normalized","H9_4day_BPK25_rep2.normalized","H9_4day_BPK25_rep3.normalized","H9_4day_DMSO_rep1.normalized","H9_4day_DMSO_rep2.normalized","H9_4day_DMSO_rep3.normalized","H9_9day_BPK25_rep1.normalized","H9_9day_BPK25_rep2.normalized","H9_9day_BPK25_rep3.normalized","H9_9day_DMSO_rep1.normalized","H9_9day_DMSO_rep2.normalized")
# head(normalized_count_df)
# # 输出结果到文件
# write.table(normalized_count_df, "0_WN_hESC_H9_BPK25_RNA_Seq_samples_count_filter18.with_normalized_counts.xls", sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)


# ##############################################################
# # 1.1 绘制基因表达热图
# ##############################################################
# dds_filter1_nc_log2_matrix <- log2(as.matrix(dds_filter1_normalized_counts)+1)
# pdf(file ="WN_hESC_H9_BPK25_RNA_Seq_samples_count_filter18_normalized_counts_log2_heatmap.pdf",width = 10, height = 12)
# pheatmap( dds_filter1_nc_log2_matrix, annotation_names_row =  FALSE, scale='row', show_rownames = F )
# dev.off()

# png(file ="WN_hESC_H9_BPK25_RNA_Seq_samples_count_filter18_normalized_counts_log2_heatmap.png",width = 1200,height = 1000)
# pheatmap( dds_filter1_nc_log2_matrix, annotation_names_row =  FALSE, scale='row', show_rownames = F )
# dev.off()

# ##############################################################
# # 1.2 绘制皮尔森相关性热图
# ##############################################################
# library(corrplot)
# dds_filter1_nc_log2_matrix <- cor(dds_filter1_normalized_counts, method = "pearson" )
# col3 <- colorRampPalette(c("blue", "white", "red"))
# # 使用pheatmap
# pdf(file ="WN_hESC_H9_BPK25_RNA_Seq_samples_count_filter18_normalized_counts_pearson_cor_matrix_heatmap.pdf",width = 10, height = 10 )
# pheatmap( dds_filter1_nc_log2_matrix )
# dev.off()

# png(file ="WN_hESC_H9_BPK25_RNA_Seq_samples_count_filter18_normalized_counts_pearson_cor_matrix_heatmap.png",width = 800, height = 800)
# pheatmap( dds_filter1_nc_log2_matrix )
# dev.off()


# ##############################################################
# # 1.3 PCA 分析及绘图
# ##############################################################
# # PCA 分析分别用 rlog 和 vst 两种方式
# # 聚类 也分别用这两种方式校准的数据
# # ref: https://yangfangs.github.io/2016/04/21/RNAseq-DEseq-analysis/
# # rlog PCA: using rlog transformed data
# rld <- rlog(dds_filter1)
# rlog_pca_result <- plotPCA(rld, intgroup=c("groups", "samples", "repeats"), returnData = TRUE)
# rlog_pca_result$groups = factor(rlog_pca_result$groups, levels = c("d0_BPK25","d0_DMSO","d4_BPK25","d4_DMSO","d9_BPK25","d9_DMSO"))
# write.table(rlog_pca_result, "WN_hESC_H9_BPK25_RNA_Seq_samples_count_filter18_PCA_result_use_rlog.txt", sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
# percentVar <- round(100 * attr(rlog_pca_result, "percentVar"))
# pdf("WN_hESC_H9_BPK25_RNA_Seq_samples_count_filter18_PCA_result_plot_use_rlog.pdf", height = 5, width = 7)
# p <- ggplot(rlog_pca_result, aes(PC1, PC2, color=groups, shape=repeats)) + geom_point(size=4) +
# # p <- ggplot(rlog_pca_result, aes(PC1, PC2, color=groups)) + geom_point(size=4) +
#   # geom_text_repel(aes(PC1, PC2, label = as.character(name)), nudge_y=-0.7) + 
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#   guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))
# p
# dev.off()

# # vsd PCA
# vsd <- vst(dds_filter1, blind = FALSE)
# vsd_pca_result <- plotPCA(vsd, intgroup=c("groups", "samples", "repeats"), returnData = TRUE)
# vsd_pca_result$groups = factor(vsd_pca_result$groups, levels = c("d0_BPK25","d0_DMSO","d4_BPK25","d4_DMSO","d9_BPK25","d9_DMSO"))
# write.table(vsd_pca_result, "WN_hESC_H9_BPK25_RNA_Seq_samples_count_filter18_PCA_result_use_vsd.txt", sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
# percentVar <- round(100 * attr(vsd_pca_result, "percentVar"))
# pdf("WN_hESC_H9_BPK25_RNA_Seq_samples_count_filter18_PCA_result_plot_use_vsd.pdf", height = 5, width = 7)
# p <- ggplot(vsd_pca_result, aes(PC1, PC2, color=groups, shape=repeats)) + geom_point(size=4) + 
# # p <- ggplot(vsd_pca_result, aes(PC1, PC2, color=groups)) + geom_point(size=4) + 
#   # geom_text_repel(aes(PC1, PC2, label = as.character(name)), nudge_y=-0.7) + 
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#   guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))
# p
# dev.off()



######################################################
# Comparison12: d4_BPK25_vs_DMSO
#####################################################
# 绘制 MA plot
library(geneplotter)
ma_for_d4_BPK25_vs_DMSO <- lfcShrink(dds_filter1_deseq, contrast=c("groups", "d4_BPK25", "d4_DMSO"), type="normal", lfcThreshold=1)
pdf("Comparison12_d4_BPK25_vs_DMSO_MA_plot.pdf", width = 7, height = 5)
plotMA(ma_for_d4_BPK25_vs_DMSO, main="MA_plot", ylim=c(-3,3), cex=.8)
abline(h=c(-1,1), col="dodgerblue", lwd=2)
dev.off()

# 差异基因分析
degs_for_d4_BPK25_vs_DMSO <- results(dds_filter1_deseq, contrast=c("groups", "d4_BPK25", "d4_DMSO"))
head(degs_for_d4_BPK25_vs_DMSO)
summary(degs_for_d4_BPK25_vs_DMSO)
# out of 22096 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 6521, 30%
# LFC < 0 (down)     : 5859, 27%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)


# stat by padj filter
table(degs_for_d4_BPK25_vs_DMSO$padj<0.05)
# FALSE  TRUE 
# 10917 11179 
table(degs_for_d4_BPK25_vs_DMSO$padj<0.05 & degs_for_d4_BPK25_vs_DMSO$log2FoldChange >= 0.585)
# FALSE  TRUE 
# 17751  4345 
table(degs_for_d4_BPK25_vs_DMSO$padj<0.05 & degs_for_d4_BPK25_vs_DMSO$log2FoldChange <= -0.585)
# FALSE  TRUE 
# 18548  3548 
table(degs_for_d4_BPK25_vs_DMSO$padj<0.05 & degs_for_d4_BPK25_vs_DMSO$log2FoldChange >= 1)
# FALSE  TRUE 
# 18986  3110 
table(degs_for_d4_BPK25_vs_DMSO$padj<0.05 & degs_for_d4_BPK25_vs_DMSO$log2FoldChange <= -1)
# FALSE  TRUE 
# 20044  2052 
table(degs_for_d4_BPK25_vs_DMSO$padj<0.01 & degs_for_d4_BPK25_vs_DMSO$log2FoldChange >= 1)
# FALSE  TRUE 
# 19352  2744 
table(degs_for_d4_BPK25_vs_DMSO$padj<0.01 & degs_for_d4_BPK25_vs_DMSO$log2FoldChange <= -1)
# FALSE  TRUE 
# 20350  1746 

# stat by p value filter
table(degs_for_d4_BPK25_vs_DMSO$pvalue<0.05)
# FALSE  TRUE 
#  9907 12189 
table(degs_for_d4_BPK25_vs_DMSO$pvalue<0.05 & degs_for_d4_BPK25_vs_DMSO$log2FoldChange >= 0.585)
# FALSE  TRUE 
# 17538  4558 
table(degs_for_d4_BPK25_vs_DMSO$pvalue<0.05 & degs_for_d4_BPK25_vs_DMSO$log2FoldChange <= -0.585)
# FALSE  TRUE 
# 18285  3811 
table(degs_for_d4_BPK25_vs_DMSO$pvalue<0.05 & degs_for_d4_BPK25_vs_DMSO$log2FoldChange >= 1)
# FALSE  TRUE 
# 18817  3279 
table(degs_for_d4_BPK25_vs_DMSO$pvalue<0.05 & degs_for_d4_BPK25_vs_DMSO$log2FoldChange <= -1)
# FALSE  TRUE 
# 19858  2238 
table(degs_for_d4_BPK25_vs_DMSO$pvalue<0.01 & degs_for_d4_BPK25_vs_DMSO$log2FoldChange >= 1)
# FALSE  TRUE 
# 19196  2900 
table(degs_for_d4_BPK25_vs_DMSO$pvalue<0.01 & degs_for_d4_BPK25_vs_DMSO$log2FoldChange <= -1)
# FALSE  TRUE 
# 20229  1867 


# 差异结果追加 经过 DESeq 标准化的 count 表达量
degs_for_d4_BPK25_vs_DMSO_add_nc <- merge(as.data.frame(degs_for_d4_BPK25_vs_DMSO), dds_filter1_normalized_counts[,c("H9_4day_BPK25_rep1","H9_4day_BPK25_rep2","H9_4day_BPK25_rep3","H9_4day_DMSO_rep1","H9_4day_DMSO_rep2","H9_4day_DMSO_rep3")], by="row.names", sort=FALSE)
head(degs_for_d4_BPK25_vs_DMSO_add_nc)
#对列名进行修改
colnames(degs_for_d4_BPK25_vs_DMSO_add_nc)<- c('gene_id',"baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","H9_4day_BPK25_rep1","H9_4day_BPK25_rep2","H9_4day_BPK25_rep3","H9_4day_DMSO_rep1","H9_4day_DMSO_rep2","H9_4day_DMSO_rep3")
# 增加基因信息
rownames(degs_for_d4_BPK25_vs_DMSO_add_nc) <- degs_for_d4_BPK25_vs_DMSO_add_nc$gene_id
head(degs_for_d4_BPK25_vs_DMSO_add_nc)
degs_all_add_geneInfo <- merge(degs_for_d4_BPK25_vs_DMSO_add_nc, exp_data[,c("gene_name","gene_type")], by="row.names", all.x = TRUE, sort = FALSE)
head(degs_all_add_geneInfo)
# 根据 padj 排序并 去掉第一列 Row.names
degs_all_add_geneInfo_sorted <- degs_all_add_geneInfo[order(degs_all_add_geneInfo$padj), -1]
head(degs_all_add_geneInfo_sorted)
dim(degs_all_add_geneInfo_sorted)
# 22096    15
# 输出全部数据
write.table(degs_all_add_geneInfo_sorted, "Comparison12_d4_BPK25_vs_DMSO_all_diff_exp_gene_info.xls", sep = '\t', col.names = TRUE, row.names=FALSE, quote = FALSE,na='')


# 追加 是否显著上下调的标签
res<-degs_for_d4_BPK25_vs_DMSO[which(degs_for_d4_BPK25_vs_DMSO$baseMean >= 0 &  degs_for_d4_BPK25_vs_DMSO$padj >= 0),]
head(res)
# 新增一列，将log2FoldChange >= 1 并且 padj < 0.01 的标注为up，log2FoldChange< = -1 并且 padj < 0.01 标准为down
# padj >= 0.01, 或  log2FoldChange 在 正负 1 之间，就标记为 non， 认为没有变化
res[which(res$log2FoldChange >= 1 & res$padj < 0.01), 'up_down'] <- 'up'
res[which((res$log2FoldChange > -1 & res$log2FoldChange < 1) | res$padj >= 0.01), 'up_down'] <- 'non'
res[which(res$log2FoldChange <= -1 & res$padj < 0.01), 'up_down'] <- 'down'
head(res)
dim(res)
# [1] 22096     7
# 差异结果追加 经过 DESeq 标准化的 count 表达量
res_add_nc <- merge(as.data.frame(res), dds_filter1_normalized_counts[,c("H9_4day_BPK25_rep1","H9_4day_BPK25_rep2","H9_4day_BPK25_rep3","H9_4day_DMSO_rep1","H9_4day_DMSO_rep2","H9_4day_DMSO_rep3")], by="row.names", all.x = TRUE, sort = FALSE)
head(res_add_nc)
dim(res_add_nc)
# [1] 22096    14
#对列名进行修改
colnames(res_add_nc)<- c('gene_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj","up_down","H9_4day_BPK25_rep1","H9_4day_BPK25_rep2","H9_4day_BPK25_rep3","H9_4day_DMSO_rep1","H9_4day_DMSO_rep2","H9_4day_DMSO_rep3")
# 增加基因信息
rownames(res_add_nc) <- res_add_nc$gene_id
head(res_add_nc)
res_add_nc_add_geneInfo <- merge(res_add_nc, exp_data[,c("gene_name","gene_type")], by="row.names", all.x = TRUE, sort = FALSE)
head(res_add_nc_add_geneInfo)
# 根据 padj 排序并 去掉第一列 Row.names
res_add_nc_add_geneInfo_sorted <- res_add_nc_add_geneInfo[order(res_add_nc_add_geneInfo$padj), -1]
head(res_add_nc_add_geneInfo_sorted)
# 保存数据
write.table(res_add_nc_add_geneInfo_sorted, "Comparison12_d4_BPK25_vs_DMSO_filtered_diff_exp_genes.baseMeanLe0_padjLe0.xls", sep = '\t', col.names = TRUE, row.names=FALSE, quote = FALSE,na='')


# 筛选显著差异表达基因，本例采用的是 abs(log2FC) >= 1 and padj <= 0.01
resSig<-res_add_nc_add_geneInfo_sorted[which(res_add_nc_add_geneInfo_sorted$padj<0.01 & abs(res_add_nc_add_geneInfo_sorted$log2FoldChange)>=1),]
head(resSig)
dim(resSig)
# [1]  4490 16
##保存数据
write.table(resSig,"Comparison12_d4_BPK25_vs_DMSO_filtered_diff_exp_genes.FC2_padj0.01.xls",sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE,na='')


######################################################
# Volcano plot for Comparison12: d4_BPK25_vs_DMSO
#####################################################
library(ggplot2)
d4_BPK25_vs_DMSO <- read.csv("Comparison12_d4_BPK25_vs_DMSO_filtered_diff_exp_genes.baseMeanLe0_padjLe0.xls",header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
d4_BPK25_vs_DMSO$up_down <- factor(d4_BPK25_vs_DMSO$up_down, levels = c("up", "non", "down"))

pdf("Comparison12_d4_BPK25_vs_DMSO_volcano_plot.pdf", width = 8, height = 6)
volcano<-ggplot(d4_BPK25_vs_DMSO,aes(x=log2FoldChange,y=-log10(padj))) 
volcano + geom_point(aes(color=up_down), size=1) +
  scale_color_manual(values=c("tomato","grey","royalblue")) +
  # scale_color_manual(values=c("grey","red")) +
  geom_point(aes(x=d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="MBD2",]$log2FoldChange, y=-log10(d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="MBD2",]$padj)), colour ="gold", size = 1.5) + 
  geom_text(aes(x=d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="MBD2",]$log2FoldChange+0.1, y=-log10(d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="MBD2",]$padj)+2, label="MBD2"), size=3, color="orange") +
  geom_point(aes(x=d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="MBD3",]$log2FoldChange, y=-log10(d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="MBD3",]$padj)), colour ="gold", size = 1.5) + 
  geom_text(aes(x=d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="MBD3",]$log2FoldChange+0.1, y=-log10(d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="MBD3",]$padj)+2, label="MBD3"), size=3, color="orange") +
  geom_point(aes(x=d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="CTCF",]$log2FoldChange, y=-log10(d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="CTCF",]$padj)), colour ="gold", size = 1.5) + 
  geom_text(aes(x=d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="CTCF",]$log2FoldChange+0.1, y=-log10(d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="CTCF",]$padj)+2, label="CTCF"), size=3, color="orange") +
  geom_point(aes(x=d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="TET1",]$log2FoldChange, y=-log10(d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="TET1",]$padj)), colour ="gold", size = 1.5) + 
  geom_text(aes(x=d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="TET1",]$log2FoldChange+0.1, y=-log10(d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="TET1",]$padj)+2, label="TET1"), size=3, color="orange") +
  geom_point(aes(x=d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="TET2",]$log2FoldChange, y=-log10(d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="TET2",]$padj)), colour ="gold", size = 1.5) + 
  geom_text(aes(x=d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="TET2",]$log2FoldChange+0.1, y=-log10(d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="TET2",]$padj)+2, label="TET2"), size=3, color="orange") +
  geom_point(aes(x=d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="TET3",]$log2FoldChange, y=-log10(d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="TET3",]$padj)), colour ="gold", size = 1.5) + 
  geom_text(aes(x=d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="TET3",]$log2FoldChange+0.1, y=-log10(d4_BPK25_vs_DMSO[d4_BPK25_vs_DMSO$gene_name=="TET3",]$padj)+2, label="TET3"), size=3, color="orange") +
  scale_x_continuous(limits = c(-8, 10), breaks=seq(-8, 10, 1)) +
  scale_y_continuous(limits = c(0, 300), breaks=seq(0, 300, 30)) +
  geom_hline(yintercept=-log(0.01, 10), colour="#990000", linetype=4) +
  geom_vline(xintercept=c(-1, 1), colour="#990000", linetype=4) +
  labs(x = "log2FC", y = "-log10(padj)", title = "d4_BPK25_vs_DMSO Volcano plot") +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank())
dev.off()



######################################################
# GSEA for  Comparison12: d4_BPK25_vs_DMSO
#####################################################
library(dplyr)
library(DOSE)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
# 重新读取 all_diff_exp_gene_info 或从现有对象选取
# 读入全部差异基因数据
# all_deseq <- read.csv("Comparison12_d4_BPK25_vs_DMSO_all_diff_exp_gene_info.xls", header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
# 提取出gene_id 和 log2FoldChange列
ensemblid_log2fdc <-dplyr::select(degs_all_add_geneInfo_sorted, gene_id, log2FoldChange)
head(ensemblid_log2fdc)
# 可以看到这时gene_id还不是整数形式，所以对这列gsub()取整，命名为 ENSEMBL 列添加到ensemblid_log2fdc中
ensemblid_log2fdc$ENSEMBL <- gsub("\\.\\d*", "", ensemblid_log2fdc[, 1])
head(ensemblid_log2fdc)
# 转换id：将ENSEMBLID转为ENTREZID
entrezid <-bitr(ensemblid_log2fdc$ENSEMBL,
                fromType = "ENSEMBL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db,
                drop = TRUE)

head(entrezid)
temp_df <- merge(ensemblid_log2fdc, entrezid, by="ENSEMBL")
head(temp_df)
#  ID 转换时 20.66% of input gene IDs are fail to map，分别输出相关数据，一方面检测合并是否正确，另一方面检测并统计 ID 转换失败的基因是什么情况
# write.table(ensemblid_log2fdc,"temp1.xls",sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE,na='')
# write.table(entrezid,"temp2.xls",sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE,na='')
# write.table(temp_df,"temp3.xls",sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE,na='')
# ID 转换失败的基因 大多是非编码（假基因、antisence、lincRNA 等）并且 abs(log2foldchang) < 1, padj 大多是显著的


# 提取entrezid 和log2fdc列
entrezid_log2fdc <- dplyr::select(temp_df, ENTREZID, log2FoldChange)
# 用arrange()按log2fdc排序，desc()表示降序
# entrezid_log2fdc_sorted <- arrange(entrezid_log2fdc, desc(log2FoldChange))
# head(entrezid_log2fdc_sorted)
# 采用如下方式排序
# https://github.com/YuLab-SMU/clusterProfiler/issues/214
entrezid_log2fdc_sorted <- entrezid_log2fdc %>% mutate(rank = rank(log2FoldChange,  ties.method = "random")) %>% arrange(desc(rank))
head(entrezid_log2fdc_sorted)
# 构建一下要GSEA分析的数据
genelist <- entrezid_log2fdc_sorted$log2FoldChange
names(genelist) <- entrezid_log2fdc_sorted[,1]
head(genelist)



# GSEA GO BP
gse_GO_BP <- gseGO(genelist,
           OrgDb = org.Hs.eg.db,
           keyType = "ENTREZID",
           nPerm = 20000,
           ont="BP",
           pvalueCutoff = 1)
head(gse_GO_BP@result, 3)
# 使用DOSE 包中的 setReadable 函数可将 富集分析结果对象中的 基因ID 转换成可读的 基因名称
gse_GO_BP_rdb <-  as.data.frame(setReadable(gse_GO_BP, 'org.Hs.eg.db', keyType = 'ENTREZID'))
head(gse_GO_BP_rdb)
# 将 GSEA GO BP 全部结果输出到表格
write.table(gse_GO_BP_rdb,"Comparison12_d4_BPK25_vs_DMSO.all_GSEA_GO_BP_reults.xls",sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE,na='')

pdf("Comparison12_d4_BPK25_vs_DMSO.all_GSEA_GO_BP_top1_plot.pdf", width = 12, height = 9)
gseaplot2(gse_GO_BP, geneSetID = 1, title = gse_GO_BP$Description[1], pvalue_table = TRUE)
dev.off()
pdf("Comparison12_d4_BPK25_vs_DMSO.all_GSEA_GO_BP_top5_plot.pdf", width = 12, height = 9)
gseaplot2(gse_GO_BP, geneSetID = 1:5, pvalue_table = TRUE)
dev.off()


# GSEA GO CC
gse_GO_CC <- gseGO(genelist,
           OrgDb = org.Hs.eg.db,
           keyType = "ENTREZID",
           nPerm = 20000,
           ont="CC",
           pvalueCutoff = 1)
head(gse_GO_CC@result, 3)
# 使用DOSE 包中的 setReadable 函数可将 富集分析结果对象中的 基因ID 转换成可读的 基因名称
gse_GO_CC_rdb <-  as.data.frame(setReadable(gse_GO_CC, 'org.Hs.eg.db', keyType = 'ENTREZID'))
head(gse_GO_CC_rdb)
# 将 GSEA GO CC 全部结果输出到表格
write.table(gse_GO_CC_rdb,"Comparison12_d4_BPK25_vs_DMSO.all_GSEA_GO_CC_reults.xls",sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE,na='')

pdf("Comparison12_d4_BPK25_vs_DMSO.all_GSEA_GO_CC_top1_plot.pdf", width = 12, height = 9)
gseaplot2(gse_GO_CC, geneSetID = 1, title = gse_GO_CC$Description[1], pvalue_table = TRUE)
dev.off()
pdf("Comparison12_d4_BPK25_vs_DMSO.all_GSEA_GO_CC_top5_plot.pdf", width = 12, height = 9)
gseaplot2(gse_GO_CC, geneSetID = 1:5, pvalue_table = TRUE)
dev.off()


# GSEA GO MF
gse_GO_MF <- gseGO(genelist,
           OrgDb = org.Hs.eg.db,
           keyType = "ENTREZID",
           nPerm = 20000,
           ont="MF",
           pvalueCutoff = 1)
head(gse_GO_MF@result, 3)
# 使用DOSE 包中的 setReadable 函数可将 富集分析结果对象中的 基因ID 转换成可读的 基因名称
gse_GO_MF_rdb <-  as.data.frame(setReadable(gse_GO_MF, 'org.Hs.eg.db', keyType = 'ENTREZID'))
head(gse_GO_MF_rdb)
# 将 GSEA GO MF 全部结果输出到表格
write.table(gse_GO_MF_rdb,"Comparison12_d4_BPK25_vs_DMSO.all_GSEA_GO_MF_reults.xls",sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE,na='')

pdf("Comparison12_d4_BPK25_vs_DMSO.all_GSEA_GO_MF_top1_plot.pdf", width = 12, height = 9)
gseaplot2(gse_GO_MF, geneSetID = 1, title = gse_GO_MF$Description[1], pvalue_table = TRUE)
dev.off()
pdf("Comparison12_d4_BPK25_vs_DMSO.all_GSEA_GO_MF_top5_plot.pdf", width = 12, height = 9)
gseaplot2(gse_GO_MF, geneSetID = 1:5, pvalue_table = TRUE)
dev.off()


# GSEA KEGG
gse_kegg <- gseKEGG(geneList = genelist, organism = "hsa", nPerm = 20000, pvalueCutoff = 1, verbose = T, use_internal_data=T)
head(gse_kegg@result,3)
gse_kegg_rdb <-  as.data.frame(setReadable(gse_kegg, 'org.Hs.eg.db', keyType = 'ENTREZID'))
head(gse_kegg_rdb)
# 将 GSEA KEGG 全部结果输出到表格
write.table(gse_kegg_rdb,"Comparison12_d4_BPK25_vs_DMSO.all_GSEA_KEGG_reults.xls",sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE,na='')

pdf("Comparison12_d4_BPK25_vs_DMSO.all_GSEA_KEGG_top1_plot.pdf", width = 12, height = 9)
gseaplot2(gse_kegg, geneSetID = 1, title = gse_kegg$Description[1], pvalue_table = TRUE)
dev.off()
pdf("Comparison12_d4_BPK25_vs_DMSO.all_GSEA_KEGG_top5_plot.pdf", width = 12, height = 9)
gseaplot2(gse_kegg, geneSetID = 1:5, pvalue_table = TRUE)
dev.off()

# 还可以使用如下方式单独对某个基因集进行富集分析或GSEA分析, gmt 文件从 GSEA 官网下载获得
# c5 <- read.gmt("c5.go.cc.v7.4.entrez.gmt")
# egmt2 <- GSEA(genelist, TERM2GENE=c5, pvalueCutoff = 1, verbose=FALSE)
# head(egmt2)



##################################################################################
# GO KEGG for Comparison12: d4_BPK25_vs_DMSO
##################################################################################
library(clusterProfiler)
library(org.Hs.eg.db)
library(topGO)
library(Rgraphviz)
deg_info <- read.table("Comparison12_d4_BPK25_vs_DMSO_filtered_diff_exp_genes.FC2_padj0.01.xls",  sep = '\t', header = T, stringsAsFactors = F)
head(deg_info)
dim(deg_info)
# 4490   16

# 后续函数需要用ENTREZID，先将 ENSEMBL ID 转换为ENTREZID
# 可以看到这时gene_id还不是整数形式，使用 gsub 获取整数部分
ensembl_id <- gsub("\\.\\d*", "", deg_info$gene_id)
head(ensembl_id)
# 转换id：将ENSEMBLID转为ENTREZID
geneid_map <-bitr(ensembl_id,
                fromType = "ENSEMBL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db,
                drop = TRUE)
# 13.81% of input gene IDs are fail to map.
head(geneid_map)
dim(geneid_map)
# 3894    2

go <- enrichGO(geneid_map$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.10, keyType = 'ENTREZID', readable = TRUE)
write.table(go,"Comparison12_d4_BPK25_vs_DMSO.DEGs_GO_all_enrichment.results.txt", sep="\t", quote=F, row.names=F)
head(go)
dim(go)
# 1007   10
dim(go[go$ONTOLOGY=='BP',])
# 763  10
dim(go[go$ONTOLOGY=='CC',])
# 130  10
dim(go[go$ONTOLOGY=='MF',])
# 114 10

# GO ALL bar plot
pdf("Comparison12_d4_BPK25_vs_DMSO.DEGs_GO_all_enrichment.barplot.pdf", width = 10, height = 6 )
barplot(go, showCategory=15, title = "GO ALL top15")
dev.off()
# GO ALL dot plot
pdf("Comparison12_d4_BPK25_vs_DMSO.DEGs_GO_all_enrichment.dotplot.pdf", width = 10, height = 6 )
dotplot(go, showCategory=15, orderBy = "GeneRatio", title = "GO ALL top15")
dev.off()


GO_BP <- enrichGO(geneid_map$ENTREZID, OrgDb = org.Hs.eg.db, ont='BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.10, keyType = 'ENTREZID', readable = TRUE)
dim(GO_BP)
# GO BP bar plot
pdf("Comparison12_d4_BPK25_vs_DMSO.DEGs_GO_BP_enrichment.barplot.pdf", width = 10, height = 6 )
barplot(GO_BP, showCategory=10, title = "GO BP top10")
dev.off()
# GO BP dot plot
pdf("Comparison12_d4_BPK25_vs_DMSO.DEGs_GO_BP_enrichment.dotplot.pdf", width = 10, height = 6 )
dotplot(GO_BP, showCategory=10, orderBy = "GeneRatio", title = "GO BP top10")
dev.off()
# GO BP network plot
pdf("Comparison12_d4_BPK25_vs_DMSO.DEGs_GO_BP_enrichment.network.pdf", width = 12, height = 12 )
plotGOgraph(GO_BP)
dev.off()


GO_CC <- enrichGO(geneid_map$ENTREZID, OrgDb = org.Hs.eg.db, ont='CC', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.10, keyType = 'ENTREZID', readable = TRUE)
dim(GO_CC)
# GO CC bar plot
pdf("Comparison12_d4_BPK25_vs_DMSO.DEGs_GO_CC_enrichment.barplot.pdf", width = 10, height = 6 )
barplot(GO_CC, showCategory=10, title = "GO CC top10")
dev.off()
# GO CC dot plot
pdf("Comparison12_d4_BPK25_vs_DMSO.DEGs_GO_CC_enrichment.dotplot.pdf", width = 10, height = 6 )
dotplot(GO_CC, showCategory=10, orderBy = "GeneRatio", title = "GO CC top10")
dev.off()
# GO CC network plot
pdf("Comparison12_d4_BPK25_vs_DMSO.DEGs_GO_CC_enrichment.network.pdf", width = 12, height = 12 )
plotGOgraph(GO_CC)
dev.off()


GO_MF <- enrichGO(geneid_map$ENTREZID, OrgDb = org.Hs.eg.db, ont='MF', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.10, keyType = 'ENTREZID', readable = TRUE)
dim(GO_MF)
# GO MF bar plot
pdf("Comparison12_d4_BPK25_vs_DMSO.DEGs_GO_MF_enrichment.barplot.pdf", width = 10, height = 6 )
barplot(GO_MF, showCategory=10, title = "GO MF top10")
dev.off()
# GO MF dot plot
pdf("Comparison12_d4_BPK25_vs_DMSO.DEGs_GO_MF_enrichment.dotplot.pdf", width = 10, height = 6 )
dotplot(GO_MF, showCategory=10, orderBy = "GeneRatio", title = "GO MF top10")
dev.off()
# GO MF network plot
pdf("Comparison12_d4_BPK25_vs_DMSO.DEGs_GO_MF_enrichment.network.pdf", width = 12, height = 12 )
plotGOgraph(GO_MF)
dev.off()


# KEGG 富集分析
kegg <- enrichKEGG(geneid_map$ENTREZID, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05, pAdjustMethod = 'BH', minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.10, use_internal_data = TRUE)
# kegg 没有 readable 参数，只能输出含 ENTREZID 的结果，所以下面要使用 setReadable 将结果转换为 基因名称的
kegg_out <-  as.data.frame(setReadable(kegg, 'org.Hs.eg.db', keyType = 'ENTREZID'))
write.table(kegg_out,"Comparison12_d4_BPK25_vs_DMSO.DEGs_KEGG_enrichme.results.txt", sep="\t", quote=F, row.names=F)
head(kegg_out)
dim(kegg_out)
# 20  9
# kegg bar plot
pdf("Comparison12_d4_BPK25_vs_DMSO.DEGs_KEGG_enrichme.barplot.pdf", width = 10, height = 6 )
barplot(kegg, showCategory=10, title = "KEGG top10")
dev.off()
# kegg dot plot
pdf("Comparison12_d4_BPK25_vs_DMSO.DEGs_KEGG_enrichme.dotplot.pdf", width = 10, height = 6 )
dotplot(kegg, showCategory=10, title = "KEGG top10")
dev.off()


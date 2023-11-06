# No-trait WGCNA
# 无表型数据的WGCNA
# Use WGCNA on a data without any trait data
# 在无表型的数据中使用WGCNA
# Author: Yuntao Zhu
# Beijing Forestry University
# Run on R 4.3
# Blog：https://yuntaobioinformatics.wordpress.com 
#######################################################


###### 0. install and library packages ######
###### 0. 安装和加载WGCNA及其所需依赖 ######
# BiocManager::install("WGCNA")
# 如果其他包缺失，也按照上述方式排查安装
library(WGCNA)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyverse)
library(gplots)
setwd("~/Documents/Rdata") # 更改为你的工作目录
getwd()
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = NULL) # 开启WGCNA多线程，请量力而行，NULL=使用所有在线的线程
# 定义导出函数解决导出表格时列名错位问题
adjusted.output <- function(df, filename) {
  gene_id <- rownames(df) # "gene_id'将是第一列的列名
  df <- cbind(gene_id, df)
  write.table(df, file = filename, sep = "\t", row.names = F, col.names = T)
}



###### 1. 数据读取与预处理 ######
dataExpr = read.table("fpkm.xls", # 这里换成基因表达量表格，FPKM或TPM，格式为tab分隔的xls或者csv或者txt都行
                      stringsAsFactors = F, 
                      header = T, 
                      sep = "\t", 
                      row.names = 1) # 第一列是基因ID，第一行是样本ID
head(dataExpr)
# 加入重要基因防止被清洗
sig.genes <- read.table("sig_genes_TX.xls", header = T, sep = "\t", row.names = 1)
# 检测输入的数据是否OK，请注意阅读Console返回的信息
gsg = goodSamplesGenes(dataExpr, verbose = 3) # 检验所有基因表达量数据是否合格可用于WGCNA分析，需要gsg$allOK返回值为TRUE
gsg$allOK
if (gsg$allOK) {
  print("Input data of all genes expression is OK!")
} else {
  print("Input data of all genes expression is not OK! Removing the offending genes...")
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(dataExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
  print("...Offending genes are removed... Checking again...")
  gsg = goodSamplesGenes(dataExpr, verbose = 3) # 再次检测
  if (gsg$allOK) {
    print("...Input data of all genes expression is OK now! ")
  } else {
    print("...Oops! Input data of all genes expression is still not OK!")
    print("Please check the input file manually!!!")
  }
}
# 单独检测重要基因的输入数据是否OK，请注意阅读Console返回的信息
gsg = goodSamplesGenes(sig.genes, verbose = 3) # 检验数据是否合格可用于WGCNA分析，需要gsg$allOK返回值为TRUE，否则进行去除缺失值步骤
gsg$allOK
if (gsg$allOK) {
  print("Input data of significant gene expression is OK!")
} else {
  print("Input data of significant gene expression is not OK! Removing the offending genes...")
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(sig.genes)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(sig.genes)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  sig.genes = sig.genes[gsg$goodSamples, gsg$goodGenes]
  print("...Offending genes are removed... Checking again...")
  gsg = goodSamplesGenes(sig.genes, verbose = 3) # 再次检测
  if (gsg$allOK) {
    print("...Input data of significant gene expression is OK now! ")
    print("WARNING: But some genes are removed to ensure the program continue to run...")
    print("...You'd better check the input file of significant genes")
  } else {
    print("...Oops! Input data of signifacant gene expression is still not OK!")
    print("Please check the input file manually!!!")
  }
}



###### 2. 数据清洗：筛选方差前x%的基因 ######
# 方法一：一步方差法（备用）
# m.vars = apply(dataExpr, 1, var)
# dataExpr.upper = dataExpr[which(m.vars > quantile(m.vars, probs = seq(0, 1, 0.25))[4]),] # 这里的是百分位数，例如0.25意味着保留前75%方差的基因
# print("Data cleaning...complete")
# reserved.gene.num <-  nrow(dataExpr.upper)
# print(paste(reserved.gene.num, "genes are reserved after data cleaning ..."))
# print("...")
# 方法二：两步法（推荐）
var.dataExpr <- as.vector(apply(as.matrix(dataExpr), 1, var, na.rm = T)) # 按行取基因的方差
nSamples <- ncol(dataExpr)
expr.threshold <- 0.0001 # 表达量筛选阈值，表达量小于等于阈值被视作无效样本，无效样本数若高于一定比例则该基因会被清洗，这是为了去除转录组的冗余噪声（大部分表达量都极低的基因）
nSample.threshold <- 0.5 # 无效样本数最大容忍比例，无效样本占比必须小于等于此比例；例如15个样本，比例0.5，则无效样本数最多为15*0.5向下取整为7个，无效样本数超过7个基因将被清洗
# 上述清洗阈值请根据自己数据实际情况判断，例如如果一半对照中基因几乎不表达，那么无效样本比例应该适当提高
sum.zero <- function(x) {
  sum(x <= expr.threshold)
}
num.missing.samples <- data.frame(samples = apply(dataExpr, 1, sum.zero)) # 统计无效样本数
sum.expr <- as.vector(apply(as.matrix(dataExpr), 1, sum))
keep.genes <- var.dataExpr > 0 & sum.expr >= 0.01 * nSamples & num.missing.samples$samples <= nSample.threshold * nSamples  # 方差大于0，表达量之和>=0.01*样本数且有效样本数超过设定比例的基因被保留
dataclean.result <- as.data.frame(table(keep.genes))
dataclean.result
num.keep.genes <- as.numeric(dataclean.result[which(dataclean.result$keep.genes == "TRUE"),2]) # 统计筛选情况
dataExpr.upper <- data.frame(dataExpr, keep.genes = keep.genes) # 将判定结果与数据结合
dataExpr.upper <- dataExpr.upper[which(dataExpr.upper$keep.genes == TRUE),] # 只保留为True的基因
dataExpr.upper <- dataExpr.upper[,-ncol(dataExpr.upper)]
print("Data cleaning step 1...complete")
print(paste(num.keep.genes, "genes are reserved after data cleaning ...")) # 第一步清洗结束
prop.mad.reserve <- 0.75 # 保留根据MAD法前百分之n的基因，例如0.7表示保留MAD前70%的基因
num.mad.reserve <- as.integer(nrow(dataExpr.upper) * prop.mad.reserve) 
dataExpr.upper <- dataExpr.upper[order(apply(dataExpr.upper, 1, mad), decreasing = T)[1:num.mad.reserve],] # MAD法第二步清洗
num.keep.genes2 <- as.numeric(nrow(dataExpr.upper))
print("Data cleaning step 2...complete")
print(paste(num.keep.genes2, "genes are reserved after data cleaning ...")) # 第二步清洗结束

# debug
# write.table(dataExpr.upper, "dataExpr.upper.xls", sep = "\t")
# debug
list1 <- row.names(sig.genes) # 重要基因
list2 <- row.names(dataExpr.upper) # 清洗后剩下的基因
list3 <- dplyr::setdiff(list1, list2) # 被清洗的重要基因
list4 <- dplyr::intersect(list1, list2) # 未被清洗的重要基因

washed.data <- anti_join(dataExpr, dataExpr.upper)
adjusted.output(washed.data, "washed_data.xls")
washed.gene.num <- nrow(washed.data)
washed.sig.genes <- dataExpr[list3,]
print(paste(washed.gene.num, "genes are washed... Gene list is saved as washed_data.xls"))
if (nrow(washed.sig.genes) > 0) {
  print(paste(nrow(washed.sig.genes), "significant genes are washed by data cleaning ..."))
  print("Here are the washed significant genes.")
  washed.sig.genes
} else {print("No significant gene is washed.")}
print("...")

used.data <- rbind(dataExpr.upper, washed.sig.genes)
adjusted.output(used.data, "used_data.xls")
used.gene.num <- nrow(used.data)
if (nrow(washed.sig.genes) > 0) {
  print(paste(nrow(washed.sig.genes), "washed significant genes are joined in the reserved genes..."))
} else {print("No significant gene is washed... Use reserved genes for WGCNA.")}
print(paste("Now totally", used.gene.num, "genes are finally used for WGCNA..."))
print("This is an union of reserved genes and significant genes...")
print("The used gene list is saved as used_data.xls, you can check your genes in this file.")
print("...")

gsg = goodSamplesGenes(used.data, verbose = 3) # 检验数据是否合格可用于WGCNA分析，需要gsg$allOK返回值为TRUE，否则进行去除缺失值步骤
gsg$allOK



###### 3. 样本聚类 ######
# 聚类前转置数据
dataExpr = as.data.frame(t(used.data))
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
print(paste("The data has", nGenes, "genes with", nSamples, "samples."))
# 样本聚类检查离群值
sampleTree = hclust(dist(dataExpr), method = 'average')
tiff(file="Fig.1.1 Sample cluster tree.tiff",width = 15,height = 15, units = "cm", compression = "lzw", res = 600)
plot1.1 <- plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
print(plot1.1)
dev.off()
save(dataExpr, file = "1.dataExpr.Rdata")
head(dataExpr)[,1:8]
# 取log10(x+1)后聚类
log.fpkm <- function(x) {log10(x+1)}
log.dataExpr <- log.fpkm(dataExpr)
log.sampleTree <- hclust(dist(log.dataExpr), method = "average")
tiff(file="Fig.1.2 log10(fpkm+1)-Sample cluster tree.tiff",width = 15,height = 15, units = "cm", compression = "lzw", res = 600)
plot1.2 <- plot(log.sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
print(plot1.2)
dev.off()



###### 4. 软阈值的选择 ######
# 设定R^2
Rsquare <- 0.9 # 一般要求R^2 =0.8，严格点可以取0.85、0.9
# 设定Power取值范围，我这里取的1-15，16-24间隔为2
# 这是因为一般不取15以上的Power
powers <- c(c(1:15), seq(from = 16, to = 24, by = 2))
sft <- pickSoftThreshold(dataExpr, 
                         powerVector = powers,
                         RsquaredCut = Rsquare)
#Warning message:executing %dopar% sequentially: no parallel backend registered 
sft$powerEstimate # WGCNA程序推荐Power，但是这个经常出问题，仅作为参考
sft$fitIndices # 人工辅助挑选软阈值，一般要求R^2 >= 0.8,或者mean.k < 1000
# Fig.2.1 软阈值
par(mfrow = c(1, 2))
tiff(file="Fig.2.1 Soft Power-1.tiff",width = 15,height = 15, units = "cm", compression = "lzw", res = 600)
plot2.1 <- plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex= 0.9 ,col = "red") # cex可以更改
abline(h = Rsquare, col = "red")
axis(2,Rsquare)
print(plot2.1)
dev.off()
# Fig.2.2 平均连接度
tiff(file="Fig.2.2 Soft Power-2.tiff",width = 15,height = 15, units = "cm", compression = "lzw", res = 600)
plot2.2 <- plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n", xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.9, col = "red")
print(plot2.2)
dev.off()
# 保存sft power挑选信息
write.table(sft$fitIndices, "Table.6 Softpower_indices.xls", sep = "\t", row.names = F)
# 根据图表辅助挑选一个最佳软阈值
power = 8



###### 5.构建共表达网络 ######
# 这里用的一步法构建，没注释的地方别随便改
net = blockwiseModules(dataExpr, 
                       power = power, 
                       maxBlockSize = 30000, # 最大BlockSize根据内存而定，这个其实不是很影响分析结果，16GB RAM 20000没啥问题，别太高把内存爆了
                       # 个人认为，仅个人认为，这个BlockSize在条件允许条件下最好略大于实际用于分析的基因数，这样就只有1个Block
                       # 其次最好应该BlockSize略大于实际分析基因数一半，这样正好划分成两个Block分析
                       # 顺便一提，个人感觉这一步分析所需时间和Blocksize大小成指数级关系
                       TOMType = "unsigned", 
                       minModuleSize = 50, # 最小模块大小，默认是30，可改，一般此值越大模块数会略微减少，否则会出现很多细碎模块
                       reassignThreshold = 0,
                       mergeCutHeight = 0.2, # 切枝高度，影响Module结果，可改，默认0.25，越大切出来模块数量越少
                       numericLabels = TRUE, 
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, 
                       saveTOMFileBase = "TOM", 
                       verbose = 3)
table(net$colors)
# 计算树之间的邻接性(邻接矩阵)
WGCNA_matrix <- dataExpr
adjacency <- adjacency(WGCNA_matrix, power = power)
# 邻接矩阵转换为拓扑重叠矩阵，以降低噪音和假相关，获得距离矩阵
TOM <- TOMsimilarity(adjacency) # 计算树之间的相似性（基因之间的相似性）
dissTOM <- 1 - TOM  # 计算树之间的距离(距离矩阵)
# 根据树之间的不相似度建立树状结构
geneTree <- hclust(as.dist(dissTOM),method = "average")



###### 6.模块可视化 ######
# grey模块的为**未分类**到模块的基因
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
tiff(file="Fig.3 Cluster Dendrogram.tiff", width = 25, height = 15, units = "cm", compression = "lzw", res = 600)
plot3 <- plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                             "Module colors",
                             dendroLabels = FALSE, hang = 0.03,
                             addGuide = TRUE, guideHang = 0.05)
print(plot3)
dev.off()
# 计算模块特征向量Mes
table(moduleColors)
MEs=net$MEs
geneTree =net$dendrograms[[1]]
# 融合相似的modules
MEDissThres <- 0.2 #整合相似度(1-n)以上的Modules，例如此值为0.2意味着整合相似度80%以上的模块;默认0.1
merge <- mergeCloseModules(WGCNA_matrix, moduleColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs #作为每个模块的基因表达趋势的展示
tiff(file = "Fig.3 MergedColors_plotDendroAndColors.tiff", width = 25, height = 15, units="cm", compression="lzw", res = 600)
p1 <- plotDendroAndColors(geneTree,
                          cbind(moduleColors,mergedColors),
                          c("DynamicTree Cut", "Merged dynamic"),
                          dendroLabels = FALSE,
                          hang = 0.03,
                          addGuide = TRUE,
                          guideHang = 0.05)
print(p1)
dev.off()
table(moduleColors)
table(mergedColors)
# 导出gene2module
gene2module <- cbind(used.data, module = moduleColors)
gene2module <- cbind(gene2module, merged.module = mergedColors)
adjusted.output(gene2module, "Table.1 Gene to module.xls")
#
write.table(as.data.frame(table(moduleColors)), "Table.2.1 Overview of moduleColors.xls", sep = "\t", row.names = F)
write.table(as.data.frame(table(mergedColors)), "Table.2.2 Overview of mergedColors.xls", sep = "\t", row.names = F)



###### 7.模块间关系可视化 ######
# 绘制模块之间相关性聚类
MEList <- moduleEigengenes(dataExpr, colors = moduleColors)
MEs <- MEList$eigengenes # 作为每个模块的基因表达趋势的展示
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
tiff(file = "Fig.4.1 METree. moduleColors.tiff",width = 18,height = 10, units = "cm", compression = "lzw", res = 600)
plot4.1 <- plot(METree,
               main = "Clustering of module eigengenes",
               xlab = "", sub = "")
print(plot4.1)
dev.off()
# 绘制merged模块之间相关性聚类
merged.MElist <- moduleEigengenes(dataExpr, colors = mergedColors)
merged.MEs <- merged.MElist$eigengenes
merged.MEDiss <- 1 - cor(merged.MEs)
merged.METree <- hclust(as.dist(merged.MEDiss), method = "average")
tiff(file = "Fig.4.2 mergedMETree. merged moduleColors.tiff",width = 18,height = 10, units = "cm", compression = "lzw", res = 600)
plot4.2 <- plot(merged.METree,
                main = "Clustering of merged module eigengenes",
                xlab = "", sub = "")
print(plot4.2)
dev.off()
# 绘制模块间相关性热图
tiff(file = "Fig.5.1 Colors Eigengeneadjacency heatmap.tiff", 
     width = 16,
     height = 24,
     units = "cm",
     compression = "lzw", 
     res = 600)
plot5.1 <- plotEigengeneNetworks(MEs,
                                 "Modules eigengene adjacencyheatmap", 
                                 marHeatmap = c(2,2,2,2), # marHeatmap 设置下、左、上、右的边距
                                 plotDendrograms = T, 
                                 xLabelsAngle = 90) 
print(plot5.1)
dev.off()
# 绘制merged模块间相关性热图
tiff(file = "Fig.5.2 mergedColors Eigengeneadjacency heatmap.tiff", 
     width = 16,
     height = 24,
     units = "cm",
     compression = "lzw", 
     res = 600)
plot5.2 <- plotEigengeneNetworks(mergedMEs,
                                 "Merged Modules eigengene adjacencyheatmap", 
                                 marHeatmap = c(2,2,2,2), # marHeatmap 设置下、左、上、右的边距
                                 plotDendrograms = T, 
                                 xLabelsAngle = 90) 
print(plot5.2)
dev.off()
# 绘制TOMplot
# 由于绘制此图消耗资源太大，且没必要，因此只绘制merge之后的
# 如果不是必要的话，我建议略过此图（我是认真的，不然运行后你肯定会后悔你为什么要运行）
tiff(file="Fig.6 TOMplot-Network heatmap plot, all genes.tiff",width = 20, height = 15, units = "cm", compression = "lzw", res = 600)
plotTOM <- dissTOM^7 # 为了更显著，用7次方
plot6 <- TOMplot(plotTOM, geneTree, mergedColors, main = "Network heatmap plot, all genes")
print(plot6)
dev.off()



###### 8. 输出模块基因 ######
# 输出模块基因kME
datME = moduleEigengenes(WGCNA_matrix,mergedColors)[[1]]
#color1=rep("pink",dim(WGCNA_matrix)[[2]]) #仅仅获得pink的gene
color1 = as.character(mergedColors)#获得全部merged module对应的基因
datKME = signedKME(WGCNA_matrix, datME)
dataExpr1 = as.data.frame(t(WGCNA_matrix));
datSummary = rownames(dataExpr1)
datout = data.frame(gene_id = datSummary, mergedModule = color1, datKME )
#write.table(datout, "pink_gene_module.xls", sep="\t", row.names=F,quote=F) #仅仅获得pink的gene
write.table(datout, "Table.3 All genes kME.xls", sep = "\t", row.names = F,quote = F)
# 自动导出所有模块信息到cytoscape
df.merged.modules <- as.data.frame(table(mergedColors))
merged.modules <- as.character(df.merged.modules[,1])
probes= colnames(WGCNA_matrix) #获得基因名
cyto.threshold <- 0.15 # 如果导出的表格太大，网络太复杂，可以适当调高此值
if (!dir.exists("./Export2Cyto")) {
  dir.create("./Export2Cyto")
} else print("Output dir already exists...")
for (merged.module in merged.modules) {
  inModule = is.finite(match(mergedColors,merged.module)) #R中的无穷大用Inf表示，要检查一个数是否为无穷，可以使用is.finite()或者is.infinite()函数，当一个数是有限大时，第一个函数将返回TRUE，第二个函数将返回FALSE
  modProbes=probes[inModule] #提取指定模块的基因名
  #Select the corresponding Topological Overlap
  modTOM =TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  cyt = exportNetworkToCytoscape(
    modTOM,
    edgeFile = paste("./Export2Cyto/Table.4.x CytoscapeInput-edges-", merged.module, "-threshold=", cyto.threshold, ".xls", sep=""),
    nodeFile = paste("./Export2Cyto/Table.4.x CytoscapeInput-nodes-", merged.module, "-threshold=", cyto.threshold, ".xls", sep=""),
    weighted = TRUE,
    threshold = cyto.threshold, 
    nodeNames = modProbes, 
    nodeAttr = mergedColors[inModule]
  )
}



###### 9. 其他工具 ######
# 单独导出某一模块信息到cytoscape
# 全自动导出时统一设定threshold可能对于不同大小的模块来说不合适，因此在这里你可以单独导出某一模块并为其单独制定threshold值
modules=c("blue") # 选择自己想研究的module
cyto.threshold <- 0.2 # 在这里更改threshold值，如果导出的表格太大，网络太复杂，可以适当调高此值
probes= colnames(WGCNA_matrix)#获得基因名
inModule = is.finite(match(mergedColors,modules)) #R中的无穷大用Inf表示，要检查一个数是否为无穷，可以使用is.finite()或者is.infinite()函数，当一个数是有限大时，第一个函数将返回TRUE，第二个函数将返回FALSE
modProbes=probes[inModule] #提取指定模块的基因名
#Select the corresponding Topological Overlap
modTOM =TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("Table.4.x CytoscapeInput-edges-", merged.module, "-threshold=", cyto.threshold, ".xls", sep=""),
  nodeFile = paste("Table.4.x CytoscapeInput-nodes-", merged.module, "-threshold=", cyto.threshold, ".xls", sep=""),
  weighted = TRUE,
  threshold = cyto.threshold, # 在这里更改threshold值
  nodeNames = modProbes, 
  nodeAttr = mergedColors[inModule]
)

# 绘制柠檬色的TOMplot
# 由于绘制此图消耗资源太大，且没必要，因此只绘制merge之后的
# 如果不是必要的话，我建议略过此图（我是认真的，不然运行后你肯定会后悔你为什么要运行）
tiff(file="Fig.6 TOMplot-Network heatmap plot, all genes-lemonchiffon.tiff",width = 20, height = 15, units = "cm", compression = "lzw", res = 600)
plotTOM <- dissTOM^7 # 为了更显著，用7次方
plot6.lemon <- TOMplot(plotTOM, geneTree, mergedColors, main = "Network heatmap plot, all genes", col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))
print(plot6.lemon)
dev.off()

# 导出重要基因的分析报告
sig.genes1 <- data.frame(gene_id = as.character(row.names(sig.genes)))
gene2module1 <- data.frame(gene_id = as.character(row.names(gene2module)), module = as.character(gene2module$merged.module))
gene2module.sig <- left_join(sig.genes1, gene2module1, by = c("gene_id" = "gene_id"))
table(gene2module.sig$module)
write.table(gene2module.sig, "Table.5.1 Sig_genes2module.xls", sep = "\t", row.names = F)
write.table(table(gene2module.sig$module), "Table.5.2 Sig_genes2module_report.xls", sep = "\t", row.names = F)

###### fin ######
# 吗的我要写死了

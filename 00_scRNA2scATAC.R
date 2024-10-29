rm(list = ls());gc()

# 加载必要的库
library(SeuratData)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)

# 安装和加载 PBMC multiome 数据集
InstallData("pbmcMultiome")

# 分别加载 RNA 和 ATAC 数据
pbmc.rna <- LoadData("pbmcMultiome", "pbmc.rna")
pbmc.atac <- LoadData("pbmcMultiome", "pbmc.atac")

# RNA 数据集加载到 Seurat 对象后，将其 RNA 数据设置为 Assay5 类型
pbmc.rna[["RNA"]] <- as(pbmc.rna[["RNA"]], Class = "Assay5")

# 去除质量不好的细胞，过滤掉标注为 "filtered" 的细胞
pbmc.rna <- subset(pbmc.rna, seurat_annotations != "filtered")
pbmc.atac <- subset(pbmc.atac, seurat_annotations != "filtered")

# 1. 对 RNA 数据进行标准分析
# 归一化处理、寻找高变异基因、数据缩放、主成分分析 (PCA) 和 UMAP 降维
pbmc.rna <- NormalizeData(pbmc.rna)
pbmc.rna <- FindVariableFeatures(pbmc.rna)
pbmc.rna <- ScaleData(pbmc.rna)
pbmc.rna <- RunPCA(pbmc.rna)
pbmc.rna <- RunUMAP(pbmc.rna, dims = 1:30)

# 2. 对 ATAC 数据进行标准分析
# 加载基因注释信息，并添加到 ATAC 数据中
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(pbmc.atac) <- annotations

# 使用 TF-IDF 标准化处理 ATAC 数据
pbmc.atac <- RunTFIDF(pbmc.atac)

# 查找高特征区域
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = "q0")

# 进行 SVD 降维，随后进行 UMAP 可视化，使用 LSI 作为降维方法
pbmc.atac <- RunSVD(pbmc.atac)
pbmc.atac <- RunUMAP(pbmc.atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# 3. 分别可视化 RNA 和 ATAC 的 UMAP 结果
p1 <- DimPlot(pbmc.rna, group.by = "seurat_annotations", label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(pbmc.atac, group.by = "orig.ident", label = FALSE) + NoLegend() + ggtitle("ATAC")
p1 + p2

# 4. 计算 ATAC 数据中的基因活性
# 使用 scATAC-seq 数据量化每个基因的转录活性
gene.activities <- GeneActivity(pbmc.atac, features = VariableFeatures(pbmc.rna))

# 将基因活性作为一个新的数据集添加到 ATAC 数据中
pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# 对基因活性数据进行归一化和缩放
DefaultAssay(pbmc.atac) <- "ACTIVITY"
pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac, features = rownames(pbmc.atac))

# 5. 识别 RNA 和 ATAC 之间的锚点，用于整合分析
# 使用 canonical correlation analysis (CCA) 识别 RNA 和 ATAC 数据之间的锚点
transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc.atac, features = VariableFeatures(object = pbmc.rna),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

# 6. 转移 RNA 注释到 ATAC 数据
# 利用之前识别的锚点将 RNA 数据中的细胞类型注释转移到 ATAC 数据中
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc.rna$seurat_annotations,
                                     weight.reduction = pbmc.atac[["lsi"]], dims = 2:30)

# 将预测的细胞类型注释结果添加到 ATAC 数据中
pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)

# 标注是否注释正确
pbmc.atac$annotation_correct <- pbmc.atac$predicted.id == pbmc.atac$seurat_annotations

# 可视化预测的注释结果和实际注释
p1 <- DimPlot(pbmc.atac, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
p2 <- DimPlot(pbmc.atac, group.by = "seurat_annotations", label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
p1 | p2

# 7. 评估注释准确率
# 计算预测与实际注释的匹配情况，并可视化
predictions <- table(pbmc.atac$seurat_annotations, pbmc.atac$predicted.id)
predictions <- predictions/rowSums(predictions)
predictions <- as.data.frame(predictions)
p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
                                                                                            low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# 对正确和错误注释的预测分数进行可视化
correct <- length(which(pbmc.atac$seurat_annotations == pbmc.atac$predicted.id))
incorrect <- length(which(pbmc.atac$seurat_annotations != pbmc.atac$predicted.id))
data <- FetchData(pbmc.atac, vars = c("prediction.score.max", "annotation_correct"))
p2 <- ggplot(data, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) +
  geom_density(alpha = 0.5) + theme_cowplot() + scale_fill_discrete(name = "Annotation Correct",
                                                                    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + scale_color_discrete(name = "Annotation Correct",
                                                                                                                                                                                  labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + xlab("Prediction Score")
p1 + p2

# 8. RNA 和 ATAC 数据联合嵌入分析
# 将 RNA 数据通过锚点“推算”到 ATAC 细胞中
genes.use <- VariableFeatures(pbmc.rna)
refdata <- GetAssayData(pbmc.rna, assay = "RNA", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc.atac[["lsi"]],
                           dims = 2:30)
pbmc.atac[["RNA"]] <- imputation

# 将 RNA 和 ATAC 数据合并成一个对象
coembed <- merge(x = pbmc.rna, y = pbmc.atac)

# 对联合数据进行 PCA 和 UMAP 降维
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

# 可视化联合嵌入的结果
DimPlot(coembed, group.by = c("orig.ident", "seurat_annotations"))

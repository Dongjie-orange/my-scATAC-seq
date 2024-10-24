rm(list = ls());gc() # 清除环境中的所有对象，释放内存

# 安装并加载所需R包，其中Signac是一个用于单细胞ATAC-seq数据分析的R包
# remotes::install_github("stuart-lab/signac", ref="develop")
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75) # 加载人类基因注释数据库
library(ggplot2)
library(patchwork)
library(AnnotationHub)

# 读取示例的ATAC-seq数据，包括计数矩阵和元数据
counts <- Read10X_h5(filename = "Data/10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "Data/10k_pbmc_ATACv2_nextgem_Chromium_Controller_singlecell.csv",
  header = TRUE,
  row.names = 1
)

# 创建ChromatinAssay对象，包含染色质开放区域的信息
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "Data/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

# 使用ChromatinAssay创建Seurat对象，以便后续的分析
pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

# 查看pbmc对象的信息和染色质开放区域的概况
pbmc
pbmc[['peaks']]
granges(pbmc)

# 过滤掉非标准染色体上的峰，保留22对常染色体和X、Y染色体上的数据
peaks.keep <- seqnames(granges(pbmc)) %in% standardChromosomes(granges(pbmc))
pbmc <- pbmc[as.vector(peaks.keep), ]

# 使用AnnotationHub获取Ensembl v98的基因注释数据
ah <- AnnotationHub()
query(ah, "EnsDb.Hsapiens.v98")
ensdb_v98 <- ah[["AH75011"]]

# 从Ensembl数据库中提取基因注释信息并添加到Seurat对象中
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations)) # 将基因组序列转换为UCSC样式
genome(annotations) <- "hg38"
Annotation(pbmc) <- annotations

# 计算每个细胞的核小体信号和TSS富集度，用于质控
pbmc <- NucleosomeSignal(object = pbmc)
pbmc <- TSSEnrichment(object = pbmc)

# 计算每个细胞在峰内的片段比例和黑名单区域比例
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- FractionCountsInRegion(
  object = pbmc, 
  assay = 'peaks',
  regions = blacklist_hg38_unified
)

# 绘制质控指标的散点图，检查细胞的测序深度和TSS富集度
DensityScatter(pbmc, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

# 根据核小体信号将细胞分组并绘制片段长度直方图
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')

# 绘制每个质控指标的小提琴图
VlnPlot(
  object = pbmc,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0,
  ncol = 5
)

# 根据质控指标过滤低质量细胞，保留高质量的细胞用于后续分析
pbmc <- subset(
  x = pbmc,
  subset = nCount_peaks > 9000 &
    nCount_peaks < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 4
)

# 对数据进行TF-IDF归一化并进行SVD降维分析
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

# 评估降维成分与测序深度的相关性，移除与测序深度强相关的成分
DepthCor(pbmc)

# 使用UMAP进行非线性降维，并使用基于图的聚类方法对细胞进行分群
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE) + NoLegend()

# 计算基因活性矩阵并将其添加到Seurat对象中
gene.activities <- GeneActivity(pbmc)
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)

# 可视化一些标志性基因的基因活性，例如B细胞、T细胞等
DefaultAssay(pbmc) <- 'RNA'
FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

# 加载预处理好的PBMC的scRNA-seq数据并进行跨模态标签转移
pbmc_rna <- readRDS("Data/pbmc_10k_v3.rds")
pbmc_rna <- UpdateSeuratObject(pbmc_rna)
transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc,
  reduction = 'cca'
)

# 将scRNA-seq的细胞类型标签转移到ATAC-seq数据中
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)

# 将预测的标签添加到pbmc对象中
pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)

# 可视化scRNA-seq和scATAC-seq的标签转移结果
plot1 <- DimPlot(
  object = pbmc_rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = pbmc,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2

# 过滤掉少于20个细胞的预测标签，保留主要的细胞类型
predicted_id_counts <- table(pbmc$predicted.id)
major_predicted_ids <- names(predicted_id_counts[predicted_id_counts > 20])
pbmc <- pbmc[, pbmc$predicted.id %in% major_predicted_ids]
Idents(pbmc) <- pbmc$predicted.id

# 重新使用峰数据进行差异峰的检测，比较不同细胞类型之间的差异
DefaultAssay(pbmc) <- 'peaks'
da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "CD4 Naive",
  ident.2 = "CD14+ Monocytes",
  test.use = 'wilcox',
  min.pct = 0.1
)

# 可视化差异峰的表达情况
plot1 <- VlnPlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("CD4 Naive","CD14+ Monocytes")
)
plot2 <- FeaturePlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

plot1 | plot2

# 找到与特定基因相关的差异峰，并绘制覆盖图
open_cd4naive <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_cd14mono <- rownames(da_peaks[da_peaks$avg_log2FC < -3, ])

closest_genes_cd4naive <- ClosestFeature(pbmc, regions = open_cd4naive)
closest_genes_cd14mono <- ClosestFeature(pbmc, regions = open_cd14mono)

head(closest_genes_cd4naive)
head(closest_genes_cd14mono)

# 对感兴趣的基因区域绘制覆盖图
pbmc <- SortIdents(pbmc)
regions_highlight <- subsetByOverlaps(StringToGRanges(open_cd4naive), LookupGeneCoords(pbmc, "CD4"))

CoveragePlot(
  object = pbmc,
  region = "CD4",
  region.highlight = regions_highlight,
  extend.upstream = 1000,
  extend.downstream = 1000
)
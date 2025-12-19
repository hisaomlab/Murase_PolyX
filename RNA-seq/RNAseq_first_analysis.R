# 必要なパッケージのインストールと読み込み
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("edgeR")
}
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx", repos = "http://cran.us.r-project.org")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "http://cran.us.r-project.org")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer", repos = "http://cran.us.r-project.org")
}
library(edgeR)
library(openxlsx)
library(ggplot2)
library(RColorBrewer)

# データのロード
file_path <- "/Users/muraseyukihiro/Desktop/PolyX/RNA-Seq/kitamura_PolyI_reanalysis/combined_TPM.csv"  # 適切なパスに置き換えてください
gene_info <- read.xlsx("/Users/muraseyukihiro/Desktop/PolyX/RNA-Seq/Gene_list_for_RNAseq.xlsx", sheet = 1)# 遺伝子IDと遺伝子名の対応表
counts <- read.csv(file_path, header = TRUE, row.names = 1)

# Gene_nameを分割してGene_idを作成
rownames_split <- strsplit(rownames(counts), "\\|")
gene_ids <- sapply(rownames_split, `[`, 1)
gene_names <- sapply(rownames_split, `[`, 2)

# Gene_idをデータフレームの行名として設定
rownames(counts) <- gene_ids

# 列名を確認してグループ名を自動取得
colnames(counts)
sample_names <- colnames(counts)
group_names <- sapply(strsplit(sample_names, "_(?=[^_]+$)", perl = TRUE), `[`, 1)  # 後ろから1つ目のアンダースコアで分割

# グループベクトルの定義
group_order <- c("PolyI", "moxYG", "mox", "EGFP", "Vector")  # 必要に応じてグループ名を指定(コントロールを後方に配置)
group <- factor(group_names)

# グループベクトルの長さが列数と一致していることを確認
n_samples <- ncol(counts)
if (length(group) != n_samples) {
  stop("The length of 'group' must equal the number of columns in 'counts'")
}

# グループ名を有効な名前に変換
group <- factor(make.names(group))

# デザインマトリックスの生成
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# DGEListオブジェクトの作成
dge <- DGEList(counts = counts, group = group)

# データの正規化
dge <- calcNormFactors(dge)

# 分散の推定
dge <- estimateDisp(dge, design = design)

# すべてのペアのDEGを抽出
all_deg_genes <- c()

# 各ペアについてglmLRTを実行し、DEGを集約
for (pair in combn(group_order, 2, simplify = FALSE)) {
  tryCatch({
    contrast <- makeContrasts(contrasts = paste(pair, collapse = "-"), levels = design)
    fit <- glmFit(dge, design)
    lrt <- glmLRT(fit, contrast = contrast)
    
    table <- topTags(lrt, n = nrow(counts))$table
    is.DEG <- table$FDR < 0.05
    DEG.names <- rownames(table)[is.DEG]
    
    all_deg_genes <- c(all_deg_genes, DEG.names)
  }, error = function(e) {
    message(sprintf("Error processing comparison %s vs %s: %s", pair[1], pair[2], e$message))
  })
}

# DEGの重複を削除
all_deg_genes <- unique(all_deg_genes)

# DEGのみを使用したデータを作成
deg_counts <- counts[all_deg_genes, ]

# 数値行列に変換
deg_counts_matrix <- as.matrix(deg_counts)

# 青赤のカラーパレットを設定
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# PCA解析
logCPM <- cpm(dge, log = TRUE)
pca <- prcomp(t(logCPM))

# 必要な数の色を生成
n_colors <- length(levels(group))
colors <- colorRampPalette(brewer.pal(min(12, n_colors), "Set3"))(n_colors)

# PCAプロット
pdf("/Users/muraseyukihiro/Desktop/PolyX/RNA-Seq/kitamura_PolyI_reanalysis/Analysis/PCA_plot.pdf", width = 10, height = 10)
plot(pca$x[,1], pca$x[,2], bg = colors[as.numeric(group)], col = "black", pch = 21, xlab = "PC1", ylab = "PC2", main = "PCA Plot")
legend("topright", legend = levels(group), col = colors, pch = 21, pt.bg = colors, cex = 0.7, bty = "n")
dev.off()

# Heatmapの作成
pdf("/Users/muraseyukihiro/Desktop/PolyX/RNA-Seq/kitamura_PolyI_reanalysis/Analysis/DEG_Heatmap.pdf", width = 10, height = 10)
heatmap(log2(deg_counts_matrix + 1), 
        col = heatmap_colors, 
        scale = "row", 
        margins = c(5, 10),
        cexRow = 0.5, 
        cexCol = 0.75,
        labCol = group,
        main = "DEG Heatmap")
dev.off()

# 結果を保存するためのExcelファイルを作成
wb <- createWorkbook()

# 各ペアについてglmLRTを実行し、結果を保存
for (pair in combn(group_order, 2, simplify = FALSE)) {
  tryCatch({
    contrast <- makeContrasts(contrasts = paste(pair, collapse = "-"), levels = design)
    fit <- glmFit(dge, design)
    lrt <- glmLRT(fit, contrast = contrast)
    
    table <- topTags(lrt, n = nrow(counts))$table
    table$Gene_id <- rownames(table)#gene_id列を追加
    table$Gene_name <- gene_info$Gene.name[match(table$Gene_id, gene_info$ID)]#gene_name列を追加
    table$Brief_description <- gene_info$Brief.description[match(table$Gene_id, gene_info$ID)]#Brief description列を対応表に沿って追加
    table <- table[, c("Gene_id", "Gene_name", colnames(table)[1:(ncol(table)-3)], "Brief_description")] # Gene_idとGene_name列を先頭に移動
    
    sheet_name <- sprintf("%s_vs_%s", pair[1], pair[2])
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, table)
    
    # Volcano Plotはpythonで書くことにした
    
    # MA Plotの作成
    plot_file_name <- sprintf("/Users/muraseyukihiro/Desktop/PolyX/RNA-Seq/kitamura_PolyI_reanalysis/Analysis/MA_plot_%s_vs_%s.pdf", pair[1], pair[2])  # 適切なパスに置き換えてください
    pdf(plot_file_name)
    plotMD(lrt, main = sprintf("MA Plot: %s vs %s", pair[1], pair[2]), ylim = c(-5, 5))
    dev.off()
    
    # Smear Plotの作成
    plot_file_name <- sprintf("/Users/muraseyukihiro/Desktop/PolyX/RNA-Seq/kitamura_PolyI_reanalysis/Analysis/Smear_plot_%s_vs_%s.pdf", pair[1], pair[2])  # 適切なパスに置き換えてください
    pdf(plot_file_name)
    is.DEG <- table$FDR < 0.05
    DEG.names <- table$Gene_id[is.DEG]
    plotSmear(lrt, de.tags = DEG.names)
    dev.off()
  }, error = function(e) {
    message(sprintf("Error processing comparison %s vs %s: %s", pair[1], pair[2], e$message))
  })
}

# Excelファイルを保存
saveWorkbook(wb, "/Users/muraseyukihiro/Desktop/PolyX/RNA-Seq/kitamura_PolyI_reanalysis/Analysis/DEG_results.xlsx", overwrite = TRUE)  # 適切なパスに置き換えてください
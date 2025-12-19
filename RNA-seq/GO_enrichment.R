# ─────────────────────────────────────────────────────
# 必要なパッケージ検証・インストール
# ─────────────────────────────────────────────────────
required_cran  <- c("openxlsx", "ggplot2", "pbapply")
required_bioc  <- c("clusterProfiler", "org.Sc.sgd.db", "enrichplot")

# CRAN
new_cran <- required_cran[!required_cran %in% rownames(installed.packages())]
if (length(new_cran)) install.packages(new_cran, repos = "http://cran.us.r-project.org")

# Bioconductor
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
new_bioc <- required_bioc[!required_bioc %in% rownames(installed.packages())]
if (length(new_bioc)) BiocManager::install(new_bioc, ask = FALSE)

lapply(c(required_cran, required_bioc), library, character.only = TRUE)


# ─────────────────────────────────────────────────────
# 共通設定
# ─────────────────────────────────────────────────────
base_dir    <- "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig5/RNA-seq/20241023_PolyE_final/Analysis"
input_file  <- file.path(base_dir, "DEG_results.xlsx")
output_dir  <- file.path(base_dir, "GO_result")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sheets      <- openxlsx::getSheetNames(input_file)
ontologies  <- c("BP", "MF", "CC")
dirs        <- list(
  up   = list(prefix = "up",   logfc =  1),
  down = list(prefix = "down", logfc = -1)
)

# ─────────────────────────────────────────────────────
# GO解析関数
# ─────────────────────────────────────────────────────
run_go_analysis <- function(res_df, direction, sheet_name, wb, out_dir) {
  # フィルタリング
  filt <- if (direction == "up") {
    res_df$FDR < 0.05 & res_df$logFC >= dirs[[direction]]$logfc
  } else {
    res_df$FDR < 0.05 & res_df$logFC <= dirs[[direction]]$logfc
  }
  subset <- res_df[filt & !is.na(res_df$Gene_id), ]
  if (nrow(subset) == 0) {
    message(sprintf("[%s][%s] 条件を満たす遺伝子がありません", sheet_name, direction))
    return(invisible(NULL))
  }
  
  # ENTREZID 取得
  entrez <- AnnotationDbi::select(
    org.Sc.sgd.db,
    keys    = subset$Gene_id,
    columns = "ENTREZID",
    keytype = "ORF"
  )$ENTREZID |> unique() |> na.omit()
  
  all_entrez <- AnnotationDbi::keys(
    org.Sc.sgd.db,
    keytype = "ENTREZID"
  )
  
  # enrichGO & simplify
  results <- lapply(ontologies, function(ont) {
    ego <- tryCatch(
      clusterProfiler::enrichGO(
        gene          = entrez,
        OrgDb         = org.Sc.sgd.db,   # 文字列ではなくオブジェクトを渡す
        keyType       = "ENTREZID",
        ont           = ont,
        universe      = all_entrez,
        pAdjustMethod = "holm",　# Holm–Bonferroni法
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05,
        readable      = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(ego) && nrow(ego@result) > 0) {
      clusterProfiler::simplify(ego, cutoff = 0.5, by = "p.adjust", measure = "Wang")
    } else {
      NULL
    }
  })
  
  # シート書き出し & dotplot保存
  for (i in seq_along(ontologies)) {
    ont <- ontologies[i]
    sheet_title <- paste0(direction, "_", ont)
    addWorksheet(wb, sheet_title)
    ego_res <- results[[i]]
    if (!is.null(ego_res) && nrow(ego_res@result) > 0) {
      writeData(wb, sheet = sheet_title, as.data.frame(ego_res))
      
      n <- nrow(ego_res@result)
      p <- dotplot(ego_res, showCategory = n, label_format = 50) + 
        scale_x_continuous(limits = c(0, 0.4)) + 
        theme(plot.margin = margin(5, 5, 5, 5))
      plot_file <- file.path(out_dir, 
                             sprintf("dotplot_%s_%s_%s.png", sheet_name, direction, ont)
      )
      ggsave(filename = plot_file, plot = p, width = 8, height = max(4, n * 0.5))
    } else {
      writeData(wb, sheet = sheet_title, sprintf("有意なGO項目がありません (%s)", ont))
    }
  }
}


# ─────────────────────────────────────────────────────
# 全シート処理（進捗バー付き）
# ─────────────────────────────────────────────────────
pbapply::pblapply(sheets, function(sheet) {
  res <- read.xlsx(input_file, sheet = sheet)
  wb  <- createWorkbook()
  for (dir in names(dirs)) {
    run_go_analysis(res, dir, sheet, wb, output_dir)
  }
  out_file <- file.path(output_dir, sprintf("GO_results_%s.xlsx", sheet))
  saveWorkbook(wb, file = out_file, overwrite = TRUE)
})

message("全シートのGO解析とdotplot作成が完了しました: ", output_dir)


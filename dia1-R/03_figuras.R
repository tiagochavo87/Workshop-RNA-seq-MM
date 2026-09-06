# =============================================================================
# Workshop RNA-seq — Mieloma Múltiplo
# 03 — FIGURAS (complementar)
#
# Material extra. Não faz parte do fluxo do workshop; serve para quem quiser
# gerar figuras publicáveis a partir da análise que já rodou.
#
# COMO USAR
#   1. Rode o 01_workshop_mieloma.R ATÉ O FIM primeiro.
#   2. Abra este arquivo e rode a seção 1 (Ctrl+Alt+T).
#   3. Daí em diante, chame as funções com os parâmetros que quiser.
#
# O QUE TEM AQUI
#   volcano()            volcano com rotulagem configurável
#   ma_plot()            MA plot, bruto ou encolhido
#   pca()                PCA colorido por qualquer coluna do metadata
#   heatmap_variancia()  heatmap NÃO circular — genes por variância
#   expressao()          boxplot de um gene ou de vários, em painéis
#   salvar_tudo()        exporta o conjunto em PNG 300 dpi e PDF vetorial
#
# Toda função aceita `arquivo=` para salvar. Sem isso, só desenha na tela.
# As figuras vão para a pasta figuras_publicacao/, criada na primeira execução —
# separada de figuras/, que é do script 01.
# =============================================================================


# 1 — Conferir o ambiente -----------------------------------------------------

.faltando <- setdiff(c("tab", "dds", "vsd", "meta", "REFERENCIA", "ALVO",
                       "PADJ", "LFC", "TESTAR_LIMIAR_NO_MODELO"),
                     ls(envir = globalenv()))
if (length(.faltando) > 0) {
  stop(sprintf(paste0(
    "Faltam objetos na memória: %s\n\n",
    "Este script continua de onde o 01_workshop_mieloma.R parou.\n",
    "Abra o 01, rode até o fim, e volte."),
    paste(.faltando, collapse = ", ")), call. = FALSE)
}

suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(ggrepel); library(pheatmap)
})

# 🔧 Pasta de saída DESTE script. Fica separada de figuras/, que é do 01_ —
#    assim as figuras publicáveis não se misturam com as da aula, e rodar este
#    script de novo não sobrescreve nada do fluxo principal.
PASTA_FIGURAS <- "figuras_publicacao"

if (!dir.exists(PASTA_FIGURAS)) {
  dir.create(PASTA_FIGURAS, recursive = TRUE, showWarnings = FALSE)
  cat("Pasta criada:", normalizePath(PASTA_FIGURAS, winslash = "/"), "\n")
} else {
  cat("Pasta de saída:", normalizePath(PASTA_FIGURAS, winslash = "/"), "\n")
}

# Paleta única, para as figuras saírem coerentes entre si
COR_REF   <- "#1b6ca8"
COR_ALVO  <- "#d1495b"
COR_NEUTRA <- "#d5d5d5"
CORES_COND <- setNames(c(COR_REF, COR_ALVO), c(REFERENCIA, ALVO))

DPI    <- 300     # 🔧 resolução dos PNG. Revista costuma pedir 300 ou 600.
TEMA   <- theme_minimal(base_size = 12) +
          theme(panel.grid.minor = element_blank(),
                plot.title = element_text(face = "bold"),
                legend.position = "top")

# Formatação de p-valor: nada de 1.00e+00
fmt_padj <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0,001")
  sub("\\.", ",", sprintf("%.3f", p))
}
estrelas <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("***"); if (p < 0.01) return("**")
  if (p < 0.05)  return("*");   "ns"
}

# O mesmo critério de significância do 01 — não reinventar aqui
.eh_sig <- function(d) {
  ok <- !is.na(d$padj) & d$padj < PADJ
  if (!TESTAR_LIMIAR_NO_MODELO) ok <- ok & abs(d$log2FoldChange) >= LFC
  ok
}

.salvar <- function(g, arquivo, larg, alt) {
  if (is.null(arquivo)) return(invisible(NULL))
  caminho <- file.path(PASTA_FIGURAS, arquivo)
  ggsave(paste0(caminho, ".png"), g, width = larg, height = alt, dpi = DPI)
  ggsave(paste0(caminho, ".pdf"), g, width = larg, height = alt)   # vetorial
  cat("  salvo:", caminho, ".png e .pdf\n", sep = "")
}

# decimal.mark explícito: sem ele o R avisa que big.mark e decimal.mark são
# ambos "." — que é justamente o caso no formato brasileiro.
cat(sprintf("Ambiente pronto. %s genes testados, %d amostras.\n",
            format(sum(!is.na(tab$padj)), big.mark = ".", decimal.mark = ",",
                   scientific = FALSE),
            nrow(meta)))
cat("Funções: volcano() ma_plot() pca() heatmap_variancia() expressao() salvar_tudo()\n\n")


# 2 — Volcano -----------------------------------------------------------------
#
#   coluna    "log2FC_bruto" (o que foi testado) ou "log2FoldChange" (encolhido)
#   rotular   vetor de símbolos para rotular. NULL = os `n_rotulos` mais significativos
#   destacar  vetor de símbolos a circular, mesmo que não sejam significativos
#   limite_x  corta o eixo x. NULL = automático
#
# ⚠️ O padrão é o log2FC BRUTO, e por um motivo: o eixo y vem do teste feito
#    sobre ele. Plotar o encolhido contra esse mesmo p mistura duas contas e
#    produz genes em posições impossíveis — alto e colado no zero.

volcano <- function(coluna = "log2FC_bruto",
                    rotular = NULL, n_rotulos = 12,
                    destacar = NULL,
                    limite_x = NULL, limite_y = NULL,
                    titulo = NULL, arquivo = NULL) {

  d <- tab[!is.na(tab$padj) & !is.na(tab[[coluna]]), ]
  d$x   <- d[[coluna]]
  d$y   <- -log10(pmax(d$padj, 1e-300))
  d$sig <- .eh_sig(d)

  rot_up   <- paste("aumentado em", ALVO)
  rot_down <- paste("reduzido em",  ALVO)
  d$grupo <- ifelse(d$sig & d$x > 0, rot_up,
             ifelse(d$sig & d$x < 0, rot_down, "ns"))

  if (is.null(rotular)) {
    alvo_rot <- head(d[d$sig, ][order(d$padj[d$sig]), "gene_name"], n_rotulos)
  } else {
    alvo_rot <- rotular
  }
  rot <- d[d$gene_name %in% alvo_rot, ]
  hl  <- if (is.null(destacar)) d[0, ] else d[d$gene_name %in% destacar, ]

  if (is.null(limite_x)) limite_x <- ceiling(quantile(abs(d$x), .999, names = FALSE)) + .5
  if (is.null(titulo))
    titulo <- sprintf("%s vs %s — %s", ALVO, REFERENCIA,
                      if (coluna == "log2FC_bruto") "log2FC bruto" else "log2FC encolhido")

  g <- ggplot(d, aes(x, y, color = grupo)) +
    geom_hline(yintercept = -log10(PADJ), linetype = "dashed", color = "grey60", linewidth = .4) +
    geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed", color = "grey60", linewidth = .4) +
    geom_point(size = 1, alpha = .6) +
    scale_color_manual(values = setNames(c(COR_ALVO, COR_REF, COR_NEUTRA),
                                         c(rot_up, rot_down, "ns")), drop = FALSE) +
    labs(x = sprintf("log2 fold change (%s / %s)", ALVO, REFERENCIA),
         y = expression(-log[10]~"(padj)"), color = NULL,
         title = titulo,
         subtitle = sprintf("%s genes · %d significativos · %s",
                            format(nrow(d), big.mark = ".", decimal.mark = ",",
                                   scientific = FALSE), sum(d$sig),
                            if (TESTAR_LIMIAR_NO_MODELO)
                              sprintf("padj < %.2f, limiar no teste", PADJ)
                            else sprintf("padj < %.2f e |LFC| >= %.1f", PADJ, LFC))) +
    coord_cartesian(xlim = c(-limite_x, limite_x), ylim = limite_y) + TEMA

  if (nrow(hl))
    g <- g + geom_point(data = hl, shape = 21, size = 3, stroke = 1,
                        color = "black", fill = NA)
  if (nrow(rot))
    g <- g + geom_text_repel(data = rot, aes(label = gene_name), color = "black",
                             size = 3.2, fontface = "bold", max.overlaps = 30,
                             min.segment.length = 0, box.padding = .4,
                             segment.color = "grey60")
  print(g)
  .salvar(g, arquivo, 8, 6.5)
  invisible(g)
}


# 3 — MA plot -----------------------------------------------------------------

ma_plot <- function(coluna = "log2FoldChange", arquivo = NULL) {
  d <- tab[!is.na(tab$padj) & !is.na(tab[[coluna]]) & tab$baseMean > 0, ]
  d$y <- d[[coluna]]; d$sig <- .eh_sig(d)

  g <- ggplot(d, aes(baseMean, y)) +
    geom_hline(yintercept = 0, color = "black", linewidth = .4) +
    geom_hline(yintercept = c(-LFC, LFC), linetype = "dashed",
               color = "grey60", linewidth = .4) +
    geom_point(data = d[!d$sig, ], color = COR_NEUTRA, size = .8, alpha = .5) +
    geom_point(data = d[d$sig, ],  color = COR_ALVO,   size = 1.4) +
    scale_x_log10(labels = scales::comma) +
    labs(x = "Expressão média normalizada (baseMean)",
         y = sprintf("log2 fold change (%s / %s)", ALVO, REFERENCIA),
         title = "MA plot",
         subtitle = if (coluna == "log2FC_bruto") "log2FC bruto"
                    else "log2FC encolhido (apeglm)") + TEMA
  print(g)
  .salvar(g, arquivo, 8, 5.5)
  invisible(g)
}


# 4 — PCA ---------------------------------------------------------------------
#
#   colorir_por  qualquer coluna do `meta`. Útil para procurar confundidor:
#                se separar por algo que não é a condição, você tem um problema.

pca <- function(colorir_por = "condition", n_genes = 2000,
                rotular_amostras = FALSE, arquivo = NULL) {

  if (!colorir_por %in% colnames(meta))
    stop("Coluna '", colorir_por, "' não existe no meta. Disponíveis: ",
         paste(colnames(meta), collapse = ", "))

  p  <- plotPCA(vsd, intgroup = colorir_por, ntop = n_genes, returnData = TRUE)
  pv <- round(100 * attr(p, "percentVar"), 1)
  p$grupo <- p[[colorir_por]]

  g <- ggplot(p, aes(PC1, PC2, color = grupo)) +
    geom_point(size = 3.5, alpha = .85) +
    labs(x = sprintf("PC1 (%.1f%% da variância)", pv[1]),
         y = sprintf("PC2 (%.1f%% da variância)", pv[2]),
         color = colorir_por,
         title = sprintf("PCA — VST, %s genes mais variáveis",
                         format(n_genes, big.mark = ".", decimal.mark = ",",
                                scientific = FALSE))) + TEMA

  if (identical(colorir_por, "condition"))
    g <- g + scale_color_manual(values = CORES_COND)
  if (rotular_amostras)
    g <- g + geom_text_repel(aes(label = substr(rownames(p), 1, 8)),
                             size = 2.6, color = "grey30", max.overlaps = 40)
  print(g)
  .salvar(g, arquivo, 7, 6)
  invisible(g)
}


# 5 — Heatmap por variância ---------------------------------------------------
#
# ⚠️ Este é o heatmap honesto. Os genes são escolhidos por VARIÂNCIA, sem olhar
#    o contraste — então os grupos separarem ou não é informação de verdade.
#    Um heatmap dos "top N por padj" separa por construção e não prova nada.

heatmap_variancia <- function(n = 50, agrupar_amostras = FALSE, arquivo = NULL) {
  m <- assay(vsd)
  variaveis <- head(order(rowVars(m), decreasing = TRUE), n)
  z <- m[variaveis, , drop = FALSE]

  simb <- if (exists("mapa_genes") && !is.null(mapa_genes))
            mapa_genes[rownames(z), "gene_name"] else rownames(z)
  rownames(z) <- ifelse(is.na(simb) | simb == "", rownames(z), simb)

  ord  <- rownames(meta)[order(meta$condition)]
  anot <- data.frame(condicao = meta[ord, "condition"], row.names = ord)

  args <- list(z[, ord], scale = "row",
               cluster_cols = agrupar_amostras, show_colnames = FALSE,
               annotation_col = anot,
               annotation_colors = list(condicao = CORES_COND),
               color = colorRampPalette(c(COR_REF, "white", COR_ALVO))(51),
               main = sprintf("Top %d genes por variância — VST, z-score", n),
               fontsize_row = if (n <= 40) 8 else 6)
  do.call(pheatmap, args)
  if (!is.null(arquivo)) {
    do.call(pheatmap, c(args, list(filename = file.path(PASTA_FIGURAS, paste0(arquivo, ".png")),
                                   width = 10, height = max(6, n * .16))))
    cat("  salvo: ", PASTA_FIGURAS, "/", arquivo, ".png\n", sep = "")
  }
  invisible(z)
}


# 6 — Expressão de qualquer gene ----------------------------------------------
#
#   expressao("MYC")                          um gene
#   expressao(c("CD38","TNFRSF17","SLAMF7"))  vários, em painéis
#   expressao(tab$gene_name[1:8])             os oito primeiros da lista
#
# Aceita QUALQUER gene testado, não só os significativos.

expressao <- function(genes, violino = FALSE, escala_livre = TRUE,
                      n_colunas = NULL, arquivo = NULL) {

  genes <- unique(genes[!is.na(genes)])
  achados <- genes[genes %in% tab$gene_name]
  if (length(achados) == 0) {
    message("Nenhum destes genes está na tabela. Procure assim:")
    message('  grep("^', substr(genes[1], 1, 3), '", tab$gene_name, value = TRUE)')
    return(invisible(NULL))
  }
  if (length(achados) < length(genes))
    message("Fora da tabela: ", paste(setdiff(genes, achados), collapse = ", "))

  partes <- lapply(achados, function(nome) {
    r <- tab[which(tab$gene_name == nome)[1], ]
    d <- plotCounts(dds, gene = r$gene_id, intgroup = "condition", returnData = TRUE)
    d$rotulo <- sprintf("%s\nlog2FC %+.2f · padj %s %s",
                        nome, r$log2FoldChange,
                        if (!is.na(r$padj) && r$padj < 0.001) "" else "=",
                        paste(fmt_padj(r$padj), estrelas(r$padj)))
    d
  })
  d <- do.call(rbind, partes)
  d$rotulo <- factor(d$rotulo, levels = unique(d$rotulo))

  g <- ggplot(d, aes(condition, count, fill = condition))
  g <- if (violino)
         g + geom_violin(alpha = .55, width = .8, trim = FALSE) +
             geom_boxplot(width = .14, outlier.shape = NA, fill = "white", alpha = .9)
       else
         g + geom_boxplot(width = .5, outlier.shape = NA, alpha = .8)

  g <- g + geom_jitter(width = .16, size = 1.6, alpha = .65) +
    scale_y_log10(labels = scales::comma) +
    scale_fill_manual(values = CORES_COND) +
    labs(x = NULL, y = "Contagem normalizada (log10)") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold", size = 9.5))

  if (length(achados) > 1) {
    if (is.null(n_colunas)) n_colunas <- min(4, length(achados))
    g <- g + facet_wrap(~ rotulo, ncol = n_colunas,
                        scales = if (escala_livre) "free_y" else "fixed")
  } else {
    g <- g + labs(title = levels(d$rotulo)[1]) +
      theme(plot.title = element_text(face = "bold", size = 11))
  }

  print(g)
  n_col <- if (length(achados) > 1) min(4, length(achados)) else 1
  n_lin <- ceiling(length(achados) / n_col)
  .salvar(g, arquivo, 2.5 * n_col + .6, 3.1 * n_lin + .4)
  invisible(g)
}


# 7 — Gerar o conjunto todo de uma vez ----------------------------------------

salvar_tudo <- function(prefixo = "fig") {
  cat("Gerando o conjunto de figuras...\n\n")
  volcano("log2FC_bruto",   arquivo = paste0(prefixo, "_volcano_bruto"))
  volcano("log2FoldChange", arquivo = paste0(prefixo, "_volcano_encolhido"))
  ma_plot(arquivo = paste0(prefixo, "_ma"))
  pca(arquivo = paste0(prefixo, "_pca"))
  heatmap_variancia(50, arquivo = paste0(prefixo, "_heatmap_variancia"))

  sig <- subset(tab, .eh_sig(tab))
  if (nrow(sig) > 0)
    expressao(head(sig$gene_name, 12), arquivo = paste0(prefixo, "_degs"))
  else
    cat("  (nenhum gene significativo — painel de DEGs não gerado)\n")

  if (exists("GENES_INTERESSE"))
    expressao(GENES_INTERESSE, arquivo = paste0(prefixo, "_genes_interesse"))

  cat(sprintf("\nPronto. PNG em 300 dpi e PDF vetorial em %s/\n", PASTA_FIGURAS))
  cat("O PDF é o que a revista quer: o texto continua sendo texto.\n")
}


# 8 — Exemplos ----------------------------------------------------------------
# Descomente e rode o que interessar.

volcano()                                        # padrão: LFC bruto
volcano("log2FoldChange")                        # com o encolhido, para comparar
volcano(rotular = c("CXCL9","CD4","IFI27"))      # rotular só estes
volcano(destacar = GENES_INTERESSE)              # circular os canônicos
volcano(limite_x = 4, n_rotulos = 20)            # eixo cortado, mais rótulos

ma_plot()
ma_plot("log2FC_bruto")

pca()
pca(colorir_por = "sample_type")                 # procurando confundidor
pca(rotular_amostras = TRUE)                     # achar o outlier

heatmap_variancia(30)
heatmap_variancia(50, agrupar_amostras = TRUE)   # veja se separa sozinho

expressao("MYC")
expressao(c("CD38","TNFRSF17","SLAMF7"), violino = TRUE)
expressao(head(subset(tab, .eh_sig(tab))$gene_name, 8))

salvar_tudo()

# -----------------------------------------------------------------------------
# Material didático complementar. Dados de acesso aberto do NCI Genomic Data
# Commons (estudo MMRF CoMMpass).
# =============================================================================


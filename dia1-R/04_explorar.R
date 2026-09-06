# =============================================================================
# Workshop RNA-seq — Mieloma Múltiplo
# 04 — EXPLORAÇÃO INTERATIVA (complementar)
#
# O 03_figuras.R gera figura estática, para publicação. Este aqui é o oposto:
# serve para VASCULHAR o resultado — passar o mouse em qualquer um dos 19 mil
# genes, buscar por nome, ordenar por qualquer coluna, folhear todos os
# boxplots.
#
# COMO USAR
#   1. Rode o 01_workshop_mieloma.R ATÉ O FIM.
#   2. Rode a seção 1 daqui (Ctrl+Alt+T).
#   3. Chame as funções. Cada uma abre no navegador e salva um HTML autônomo —
#      o arquivo continua interativo depois, sem precisar do R.
#
# O QUE TEM AQUI
#   volcano_web()      volcano com TODOS os genes no hover
#   ma_web()           MA plot interativo
#   pca_web()          PCA com identificação de amostra no hover
#   heatmap_web()      heatmap com gene, amostra e z-score no hover
#   gene_web()         boxplot interativo de um gene
#   tabela_web()       tabela pesquisável e ordenável, todos os genes
#   pdf_todos_genes()  PDF paginado com o boxplot de cada gene
#   painel_web()       gera tudo de uma vez, com um índice HTML
# =============================================================================


# 1 — Ambiente ----------------------------------------------------------------

.faltando <- setdiff(c("tab", "dds", "vsd", "meta", "REFERENCIA", "ALVO",
                       "PADJ", "LFC", "TESTAR_LIMIAR_NO_MODELO"),
                     ls(envir = globalenv()))
if (length(.faltando) > 0) {
  stop(sprintf(paste0(
    "Faltam objetos na memória: %s\n\n",
    "Este script continua de onde o 01_workshop_mieloma.R parou."),
    paste(.faltando, collapse = ", ")), call. = FALSE)
}

suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })

# plotly não vem no 00_instalar_pacotes.R — é dependência só deste script.
if (!requireNamespace("plotly", quietly = TRUE)) {
  cat("O plotly não está instalado. Ele é necessário só para este script.\n")
  cat("Instalando agora (1 a 3 minutos)...\n\n")
  install.packages("plotly")
}
if (!requireNamespace("plotly", quietly = TRUE))
  stop("Não consegui instalar o plotly. Rode install.packages(\"plotly\") à mão.",
       call. = FALSE)
suppressPackageStartupMessages(library(plotly))

# 🔧 Pasta de saída deste script
PASTA_WEB <- "exploracao"
if (!dir.exists(PASTA_WEB)) {
  dir.create(PASTA_WEB, recursive = TRUE, showWarnings = FALSE)
  cat("Pasta criada:", normalizePath(PASTA_WEB, winslash = "/"), "\n")
} else {
  cat("Pasta de saída:", normalizePath(PASTA_WEB, winslash = "/"), "\n")
}

COR_REF <- "#1b6ca8"; COR_ALVO <- "#d1495b"; COR_NEUTRA <- "#d5d5d5"
CORES_COND <- setNames(c(COR_REF, COR_ALVO), c(REFERENCIA, ALVO))

.fmt_padj <- function(p) ifelse(is.na(p), "NA",
                         ifelse(p < 0.001, "< 0,001",
                                sub("\\.", ",", sprintf("%.3f", p))))
.estrelas <- function(p) ifelse(is.na(p), "ns",
                         ifelse(p < 0.001, "***",
                         ifelse(p < 0.01,  "**",
                         ifelse(p < 0.05,  "*", "ns"))))
.eh_sig <- function(d) {
  ok <- !is.na(d$padj) & d$padj < PADJ
  if (!TESTAR_LIMIAR_NO_MODELO) ok <- ok & abs(d$log2FoldChange) >= LFC
  ok
}
.mil <- function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE)

.salvar_web <- function(p, arquivo) {
  if (is.null(arquivo)) return(invisible(NULL))
  caminho <- file.path(PASTA_WEB, paste0(arquivo, ".html"))
  htmlwidgets::saveWidget(plotly::as_widget(p), file = normalizePath(caminho, mustWork = FALSE),
                          selfcontained = TRUE, title = arquivo)
  cat("  salvo:", caminho, " (abre no navegador, continua interativo)\n")
}

cat(sprintf("Ambiente pronto. %s genes, %d amostras.\n",
            .mil(sum(!is.na(tab$padj))), nrow(meta)))
cat("Funções: volcano_web() ma_web() pca_web() heatmap_web() gene_web()\n")
cat("         tabela_web() pdf_todos_genes() painel_web()\n\n")


# 2 — Volcano interativo ------------------------------------------------------
#
# Todos os genes testados estão aqui. Passe o mouse em qualquer ponto e veja
# nome, baseMean, os dois fold changes e o padj. Use a lupa para dar zoom,
# clique na legenda para ligar e desligar grupos, e o ícone da câmera para
# exportar PNG.

volcano_web <- function(coluna = "log2FC_bruto", destacar = NULL,
                        limite_x = NULL, arquivo = "volcano_interativo") {

  d <- tab[!is.na(tab$padj) & !is.na(tab[[coluna]]), ]
  d$x <- d[[coluna]]
  d$y <- -log10(pmax(d$padj, 1e-300))
  d$sig <- .eh_sig(d)
  d$situacao <- ifelse(d$sig & d$x > 0, paste("aumentado em", ALVO),
                ifelse(d$sig & d$x < 0, paste("reduzido em", ALVO),
                       "não significativo"))
  d$dica <- sprintf(
    paste0("<b>%s</b><br>baseMean: %s<br>log2FC bruto: %.3f",
           "<br>log2FC encolhido: %.3f<br>padj: %s %s"),
    d$gene_name, .mil(round(d$baseMean)), d$log2FC_bruto,
    d$log2FoldChange, .fmt_padj(d$padj), .estrelas(d$padj))

  if (is.null(limite_x))
    limite_x <- ceiling(quantile(abs(d$x), .999, names = FALSE)) + .5

  cores <- setNames(c(COR_ALVO, COR_REF, COR_NEUTRA),
                    c(paste("aumentado em", ALVO), paste("reduzido em", ALVO),
                      "não significativo"))

  p <- plot_ly(d, x = ~x, y = ~y, color = ~situacao, colors = cores,
               type = "scatter", mode = "markers",
               marker = list(size = 5, line = list(width = 0)),
               text = ~dica, hoverinfo = "text") %>%
    layout(
      title = list(text = sprintf("<b>%s vs %s</b> — %s · %d significativos de %s",
                                  ALVO, REFERENCIA,
                                  if (coluna == "log2FC_bruto") "log2FC bruto"
                                  else "log2FC encolhido",
                                  sum(d$sig), .mil(nrow(d))),
                   x = 0, font = list(size = 15)),
      xaxis = list(title = sprintf("log2 fold change (%s / %s)", ALVO, REFERENCIA),
                   range = c(-limite_x, limite_x), gridcolor = "#f0f0f0",
                   zeroline = FALSE),
      yaxis = list(title = "-log10 (padj)", gridcolor = "#f0f0f0", zeroline = FALSE),
      shapes = list(
        list(type = "line", x0 = -limite_x, x1 = limite_x,
             y0 = -log10(PADJ), y1 = -log10(PADJ),
             line = list(dash = "dash", color = "grey", width = 1)),
        list(type = "line", x0 = -LFC, x1 = -LFC, y0 = 0, y1 = max(d$y) * 1.05,
             line = list(dash = "dash", color = "grey", width = 1)),
        list(type = "line", x0 =  LFC, x1 =  LFC, y0 = 0, y1 = max(d$y) * 1.05,
             line = list(dash = "dash", color = "grey", width = 1))),
      legend = list(orientation = "h", y = 1.08, title = list(text = "")),
      plot_bgcolor = "white", margin = list(t = 70))

  if (!is.null(destacar)) {
    hl <- d[d$gene_name %in% destacar, ]
    if (nrow(hl))
      p <- p %>% add_trace(data = hl, x = ~x, y = ~y, type = "scatter",
                           mode = "markers+text", text = ~gene_name,
                           textposition = "top center",
                           textfont = list(size = 10, color = "black"),
                           marker = list(size = 11, color = "rgba(0,0,0,0)",
                                         line = list(width = 1.6, color = "black")),
                           name = "destacados", hoverinfo = "skip",
                           inherit = FALSE)
  }
  .salvar_web(p, arquivo)
  p
}


# 3 — MA interativo -----------------------------------------------------------

ma_web <- function(coluna = "log2FoldChange", arquivo = "ma_interativo") {
  d <- tab[!is.na(tab$padj) & !is.na(tab[[coluna]]) & tab$baseMean > 0, ]
  d$y <- d[[coluna]]; d$sig <- .eh_sig(d)
  d$situacao <- ifelse(d$sig, "significativo", "não significativo")
  d$dica <- sprintf("<b>%s</b><br>baseMean: %s<br>log2FC: %.3f<br>padj: %s",
                    d$gene_name, .mil(round(d$baseMean)), d$y, .fmt_padj(d$padj))

  p <- plot_ly(d, x = ~baseMean, y = ~y, color = ~situacao,
               colors = setNames(c(COR_ALVO, COR_NEUTRA),
                                 c("significativo", "não significativo")),
               type = "scatter", mode = "markers",
               marker = list(size = 5), text = ~dica, hoverinfo = "text") %>%
    layout(title = list(text = "<b>MA plot</b>", x = 0),
           xaxis = list(title = "baseMean (escala log)", type = "log",
                        gridcolor = "#f0f0f0"),
           yaxis = list(title = sprintf("log2 fold change (%s / %s)", ALVO, REFERENCIA),
                        gridcolor = "#f0f0f0"),
           legend = list(orientation = "h", y = 1.08, title = list(text = "")),
           plot_bgcolor = "white", margin = list(t = 60))
  .salvar_web(p, arquivo)
  p
}


# 4 — PCA interativo ----------------------------------------------------------
#
# O hover diz QUAL amostra é cada ponto. É assim que se identifica o outlier
# sem ficar contando pontos no gráfico estático.

pca_web <- function(colorir_por = "condition", n_genes = 2000,
                    arquivo = "pca_interativo") {
  if (!colorir_por %in% colnames(meta))
    stop("Coluna '", colorir_por, "' não existe. Disponíveis: ",
         paste(colnames(meta), collapse = ", "))

  d  <- plotPCA(vsd, intgroup = colorir_por, ntop = n_genes, returnData = TRUE)
  pv <- round(100 * attr(d, "percentVar"), 1)
  d$grupo <- as.character(d[[colorir_por]])
  d$amostra <- rownames(d)
  extras <- setdiff(intersect(c("case_id", "condition", "sample_type"),
                              colnames(meta)), colorir_por)
  d$dica <- sprintf("<b>%s</b><br>%s: %s%s",
    substr(d$amostra, 1, 12), colorir_por, d$grupo,
    if (length(extras))
      paste0("<br>", paste(sprintf("%s: %s", extras,
             sapply(extras, function(k) as.character(meta[d$amostra, k]))),
             collapse = "<br>")) else "")

  cores <- if (identical(colorir_por, "condition")) CORES_COND else NULL
  p <- plot_ly(d, x = ~PC1, y = ~PC2, color = ~grupo, colors = cores,
               type = "scatter", mode = "markers",
               marker = list(size = 11, line = list(width = 1, color = "white")),
               text = ~dica, hoverinfo = "text") %>%
    layout(title = list(text = sprintf("<b>PCA</b> — VST, %s genes mais variáveis",
                                       .mil(n_genes)), x = 0),
           xaxis = list(title = sprintf("PC1 (%.1f%%)", pv[1]), gridcolor = "#f0f0f0"),
           yaxis = list(title = sprintf("PC2 (%.1f%%)", pv[2]), gridcolor = "#f0f0f0"),
           legend = list(orientation = "h", y = 1.08, title = list(text = "")),
           plot_bgcolor = "white", margin = list(t = 60))
  .salvar_web(p, arquivo)
  p
}


# 5 — Heatmap interativo ------------------------------------------------------

heatmap_web <- function(n = 60, por = c("variancia", "padj"),
                        arquivo = "heatmap_interativo") {
  por <- match.arg(por)
  m <- assay(vsd)
  if (por == "variancia") {
    idx <- head(order(rowVars(m), decreasing = TRUE), n)
    titulo <- sprintf("Top %d genes por VARIÂNCIA — independente do contraste", n)
  } else {
    ids <- head(tab$gene_id[!is.na(tab$padj)], n)
    idx <- match(ids, rownames(m)); idx <- idx[!is.na(idx)]
    titulo <- sprintf("Top %d genes por padj — ⚠️ separação circular", n)
  }
  z <- m[idx, , drop = FALSE]
  z <- t(scale(t(z)))

  simb <- if (exists("mapa_genes") && !is.null(mapa_genes))
            mapa_genes[rownames(z), "gene_name"] else rownames(z)
  rownames(z) <- ifelse(is.na(simb) | simb == "", rownames(z), simb)
  ord <- rownames(meta)[order(meta$condition)]
  z <- z[, ord, drop = FALSE]

  dica <- outer(rownames(z), colnames(z),
                function(g, a) sprintf("<b>%s</b><br>amostra: %s<br>grupo: %s",
                                       g, substr(a, 1, 12),
                                       as.character(meta[a, "condition"])))
  dica <- matrix(dica, nrow = nrow(z))

  p <- plot_ly(x = seq_len(ncol(z)), y = rownames(z), z = z, type = "heatmap",
               colors = colorRamp(c(COR_REF, "white", COR_ALVO)),
               zmid = 0, zmin = -2.5, zmax = 2.5,
               text = dica, hoverinfo = "text+z",
               colorbar = list(title = "z-score")) %>%
    layout(title = list(text = paste0("<b>", titulo, "</b>"), x = 0),
           xaxis = list(title = sprintf("amostras (%s à esquerda, %s à direita)",
                                        REFERENCIA, ALVO),
                        showticklabels = FALSE),
           yaxis = list(title = "", tickfont = list(size = 8)),
           margin = list(t = 60, l = 110))
  .salvar_web(p, arquivo)
  p
}


# 6 — Boxplot interativo de um gene -------------------------------------------
#
# O hover diz qual amostra é cada ponto — útil para descobrir se o mesmo
# paciente extremo aparece em vários genes.

gene_web <- function(nome, arquivo = NULL) {
  linhas <- which(tab$gene_name == nome)
  if (length(linhas) == 0) {
    message("'", nome, "' não está entre os genes testados.")
    return(invisible(NULL))
  }
  r <- tab[linhas[1], ]
  d <- plotCounts(dds, gene = r$gene_id, intgroup = "condition", returnData = TRUE)
  d$amostra <- rownames(d)
  d$dica <- sprintf("<b>%s</b><br>%s<br>contagem: %s",
                    substr(d$amostra, 1, 12), d$condition, .mil(round(d$count)))

  p <- plot_ly(d, x = ~condition, y = ~count, color = ~condition,
               colors = CORES_COND, type = "box", boxpoints = "all",
               jitter = .4, pointpos = 0, text = ~dica, hoverinfo = "text",
               marker = list(size = 6, opacity = .7)) %>%
    layout(title = list(text = sprintf("<b>%s</b> — log2FC %+.2f · padj %s %s",
                                       nome, r$log2FoldChange,
                                       .fmt_padj(r$padj), .estrelas(r$padj)), x = 0),
           yaxis = list(title = "Contagem normalizada (log)", type = "log",
                        gridcolor = "#f0f0f0"),
           xaxis = list(title = ""), showlegend = FALSE,
           plot_bgcolor = "white", margin = list(t = 60))
  .salvar_web(p, if (is.null(arquivo)) NULL else arquivo)
  p
}


# 7 — Tabela pesquisável, TODOS os genes --------------------------------------
#
# HTML autônomo, sem dependência de pacote. Digite parte do nome de um gene
# para filtrar; clique no cabeçalho para ordenar. Funciona com os 19 mil.

tabela_web <- function(arquivo = "tabela_completa", ordenar_por = "padj") {
  d <- tab[order(tab[[ordenar_por]], na.last = TRUE), ]
  d$sig <- .eh_sig(d)

  esc <- function(x) {
    x <- as.character(x); x[is.na(x)] <- ""
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;",  x, fixed = TRUE)
    gsub(">", "&gt;", x, fixed = TRUE)
  }

  linhas <- sprintf(
    '<tr class="%s"><td>%s</td><td class="n">%s</td><td class="n">%+.3f</td><td class="n">%+.3f</td><td class="n">%s</td><td class="c">%s</td><td class="g">%s</td></tr>',
    ifelse(d$sig, "sig", ""), esc(d$gene_name), .mil(round(d$baseMean)),
    d$log2FC_bruto, d$log2FoldChange, .fmt_padj(d$padj), .estrelas(d$padj),
    esc(d$gene_id))

  html <- c(
'<!DOCTYPE html><html lang="pt-BR"><head><meta charset="utf-8">',
sprintf('<title>Genes testados — %s vs %s</title>', ALVO, REFERENCIA),
'<style>
 body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;margin:0;padding:24px 28px;color:#1a1a1a;background:#fafafa}
 h1{font-size:20px;margin:0 0 4px} .sub{color:#666;font-size:13px;margin-bottom:18px}
 .barra{position:sticky;top:0;background:#fafafa;padding:10px 0 14px;z-index:5}
 input{font-size:15px;padding:9px 13px;width:340px;border:1px solid #ccc;border-radius:6px}
 button{font-size:13px;padding:9px 14px;margin-left:8px;border:1px solid #ccc;
        border-radius:6px;background:#fff;cursor:pointer}
 button.on{background:#1f4e79;color:#fff;border-color:#1f4e79}
 #n{margin-left:12px;color:#666;font-size:13px}
 table{border-collapse:collapse;width:100%;background:#fff;font-size:13px}
 th{background:#e8edf3;color:#1f4e79;text-align:left;padding:9px 11px;
    position:sticky;top:62px;cursor:pointer;user-select:none;white-space:nowrap}
 th:hover{background:#dbe4ee} td{padding:7px 11px;border-top:1px solid #eee}
 td.n{text-align:right;font-variant-numeric:tabular-nums}
 td.c{text-align:center;font-weight:700;color:#9b2226}
 td.g{color:#999;font-size:11px;font-family:ui-monospace,Consolas,monospace}
 tr.sig{background:#fff6f7} tr:hover{background:#eef4fb}
</style></head><body>',
sprintf('<h1>Genes testados — %s vs %s</h1>', ALVO, REFERENCIA),
sprintf('<div class="sub">%s genes · %d significativos · critério: %s</div>',
        .mil(nrow(d)), sum(d$sig),
        if (TESTAR_LIMIAR_NO_MODELO)
          sprintf("padj &lt; %.2f, limiar de %.1f dentro do teste", PADJ, LFC)
        else sprintf("padj &lt; %.2f e |log2FC| &ge; %.1f", PADJ, LFC)),
'<div class="barra"><input id="q" placeholder="Digite parte do nome de um gene…" autofocus>',
'<button id="so" onclick="alterna()">só significativos</button><span id="n"></span></div>',
'<table id="t"><thead><tr>',
'<th onclick="ordena(0,0)">gene</th><th onclick="ordena(1,1)">baseMean</th>',
'<th onclick="ordena(2,1)">log2FC bruto</th><th onclick="ordena(3,1)">log2FC encolhido</th>',
'<th onclick="ordena(4,1)">padj</th><th onclick="ordena(5,0)">&nbsp;</th>',
'<th onclick="ordena(6,0)">gene_id</th></tr></thead><tbody>',
linhas,
'</tbody></table>
<script>
 const corpo=document.getElementById("t").tBodies[0];
 const todas=Array.from(corpo.rows); let soSig=false;
 function num(s){s=s.replace(/\\./g,"").replace(",",".").replace("< 0","0");
                 const v=parseFloat(s); return isNaN(v)?Infinity:v;}
 function filtra(){
   const q=document.getElementById("q").value.trim().toUpperCase();
   let n=0; corpo.replaceChildren();
   for(const tr of todas){
     if(soSig && !tr.classList.contains("sig")) continue;
     if(q && !tr.cells[0].textContent.toUpperCase().includes(q)) continue;
     corpo.appendChild(tr); n++;
     if(n>3000) break;
   }
   document.getElementById("n").textContent =
     n>3000 ? "mostrando 3.000 — refine a busca" : n+" gene(s)";
 }
 function alterna(){ soSig=!soSig;
   document.getElementById("so").classList.toggle("on",soSig); filtra(); }
 let inv=false;
 function ordena(col,numerico){
   inv=!inv;
   todas.sort((a,b)=>{const x=a.cells[col].textContent,y=b.cells[col].textContent;
     const r=numerico?num(x)-num(y):x.localeCompare(y,"pt");return inv?-r:r;});
   filtra();
 }
 document.getElementById("q").addEventListener("input",filtra);
 filtra();
</script></body></html>')

  caminho <- file.path(PASTA_WEB, paste0(arquivo, ".html"))
  writeLines(html, caminho, useBytes = TRUE)
  cat("  salvo:", caminho, "\n")
  cat("  Abra no navegador: busca por nome, ordena por qualquer coluna.\n")
  if (interactive()) try(utils::browseURL(normalizePath(caminho)), silent = TRUE)
  invisible(caminho)
}


# 8 — PDF paginado com o boxplot de TODOS os genes ----------------------------
#
# 12 painéis por página. Com 19 mil genes são ~1.600 páginas e uns 20 minutos.
# O padrão é limitar aos genes com padj abaixo de `ate_padj`; passe todos = TRUE
# para gerar o transcriptoma inteiro.

pdf_todos_genes <- function(ate_padj = 0.20, todos = FALSE,
                            por_pagina = 12, arquivo = "todos_os_genes") {

  d <- tab[!is.na(tab$padj), ]
  d <- d[order(d$padj), ]
  if (!todos) {
    d <- d[d$padj <= ate_padj, ]
    if (nrow(d) == 0) {
      cat("Nenhum gene com padj <=", ate_padj, "— aumente o corte ou use todos = TRUE.\n")
      return(invisible(NULL))
    }
  }

  n_pag <- ceiling(nrow(d) / por_pagina)
  cat(sprintf("Gerando %s gene(s) em %s página(s)...\n", .mil(nrow(d)), .mil(n_pag)))
  if (n_pag > 200)
    cat("  Isso vai demorar. Estimativa:", round(n_pag * 0.6 / 60), "a",
        round(n_pag * 1.2 / 60), "minutos.\n")

  caminho <- file.path(PASTA_WEB, paste0(arquivo, ".pdf"))
  grDevices::pdf(caminho, width = 11, height = 8.5)
  on.exit({ if (grDevices::dev.cur() > 1) grDevices::dev.off() }, add = TRUE)

  t0 <- Sys.time()
  for (pg in seq_len(n_pag)) {
    bloco <- d[((pg - 1) * por_pagina + 1):min(pg * por_pagina, nrow(d)), ]
    partes <- lapply(seq_len(nrow(bloco)), function(k) {
      r <- bloco[k, ]
      x <- plotCounts(dds, gene = r$gene_id, intgroup = "condition",
                      returnData = TRUE)
      x$rotulo <- sprintf("%s\nlog2FC %+.2f · padj %s %s",
                          r$gene_name, r$log2FoldChange,
                          .fmt_padj(r$padj), .estrelas(r$padj))
      x
    })
    dd <- do.call(rbind, partes)
    dd$rotulo <- factor(dd$rotulo, levels = unique(dd$rotulo))

    g <- ggplot(dd, aes(condition, count, fill = condition)) +
      geom_boxplot(width = .5, outlier.shape = NA, alpha = .8) +
      geom_jitter(width = .16, size = 1.1, alpha = .6) +
      scale_y_log10() +
      scale_fill_manual(values = CORES_COND) +
      facet_wrap(~ rotulo, ncol = 4, scales = "free_y") +
      labs(x = NULL, y = "Contagem normalizada (log10)",
           caption = sprintf("página %d de %d", pg, n_pag)) +
      theme_minimal(base_size = 9) +
      theme(legend.position = "none", panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold", size = 7.5))
    print(g)

    if (pg %% 25 == 0 || pg == n_pag) {
      decorrido <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      cat(sprintf("   %s/%s páginas (%.1f min)\n", .mil(pg), .mil(n_pag), decorrido))
    }
  }
  grDevices::dev.off()
  cat("  salvo:", caminho, "\n")
  invisible(caminho)
}


# 9 — Gerar tudo, com índice --------------------------------------------------

painel_web <- function(incluir_pdf = FALSE) {
  cat("Gerando o painel de exploração...\n\n")
  volcano_web("log2FC_bruto",   arquivo = "volcano_bruto")
  volcano_web("log2FoldChange", arquivo = "volcano_encolhido")
  ma_web();  pca_web();  heatmap_web(60)
  tabela_web()

  sig <- subset(tab, .eh_sig(tab))
  if (nrow(sig) > 0)
    for (g in head(sig$gene_name, 6))
      gene_web(g, arquivo = paste0("gene_", gsub("[^A-Za-z0-9]", "_", g)))

  if (incluir_pdf) pdf_todos_genes()

  # índice
  arqs <- sort(list.files(PASTA_WEB, pattern = "\\.(html|pdf)$"))
  arqs <- setdiff(arqs, "index.html")
  itens <- sprintf('<li><a href="%s">%s</a></li>', arqs, sub("\\.[a-z]+$", "", arqs))
  writeLines(c(
    '<!DOCTYPE html><html lang="pt-BR"><head><meta charset="utf-8">',
    '<title>Exploração — RNA-seq mieloma</title><style>',
    'body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;max-width:640px;',
    'margin:60px auto;color:#1a1a1a}h1{font-size:22px}li{margin:8px 0;font-size:15px}',
    'a{color:#1f4e79}</style></head><body>',
    sprintf('<h1>Exploração — %s vs %s</h1>', ALVO, REFERENCIA),
    '<p style="color:#666">Cada arquivo é autônomo: abre no navegador sem precisar do R.</p>',
    '<ul>', itens, '</ul></body></html>'),
    file.path(PASTA_WEB, "index.html"), useBytes = TRUE)

  cat(sprintf("\nPronto. %d arquivos em %s/\n", length(arqs), PASTA_WEB))
  cat("Comece pelo index.html\n")
  if (interactive())
    try(utils::browseURL(normalizePath(file.path(PASTA_WEB, "index.html"))),
        silent = TRUE)
}


# 10 — Exemplos ---------------------------------------------------------------

volcano_web()                                  # todos os genes no hover
volcano_web(destacar = GENES_INTERESSE)
volcano_web("log2FoldChange", limite_x = 5)

ma_web()
pca_web()
pca_web(colorir_por = "sample_type")

heatmap_web(60)                                # por variância
heatmap_web(40, por = "padj")                  # o circular, para comparar

gene_web("MYC")
tabela_web()                                   # busca em todos os genes

pdf_todos_genes()                              # padj <= 0,20
pdf_todos_genes(todos = TRUE)                  # o transcriptoma inteiro

painel_web()                                   # tudo, com índice

# -----------------------------------------------------------------------------
# Material didático complementar. Dados de acesso aberto do NCI Genomic Data
# Commons (estudo MMRF CoMMpass).
# =============================================================================

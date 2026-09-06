# =============================================================================
# Workshop RNA-seq — Mieloma Múltiplo
# 02 — EXPLORAR GENES, UM DE CADA VEZ
#
# O gráfico que desmente. Uma lista de DEGs é uma afirmação estatística; o
# boxplot das contagens é o que mostra se ela se sustenta. Se o efeito
# "significativo" vem de duas amostras extremas, aparece aqui na hora.
#
# COMO USAR
#   1. Rode o 01_workshop_mieloma.R ATÉ O FIM primeiro.
#      Este script usa os objetos que ele deixou na memória.
#   2. Abra este arquivo e rode a seção 1 (Ctrl+Alt+T).
#   3. Daí em diante é mão na massa: troque o nome do gene e rode de novo.
#
# ⚠️ Não reinicie o R entre um script e outro. Se reiniciar, rode o 01 de novo.
# =============================================================================


# 1 — Conferir que o ambiente está pronto -------------------------------------

.faltando <- setdiff(
  c("tab", "dds", "meta", "REFERENCIA", "ALVO", "PADJ", "LFC"),
  ls(envir = globalenv())
)

if (length(.faltando) > 0) {
  stop(sprintf(paste0(
    "Faltam objetos na memória: %s\n\n",
    "Este script continua de onde o 01_workshop_mieloma.R parou.\n",
    "Abra o 01, rode até o fim (Ctrl+Shift+S faz o arquivo inteiro), e volte."),
    paste(.faltando, collapse = ", ")), call. = FALSE)
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

if (!dir.exists("figuras")) dir.create("figuras", showWarnings = FALSE)

cat("Ambiente pronto.\n")
cat(sprintf("  genes testados : %s\n", format(sum(!is.na(tab$padj)), big.mark = ".")))
cat(sprintf("  amostras       : %d (%s)\n", nrow(meta),
            paste(names(table(meta$condition)), table(meta$condition),
                  sep = "=", collapse = ", ")))
cat(sprintf("  critério de DEG: padj < %.2f\n\n", PADJ))


# 2 — As duas funções ---------------------------------------------------------

# ver_gene()   -> os números, em texto
# plotar_gene() -> o boxplot
#
# Use as duas juntas. O número diz o que o modelo achou; o gráfico diz se você
# acredita nele.

# Formata o padj do jeito que se reporta num artigo, não em notação científica.
#   >= 0,001  ->  o valor, com três casas         (padj = 0,039)
#   <  0,001  ->  o limite                        (padj < 0,001)
# As estrelas seguem a convenção usual: *** < 0,001, ** < 0,01, * < 0,05.
fmt_padj <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0,001")
  sub("\\.", ",", sprintf("%.3f", p))
}

estrelas <- function(p) {
  if (is.na(p))    return("ns")
  if (p < 0.001)   return("***")
  if (p < 0.01)    return("**")
  if (p < 0.05)    return("*")
  "ns"
}


ver_gene <- function(nome) {
  linhas <- which(tab$gene_name == nome)
  if (length(linhas) == 0) {
    message(sprintf("'%s' não está entre os genes testados.", nome))
    message("  Pode ter caído no filtro do M3, ou o símbolo pode estar escrito diferente.")
    message("  Tente:  grep(\"", substr(nome, 1, 3), "\", tab$gene_name, value = TRUE)")
    return(invisible(NULL))
  }
  r <- tab[linhas[1], ]
  cts <- counts(dds, normalized = TRUE)[r$gene_id, ]
  por_grupo <- tapply(cts, meta$condition, function(x)
    c(min = min(x), mediana = median(x), max = max(x)))

  cat(sprintf("\n%s  (%s)\n", nome, r$gene_id))
  cat(strrep("-", 58), "\n")
  cat(sprintf("  baseMean       : %s\n", format(round(r$baseMean), big.mark = ".")))
  cat(sprintf("  log2FC bruto   : %+.3f\n", r$log2FC_bruto))
  cat(sprintf("  log2FC encolh. : %+.3f\n", r$log2FoldChange))
  cat(sprintf("  padj           : %s  %s\n", fmt_padj(r$padj), estrelas(r$padj)))
  cat(sprintf("  significativo  : %s\n",
              if (!is.na(r$padj) && r$padj < PADJ) "SIM" else "não"))
  cat("\n  contagens normalizadas por grupo:\n")
  for (g in names(por_grupo)) {
    v <- por_grupo[[g]]
    cat(sprintf("    %-12s min %9.0f | mediana %9.0f | máx %9.0f\n",
                g, v["min"], v["mediana"], v["max"]))
  }

  # A razão máx/mediana é o alarme: se for enorme, o efeito depende de poucas
  # amostras, e o boxplot vai mostrar isso.
  razoes <- vapply(por_grupo, function(v)
    v["max"] / max(v["mediana"], 1), numeric(1))
  if (max(razoes) > 50) {
    cat(sprintf("\n  ⚠️ máx/mediana = %.0fx em '%s'. Efeito provavelmente sustentado\n",
                max(razoes), names(which.max(razoes))))
    cat("     por poucas amostras. Rode plotar_gene() e olhe os pontos.\n")
  }
  cat("\n")
  invisible(r)
}


plotar_gene <- function(nome, salvar = TRUE) {
  linhas <- which(tab$gene_name == nome)
  if (length(linhas) == 0) {
    message(sprintf("'%s' não está entre os genes testados.", nome))
    return(invisible(NULL))
  }
  r <- tab[linhas[1], ]
  d <- plotCounts(dds, gene = r$gene_id, intgroup = "condition", returnData = TRUE)

  g <- ggplot(d, aes(condition, count, fill = condition)) +
    geom_boxplot(width = .5, outlier.shape = NA, alpha = .8) +
    geom_jitter(width = .18, size = 2, alpha = .7) +
    scale_y_log10() +
    scale_fill_manual(values = setNames(c("#1b6ca8", "#d1495b"), c(REFERENCIA, ALVO))) +
    labs(title = sprintf("%s   —   log2FC = %+.2f   ·   padj %s %s",
                         nome, r$log2FoldChange,
                         if (!is.na(r$padj) && r$padj < 0.001) "" else "=",
                         paste(fmt_padj(r$padj), estrelas(r$padj))),
         subtitle = sprintf("baseMean = %s   ·   log2FC bruto = %+.2f   ·   %s",
                            format(round(r$baseMean), big.mark = "."), r$log2FC_bruto,
                            if (!is.na(r$padj) && r$padj < PADJ)
                              "significativo" else "não significativo"),
         x = NULL, y = "Contagem normalizada (log10)") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none", panel.grid.minor = element_blank())

  # Se o painel Plots estiver pequeno, avisa em vez de parar o script.
  tryCatch(print(g), error = function(e) {
    cat("(figura não coube no painel Plots. Aumente o painel ou clique em Zoom.)\n")
  })
  if (salvar) {
    caminho <- sprintf("figuras/gene_%s.png", nome)
    ggsave(caminho, g, width = 5, height = 5, dpi = 200)
    cat("  salvo:", caminho, "\n")
  }
  invisible(g)
}


# Vários genes de uma vez, para comparar lado a lado.
plotar_varios <- function(nomes, salvar = TRUE, arquivo = "genes_comparados") {
  achados <- nomes[nomes %in% tab$gene_name]
  if (length(achados) == 0) {
    message("Nenhum destes genes está na tabela.")
    return(invisible(NULL))
  }
  if (length(achados) < length(nomes)) {
    message("Fora da tabela: ", paste(setdiff(nomes, achados), collapse = ", "))
  }

  d <- do.call(rbind, lapply(achados, function(n) {
    r <- tab[which(tab$gene_name == n)[1], ]
    x <- plotCounts(dds, gene = r$gene_id, intgroup = "condition", returnData = TRUE)
    x$gene <- sprintf("%s\npadj %s %s  %s", n,
                      if (!is.na(r$padj) && r$padj < 0.001) "" else "=",
                      fmt_padj(r$padj), estrelas(r$padj))
    x
  }))
  d$gene <- factor(d$gene, levels = unique(d$gene))

  g <- ggplot(d, aes(condition, count, fill = condition)) +
    geom_boxplot(width = .55, outlier.shape = NA, alpha = .8) +
    geom_jitter(width = .18, size = 1.4, alpha = .6) +
    scale_y_log10() +
    scale_fill_manual(values = setNames(c("#1b6ca8", "#d1495b"), c(REFERENCIA, ALVO))) +
    facet_wrap(~ gene, scales = "free_y", nrow = 1) +
    labs(x = NULL, y = "Contagem normalizada (log10)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold", size = 10))

  tryCatch(print(g), error = function(e) {
    cat("(figura não coube no painel Plots. Aumente o painel ou clique em Zoom.)\n")
  })
  if (salvar) {
    caminho <- sprintf("figuras/%s.png", arquivo)
    ggsave(caminho, g, width = 2.6 * length(achados) + 1, height = 4.5, dpi = 200)
    cat("  salvo:", caminho, "\n")
  }
  invisible(g)
}

cat("Funções carregadas: ver_gene(), plotar_gene(), plotar_varios()\n")
cat("Cada figura é salva em figuras/ automaticamente (salvar = FALSE desliga).\n")
cat("Estrelas: *** padj < 0,001   ** < 0,01   * < 0,05   ns = não significativo\n\n")


# 3 — Exercício 1: o gene do topo da sua lista --------------------------------
# Este é o gene com menor padj. Antes de olhar o gráfico, faça a aposta: você
# espera duas nuvens de pontos bem separadas?

.topo <- tab$gene_name[1]
cat("Gene do topo da lista:", .topo, "\n")

ver_gene(.topo)
plotar_gene(.topo)

# 💬 A separação é convincente? Os grupos se sobrepõem muito? Um punhado de
#    pontos está carregando a diferença sozinho?


# 4 — Exercício 2: um gene NÃO significativo ----------------------------------
# Para calibrar o olho, você precisa ver como é um gene que não mudou. Pegamos
# um do meio da tabela, com expressão parecida com a do gene do topo.

.nao_sig <- {
  candidatos <- subset(tab, !is.na(padj) & padj > 0.5 &
                            baseMean > tab$baseMean[1] * 0.5 &
                            baseMean < tab$baseMean[1] * 2)
  if (nrow(candidatos) > 0) candidatos$gene_name[1] else tab$gene_name[nrow(tab)]
}
cat("Gene não significativo, de expressão parecida:", .nao_sig, "\n")

ver_gene(.nao_sig)
plotar_gene(.nao_sig)

# 💬 Compare os dois gráficos lado a lado no painel Plots (setas ← →).
#    A diferença entre "significativo" e "não significativo" é essa.


# 5 — Exercício 3: o gene que o encolhimento derrubou -------------------------
# Procura o gene onde a diferença entre o fold change bruto e o encolhido foi
# maior. É a linha da tabela que parece se contradizer: efeito quase zero ao
# lado de um p pequeno.

.encolhidos <- tab[!is.na(tab$log2FC_bruto) & !is.na(tab$log2FoldChange), ]
.encolhidos$queda <- abs(.encolhidos$log2FC_bruto) - abs(.encolhidos$log2FoldChange)
.encolhidos <- .encolhidos[order(-.encolhidos$queda), ]

cat("\nOnde o encolhimento foi mais agressivo:\n")
print(head(.encolhidos[, c("gene_name", "baseMean", "log2FC_bruto",
                           "log2FoldChange", "padj")], 5),
      row.names = FALSE, digits = 3)

ver_gene(.encolhidos$gene_name[1])
plotar_gene(.encolhidos$gene_name[1])

# 💬 Olhe a escala do eixo y. A mediana está em dezenas e alguns pontos em
#    centenas de milhares. O apeglm viu essa incerteza e puxou a estimativa
#    para perto de zero. Foi a decisão certa.


# 6 — Exercício 4: os alvos terapêuticos --------------------------------------
# Aqui a pergunta deixa de ser estatística e passa a ser clínica.
#
#   TNFRSF17 = BCMA — alvo de CAR-T e de anticorpos biespecíficos
#   CD38     = alvo do daratumumabe
#   SLAMF7   = alvo do elotuzumabe
#
# Perda de antígeno é um mecanismo conhecido de escape terapêutico: o tumor
# deixa de expressar o alvo e a droga para de funcionar. Se isso acontecesse
# de forma generalizada entre o diagnóstico e a recidiva, apareceria aqui.

ALVOS_TERAPEUTICOS <- c("TNFRSF17", "CD38", "SLAMF7")

for (g in ALVOS_TERAPEUTICOS) ver_gene(g)

plotar_varios(ALVOS_TERAPEUTICOS, arquivo = "alvos_terapeuticos")

# 💬 PARE E DISCUTA
#
#  - Algum dos três mudou? (Nesta coorte, não.)
#
#  - O que isso significa: os alvos continuam expressos na recidiva. A perda de
#    antígeno existe e está descrita, mas é um fenômeno de subgrupo — ocorre em
#    pacientes específicos, sobretudo depois da exposição à droga. Uma média de
#    30 pacientes não captura isso.
#
#  - O que isso NÃO significa: que a perda de antígeno não importa. Significa
#    que ESTE desenho não consegue vê-la. Para isso seria preciso comparar cada
#    paciente antes e depois do tratamento com aquele alvo específico.
#
#  ⚠️ Repare no baseMean dos três. São dezenas de milhares. Ausência de mudança
#     não é ausência de expressão — é a frase mais importante do workshop, e
#     estes três genes são a demonstração dela.


# 7 — Agora é com você --------------------------------------------------------
# Troque o nome e rode a linha (Ctrl+Enter). Alguns pontos de partida:
#
#   ver_gene("MYC")        translocação envolvendo MYC é evento tardio no mieloma
#   ver_gene("NSD2")       o gene da t(4;14), junto com FGFR3
#   ver_gene("CCND1")      t(11;14) — o subtipo mais comum
#   ver_gene("XIST")       controle positivo do M8; separa por sexo, não por condição
#
# E se o símbolo não for encontrado, procure assim:
#
#   grep("^HLA", tab$gene_name, value = TRUE)
#   subset(tab, gene_name %in% c("CD38", "TNFRSF17"))

plotar_gene("MYC")


# 8 — O exercício que fecha ---------------------------------------------------
# Pegue os genes significativos da SUA rodada e olhe um por um. Se algum deles
# tiver o padrão do exercício 5 — mediana baixa e um punhado de pontos altos —
# ele não deveria estar no seu resumo, mesmo tendo passado no teste.

.sig <- subset(tab, !is.na(padj) & padj < PADJ)

if (nrow(.sig) == 0) {
  cat("\nNenhum gene significativo nesta configuração.\n")
  cat("Isso também é um resultado — e é o assunto do M4.\n")
} else {
  cat(sprintf("\nSeus %d gene(s) significativo(s): %s\n",
              nrow(.sig), paste(.sig$gene_name, collapse = ", ")))
  cat("Rode plotar_varios() abaixo e olhe cada um antes de reportar.\n\n")
  plotar_varios(head(.sig$gene_name, 6), arquivo = "meus_degs")
}

# 💬 A pergunta final: você assinaria embaixo de cada um destes?
#
#    Um resultado que você não olhou é um resultado que você não conferiu.
#    O teste estatístico é uma triagem, não um veredito.

# -----------------------------------------------------------------------------
# Material didático. Dados de acesso aberto do NCI Genomic Data Commons
# (estudo MMRF CoMMpass).
# =============================================================================

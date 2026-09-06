# =============================================================================
# Workshop RNA-seq — Mieloma Múltiplo
# 00 — INSTALAÇÃO DOS PACOTES
#
# >>> RODE ESTE ARQUIVO ANTES DO WORKSHOP. NÃO DEIXE PARA O DIA. <<<
#
# Como usar:
#   1. Abra este arquivo no RStudio
#   2. Clique em "Source" (canto superior direito) ou aperte Ctrl+Shift+S
#   3. Espere. A primeira instalação leva de 3 a 10 minutos.
#   4. No fim, ele imprime um relatório. Se aparecer qualquer ❌,
#      copie o relatório inteiro e mande para o instrutor.
#
# Você só precisa rodar isto UMA VEZ.
# =============================================================================


# 1. Diagnóstico do ambiente --------------------------------------------------

cat("\n=== AMBIENTE ===\n")
cat("R           :", R.version.string, "\n")
cat("Plataforma  :", R.version$platform, "\n")
cat("Sistema     :", Sys.info()[["sysname"]], Sys.info()[["release"]], "\n")
cat("Biblioteca  :", .libPaths()[1], "\n")

# getRversion() compara versões componente a componente. Não use
# as.numeric(paste0(major, ".", minor)): "4.10" vira 4.1, e um R 4.10 seria
# diagnosticado como anterior ao 4.3.
if (getRversion() < "4.3.0") {
  cat("\n⚠️  ATENÇÃO: seu R é anterior à versão 4.3.\n")
  cat("   O Bioconductor atual exige uma versão recente. Atualize antes de continuar:\n")
  cat("   https://cran.r-project.org/\n\n")
  cat("   No Windows, o pacote 'installr' facilita: install.packages('installr')\n\n")
}

# Se o script rodar fora do RStudio (Rscript), o repositório do CRAN pode estar
# como "@CRAN@" e a instalação trava esperando a escolha de um espelho.
.repos <- getOption("repos")
if (is.null(.repos[["CRAN"]]) || .repos[["CRAN"]] %in% c("", "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  cat("Espelho do CRAN definido: cloud.r-project.org\n")
}


# 2. Biblioteca gravável — e que continue existindo depois --------------------
# Causa nº 1 de falha em computador corporativo: a biblioteca padrão é
# somente leitura.

# A documentação do próprio R descreve file.access(mode = 2) como um guia
# aproximado no Windows, onde costuma dizer que a pasta é gravável quando não é.
# A única prova é tentar escrever.
eh_gravavel <- function(caminho) {
  if (!dir.exists(caminho)) return(FALSE)
  teste <- file.path(caminho, paste0(".teste_escrita_", Sys.getpid()))
  criou <- suppressWarnings(tryCatch(file.create(teste), error = function(e) FALSE))
  if (isTRUE(criou)) unlink(teste)
  isTRUE(criou)
}

# Grava R_LIBS_USER no ~/.Renviron. Sem isto, o .libPaths() abaixo vale só para
# ESTA sessão: no dia do workshop o aluno abre o RStudio, a biblioteca pessoal
# não está no caminho, e o library(DESeq2) falha como se nada tivesse sido
# instalado.
fixar_biblioteca <- function(caminho) {
  renviron <- path.expand("~/.Renviron")
  linha <- sprintf('R_LIBS_USER="%s"',
                   normalizePath(caminho, winslash = "/", mustWork = FALSE))
  anteriores <- if (file.exists(renviron)) readLines(renviron, warn = FALSE) else character(0)

  if (any(trimws(anteriores) == linha)) return(TRUE)   # já está lá

  # remove qualquer R_LIBS_USER antigo, para não ficarem duas definições
  anteriores <- anteriores[!grepl("^\\s*R_LIBS_USER\\s*=", anteriores)]
  tryCatch({ writeLines(c(anteriores, linha), renviron); TRUE },
           error = function(e) FALSE)
}

lib_fixada <- NA
if (!eh_gravavel(.libPaths()[1])) {
  cat("\n⚠️  Sua biblioteca de pacotes não é gravável.\n")
  cat("   Vamos criar uma biblioteca pessoal na sua pasta de usuário.\n")

  lib_pessoal <- path.expand("~/R/biblioteca")
  dir.create(lib_pessoal, recursive = TRUE, showWarnings = FALSE)

  if (!eh_gravavel(lib_pessoal)) {
    cat("\n❌ Também não consigo escrever em", lib_pessoal, "\n")
    cat("   Isto precisa do instrutor ou do suporte de TI. Pare por aqui.\n")
    stop("Sem biblioteca gravável.", call. = FALSE)
  }

  .libPaths(c(lib_pessoal, .libPaths()))
  cat("   Nova biblioteca:", lib_pessoal, "\n")

  lib_fixada <- fixar_biblioteca(lib_pessoal)
  if (isTRUE(lib_fixada)) {
    cat("   ✅ Gravado em ~/.Renviron — vale para as próximas sessões também.\n")
    cat("   ⚠️  FECHE E REABRA O RSTUDIO quando este script terminar.\n\n")
  } else {
    cat("   ⚠️  Não consegui gravar o ~/.Renviron. A biblioteca vale só nesta\n")
    cat("      sessão. No dia do workshop, rode ANTES de qualquer outra coisa:\n")
    cat(sprintf('        .libPaths(c("%s", .libPaths()))\n\n', lib_pessoal))
  }
}


# 3. Instalação ---------------------------------------------------------------

t0 <- Sys.time()

# --- CRAN ---
pacotes_cran <- c("BiocManager", "ggplot2", "ggrepel", "pheatmap", "jsonlite")
faltando <- pacotes_cran[!sapply(pacotes_cran, requireNamespace, quietly = TRUE)]

if (length(faltando)) {
  cat("\nInstalando do CRAN:", paste(faltando, collapse = ", "), "\n")
  install.packages(faltando, quiet = TRUE)
} else {
  cat("\nPacotes do CRAN: já instalados.\n")
}

# --- Bioconductor ---
# No Windows e no macOS isto baixa binários prontos e leva poucos minutos.
# No Linux pode compilar do código-fonte e demorar bem mais.
pacotes_bioc <- c("DESeq2", "apeglm")
faltando_bioc <- pacotes_bioc[!sapply(pacotes_bioc, requireNamespace, quietly = TRUE)]

if (length(faltando_bioc)) {
  cat("Instalando do Bioconductor:", paste(faltando_bioc, collapse = ", "), "\n")
  cat("(esta é a parte demorada — não interrompa)\n\n")

  if (Sys.info()[["sysname"]] == "Linux") {
    # Binários de CRAN para Linux, quando disponíveis para a distribuição
    codinome <- tryCatch({
      os <- readLines("/etc/os-release", warn = FALSE)
      linha <- grep("^VERSION_CODENAME=", os, value = TRUE)[1]
      if (is.na(linha)) NA_character_ else gsub('^VERSION_CODENAME=|"', "", linha)
    }, error = function(e) NA_character_)

    if (!is.na(codinome) && nzchar(codinome)) {
      # Mantém o CRAN na lista: o P3M sozinho deixaria o BiocManager sem
      # espelho de fallback se a distribuição não tiver binário do pacote.
      options(repos = c(
        P3M  = sprintf("https://packagemanager.posit.co/cran/__linux__/%s/latest",
                       codinome),
        CRAN = "https://cloud.r-project.org"))
      cat("Linux detectado — usando binários do P3M para:", codinome, "\n\n")
    }
  }

  BiocManager::install(faltando_bioc, ask = FALSE, update = FALSE)
} else {
  cat("Pacotes do Bioconductor: já instalados.\n")
}

tempo <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)


# 4. Relatório final ----------------------------------------------------------
# Se algo der errado, é este bloco que você manda para o instrutor.

cat("\n\n")
cat("=========================================================\n")
cat("  RELATÓRIO — copie tudo isto se precisar pedir ajuda\n")
cat("=========================================================\n")
cat("R          :", R.version.string, "\n")
cat("Plataforma :", R.version$platform, "\n")
cat("Biblioteca :", .libPaths()[1], "\n")
if (!is.na(lib_fixada)) {
  cat("~/.Renviron:", if (isTRUE(lib_fixada)) "gravado" else "FALHOU", "\n")
}
cat("Tempo      :", tempo, "min\n\n")

todos <- c(pacotes_cran, pacotes_bioc)
falhas <- character(0)

for (p in todos) {
  ok <- requireNamespace(p, quietly = TRUE)
  if (ok) {
    cat(sprintf("  %-14s ✅  %s\n", p, as.character(packageVersion(p))))
  } else {
    cat(sprintf("  %-14s ❌  NÃO INSTALADO\n", p))
    falhas <- c(falhas, p)
  }
}

cat("\n")
if (length(falhas) == 0) {
  # Teste real: uma análise minúscula do começo ao fim, INCLUINDO o lfcShrink
  # com apeglm — que é o que o M4 faz. Sem esta parte, um apeglm quebrado só
  # apareceria no meio da aula.
  teste <- tryCatch({
    suppressPackageStartupMessages(library(DESeq2))
    set.seed(1)
    cts <- matrix(rnbinom(600, mu = 100, size = 10), ncol = 6)
    rownames(cts) <- paste0("g", 1:100); colnames(cts) <- paste0("s", 1:6)
    cd <- data.frame(grupo = factor(rep(c("a", "b"), each = 3)),
                     row.names = colnames(cts))
    d <- DESeqDataSetFromMatrix(cts, cd, ~ grupo)
    d <- DESeq(d, quiet = TRUE)
    r <- results(d)

    coef_alvo <- resultsNames(d)[2]                      # a mesma chamada do M4
    rs <- lfcShrink(d, coef = coef_alvo, type = "apeglm", quiet = TRUE)

    nrow(r) == 100 && nrow(rs) == 100 && !all(is.na(rs$log2FoldChange))
  }, error = function(e) {
    cat("Erro no teste:", conditionMessage(e), "\n"); FALSE
  })

  if (isTRUE(teste)) {
    cat("✅ TUDO PRONTO. DESeq2 e apeglm rodaram uma análise de teste.\n")
    cat("   Você está preparado para o workshop. Até lá!\n")
    if (isTRUE(lib_fixada)) {
      cat("\n   ⚠️  Lembre: feche e reabra o RStudio antes do workshop.\n")
    }
  } else {
    cat("⚠️  Os pacotes instalaram, mas o teste falhou.\n")
    cat("   Mande este relatório para o instrutor.\n")
  }
} else {
  cat("❌ Falharam:", paste(falhas, collapse = ", "), "\n")
  cat("   Mande este relatório inteiro para o instrutor.\n")
  cat("   Não tente resolver sozinho no dia do workshop.\n")
}
cat("=========================================================\n")

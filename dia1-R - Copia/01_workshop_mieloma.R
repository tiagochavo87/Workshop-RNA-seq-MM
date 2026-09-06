# =============================================================================
# Workshop RNA-seq — Expressão diferencial em Mieloma Múltiplo
# Trilha R · DESeq2 (Bioconductor)
#
# Dado      : MMRF-CoMMpass (NCI Genomic Data Commons), STAR - Counts
# Contraste : mieloma ao diagnóstico  ×  mieloma na recidiva
#
# -----------------------------------------------------------------------------
# COMO USAR ESTE SCRIPT
#
#   Abra o arquivo workshop-mieloma.Rproj ANTES de abrir este script.
#   Isso garante que o R esteja na pasta certa. Se não abriu pelo .Rproj:
#   Session → Set Working Directory → To Source File Location
#
#   NÃO clique em "Source". Vamos executar bloco por bloco.
#
#   Atalhos do RStudio que valem a pena decorar hoje:
#     Ctrl+Enter        roda a linha (ou a seleção) onde está o cursor
#     Ctrl+Alt+T        roda a SEÇÃO inteira em que o cursor está  ← o principal
#     Ctrl+Shift+O      abre o índice de seções (use como roteiro da aula)
#     Alt+-             escreve  <-
#     F1 sobre a função abre a ajuda no painel da direita
#
#   Onde for para você mexer, está marcado com 🔧.
#
# ⚠️  Material didático. Demonstra o método; não é um estudo com rigor
#     metodológico. As limitações estão na seção M7 e devem ser lidas.
# =============================================================================


# M0 — Preparar o ambiente ----------------------------------------------------

# Se algum destes der erro, você pulou o 00_instalar_pacotes.R. Rode-o primeiro.
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(jsonlite)
})

# Formata p-valor sem notação científica: >= 0,001 mostra o valor, abaixo mostra
# o limite. Usada nas mensagens do M8.
fmt_padj <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0,001")
  sub("\\.", ",", sprintf("%.3f", p))
}

# Fecha qualquer device gráfico que tenha ficado aberto de uma execução
# anterior. Sem isto, um png() sem dev.off() faz a figura seguinte ser gravada
# no arquivo errado — e o dev.off() devolve "pdf 4" em vez de "null device 1".
graphics.off()

# Redesenha uma figura no painel Plots do RStudio. Se o painel estiver pequeno
# demais, o R devolve "figure margins too large" e PARA o script — por causa de
# uma cópia na tela, com o PNG já salvo em disco. Aqui isso vira um aviso.
na_tela <- function(expr) {
  tryCatch(expr, error = function(e) {
    cat("(figura não coube no painel Plots:", conditionMessage(e), ")\n")
    cat(" Aumente o painel arrastando a divisória, ou clique em Zoom.",
        "O PNG já foi salvo em figuras/.\n")
    try(while (dev.cur() > 1) dev.off(), silent = TRUE)   # não deixa device meio aberto
    invisible(NULL)
  })
}

# Pastas de saída, na estrutura do pipeline original
for (d in c("dados", "figuras", "resultados")) dir.create(d, showWarnings = FALSE)

cat("R      :", R.version.string, "\n")
cat("DESeq2 :", as.character(packageVersion("DESeq2")), "\n")
cat("apeglm :", as.character(packageVersion("apeglm")), "\n\n")
cat("Anote estas versões — elas entram na seção de métodos.\n")

# 👉 Olhe o painel Environment (canto superior direito). Está vazio.
#    Ao longo do dia ele vai encher de objetos, e você vai poder clicar em
#    cada um para ver o que tem dentro. Essa é a vantagem do RStudio sobre
#    um notebook: o estado da análise fica visível o tempo todo.


# 🔧 PAINEL DE CONTROLE -------------------------------------------------------
# A única seção que você mexe durante o workshop.
# Ao mudar um valor, rode desta seção para baixo (Ctrl+Alt+T em cada uma).

# --- Amostragem ---
N_POR_GRUPO <- 30        # amostras por grupo
SEED        <- 42        # semente da subamostragem (reprodutibilidade)

# --- Filtro de genes ---
MIN_CONTAGENS <- 10      # contagens mínimas...
MIN_AMOSTRAS  <- NA      # ...em quantas amostras? NA = tamanho do menor grupo

REMOVER_IG <- TRUE       # remove segmentos V(D)J de imunoglobulina/TCR.
                         # Cada clone de plasmócito tem seu próprio rearranjo:
                         # IGKV/IGHV/IGLV têm expressão idiossincrática por
                         # paciente e dominam o contraste por acaso. Em mieloma
                         # isto é praxe. Deixe FALSE uma vez para VER o efeito.

# --- Download ---
TENTATIVAS_DOWNLOAD <- 4    # o GDC corta conexões longas em Wi-Fi de auditório
TIMEOUT_DOWNLOAD    <- 600  # segundos por arquivo

# --- Significância ---
PADJ <- 0.05             # FDR - (Adjusted p-value) - Taxa de Falsa Descoberta (False Discovery Rate)
LFC  <- 1.0              # |log2 fold change| mínimo (1 = o dobro)
# Limiar de Fold Change em escala log2:
# 1.0  -> O gene precisa, no mínimo, duplicar (+2x) ou cair pela metade (0.5x)
# O uso de valor absoluto (| |) aplica o corte tanto para genes ativados quanto reprimidos
TESTAR_LIMIAR_NO_MODELO <- TRUE   # TRUE  -> results(lfcThreshold = LFC) : O teste estatístico prova 
                         #que a variação foi MAIOR que o limiar (LFC). É o método mais rigoroso.
                                  # FALSE -> testa contra zero e filtra depois :O teste prova 
#apenas que a variação é DIFERENTE DE ZERO. O corte biológico (LFC) é feito mecanicamente depois.

# --- Genes para destacar nos gráficos ---
GENES_INTERESSE <- c("MYC", "CCND1", "CCND2", "NSD2", "FGFR3", "MAF",
                     "TP53", "KRAS", "NRAS", "TNFRSF17", "SLAMF7", "CD38")

# --- Rótulos do contraste (confirme com a saída do M1) ---
MAPA_SAMPLE_TYPE <- c(
  "Primary Blood Derived Cancer - Bone Marrow"   = "diagnostico",
  "Recurrent Blood Derived Cancer - Bone Marrow" = "recidiva"
)
REFERENCIA <- "diagnostico"   # grupo de referência
ALVO       <- "recidiva"      # numerador do fold change

# --- 🔧 Release do GDC: PREENCHA antes de rodar ---
# Vai direto para o registro do M7. O GDC reprocessa os dados entre releases:
# sem esta linha, "baixado do GDC" não identifica nada.
# Veja em: gdc.cancer.gov/about-data/data-release
RELEASE_GDC <- "Data Release 46.0 - August 10, 2026"

# --- Origem dos dados ---
# TRUE  -> usa a matriz congelada do Zenodo, mesmo com o GDC no ar.
# FALSE -> tenta o GDC ao vivo e cai para o Zenodo se falhar.
FORCAR_CONGELADA <- TRUE

# --- Plano B: matriz congelada pelo instrutor ---
# Publicada no Zenodo com DOI. O "?download=1" é obrigatório: sem ele o
# servidor devolve a página HTML de preview, não o arquivo.
DOI_DADOS    <- "10.5281/zenodo.22399365"
URL_COUNTS   <- "https://zenodo.org/records/22399365/files/mm_counts.csv.gz?download=1"
URL_METADATA <- "https://zenodo.org/records/22399365/files/mm_metadata.csv?download=1"
URL_GENES    <- "https://zenodo.org/records/22399365/files/gene_names.csv.gz?download=1"

set.seed(SEED)
options(timeout = TIMEOUT_DOWNLOAD)

# Por que usar a matriz congelada?
# 1. Reprodutibilidade: O GDC atualiza dados, anotações de genes e versões de pipelines com frequência. 
#    A matriz do Zenodo garante que você e qualquer outra pessoa rodem exatamente com os mesmos números, hoje ou daqui a 5 anos.
# 2. Resiliência: Servidores públicos como o GDC caem, mudam de API ou sofrem lentidão. 
#    O arquivo congelado funciona offline e elimina gargalos na execução do script.


# Formata número no padrão brasileiro. Sem o decimal.mark explícito, o R avisa
# que big.mark e decimal.mark são ambos "." — que é justamente o caso aqui.
fmt <- function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE)

cat("Configuração carregada.\n")


# M1 — Onde mora o dado: a API do GDC -----------------------------------------
# O GDC guarda dados genômicos de câncer de vários projetos, processados por um
# pipeline único — e é essa padronização que torna as amostras comparáveis.
#
# Nossa consulta é um filtro estruturado:
#   projeto = MMRF-COMMPASS, tipo = expressão gênica, workflow = STAR-Counts,
#   acesso = aberto.

GDC_FILES <- "https://api.gdc.cancer.gov/files"
GDC_DATA  <- "https://api.gdc.cancer.gov/data"

FILTROS <- list(
  op = "and",
  content = list(
    list(op = "in", content = list(field = "cases.project.project_id",
                                   value = list("MMRF-COMMPASS"))),
    list(op = "in", content = list(field = "data_type",
                                   value = list("Gene Expression Quantification"))),
    list(op = "in", content = list(field = "analysis.workflow_type",
                                   value = list("STAR - Counts"))),
    list(op = "in", content = list(field = "access", value = list("open")))
  )
)
FILTROS_JSON <- toJSON(FILTROS, auto_unbox = TRUE)

gdc_get <- function(params, endpoint = GDC_FILES) {
  url <- paste0(endpoint, "?",
                paste(sprintf("%s=%s", names(params),
                              vapply(params, URLencode, character(1), reserved = TRUE)),
                      collapse = "&"))
  fromJSON(url, simplifyVector = FALSE)
}

# --- Que grupos existem neste dado? ---
# Antes de baixar qualquer coisa, olhe os metadados. Esta é a diferença entre
# fazer análise e apertar botões.

ONLINE <- !FORCAR_CONGELADA
# Controla a fonte dos dados: se FORCAR_CONGELADA for TRUE, ONLINE vira FALSE (usa o arquivo local).
tabela <- if (!ONLINE) NULL else tryCatch({
# Tenta buscar os dados direto do servidor do GDC via API
  r <- gdc_get(list(filters = FILTROS_JSON,
                    facets  = "cases.samples.sample_type",
                    size    = "0", format = "JSON"))
  b <- r$data$aggregations$`cases.samples.sample_type`$buckets
  # Usa unlist (e não vapply) porque o servidor do GDC às vezes devolve 
  # os números como inteiros e outras vezes como decimais (double).
  # O unlist aceita essa variação de tipo sem travar o código.
  data.frame(sample_type = unlist(lapply(b, `[[`, "key")),
             n           = as.integer(unlist(lapply(b, `[[`, "doc_count"))),
             stringsAsFactors = FALSE)
}, error = function(e) { ONLINE <<- FALSE; NULL })

if (ONLINE) {
  cat("Tipos de amostra disponíveis (arquivos de acesso aberto):\n\n")
  print(tabela, row.names = FALSE)
} else if (FORCAR_CONGELADA) {
  cat("FORCAR_CONGELADA = TRUE — usando a matriz congelada do Zenodo.\n")
  cat("Para consultar o GDC ao vivo, mude para FALSE no painel de controle.\n")
} else {
  cat("Não consegui falar com o GDC agora.\n")
  cat("Sem problema — vamos usar a matriz congelada pelo instrutor.\n")
}

# 💬 PARE E DISCUTA
#
#  - Existe algum grupo "normal" ou "tecido saudável"?  NÃO.
#    Todas as amostras são de pacientes com mieloma.
#
#  - Isso muda a pergunta possível. Não dá para perguntar "o que muda no
#    mieloma em relação ao plasmócito normal?". Dá para perguntar "o que muda
#    entre o diagnóstico e a recidiva?".
#
#  - Os nomes na coluna sample_type são os que devem estar em MAPA_SAMPLE_TYPE.
#    Se não baterem, corrija no PAINEL DE CONTROLE e rode de novo.
#
#  ⚠️ O CoMMpass é longitudinal: o mesmo paciente pode ter amostra nos dois
#     grupos. Amostras do mesmo paciente NÃO são independentes — isso viola uma
#     premissa do modelo. A seção seguinte mantém uma amostra por paciente.


# M1b — Selecionar as amostras ------------------------------------------------
# O CoMMpass é longitudinal: o mesmo paciente pode ter amostra no diagnóstico E
# na recidiva. Usar as duas violaria a premissa de independência do modelo. Por
# isso: uma amostra por paciente, e subamostragem balanceada entre os grupos.

# listar_arquivos(): pergunta ao GDC quais arquivos existem e devolve um por
# linha, já com a condição (diagnóstico/recidiva) traduzida.
listar_arquivos <- function() {
  # 1. Monta a lista de campos que queremos de volta na resposta da API.
  #    file_id identifica o arquivo; case_id/submitter_id identificam o paciente;
  #    sample_type é o rótulo que vira a condição do contraste.
  campos <- paste("file_id", "file_name", "cases.case_id", "cases.submitter_id",
                  "cases.samples.sample_type", "cases.samples.submitter_id",
                  sep = ",")

  # 2. Faz a chamada. size = "20000" garante que a resposta não venha paginada:
  #    pedimos todos os hits de uma vez.
  r <- gdc_get(list(filters = FILTROS_JSON, fields = campos,
                    format = "JSON", size = "20000"))

  # 3. Percorre cada arquivo retornado e extrai só o que interessa.
  linhas <- lapply(r$data$hits, function(h) {
    caso <- h$cases[[1]]; am <- caso$samples[[1]]
    # traduz o sample_type do GDC para "diagnostico"/"recidiva" via MAPA_SAMPLE_TYPE
    cond <- unname(MAPA_SAMPLE_TYPE[am$sample_type])
    # tipos fora do contraste (ex.: sangue periférico) viram NA e são descartados
    if (is.na(cond)) return(NULL)
    data.frame(
      file_id     = h$file_id,
      case_id     = if (!is.null(caso$submitter_id)) caso$submitter_id else caso$case_id,
      sample_type = am$sample_type,
      condition   = cond,
      stringsAsFactors = FALSE)
  })
  # 4. Empilha todas as linhas num único data.frame.
  do.call(rbind, linhas)
}

# selecionar(): do conjunto completo, escolhe as amostras que vão para a análise.
selecionar <- function(meta) {
  # 1. Embaralha as linhas. Sem isto, o passo seguinte manteria sempre a primeira
  #    amostra de cada paciente na ordem em que o GDC devolveu — arbitrário, não
  #    aleatório. É o embaralhamento que torna a escolha um sorteio de verdade.
  meta <- meta[sample(nrow(meta)), ]
  # 2. Mantém uma linha por paciente. Como já embaralhamos, qual amostra fica
  #    é sorteio.
  meta <- meta[!duplicated(meta$case_id), ]
  # 3. Dentro de cada condição, pega as N primeiras. O min() protege o caso de um
  #    grupo ter menos amostras que o pedido (a recidiva tem só 80 no total).
  partes <- lapply(split(meta, meta$condition),
                   function(b) head(b, min(N_POR_GRUPO, nrow(b))))
  # 4. Junta os grupos de volta num só data.frame.
  do.call(rbind, partes)
}

# Executa a seleção — só quando estamos no modo online (consultando o GDC).
# No modo congelado, as amostras já vêm definidas na matriz do Zenodo.
if (ONLINE) {
  todos <- listar_arquivos()
  cat("Arquivos elegíveis:", nrow(todos), "\n")
  print(table(todos$condition))          # quantos de cada condição existem no total

  meta <- selecionar(todos)
  rownames(meta) <- meta$file_id         # o file_id passa a identificar cada linha
  cat("\nApós 1 amostra/paciente e subamostragem:\n")
  print(table(meta$condition))           # confirma o balanceamento final
}


# M1c — Baixar os arquivos ----------------------------------------------------
# ~4 MB por arquivo. Com 60 amostras são poucos minutos.
#
# ⚠️ A ARMADILHA QUE ESTA SEÇÃO EVITA
# O GDC corta conexões longas — em Wi-Fi de auditório, com alunos baixando ao
# mesmo tempo, isso é regra e não exceção. Um download interrompido deixa um
# arquivo PARCIAL no disco. Se a função apenas perguntar "o arquivo existe?",
# ela reaproveita o lixo na execução seguinte: a matriz sai com NA, o colSums
# vira NA, o filtro de genes devolve um número absurdo, e o erro só aparece
# lá no DESeq2 — "NA values are not allowed in the count matrix" — a 200 linhas
# de distância da causa.
#
# Três defesas: baixa para arquivo temporário e só renomeia se estiver íntegro;
# valida pelo número de linhas; tenta de novo com espera crescente.

LINHAS_ESPERADAS <- 60665L  # GENCODE v36: comentário + cabeçalho + 4 N_* + 60.660 genes
TOLERANCIA_LINHAS <- 10L    # margem para variação entre releases do GDC

# arquivo_ok(): responde "este arquivo está íntegro?". Três checagens em cascata,
# da mais barata para a mais cara — só chega a contar linhas se passar nas duas
# primeiras.
arquivo_ok <- function(caminho) {
  if (!file.exists(caminho))     return(FALSE)   # 1. existe?
  if (file.size(caminho) < 3e6)  return(FALSE)   # 2. tem tamanho plausível (>3 MB)?
  # 3. tem o número esperado de linhas? Arquivo truncado falha aqui.
  n <- tryCatch(length(readLines(caminho, warn = FALSE)), error = function(e) 0L)
  abs(n - LINHAS_ESPERADAS) <= TOLERANCIA_LINHAS
}

# baixar(): baixa UM arquivo de contagem, com cache e retentativa.
baixar <- function(fid) {
  destino <- file.path("dados", paste0(fid, ".tsv"))
  # 1. Se já existe um arquivo íntegro em disco, reaproveita e sai. É o que torna
  #    o script retomável: rodar de novo não rebaixa o que já veio.
  if (arquivo_ok(destino)) return(TRUE)
  unlink(destino)                         # se existe mas está corrompido, descarta

  # 2. Tenta baixar até TENTATIVAS_DOWNLOAD vezes.
  for (k in seq_len(TENTATIVAS_DOWNLOAD)) {
    tmp <- paste0(destino, ".parcial")    # baixa para um nome temporário
    ok <- tryCatch({
      # mode = "wb" é obrigatório no Windows: sem ele o R corrompe o binário.
      download.file(paste0(GDC_DATA, "/", fid), destfile = tmp,
                    quiet = TRUE, mode = "wb", method = "libcurl")
      arquivo_ok(tmp)                      # o que baixou está íntegro?
    }, error = function(e) FALSE, warning = function(w) FALSE)

    # 3. Só renomeia para o nome definitivo se passou na validação. Assim o
    #    arquivo com nome final é SEMPRE íntegro — nunca um download pela metade.
    if (isTRUE(ok)) { file.rename(tmp, destino); return(TRUE) }
    unlink(tmp)                            # falhou: limpa o parcial
    # 4. Espera antes de tentar de novo: 2s, 4s, 8s. Bater de novo na hora num
    #    servidor que acabou de recusar raramente ajuda.
    if (k < TENTATIVAS_DOWNLOAD) Sys.sleep(2^k)
  }
  FALSE                                    # esgotou as tentativas
}

# Baixa todas as amostras, com uma barra de progresso, e junta as que falharem.
if (ONLINE) {
  ids    <- meta$file_id
  falhas <- character(0)                   # acumula os file_id que não baixaram
  pb <- txtProgressBar(min = 0, max = length(ids), style = 3)
  for (i in seq_along(ids)) {
    if (!baixar(ids[i])) falhas <- c(falhas, ids[i])
    setTxtProgressBar(pb, i)
  }
  close(pb); cat("\n")

  if (length(falhas) > 0) {
    cat(sprintf("⚠️  %d de %d arquivos falharam após %d tentativas.\n",
                length(falhas), length(ids), TENTATIVAS_DOWNLOAD))
    cat("Descartando essas amostras e seguindo com o restante.\n")
    meta <- meta[!meta$file_id %in% falhas, , drop = FALSE]
    rownames(meta) <- meta$file_id
    print(table(meta$condition))

    if (nrow(meta) == 0) {
      ONLINE <- FALSE
      cat("Nenhum arquivo baixado. Caindo para a matriz congelada.\n")
    } else if (min(table(meta$condition)) < 5) {
      stop("Menos de 5 amostras em um dos grupos. Rode de novo ou use a matriz congelada.")
    }
  } else {
    cat("Download concluído —", length(ids), "arquivos validados.\n")
  }
}


# M2 — Da contagem à matriz ---------------------------------------------------
# Antes de montar a matriz, ABRA UM ARQUIVO. Muita gente analisa RNA-seq a vida
# inteira sem nunca ter olhado o arquivo bruto.

# Dica rápida: Como abrir e olhar o arquivo bruto (ex: mm_counts.csv.gz ou um arquivo .tsv/.txt)
# Como ele está compactado (.gz), você pode abri-lo direto no R ou usar um editor de texto se for leve:

# Opção 1: No R (lê as primeiras 6 linhas sem carregar o arquivo inteiro na memória)
# Baixa temporariamente só para ver, sem salvar permanentemente na pasta dados/
temp_test <- tempfile(fileext = ".csv.gz")
download.file("https://zenodo.org/records/22399365/files/mm_counts.csv.gz?download=1", temp_test, mode = "wb")

# Lê as primeiras 6 linhas direto do arquivo temporário compactado
head(read.csv(gzfile(temp_test), nrows = 6))

# Opção 2: Inspecionando o texto bruto direto da URL (só é útil com a matriz congelada, 
# já que o if(ONLINE) original era exclusivo para os arquivos .tsv individuais do GDC)
# Opcional: Se quiser inspecionar o texto bruto diretamente da URL (sem pastas locais)
if (ONLINE) {
  url_congelada <- "https://zenodo.org/records/22399365/files/mm_counts.csv.gz?download=1"
  con <- gzcon(url(url_congelada, "rb"))
  writeLines(readLines(con, n = 10, warn = FALSE))
  close(con)
}

# ⚠️ ⚠️ ATENÇÃO: Dependendo de como você está rodando, existem duas realidades:
# 1. Se estivesse usando os arquivos .tsv individuais brutos do GDC: 
#    Eles contêm cabeçalhos de comentários do STAR, métricas de controle de 
#    qualidade no final (como N_unmapped, N_ambiguous) e várias colunas de 
#    estrandamento (unstranded, stranded_first, stranded_second).
# 
# 2. Como estamos usando a MATRIZ CONGELADA (mm_counts.csv.gz):
#    Ela já vem limpa e integrada! O arquivo contém apenas uma tabela pronta 
#    onde a primeira coluna é o 'gene_id' (com a versão do Ensembl) e as 
#    demais colunas representam diretamente os IDs de cada amostra/paciente, 
#    tendo apenas os valores numéricos de contagem.

# Qual coluna de contagem usar? Depende do kit da biblioteca. Regra prática do
# GDC: compare N_ambiguous entre as três; se uma das colunas stranded tem
# N_ambiguous muito menor, a biblioteca é stranded.
# ⚠️ Dados de strandness diferente NÃO se comparam diretamente.

# -----------------------------------------------------------------------------
# ATENÇÃO: O bloco abaixo (funções ler_star e montar_matriz) SÓ É USADO se 
# ONLINE = TRUE. Se você estiver usando a matriz congelada (ONLINE = FALSE / 
# FORCAR_CONGELADA = TRUE), o script pula tudo isso e baixa a matriz pronta direto.
# -----------------------------------------------------------------------------

USAR_COLUNA <- "unstranded"

ler_star <- function(caminho) {
  linhas <- readLines(caminho, warn = FALSE)
  idx <- grep("^gene_id", linhas)[1]
  if (is.na(idx)) stop("Cabeçalho 'gene_id' não encontrado em ", caminho)
  read.delim(text = paste(linhas[idx:length(linhas)], collapse = "\n"),
             stringsAsFactors = FALSE, check.names = FALSE)
}

montar_matriz <- function(file_ids) {
  mapa <- NULL; cols <- list()
  pb <- txtProgressBar(min = 0, max = length(file_ids), style = 3)
  for (i in seq_along(file_ids)) {
    df <- ler_star(file.path("dados", paste0(file_ids[i], ".tsv")))
    df <- df[!startsWith(as.character(df$gene_id), "N_"), ]   # fora as linhas N_*
    if (is.null(mapa)) {
      mapa <- unique(df[, c("gene_id", "gene_name", "gene_type")])
      rownames(mapa) <- mapa$gene_id
    }
    v <- as.integer(df[[USAR_COLUNA]]); names(v) <- df$gene_id
    cols[[file_ids[i]]] <- v[mapa$gene_id]
    setTxtProgressBar(pb, i)
  }
  close(pb)
  list(counts = do.call(cbind, cols), mapa = mapa)
}

if (ONLINE) {
  .m <- montar_matriz(meta$file_id)
  counts <- .m$counts
  mapa_genes <- .m$mapa
  rm(.m)
} else {
  if (!nzchar(URL_COUNTS)) stop("Sem GDC e sem URL de backup. Avise o instrutor.")

  # read.csv() NÃO descompacta .gz vindo de URL: o "?download=1" esconde a
  # extensão e o servidor manda os bytes crus. Sem isto, o erro é um
  # "invalid multibyte string" que não diz nada. Baixe e leia com gzfile().
  baixar_csv <- function(url, comprimido = grepl("\\.gz", url)) {
    tmp <- tempfile(fileext = if (comprimido) ".csv.gz" else ".csv")
    ok <- FALSE
    for (tentativa in seq_len(TENTATIVAS_DOWNLOAD)) {
      ok <- tryCatch({
        download.file(url, tmp, mode = "wb", quiet = TRUE)   # mode="wb": Windows
        file.exists(tmp) && file.size(tmp) > 0
      }, error = function(e) FALSE)
      if (ok) break
      Sys.sleep(2 * tentativa)
    }
    if (!ok) stop("Não consegui baixar: ", url)
    read.csv(if (comprimido) gzfile(tmp) else tmp,
             row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  }

  counts <- as.matrix(baixar_csv(URL_COUNTS))
  meta   <- baixar_csv(URL_METADATA)

  # O mapa de símbolos vem junto. Sem ele o REMOVER_IG não roda, os gráficos
  # ficam sem nome de gene e o controle positivo do M8 fica indisponível —
  # e o Módulo 9 do Dia 2 passa a comparar universos de teste diferentes.
  mapa_genes <- if (nzchar(URL_GENES)) baixar_csv(URL_GENES) else NULL
  if (!is.null(mapa_genes)) {
    mapa_genes$gene_id <- rownames(mapa_genes)
    faltando <- setdiff(rownames(counts), rownames(mapa_genes))
    if (length(faltando) > 0)
      warning(sprintf("%d gene(s) da matriz sem símbolo no mapa.", length(faltando)))
  }

  cat(sprintf("Usando a matriz congelada pelo instrutor (DOI %s).\n", DOI_DADOS))
  if (is.null(mapa_genes))
    cat("⚠️ Sem URL_GENES: símbolos indisponíveis. Avise o instrutor.\n")
}

# A checagem mais importante do script inteiro: amostras alinhadas?
comuns <- intersect(colnames(counts), rownames(meta))
counts <- counts[, comuns, drop = FALSE]
meta   <- meta[comuns, , drop = FALSE]
# Se passar em silêncio, deu tudo certo: as colunas da matriz batem exatamente com as linhas dos metadados.
# Se houver qualquer divergência (ex: amostras fora de ordem ou faltando), o stopifnot() para 
# a execução na hora, evitando erros silenciosos lá na frente na análise diferencial.
stopifnot(identical(colnames(counts), rownames(meta)))

n_na <- sum(is.na(counts))
if (n_na > 0) {
  ruins <- names(which(colSums(is.na(counts)) > 0))
  stop(sprintf(paste0("%s valores NA na matriz, em %d amostra(s).\n",
                      "Arquivos suspeitos (truncados no download):\n  %s\n",
                      "Apague-os da pasta 'dados/' e rode a seção M1c de novo."),
               fmt(n_na), length(ruins),
               paste(file.path("dados", paste0(ruins, ".tsv")), collapse = "\n  ")))
}

cat(sprintf("\nMatriz: %s genes x %d amostras — zero NA\n",
            fmt(nrow(counts)), ncol(counts)))
print(table(meta$condition))

# 👉 Rode  if (interactive()) View(meta)   # fora do RStudio, View() não existe  para abrir a tabela numa aba, com filtro e ordenação.
#    É assim que se olha um data.frame no RStudio — não com print().

View(meta)

# 💬 A matriz tem ~60 mil linhas, mas o genoma humano tem ~20 mil genes
#    codificantes. Por quê? Porque a anotação inclui pseudogenes, lncRNAs,
#    miRNAs e transcritos processados — veja a coluna gene_type. A maioria terá
#    contagem zero ou quase, e vai cair no filtro da próxima seção.


# M3 — Olhe o dado antes de testar o dado -------------------------------------
# A regra mais importante do workshop. Um outlier ou um efeito de lote
# encontrado aqui vale mais que qualquer refinamento estatístico depois.

# --- 3.1 Tamanho de biblioteca ---
# 1. Calcula o total de reads por amostra (soma de todas as colunas da matriz) 
# e divide por 1 milhão para facilitar a visualização no gráfico (milhões de reads).
bibl  <- colSums(counts) / 1e6

# 2. Organiza a ordem em que as amostras vão aparecer no eixo X do gráfico,
# agrupando-as por condição biológica (ex: controle vs. tratado) para facilitar a comparação.
ordem <- rownames(meta)[order(meta$condition)]

# 3. Define as cores de cada barra com base na condição da amostra:
# Usa uma cor para a condição de referência e outra para a condição alvo.
cores <- ifelse(meta[ordem, "condition"] == REFERENCIA, "#1b6ca8", "#d1495b")

# 4. Abre o dispositivo gráfico PNG para salvar a imagem gerada automaticamente na pasta 'figuras/'.
png("figuras/01_tamanho_biblioteca.png", width = 1400, height = 600, res = 130)

# 5. Plasta o gráfico de barras (barplot):
# - 'bibl[ordem]': plota os tamanhos na ordem agrupada pelas condições.
# - 'col = cores': pinta cada barra de acordo com o grupo biológico.
# - 'names.arg = rep("", length(ordem))': esconde os nomes das amostras no eixo X 
#   para evitar poluição visual (já que são muitas colunas/amostras juntas).
barplot(bibl[ordem], col = cores, border = NA, names.arg = rep("", length(ordem)),
        ylab = "Milhões de reads", xlab = "Amostras (agrupadas por condição)",
        main = "Tamanho da biblioteca por amostra")

# 6. Adiciona a legenda no canto superior direito para identificar o que cada cor representa,
# removendo a borda (border=NA) e o fundo da caixa (bty="n") para um visual mais limpo.
legend("topright", legend = c(REFERENCIA, ALVO), fill = c("#1b6ca8", "#d1495b"),
       border = NA, bty = "n")

# 7. Fecha o arquivo PNG, salvando o gráfico definitivamente no disco.
dev.off()

# Repete na tela (o painel Plots guarda o histórico — use as setas ← →)
na_tela({
  barplot(bibl[ordem], col = cores, border = NA, names.arg = rep("", length(ordem)),
          ylab = "Milhões de reads", xlab = "Amostras (agrupadas por condição)",
          main = "Tamanho da biblioteca por amostra")
  legend("topright", legend = c(REFERENCIA, ALVO), fill = c("#1b6ca8", "#d1495b"),
         border = NA, bty = "n")
})

# --- 3.2 Relatório no Console e Validação do Tamanho de Biblioteca ---
# 1. Imprime estatísticas resumidas das profundidades de sequenciamento:
#    - Mostra o tamanho da menor biblioteca (em milhões de reads)
#    - Mostra o tamanho da maior biblioteca
#    - Calcula a razão máx/mín para ver a amplitude de variação entre as amostras
cat(sprintf("\nMenor: %.1f M   Maior: %.1f M   Razão máx/mín: %.1fx\n",
            min(bibl), max(bibl), max(bibl) / min(bibl)))

# 2. Trava de segurança contra dados corrompidos ou valores ausentes:
#    Se qualquer amostra resultou em NA na soma, o script para imediatamente 
#    com uma mensagem orientando a voltar ao passo de montagem da matriz.
if (any(is.na(bibl))) stop("NA no tamanho de biblioteca — a matriz tem buracos. Volte ao M1c.")

# 3. Orientação prática de controle de qualidade (QC):
#    Até 3x de diferença entre a menor e a maior biblioteca é normal/confortável.
#    Diferenças extremas (razões muito altas) podem introduzir viés de contagem 
#    e exigir atenção extra na normalização do DESeq2.
cat("Regra prática: até ~3x é confortável. Muito acima disso, desconfie.\n")

# M3b — 🔧 Filtro de genes: a primeira decisão real ---------------------------
# Genes com quase nenhuma contagem não têm poder para detectar nada, e ainda
# aumentam o peso da correção para testes múltiplos.

n_menor      <- min(table(meta$condition))
min_amostras <- if (is.na(MIN_AMOSTRAS)) n_menor else MIN_AMOSTRAS

keep_A <- rowSums(counts) >= MIN_CONTAGENS                    # soma total
keep_B <- rowSums(counts >= MIN_CONTAGENS) >= min_amostras    # presença consistente

cat(sprintf("Genes totais                          : %s\n", fmt(nrow(counts))))
cat(sprintf("A) soma total >= %-3d                  : %s\n", MIN_CONTAGENS, fmt(sum(keep_A))))
cat(sprintf("B) >= %d contagens em >= %d amostras   : %s\n", MIN_CONTAGENS, min_amostras, fmt(sum(keep_B))))
cat(sprintf("\nDiferença: %s genes que o critério A deixa passar\n",
            fmt(sum(keep_A) - sum(keep_B))))
cat("(um gene com 10 reads numa única amostra passa em A, não passa em B)\n")

# B é o critério do filterByExpr (edgeR): exige presença consistente, não um
# pico isolado. Não existe filtro universalmente certo — existe filtro
# DECLARADO, com o número de genes antes e depois escrito no método.
counts_f <- counts[keep_B, , drop = FALSE]
cat(sprintf("\nSeguindo com o critério B: %s genes.\n", fmt(nrow(counts_f))))


# M3b2 — Segmentos V(D)J de imunoglobulina ----
# Cada paciente de mieloma tem um clone de plasmócitos com rearranjo V(D)J
# ÚNICO. Os segmentos IGKV/IGHV/IGLV ficam com expressão altíssima e totalmente
# idiossincrática por amostra. Com n=30 por grupo, eles dominam o topo da lista
# por acaso — não por biologia de recidiva.
#
# 💬 Rode uma vez com REMOVER_IG <- FALSE e olhe o volcano. Depois volte para
#    TRUE. A diferença entre os dois gráficos é a aula inteira.

IG_PADRAO <- "^(IGKV|IGHV|IGLV|IGKJ|IGHJ|IGLJ|IGKC|IGHG|IGHA|IGHM|IGHD|IGLC|TRBV|TRAV|TRGV|TRDV)"

if (REMOVER_IG && !is.null(mapa_genes)) {
  simbolos <- mapa_genes[rownames(counts_f), "gene_name"]
  eh_ig    <- !is.na(simbolos) & grepl(IG_PADRAO, simbolos)
  cat(sprintf("Segmentos V(D)J removidos: %s genes\n", fmt(sum(eh_ig))))
  counts_f <- counts_f[!eh_ig, , drop = FALSE]
  cat(sprintf("Matriz final para o teste: %s genes.\n", fmt(nrow(counts_f))))
} else if (REMOVER_IG) {
  cat("⚠️ REMOVER_IG é TRUE mas mapa_genes não existe (matriz congelada sem símbolos).\n")
  cat("   Os segmentos de Ig vão para o teste — espere IGKV/IGHV no topo da lista.\n")
}


# M3c — O objeto DESeqDataSet -------------------------------------------------
# Guarda contagens, metadados e desenho experimental juntos. É a peça central
# da análise, e a razão de o DESeq2 ser difícil de usar errado.
#
# A ordem dos níveis do fator importa: o PRIMEIRO nível é a referência.

# 1. Transforma condition em fator com ordem EXPLÍCITA. O primeiro nível é a
#    referência (o denominador do fold change). Inverter aqui inverte o sinal de
#    todos os log2FC — e nada no resultado avisaria que isso aconteceu.
meta$condition <- factor(meta$condition, levels = c(REFERENCIA, ALVO))

# 2. Cria o objeto que amarra contagens + metadados + desenho. design = ~ condition
#    diz ao DESeq2 que queremos modelar a expressão em função da condição.
dds <- DESeqDataSetFromMatrix(countData = counts_f,
                              colData   = meta,
                              design    = ~ condition)

# 3. Confirma qual grupo ficou como referência. Vale conferir sempre.
cat("Referência do contraste:", levels(dds$condition)[1], "\n")
dds

# 👉 Veja as contagens dentro do objeto:  head(counts(dds), 10)
#    Ou clique na setinha azul ao lado de 'dds' no painel Environment para
#    navegar a estrutura interna completa.
head(counts(dds), 10)


# M3d — Transformação e PCA ---------------------------------------------------
# Contagem bruta tem variância dependente da média: gene muito expresso varia
# muito em valor absoluto. PCA e heatmap sobre o dado cru ficam dominados pelos
# genes mais expressos. A VST corrige isso — e NÃO é log2(x + 1): ela é
# calculada a partir da relação média-variância do próprio dado.

# PARA QUE SERVE O PCA AQUI
# O PCA pega as ~19.000 dimensões (uma por gene) e resume cada amostra em dois
# números — PC1 e PC2 — que capturam a maior parte da variação entre elas. Serve
# para responder, de um olhar, três perguntas ANTES de qualquer teste:
#
#   1. As amostras do mesmo grupo se parecem entre si? Se diagnóstico e recidiva
#      formam duas nuvens separadas, há sinal forte. Se estão misturadas, o
#      efeito é sutil (é o caso deste dado — veja a nota no fim da seção).
#
#   2. Existe algum OUTLIER? Uma amostra isolada, longe de todas as outras,
#      costuma ser problema técnico — e é melhor descobrir agora do que depois.
#
#   3. O que MAIS separa as amostras é a condição, ou é outra coisa? Se o PC1
#      separar por sexo, por lote ou por data de processamento, esse fator é um
#      confundidor e precisa entrar no modelo como covariável.
#
# O PCA não testa hipótese nenhuma — ele é o mapa que se olha antes de viajar.

# 1. Aplica a VST. blind = FALSE permite que a transformação use o desenho na
#    estimativa — adequado para exploração. (blind = TRUE seria para QC cego ao
#    desenho.)
vsd <- vst(dds, blind = FALSE)

# 2. Calcula o PCA sobre os 2000 genes mais variáveis. returnData = TRUE devolve
#    as coordenadas para plotarmos com o ggplot.
p  <- plotPCA(vsd, intgroup = "condition", ntop = 2000, returnData = TRUE)
# 3. Extrai a % da variância explicada por PC1 e PC2, para os rótulos dos eixos.
pv <- round(100 * attr(p, "percentVar"))

g_pca <- ggplot(p, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4, alpha = .85) +
  scale_color_manual(values = setNames(c("#1b6ca8", "#d1495b"), c(REFERENCIA, ALVO))) +
  labs(x = paste0("PC1 (", pv[1], "% da variância)"),
       y = paste0("PC2 (", pv[2], "% da variância)"),
       title = "PCA — VST, top 2000 genes variáveis", color = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

print(g_pca)
ggsave("figuras/02_pca.png", g_pca, width = 7, height = 6, dpi = 200)

# 💬 PARE E DISCUTA — e não se assuste
#
#  Os grupos separam? Provavelmente NÃO, ou pouco. E está tudo bem.
#
#  Mieloma tem subtipos citogenéticos muito distintos — t(4;14), t(11;14),
#  hiperdiploidia. A variação ENTRE PACIENTES é maior que a variação entre
#  diagnóstico e recidiva, e o PC1 quase sempre capta o subtipo, não o momento
#  da doença.
#
#  Duas lições:
#   1. PCA que não separa ≠ ausência de genes diferenciais. Significa que o
#      efeito é menor que a variação entre indivíduos. O modelo ainda pode
#      encontrá-lo, porque testa gene a gene em vez de olhar a variância global.
#   2. Se o PC1 separasse por SEXO ou por data de processamento, seria
#      confundidor, e entraria no modelo como covariável (~ sexo + condition).


# M4 — Expressão diferencial --------------------------------------------------
# O que DESeq() faz, em quatro etapas:
#
#  1. NORMALIZAÇÃO POR MEDIANA DE RAZÕES  (estimateSizeFactors)
#     Não é dividir pelo total de reads. Para cada gene, calcula-se a razão em
#     relação à média geométrica entre amostras; o fator de tamanho da amostra
#     é a MEDIANA dessas razões.
#     Isso importa MUITO em mieloma: plasmócito é uma fábrica de anticorpo, e
#     IGHG1, IGKC e IGLC1 podem responder sozinhos por uma fatia enorme da
#     biblioteca — fatia que varia conforme o isotipo do mieloma. Com
#     normalização por soma total, esses poucos genes sequestrariam a escala de
#     todas as amostras. A mediana é imune a eles.
#
#  2. DISPERSÃO COM SHRINKAGE  (estimateDispersions)
#     Genes de expressão parecida emprestam informação uns aos outros. É o que
#     torna o método utilizável com poucas réplicas.
#
#  3. TESTE DE WALD  (nbinomWaldTest), gene a gene.
#
#  4. CORREÇÃO BH + FILTRAGEM INDEPENDENTE, dentro de results().

dds <- DESeq(dds)

sf <- sizeFactors(dds)
cat(sprintf("Fatores de tamanho: min %.2f, mediana %.2f, máx %.2f\n\n",
            min(sf), median(sf), max(sf)))

png("figuras/03_dispersao.png", width = 1200, height = 900, res = 130)
plotDispEsts(dds, main = "Estimativa de dispersão")
dev.off()
na_tela(plotDispEsts(dds, main = "Estimativa de dispersão"))

cat("Como ler o gráfico de dispersão:\n")
cat("  pontos pretos = estimativa por gene\n")
cat("  linha vermelha = tendência ajustada\n")
cat("  pontos azuis  = estimativa final, puxada em direção à tendência\n")
cat("Os azuis são o shrinkage da etapa 2 acontecendo.\n")


# M4b — O teste, e onde entra o limiar de fold change -------------------------
# Filtrar padj < 0.05 & |LFC| >= 1 DEPOIS do teste não é o mesmo que testar a
# hipótese "o efeito é maior que 1". O primeiro não controla o FDR para a
# pergunta que você está fazendo. O DESeq2 permite declarar o limiar DENTRO do
# teste, com lfcThreshold.

if (TESTAR_LIMIAR_NO_MODELO) {
  # Caminho recomendado: o limiar entra DENTRO do teste. A hipótese nula passa a
  # ser "o efeito é <= LFC", e o padj já responde à pergunta que interessa.
  # contrast = c("condition", ALVO, REFERENCIA) fixa a direção: ALVO / REFERENCIA.
  res  <- results(dds, contrast = c("condition", ALVO, REFERENCIA),
                  alpha = PADJ, lfcThreshold = LFC, altHypothesis = "greaterAbs")
  modo <- sprintf("H0: |log2FC| <= %.1f (limiar dentro do teste)", LFC)
} else {
  # Caminho comum (e pior): testa contra zero e filtra por |LFC| depois. O padj
  # deixa de controlar o FDR da afirmação que você acaba fazendo.
  res  <- results(dds, contrast = c("condition", ALVO, REFERENCIA), alpha = PADJ)
  modo <- "H0: log2FC = 0 (limiar aplicado depois — o jeito comum)"
}

cat("Modo de teste:", modo, "\n\n")
summary(res)

# O summary(res) acima termina com duas notas de rodapé, [1] e [2]. Não são erro
# — explicam as duas causas de um gene receber padj = NA:
#
#   [1] cooksCutoff        — o gene tem um OUTLIER. A distância de Cook detectou
#                            uma amostra com contagem tão extrema que distorceria
#                            o resultado. O gene é retirado do teste.
#   [2] independentFiltering — o gene tem EXPRESSÃO BAIXA DEMAIS. Sem poder
#                            estatístico, ele só pesaria na correção múltipla, e
#                            é descartado automaticamente.
#
# São esses dois mecanismos que explicam a diferença entre "genes na matriz" e
# "genes testados". Para ajustá-los: ?results, argumentos cooksCutoff e
# independentFiltering.


# 👉 Aperte F1 sobre  results  para ver todos os argumentos. O painel de ajuda
#    do RStudio é onde mora a documentação real do DESeq2.


# M4c — O encolhimento do fold change -----------------------------------------
# Um gene com 12 reads pode ter log2FoldChange de +8. É ruído. lfcShrink() com
# apeglm encolhe o efeito na proporção da incerteza da estimativa.
#
# REGRA: teste com o LFC bruto; ranqueie e plote com o LFC encolhido.
#
# Nota: apeglm trabalha sobre um COEFICIENTE do modelo, não sobre um contrast
# arbitrário — por isso usamos resultsNames() para achar o nome certo.

# 1. Lista os coeficientes do modelo. O apeglm precisa do NOME de um coeficiente,
#    não de um contraste — por isso não reaproveitamos o contrast do M4b.
nomes <- resultsNames(dds)
cat("Coeficientes disponíveis:\n"); print(nomes)

# 2. Acha, na lista acima, o nome do coeficiente que corresponde ao nosso
#    contraste — é o que contém ALVO no nome (algo como "condition_recidiva_vs_
#    diagnostico"). Buscar com grep em vez de escrever o nome à mão evita que o
#    script quebre quando o DESeq2 muda esses nomes entre versões.
#    O [1] pega o primeiro resultado, caso haja mais de um.
coef_alvo <- grep(ALVO, nomes, value = TRUE)[1]
cat("\nUsando o coeficiente:", coef_alvo, "\n")

# 3. Aplica o encolhimento (shrinkage) a esse coeficiente. O apeglm reavalia
#    cada fold change à luz da incerteza da estimativa: gene bem medido quase
#    não muda; gene com poucas contagens (estimativa incerta) é puxado para
#    perto de zero. quiet = TRUE só silencia as mensagens de progresso.
res_shrunk <- lfcShrink(dds, coef = coef_alvo, type = "apeglm", quiet = TRUE)

# 4. Monta uma tabelinha para VER o encolhimento acontecendo. Junta, para cada
#    gene, a expressão média (baseMean), o fold change bruto e o encolhido.
comp <- data.frame(baseMean      = res$baseMean,
                   LFC_bruto     = res$log2FoldChange,
                   LFC_encolhido = res_shrunk$log2FoldChange)

# 5. Ordena por baseMean (do menor para o maior) e mostra os 10 genes de MENOR
#    expressão. É neles que o efeito é mais dramático: compare as colunas
#    LFC_bruto e LFC_encolhido — o bruto tem valores grandes e instáveis, o
#    encolhido está quase todo perto de zero. É o apeglm dizendo "não confio
#    nesse efeito, há contagem de menos para sustentá-lo".
comp <- comp[order(comp$baseMean), ]
cat("\nEfeito do shrinkage nos genes de MENOR expressão:\n")
print(round(head(comp, 10), 2))

# M4d — A tabela de resultados ------------------------------------------------
# padj vem do teste (com limiar); log2FoldChange vem do objeto encolhido.

# Monta a tabela final juntando as duas análises. As duas colunas de fold change
# ficam lado a lado DE PROPÓSITO: o p vem do bruto (que foi testado), o valor que
# se reporta e se plota é o encolhido.
tab <- data.frame(
  gene_id        = rownames(res),
  baseMean       = res$baseMean,                # expressão média do gene
  log2FC_bruto   = res$log2FoldChange,          # o que FOI TESTADO
  log2FoldChange = res_shrunk$log2FoldChange,   # o que se REPORTA e se plota
  lfcSE          = res_shrunk$lfcSE,            # erro-padrão do LFC encolhido
  pvalue         = res$pvalue,
  padj           = res$padj,                    # p-valor corrigido (BH)
  stringsAsFactors = FALSE
)

# Anexa símbolo e tipo de cada gene, se o mapa estiver disponível. Sem ele,
# usamos o próprio Ensembl ID como nome.
if (!is.null(mapa_genes)) {
  tab$gene_name <- mapa_genes[tab$gene_id, "gene_name"]
  tab$gene_type <- mapa_genes[tab$gene_id, "gene_type"]
} else {
  tab$gene_name <- tab$gene_id
}

# Ordena por padj (mais significativos no topo) e limpa os nomes de linha.
tab <- tab[order(tab$padj), ]
rownames(tab) <- NULL

# Espia rápida: as 6 primeiras linhas e as dimensões, só para confirmar que a
# tabela foi montada como esperado. A visão completa vem logo abaixo, com o
# resumo dos DEGs e o View(tab).
cat(sprintf("Tabela montada: %s genes x %d colunas\n", fmt(nrow(tab)), ncol(tab)))
head(tab)

# --- O QUE CONTA COMO SIGNIFICATIVO ---
# Isto depende de COMO o limiar foi aplicado, e a diferença não é cosmética.
#
#   limiar DENTRO do teste  -> o padj já responde "|LFC| > limiar?".
#                              Filtrar por |LFC| de novo é redundante, e ainda
#                              por cima usaria o LFC ENCOLHIDO, que não foi o
#                              testado.
#   limiar DEPOIS do teste  -> o padj responde "LFC ≠ 0?". Aí o filtro por
#                              tamanho de efeito é necessário — mas o padj
#                              deixa de controlar o FDR da afirmação que você
#                              está fazendo, porque foi calculado para outra
#                              hipótese. É o jeito comum e é o pior.
if (TESTAR_LIMIAR_NO_MODELO) {
  sig <- subset(tab, !is.na(padj) & padj < PADJ)
  criterio <- sprintf("padj < %.2f (o limiar de %.1f já está no teste)", PADJ, LFC)
} else {
  sig <- subset(tab, !is.na(padj) & padj < PADJ & abs(log2FoldChange) >= LFC)
  criterio <- sprintf("padj < %.2f E |log2FC| >= %.1f (filtro post-hoc)", PADJ, LFC)
}
testados <- sum(!is.na(tab$padj))

cat(sprintf("Genes testados                  : %s\n", fmt(testados)))
cat(sprintf("Genes com padj = NA (filtrados) : %s\n", fmt(sum(is.na(tab$padj)))))
cat(sprintf("Critério                        : %s\n", criterio))
cat(sprintf("DEGs                            : %s\n", fmt(nrow(sig))))
cat(sprintf("   aumentados em '%s' : %d\n", ALVO, sum(sig$log2FoldChange > 0)))
cat(sprintf("   reduzidos  em '%s' : %d\n", ALVO, sum(sig$log2FoldChange < 0)))
cat(sprintf("\nProporção do transcriptoma testado: %.3f%%  (%d de %s)\n",
            100 * nrow(sig) / max(testados, 1), nrow(sig), fmt(testados)))

# --- Confere com o summary(res) lá de cima? ---
# Com o limiar dentro do teste, tem que bater: os dois usam só o padj.
post_hoc <- subset(tab, !is.na(padj) & padj < PADJ & abs(log2FoldChange) >= LFC)
perdidos <- setdiff(sig$gene_id, post_hoc$gene_id)

cat(sprintf("\nDEGs pelo critério declarado        : %d\n", nrow(sig)))
cat(sprintf("DEGs se filtrássemos |LFC| de novo   : %d\n", nrow(post_hoc)))
if (length(perdidos) > 0) {
  cat(sprintf("\n  %d gene(s) sumiriam no filtro redundante:\n", length(perdidos)))
  print(tab[tab$gene_id %in% perdidos,
            c("gene_name", "baseMean", "log2FC_bruto", "log2FoldChange", "padj")],
        row.names = FALSE, digits = 3)
  cat("\n  💬 O teste disse que o efeito destes genes excede o limiar. O\n")
  cat("     encolhimento puxou a ESTIMATIVA para perto de zero, mas isso é\n")
  cat("     outra pergunta. Filtrar aqui seria descartar o resultado do teste\n")
  cat("     que você mesmo escolheu fazer.\n")
}

cat("\n--- Top 15 por padj ---\n")
print(head(tab[, c("gene_name", "baseMean", "log2FC_bruto", "log2FoldChange", "padj")], 15),
      digits = 3)

# 👉 Agora rode:  View(tab)
#    Abre a tabela numa aba, com busca e ordenação por coluna. É a melhor
#    forma de explorar a lista — clique no cabeçalho de padj, depois no de
#    log2FoldChange, e veja como a resposta muda conforme o critério.

View(tab)

# 💬 O TESTE DE SANIDADE MAIS ÚTIL QUE EXISTE
#
#  Olhe a PROPORÇÃO do transcriptoma que saiu significativa.
#    < 5%   plausível para um contraste sutil como este
#    > 25%  bandeira vermelha
#
#  Quando um quarto do transcriptoma é "diferencial", o contraste quase sempre
#  está capturando diferença de TECIDO ou de COMPOSIÇÃO CELULAR, não regulação
#  gênica fina. (No pipeline de colangiocarcinoma, tumor × ducto biliar normal
#  deu 10.774 DEGs em 38.846 genes — 28%. É um resultado real, mas o que ele
#  mede sobretudo é que fígado tumoral não é ducto biliar.)
#
#  E os genes com padj = NA não são erro: ou foram removidos pela filtragem
#  independente (expressão baixa demais para ter poder), ou marcados como
#  outlier pela distância de Cook. Quando reportar "X genes testados", o número
#  certo é o de genes com padj não nulo — não o de linhas da matriz.

# 🔧 EXPERIMENTE AGORA
#   1. Volte ao PAINEL DE CONTROLE, mude PADJ de 0.05 para 0.01 e rode do M4
#      para baixo. Quantos DEGs sobraram?
#   2. Mude TESTAR_LIMIAR_NO_MODELO para FALSE. O número mudou muito?


# M5 — Visualização honesta ---------------------------------------------------

# --- 5.1 MA plot: o mais informativo e o menos usado ---
# O MA plot mostra o fold change (eixo y) contra a expressão média (eixo x).
# É onde o efeito do encolhimento fica visível: comparamos os dois lado a lado.

# 1. Abre o PNG e divide a área em dois painéis (1 linha, 2 colunas).
png("figuras/04_ma_plot.png", width = 1600, height = 700, res = 130)
par(mfrow = c(1, 2))
# 2. Painel da esquerda: fold change BRUTO. Repare na nuvem larga à esquerda
#    (genes de baixa expressão com efeitos grandes e sem sustentação).
plotMA(res,        ylim = c(-6, 6), main = "LFC bruto")
# 3. Painel da direita: fold change ENCOLHIDO. A nuvem da esquerda colapsa em
#    direção a zero — é o apeglm corrigindo as estimativas mal medidas.
plotMA(res_shrunk, ylim = c(-6, 6), main = "LFC encolhido (apeglm)")
# 4. Restaura o layout de painel único e fecha o arquivo.
par(mfrow = c(1, 1))
dev.off()

na_tela({
  par(mfrow = c(1, 2))
  plotMA(res,        ylim = c(-6, 6), main = "LFC bruto")
  plotMA(res_shrunk, ylim = c(-6, 6), main = "LFC encolhido (apeglm)")
  par(mfrow = c(1, 1))
})

cat("Compare os dois painéis. À esquerda (baixa expressão), o bruto tem uma\n")
cat("nuvem larga de efeitos grandes e sem sustentação. Depois do shrinkage ela\n")
cat("colapsa em direção a zero. A nuvem deve estar centrada em zero — é o sinal\n")
cat("de que a normalização funcionou.\n")


# M5b — Volcano ---------------------------------------------------------------
#
# ⚠️ O ERRO MAIS COMUM EM VOLCANO, E QUASE NINGUÉM PERCEBE:
#    plotar o log2FC ENCOLHIDO no eixo x contra o padj que veio do teste do
#    log2FC BRUTO. São duas grandezas de contas diferentes. O resultado é um
#    gráfico onde um gene pode aparecer altíssimo (padj minúsculo) e colado no
#    zero (LFC encolhido) — uma posição que nenhum gene deveria poder ocupar.
#
#    Aqui isso acontece com o FGFR3. Em vez de esconder, mostramos os dois
#    painéis lado a lado: é a melhor figura da aula.

# volcano(): desenha um volcano. Recebe qual coluna de fold change usar no eixo x
# (bruto ou encolhido), o título e o subtítulo. Devolve o objeto ggplot.
volcano <- function(coluna_lfc, titulo, subtitulo) {
  # 1. Descarta genes sem padj ou sem fold change (não têm lugar no gráfico).
  d <- tab[!is.na(tab$padj) & !is.na(tab[[coluna_lfc]]), ]
  # 2. Eixo x = fold change escolhido. Eixo y = -log10(padj): quanto menor o
  #    padj, mais alto o ponto. O pmax evita log10(0) = infinito.
  d$x <- d[[coluna_lfc]]
  d$y <- -log10(pmax(d$padj, 1e-300))

  rot_up   <- paste("aumentado em", ALVO)
  rot_down <- paste("reduzido em",  ALVO)
  # A cor segue o MESMO critério de significância da tabela. Com o limiar no
  # teste, o padj basta e a direção vem do sinal do LFC.
  sig_aqui <- if (TESTAR_LIMIAR_NO_MODELO) d$padj < PADJ
              else d$padj < PADJ & abs(d$x) >= LFC
  d$grupo  <- ifelse(sig_aqui & d$x > 0, rot_up,
              ifelse(sig_aqui & d$x < 0, rot_down, "ns"))

  # 3. Separa os genes de interesse, que ganharão rótulo e círculo preto.
  hl <- d[d$gene_name %in% GENES_INTERESSE, ]

  # 4. Monta o gráfico camada por camada.
  ggplot(d, aes(x, y, color = grupo)) +
    geom_point(size = .9, alpha = .6) +                    # todos os genes
    scale_color_manual(values = setNames(c("#d1495b", "#1b6ca8", "#d5d5d5"),
                                         c(rot_up, rot_down, "ns")),
                       drop = FALSE) +
    geom_point(data = hl, shape = 21, size = 3.2, stroke = 1.1,
               color = "black", fill = NA) +
    geom_text_repel(data = hl, aes(label = gene_name), color = "black",
                    fontface = "bold", size = 3.2, max.overlaps = 20,
                    min.segment.length = 0, box.padding = .45) +
    geom_hline(yintercept = -log10(PADJ), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed", color = "grey50") +
    labs(x = sprintf("log2 Fold Change (%s / %s)", ALVO, REFERENCIA),
         y = "-log10 (padj)", color = NULL,
         title = titulo, subtitle = subtitulo) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(), legend.position = "top")
}

# --- O volcano correto: eixo x e eixo y vindos do MESMO cálculo ---
g_volc <- volcano("log2FC_bruto",
                  sprintf("Mieloma múltiplo — %s vs %s (MMRF-CoMMpass)", ALVO, REFERENCIA),
                  "log2FC bruto — a mesma grandeza que o teste avaliou")
print(g_volc)
ggsave("figuras/05_volcano.png", g_volc, width = 9, height = 7.5, dpi = 200)

# --- O mesmo gráfico com o LFC encolhido, para comparar ---
g_volc_shrunk <- volcano("log2FoldChange",
                         "O mesmo contraste, com o log2FC encolhido",
                         "eixo x encolhido, eixo y do teste bruto — repare no FGFR3")
print(g_volc_shrunk)
ggsave("figuras/05b_volcano_encolhido.png", g_volc_shrunk, width = 9, height = 7.5, dpi = 200)

cat("Compare 05_volcano.png com 05b_volcano_encolhido.png.\n")
cat("No primeiro, o FGFR3 aparece longe do zero: é o efeito que o teste viu.\n")
cat("No segundo, ele desliza para x = 0 e fica pendurado no alto, sozinho —\n")
cat("porque o eixo y continua vindo do teste bruto. Essa posição impossível é\n")
cat("a assinatura visual de um gene cujo efeito não se sustenta.\n\n")
cat("💬 Qual dos dois você poria no artigo? E qual você poria no suplementar?\n")


# M5c — Heatmap, e o alerta mais importante da aula ---------------------------
#
# ⚠️ Um heatmap dos "top N genes por padj" com clustering das amostras separa os
#    grupos POR CONSTRUÇÃO. Você selecionou os genes PORQUE eles separam os
#    grupos, e depois mostrou que eles separam os grupos. É raciocínio circular
#    — e é o gráfico mais comum e mais enganoso em artigo de expressão.
#
# O heatmap abaixo é demonstração. Leia o comentário no fim da seção.

top <- head(subset(tab, !is.na(padj)), 30)
mat <- assay(vsd)[top$gene_id, , drop = FALSE]
rownames(mat) <- ifelse(is.na(top$gene_name) | top$gene_name == "",
                        top$gene_id, top$gene_name)

ord  <- rownames(meta)[order(meta$condition)]
anot <- data.frame(condicao = meta[ord, "condition"], row.names = ord)

pheatmap(mat[, ord],
         scale = "row", cluster_cols = FALSE, show_colnames = FALSE,
         annotation_col    = anot,
         annotation_colors = list(condicao = setNames(c("#1b6ca8", "#d1495b"),
                                                      c(REFERENCIA, ALVO))),
         color = colorRampPalette(c("#1b6ca8", "white", "#d1495b"))(51),
         main  = "Top 30 genes por padj — VST, z-score por gene",
         filename = "figuras/06_heatmap.png", width = 11, height = 9)

na_tela(pheatmap(mat[, ord],
         scale = "row", cluster_cols = FALSE, show_colnames = FALSE,
         annotation_col    = anot,
         annotation_colors = list(condicao = setNames(c("#1b6ca8", "#d1495b"),
                                                      c(REFERENCIA, ALVO))),
         color = colorRampPalette(c("#1b6ca8", "white", "#d1495b"))(51),
         main  = "Top 30 genes por padj — VST, z-score por gene"))

# Para um heatmap que signifique alguma coisa:
#   - selecione genes por VARIÂNCIA (independente do contraste), ou
#   - use um conjunto de genes definido A PRIORI (uma assinatura publicada).


# M5d — O gráfico que desmente: contagem de um gene só ------------------------
# Aqui a mentira morre. Se o efeito "significativo" vem de duas amostras
# extremas, o gráfico mostra na hora.

# plotar_gene(): boxplot das contagens de UM gene, com os pontos individuais por
# cima. É a ferramenta que revela se um efeito "significativo" se apoia em poucas
# amostras extremas.
plotar_gene <- function(nome, salvar = FALSE) {
  # 1. Localiza o gene pelo símbolo. Se não estiver na tabela, avisa e sai.
  linhas <- which(tab$gene_name == nome)
  if (length(linhas) == 0) {
    message(sprintf("'%s' não está entre os genes testados (pode ter caído no filtro).", nome))
    return(invisible(NULL))
  }
  r <- tab[linhas[1], ]
  # 2. plotCounts() extrai as contagens normalizadas do gene, por amostra, já
  #    com a condição de cada uma. returnData = TRUE devolve os dados em vez de
  #    plotar direto — assim montamos o gráfico do nosso jeito com o ggplot.
  d <- plotCounts(dds, gene = r$gene_id, intgroup = "condition", returnData = TRUE)
  padj_txt <- if (is.na(r$padj)) "NA" else format(r$padj, digits = 3, scientific = TRUE)

  g <- ggplot(d, aes(condition, count, fill = condition)) +
    geom_boxplot(width = .5, outlier.shape = NA, alpha = .8) +
    geom_jitter(width = .18, size = 2, alpha = .7) +
    scale_y_log10() +
    scale_fill_manual(values = setNames(c("#1b6ca8", "#d1495b"), c(REFERENCIA, ALVO))) +
    labs(title = sprintf("%s   —   log2FC = %+.2f,  padj = %s",
                         nome, r$log2FoldChange, padj_txt),
         x = NULL, y = "Contagem normalizada (log10)") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none", panel.grid.minor = element_blank())

  print(g)
  if (salvar) ggsave(sprintf("figuras/gene_%s.png", nome), g,
                     width = 5, height = 5, dpi = 200)
  invisible(g)
}

# 🔧 Troque o nome do gene e rode de novo (Ctrl+Enter na linha)
plotar_gene("MYC")

# EXPERIMENTE:
#   1. plotar_gene(tab$gene_name[1])   — o gene do topo da sua lista
#   2. um gene com padj alto           — veja a diferença
#   3. plotar_gene("TNFRSF17")         — BCMA, alvo de CAR-T e biespecíficos
#   4. plotar_gene("CD38")             — alvo do daratumumabe
#
# As duas últimas têm relevância clínica direta: perda de antígeno é um
# mecanismo conhecido de escape terapêutico.


# M6 — Da lista ao significado biológico --------------------------------------

# Filtra a tabela para os genes canônicos de mieloma (definidos no painel de
# controle, em GENES_INTERESSE) e mostra como cada um se comportou no contraste.
alvo_tab <- subset(tab, gene_name %in% GENES_INTERESSE,
                   select = c(gene_name, baseMean, log2FoldChange, padj))

if (nrow(alvo_tab) > 0) {
  alvo_tab$status <- ifelse(is.na(alvo_tab$padj), "filtrado / NA",
                     ifelse(alvo_tab$gene_name %in% sig$gene_name,
                            "significativo", "não significativo"))
  print(alvo_tab[order(alvo_tab$padj), ], row.names = FALSE, digits = 3)
} else {
  cat("Nenhum dos genes de interesse sobreviveu ao filtro.\n")
}

# 💬 Se MYC ou os genes de translocação (NSD2, FGFR3, MAF, CCND1) não aparecerem
#    como diferenciais, isso NÃO significa que não importam no mieloma.
#    Significa que a relevância deles é ESTRUTURAL — translocação, amplificação,
#    posicionamento junto ao enhancer de IGH — e não necessariamente
#    transcricional entre estes dois momentos da doença.


# M6b — Exportar --------------------------------------------------------------
# Quatro arquivos, e dois deles são os que quase ninguém exporta:
#
#  - genes_de_fundo.txt : o UNIVERSO correto para enriquecimento são os genes que
#    entraram no teste, não "todos os genes humanos". Fundo errado é a causa nº 1
#    de enriquecimento inflado.
#
#  - amostras_R.csv : a lista exata de amostras sorteadas. Ver nota abaixo.

write.csv(tab, "resultados/deseq_results_R.csv", row.names = FALSE)

# --- A lista de amostras: por que ela precisa existir ---
# A seleção das 30+30 amostras é ALEATÓRIA. O modelo estatístico não é — o
# DESeq2 é determinístico do começo ao fim — mas a escolha de QUAIS pacientes
# entram nele é um sorteio.
#
# E set.seed(42) no R NÃO produz o mesmo sorteio que random_state=42 no pandas.
# São geradores pseudoaleatórios diferentes, com algoritmos de permutação
# diferentes. O número é o mesmo; a sequência que ele gera, não.
#
# Sem este arquivo, a trilha Python sorteia OUTROS 60 pacientes, e o Módulo 9
# acaba medindo o sorteio em vez de medir a implementação. Já aconteceu neste
# workshop: correlação de fold change de 0,37 e concordância de sinal de 61%
# (o acaso puro é 50%).
#
# A solução não é "usar a mesma semente nos dois". É NÃO SORTEAR DUAS VEZES:
# um sorteia e grava; o outro lê.
write.csv(data.frame(file_id   = rownames(meta),
                     case_id   = meta$case_id,
                     condition = meta$condition),
          "resultados/amostras_R.csv", row.names = FALSE)

lista_degs <- unique(sig$gene_name[!is.na(sig$gene_name) & sig$gene_name != ""])
fundo      <- unique(tab$gene_name[!is.na(tab$padj) & !is.na(tab$gene_name)])

writeLines(lista_degs, "resultados/degs_para_enriquecimento.txt")
writeLines(fundo,      "resultados/genes_de_fundo.txt")

# Matriz por símbolo, para quem quiser conferir no iDEP
if (!is.null(mapa_genes)) {
  sym <- mapa_genes[rownames(counts_f), "gene_name"]
  ok  <- !is.na(sym) & sym != ""
  agg <- rowsum(counts_f[ok, , drop = FALSE], group = sym[ok])
  write.csv(agg, "resultados/mm_counts_symbols.csv")
  write.csv(data.frame(condition = meta$condition, row.names = rownames(meta)),
            "resultados/mm_metadata_idep.csv")
}

cat(sprintf("DEGs exportados : %s\n", fmt(length(lista_degs))))
cat(sprintf("Genes de fundo  : %s\n", fmt(length(fundo))))
cat(sprintf("Amostras salvas : %d (resultados/amostras_R.csv)\n", nrow(meta)))
cat("\nArquivos em resultados/ — veja no painel Files (canto inferior direito).\n")

# ⚠️ GUARDE DOIS ARQUIVOS para o dia 2:
#      resultados/deseq_results_R.csv  — a tabela de resultados
#      resultados/amostras_R.csv       — a lista de amostras
#    Os dois vão para o Colab, no Módulo 9. Sem o segundo, a comparação não vale.

# ENRIQUECIMENTO — ShinyGO: https://bioinformatics.sdstate.edu/go/
#   1. Cole o conteúdo de degs_para_enriquecimento.txt
#   2. Em "Custom background", cole genes_de_fundo.txt   ← não pule este passo
#   3. Espécie: Human. Rode.
#
# ⚠️ Cuidado com termos genéricos: quase toda lista grande de DEGs em câncer
#    enriquece "ciclo celular". Conjuntos com mais de 500 genes enriquecem por
#    inércia. Olhe o fold enrichment e o tamanho do conjunto, não só o p-valor.


# M7 — Reprodutibilidade ------------------------------------------------------
# Rode e GUARDE a saída. Sem isso a análise não é reproduzível: o GDC
# reprocessa os dados periodicamente, e "baixado do GDC" não identifica nada.

registro <- sprintf("
=== REGISTRO DA ANÁLISE (trilha R) ===
Data de execução      : %s
Fonte                 : NCI GDC, projeto MMRF-COMMPASS
Workflow              : STAR - Counts, coluna '%s'
Release do GDC        : %s

Desenho
  Contraste           : %s vs %s (referência = %s)
  n por grupo         : %s
  Seleção             : 1 amostra por paciente; subamostragem aleatória
  Semente             : %d

Processamento
  Genes antes filtro  : %s
  Filtro              : >= %d contagens em >= %d amostras
  Genes após filtro   : %s
  Transformação (EDA) : VST

Estatística
  Modelo              : ~ condition
  Teste               : %s
  Shrinkage do LFC    : lfcShrink(type = 'apeglm')
  Correção múltipla   : Benjamini-Hochberg
  Critério de DEG     : %s
  DEGs                : %s
",
Sys.time(), USAR_COLUNA,
if (nzchar(RELEASE_GDC)) RELEASE_GDC else "<<< PREENCHA: gdc.cancer.gov/about-data/data-release >>>",
ALVO, REFERENCIA, REFERENCIA,
paste(names(table(meta$condition)), table(meta$condition), collapse = ", "),
SEED,
fmt(nrow(counts)), MIN_CONTAGENS, min_amostras,
fmt(nrow(counts_f)), modo, criterio,
fmt(nrow(sig)))

# Imprime o registro na tela e, em seguida, grava tudo em arquivo.
cat(registro)

# sink() redireciona a saída do console para o arquivo. Entre o sink() de abrir e
# o sink() de fechar, nada aparece na tela — vai tudo para o REGISTRO.txt.
# sessionInfo() anexa a versão do R e de cada pacote: é o que permite reproduzir
# a análise no futuro.
sink("resultados/REGISTRO.txt")
cat(registro)
cat("\n=== VERSÕES ===\n")
print(sessionInfo())
sink()   # ⚠️ este segundo sink() é obrigatório: sem ele o console fica mudo

cat("\nRegistro salvo em resultados/REGISTRO.txt\n")
print(sessionInfo())

# LIMITAÇÕES DESTE DESENHO — precisam estar escritas
#
#  1. Não há plasmócito normal. Nada aqui permite afirmar o que difere entre
#     mieloma e tecido normal. O contraste é entre dois momentos da doença.
#  2. O tratamento é confundidor não controlado. Entre o diagnóstico e a
#     recidiva o paciente recebeu terapia; biologia da progressão e resposta ao
#     tratamento não são separáveis neste desenho.
#  3. Pureza variável: a fração de plasmócitos tumorais na seleção CD138+ difere
#     entre pacientes. Diferença de expressão pode ser diferença de composição
#     celular.
#  4. Sem correção de lote nem de covariáveis clínicas.
#  5. Sem coorte de validação independente.
#  6. Subamostragem didática de uma coorte maior.
#
# ✍️ EXERCÍCIO: escreva em cinco linhas a seção de métodos desta análise.
#    Passe para o colega ao lado. Ele tenta apontar o que falta para reproduzir
#    o resultado sem falar com você.


# M8 — Desafio final: o controle positivo -------------------------------------
#
# PROBLEMA: como saber se todo o pipeline — download, matriz, filtro, modelo —
# está correto? Um resultado plausível pode estar errado sem dar erro nenhum.
#
# SOLUÇÃO: rodar um contraste cuja resposta já conhecemos de antemão. Se
# comparássemos amostras de homens contra mulheres, o topo da lista TERIA que
# ser XIST (ligado em mulheres) e os genes do cromossomo Y (ligados em homens).
# Isso é biologia conhecida. Se o pipeline devolver exatamente isso, ele
# funciona. É um "controle positivo": um teste com gabarito.
#
# O OBSTÁCULO: para separar homens de mulheres, precisaríamos do sexo de cada
# paciente. Mas o MMRF-CoMMpass NÃO informa o sexo no tier de acesso aberto —
# o campo demographic.gender vem vazio para todos os casos (testado nos dois
# endpoints do GDC, /cases e /files).
#
#   💡 Lição à parte: o nível de acesso ao METADADO é uma barreira de
#      reprodutibilidade tão real quanto código fechado. Volta à ideia do M1.
#
# A SAÍDA: inferir o sexo a partir da própria expressão. Os genes do Y só são
# expressos em homens — então dá para deduzir o sexo olhando para eles. Mas em
# vez de confiar num gene só, usamos vários e os fazemos VOTAR: cada gene do Y
# classifica as amostras por conta própria (k-means em 2 grupos), e depois
# vemos se concordam entre si. Genes independentes do mesmo cromossomo TÊM que
# chegar ao mesmo resultado — se concordam, o rótulo é confiável.
#
# ⚠️ A ARMADILHA FINAL: o rótulo de sexo saiu da expressão dos genes do Y.
# Então encontrar os genes do Y no topo do contraste é TAUTOLOGIA — eles
# definiram os grupos, é óbvio que vão separar. Isso se chama controle
# degradado. Quem realmente valida o pipeline é o XIST: ele está no cromossomo
# X, NÃO participou da rotulagem, e mesmo assim deve aparecer no topo. Se o XIST
# aparecer, é prova independente de que tudo funciona.


GENES_Y <- c("RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "USP9Y", "TXLNGY", "ZFY", "NLGN4Y")
CONCORDANCIA_MIN <- 0.80

if (is.null(mapa_genes)) {
  cat("Sem mapa_genes (matriz congelada sem símbolos) — controle positivo indisponível.\n")
} else {

  # --- localiza os genes do Y na matriz VST ---
  # 1. Pega o símbolo de cada linha da matriz e monta uma função que devolve a
  #    posição de um gene pelo nome (ou NA se ele não sobreviveu ao filtro do M3).
  simb  <- mapa_genes[rownames(assay(vsd)), "gene_name"]
  achar <- function(nome) { k <- which(!is.na(simb) & simb == nome); if (length(k)) k[1] else NA_integer_ }
  # 2. Aplica a busca a todos os genes do Y e descarta os que não foram achados.
  idx   <- vapply(GENES_Y, achar, integer(1))
  idx   <- idx[!is.na(idx)]

  cat(sprintf("Genes do Y disponíveis (%d/%d): %s\n",
              length(idx), length(GENES_Y), paste(names(idx), collapse = ", ")))

  if (length(idx) < 3) {
    cat("⚠️ Menos de 3 genes do Y sobreviveram ao filtro do M3.\n")
    cat("   Baixe MIN_CONTAGENS/MIN_AMOSTRAS no painel e rode de novo a partir do M3.\n")
  } else {

    # 3. Extrai a submatriz só com os genes do Y (linhas) x amostras (colunas).
    V <- assay(vsd)[idx, , drop = FALSE]
    rownames(V) <- names(idx)

    # --- cada gene classifica as amostras sozinho ---
    # 4. Para cada gene do Y, roda um k-means de 2 grupos na expressão dele. O
    #    cluster de MAIOR expressão é o masculino. O gene "vota" 1 (male) ou 0
    #    (female) em cada amostra. A matriz votos tem genes nas linhas, amostras
    #    nas colunas. set.seed fixa o resultado do k-means.
    set.seed(SEED)
    votos <- t(apply(V, 1, function(v) {
      km    <- kmeans(matrix(v, ncol = 1), centers = 2, nstart = 25)
      alto  <- which.max(km$centers)          # cluster de maior expressão = male
      as.integer(km$cluster == alto)
    }))
    colnames(votos) <- colnames(V)

    # 5. Consolida os votos. consenso = fração de genes que chamou a amostra de
    #    "male". Acima de 0,5, o rótulo é male. concord mede o quanto os genes
    #    concordaram entre si naquela amostra (1,0 = unânime) — é a nossa medida
    #    de confiança no rótulo.
    consenso <- colMeans(votos)
    rotulo   <- ifelse(consenso > 0.5, "male", "female")
    concord  <- pmax(consenso, 1 - consenso)

    # --- diagnóstico por gene: quem destoa? ---
    separacao <- apply(V, 1, function(v) {
      km <- kmeans(matrix(v, ncol = 1), centers = 2, nstart = 25)
      abs(diff(as.vector(km$centers)))
    })
    diag_y <- data.frame(
      separacao_clusters  = round(separacao, 2),
      n_male              = rowSums(votos),
      concorda_c_consenso = round(rowMeans(votos == matrix(as.integer(rotulo == "male"),
                                            nrow = nrow(votos), ncol = ncol(votos),
                                            byrow = TRUE)), 3)
    )
    cat("\n--- diagnóstico por gene ---\n")
    print(diag_y[order(-diag_y$concorda_c_consenso), ])

    cat(sprintf("\nrótulo final: %s\n", paste(names(table(rotulo)), table(rotulo),
                                              sep = "=", collapse = "  ")))
    cat(sprintf("concordância média: %.3f   mínima: %.3f\n", mean(concord), min(concord)))

    # --- inspeção visual ---
    escore <- colMeans(V)
    png("figuras/09_controle_sexo.png", width = 1500, height = 500, res = 130)
    par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
    plot(escore, jitter(rep(0, length(escore)), amount = .08),
         col = ifelse(rotulo == "male", "#1b6ca8", "#d1495b"), pch = 19, cex = 1.2,
         yaxt = "n", ylab = "", xlab = sprintf("escore Y — média VST de %d genes", nrow(V)),
         main = "Escore agregado do cromossomo Y")
    legend("topleft", c("female", "male"), col = c("#d1495b", "#1b6ca8"),
           pch = 19, bty = "n")
    hist(concord, breaks = seq(.5, 1, by = .05), col = "#3d5a80", border = "white",
         xlab = "concordância entre genes, por amostra", ylab = "n amostras",
         main = "Votos alinhados = rótulo confiável")
    abline(v = CONCORDANCIA_MIN, lty = 2, col = "#d1495b", lwd = 2)
    dev.off()
    cat("Figura salva: figuras/09_controle_sexo.png\n")

    # --- checagens ---
    if (mean(concord) < CONCORDANCIA_MIN) {
      cat(sprintf("\n⚠️ Os genes do Y não concordam entre si (média %.2f).\n", mean(concord)))
      cat("   Olhe a tabela acima: se UM gene destoa, tire-o de GENES_Y. Se todos\n")
      cat("   destoam, o sinal de Y não está nesta amostragem.\n")
    } else if (length(unique(rotulo)) < 2 || min(table(rotulo)) < 3) {
      cat(sprintf("\n⚠️ Grupos insuficientes: %s. Aumente N_POR_GRUPO.\n",
                  paste(names(table(rotulo)), table(rotulo), sep = "=", collapse = ", ")))
    } else{
      
      # --- modelo ---
      # 6. Agora repete a análise de expressão diferencial, mas com SEXO no lugar
      #    de condição. Anexa o rótulo inferido ao metadata e confere o
      #    alinhamento amostras x colunas antes de ajustar (mesma trava do M2).
      meta_s <- meta
      meta_s$sexo <- factor(rotulo[rownames(meta_s)], levels = c("female", "male"))
      stopifnot(identical(colnames(counts_f), rownames(meta_s)))

      # 7. Ajusta o modelo ~ sexo e testa male vs female. É o mesmo DESeq() do M4,
      #    só que num contraste cuja resposta certa nós já conhecemos de antemão.
      dds_s <- DESeqDataSetFromMatrix(counts_f, meta_s, ~ sexo)
      dds_s <- DESeq(dds_s, quiet = TRUE)
      rs <- as.data.frame(results(dds_s, contrast = c("sexo", "male", "female")))
      rs$gene_name <- mapa_genes[rownames(rs), "gene_name"]
      rs <- rs[order(rs$padj), ]
      # Mostra os 15 genes de menor padj no contraste por sexo. O esperado é ver
      # XIST e os genes do Y (RPS4Y1, DDX3Y, UTY...) no topo — são eles que
      # diferenciam homens de mulheres.
      cat("\n=== TOP 15 — contraste por sexo ===\n")
      print(head(rs[, c("gene_name", "log2FoldChange", "padj")], 15), digits = 3)
      
      # --- veredito: só o XIST vale ---
      # ATENÇÃO à lógica: o rótulo de sexo foi inferido dos genes do Y. Então
      # achar os genes do Y no topo é tautologia — eles definiram os grupos, é
      # óbvio que separam. Quem valida de verdade é o XIST, que está no
      # cromossomo X e NÃO participou da rotulagem.
      
      # 1. Procura o XIST na tabela de resultados.
      xist <- rs[which(rs$gene_name == "XIST"), ]
      cat("\n", strrep("=", 62), "\n", sep = "")
      
      if (nrow(xist) > 0) {
        # 2. Extrai o fold change e o padj do XIST.
        lfc <- xist$log2FoldChange[1]; pa <- xist$padj[1]
        cat(sprintf("XIST — o teste que vale: log2FC = %+.2f | padj = %s\n",
                    lfc, fmt_padj(pa)))
        # traduz o log2FC para "quantas vezes", que se lê sem fazer conta
        cat(sprintf("  (log2FC %+.2f = XIST ~%.0fx menos expresso no grupo 'male')\n",
                    lfc, 2^abs(lfc)))
        
        # 3. O teste do gabarito: o XIST tem que estar REPRIMIDO nos homens
        #    (lfc < -1) e com padj significativo. Se sim, o pipeline inteiro —
        #    do download ao teste — funcionou, provado por um gene que não
        #    entrou na rotulagem.
        if (!is.na(pa) && lfc < -1 && pa < 0.05) {
          cat("\n✅ Pipeline validado. O XIST saiu reprimido nos homens, como a\n")
          cat("   biologia manda. Como ele está no cromossomo X e não participou\n")
          cat("   da inferência de sexo, isto é prova INDEPENDENTE de que cada\n")
          cat("   etapa (matriz, alinhamento amostra-metadado, normalização,\n")
          cat("   teste) está correta de ponta a ponta.\n")
        } else {
          # 4. Se o XIST NÃO se comportou como esperado, algo quebrou. As três
          #    causas mais prováveis, em ordem:
          cat("\n❌ XIST não se comportou como esperado. Revise:\n")
          cat("   (a) o alinhamento counts x meta no M2 (amostra trocada de rótulo)\n")
          cat("   (b) a direção do contraste (male vs female invertida)\n")
          cat("   (c) a tabela de diagnóstico acima — os grupos podem estar trocados\n")
        }
      } else {
        # 5. Se o XIST não passou no filtro do M3, ficamos sem a testemunha
        #    independente. Os genes do Y no topo, sozinhos, não provam nada.
        cat("⚠️ XIST não sobreviveu ao filtro do M3 — sem evidência independente.\n")
        cat("   Os genes do Y no topo são CIRCULARES e não validam nada.\n")
      }
      cat(strrep("=", 62), "\n")
      
      # 6. Salva a tabela completa do contraste por sexo, para consulta.
      write.csv(rs, "resultados/controle_positivo_sexo.csv", row.names = TRUE)
    }
  }
}

      
# FIM DA TRILHA R -------------------------------------------------------------
#
# Você percorreu o caminho completo em DESeq2: API pública → matriz → QC →
# modelo → lista de genes → via biológica → registro reprodutível.
#
# O QUE VEM AGORA
#
# O mesmo dado, o mesmo contraste, os mesmos parâmetros — em Python, com
# PyDESeq2. Não é repetição, é a pergunta que fecha o workshop:
#
#   Duas implementações independentes do mesmo método estatístico chegam ao
#   mesmo resultado?
#
# Se chegam, o achado é da biologia, não da ferramenta. Se divergem em algum
# ponto, descobrir onde e por quê ensina mais do que qualquer slide.
#
# ⚠️ ANTES DE FECHAR: guarde DOIS arquivos da pasta resultados/
#
#      deseq_results_R.csv  — a tabela de resultados
#      amostras_R.csv       — a lista das amostras sorteadas
#
#    Os dois vão para o Colab. O segundo é o que garante que a trilha Python
#    analise OS MESMOS pacientes: sem ele, o Python sorteia outro conjunto e o
#    Módulo 9 mede o sorteio, não a implementação.
#
#    No notebook, confirme que USAR_AMOSTRAS_DO_R = True no painel de controle.
#
# -----------------------------------------------------------------------------
# Material didático. Dados de acesso aberto do NCI Genomic Data Commons
# (estudo MMRF CoMMpass). Adaptado do pipeline chol-expression-portfolio.
# =============================================================================


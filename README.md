<p align="center">
  <img src="docs/logo.svg" alt="Workshop de RNA-seq — expressão diferencial em mieloma múltiplo" width="720">
</p>

<p align="center">
  <a href="https://doi.org/10.5281/zenodo.22399365"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.22399365.svg" alt="DOI"></a>
</p>

Material de um workshop de dois dias sobre análise de expressão diferencial em
RNA-seq, usando dados reais de mieloma múltiplo do **MMRF-CoMMpass** (camada
aberta do NCI Genomic Data Commons). O contraste é diagnóstico contra recidiva,
com cerca de 30 amostras por grupo.

A mesma análise é feita duas vezes: no primeiro dia em R com DESeq2, no segundo
em Python com PyDESeq2. Comparar os dois resultados não é um apêndice — é o
exercício central, e é onde aparece por que reprodutibilidade dá trabalho.

## Estrutura

| Pasta | Conteúdo |
|---|---|
| `dia0-tratamento/` | Módulo 0.5 — do FASTQ à matriz de contagens: FastQC, fastp, HISAT2, samtools, featureCounts e MultiQC. Roda no Colab, com dados de teste de levedura. |
| `dia1-R/` | Instalação de pacotes, análise principal com DESeq2, exploração por gene, figuras para publicação e um explorador em HTML. |
| `dia2-python/` | O mesmo fluxo em PyDESeq2, mais o módulo de concordância entre as duas linguagens. |
| `docs/` | Guia de código, roteiro do instrutor e material de apoio. |

## Dados

A matriz congelada está no Zenodo, não no GDC. A razão é prática: trinta pessoas
consultando a API do GDC ao mesmo tempo, no wi-fi de um auditório, não terminam a
aula. E há uma razão melhor: um DOI fixa a versão do dado, coisa que uma consulta
à API não faz.

```
https://doi.org/10.5281/zenodo.22399365
```

Os scripts baixam os arquivos sozinhos. Não é preciso conta no GDC nem acesso
controlado ao dbGaP.

## Como usar

**Antes do Dia 1 — `dia0-tratamento/`**
Abra o notebook no Google Colab e execute do começo ao fim. Leva de 15 a 20
minutos e mostra de onde vem uma matriz de contagens. Nenhuma instalação local.

**Dia 1 — `dia1-R/`**
No RStudio Desktop, rode `00_instalar_pacotes.R` uma vez e depois
`01_workshop_mieloma.R`. O script está dividido em módulos M0 a M8; use
`Ctrl+Shift+O` para navegar entre eles e `Ctrl+Alt+T` para executar a seção onde
o cursor estiver.

**Dia 2 — `dia2-python/`**
No Colab, abra o notebook do PyDESeq2. O módulo 9 refaz a comparação com o
resultado do Dia 1 e mede a concordância gene a gene.

## O que esperar do resultado

O contraste diagnóstico contra recidiva devolve **poucos genes diferenciais** — e
isso está certo. Ausência de mudança não é ausência de expressão: as células ainda
são plasmócitos malignos nos dois momentos, e o que separa recidiva de diagnóstico
é mais subclonal do que transcricional em larga escala.

Os dois pipelines chegam a listas de tamanhos diferentes, com correlação de
log2FoldChange acima de 0,999 entre eles. A discussão de por que isso acontece —
filtro independente, convenção de p-valor composto, denominador do
Benjamini-Hochberg — é o assunto do módulo de concordância.

Duas decisões metodológicas que valem por metade do workshop:

- **Segmentos de imunoglobulina são removidos antes do ajuste do modelo.** Em
  mieloma, o rearranjo V(D)J é específico do clone de cada paciente. Se ficarem na
  matriz, dominam o contraste por acaso.
- **A identidade das amostras é travada entre as duas linguagens.**
  `set.seed(42)` em R e `random_state=42` em pandas não sorteiam os mesmos
  pacientes. O arquivo `amostras_R.csv` existe para resolver isso.

## Requisitos

R 4.3 ou superior com DESeq2 e apeglm, instalados pelo script do Dia 1. Para os
notebooks, só um navegador — o Colab cuida do resto.

## Licença

Código sob licença MIT. Material didático sob CC BY 4.0. Os dados de expressão são
do MMRF-CoMMpass, distribuídos pelo NCI Genomic Data Commons sob os termos do
consórcio; cite a fonte original ao reutilizá-los.

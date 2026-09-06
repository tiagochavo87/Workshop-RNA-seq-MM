<div align="center">

# Workshop de RNA-seq
### Expressão diferencial em mieloma múltiplo

**A mesma análise, duas vezes: R e Python.**
Não pela redundância — pela calibragem do que se deve levar a sério.

<br>

[![Dados](https://img.shields.io/badge/dados-10.5281%2Fzenodo.22399365-1f4e79?style=flat-square)](https://doi.org/10.5281/zenodo.22399365)
[![R](https://img.shields.io/badge/R-%E2%89%A5%204.3-276DC3?style=flat-square&logo=r&logoColor=white)](https://cran.r-project.org/)
[![DESeq2](https://img.shields.io/badge/DESeq2-1.50-9b2226?style=flat-square)](https://bioconductor.org/packages/DESeq2/)
[![PyDESeq2](https://img.shields.io/badge/PyDESeq2-0.5-3776AB?style=flat-square&logo=python&logoColor=white)](https://pydeseq2.readthedocs.io/)
[![Licença](https://img.shields.io/badge/licen%C3%A7a-CC%20BY%204.0-0b7a75?style=flat-square)](#licença)

</div>

---

## O que você vai fazer aqui

Percorrer o caminho completo de uma análise de expressão diferencial — **da API pública
até a lista de genes** — e depois refazer tudo em outra linguagem para descobrir o quanto
o resultado depende da ferramenta.

O dado é real: RNA-seq de medula óssea de pacientes com mieloma múltiplo, do estudo
[MMRF-CoMMpass](https://portal.gdc.cancer.gov/projects/MMRF-COMMPASS), na camada de
acesso aberto do NCI Genomic Data Commons.

O contraste é **diagnóstico × recidiva**, com 30 pacientes de cada lado.

> [!NOTE]
> Não existe tecido normal nesta coorte. Isso não é uma limitação a contornar — é o que
> define a pergunta que se pode fazer. A primeira lição do workshop está aí.

---

## Comece por aqui

### Antes do workshop — 10 minutos

Abra o RStudio e rode o instalador. Ele detecta problemas de permissão, cria uma
biblioteca pessoal se precisar, e testa o DESeq2 com uma análise de brinquedo.

```r
source("dia1-R/00_instalar_pacotes.R")
```

Se aparecer qualquer ❌ no relatório final, **copie o relatório inteiro e mande para o
instrutor.** Não tente resolver no dia.

<br>

### Dia 1 — R e DESeq2

1. Duplo clique em **`dia1-R/workshop-mieloma.Rproj`**
   Isso abre o RStudio já apontado para a pasta certa. Sem isso, as figuras vão parar em Documentos.
2. Abra o `01_workshop_mieloma.R`
3. **`Ctrl+Shift+O`** — abre o painel *Outline*, que é a agenda da aula
4. **`Ctrl+Alt+T`** — executa a seção onde o cursor está

<br>

### Dia 2 — Python e PyDESeq2

[![Abrir no Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/tiagochavo87/Workshop-RNA-seq-MM/blob/main/dia2-python/workshop_rnaseq_mieloma.ipynb)

Clique no botão, faça **Arquivo → Salvar uma cópia no Drive**, e trabalhe na sua cópia.

Você vai precisar de dois arquivos gerados no Dia 1, da pasta `resultados/`:

| Arquivo | Para quê |
|---|---|
| `deseq_results_R.csv` | a tabela que o Módulo 9 compara |
| `amostras_R.csv` | trava as mesmas amostras nas duas trilhas |

> [!IMPORTANT]
> Sem o `amostras_R.csv`, as duas trilhas sorteiam pacientes diferentes.
> `set.seed(42)` no R e `random_state=42` no pandas **não** produzem o mesmo sorteio —
> são geradores distintos. O Módulo 9 acabaria medindo o sorteio em vez da implementação.

---

## O que tem em cada arquivo

```
Workshop-RNA-seq-MM/
│
├── dia1-R/
│   ├── workshop-mieloma.Rproj      abra este primeiro
│   ├── 00_instalar_pacotes.R       rode antes do workshop
│   ├── 01_workshop_mieloma.R       a análise completa, M0 a M8
│   └── 02_explorar_genes.R         gene a gene, depois do 01
│
├── dia2-python/
│   └── workshop_rnaseq_mieloma.ipynb    a mesma análise + Módulo 9
│
└── docs/
    ├── guia_codigo_R.docx          o que cada seção do código faz, e por quê
    └── slides_workshop.pdf         os slides da aula
```

Os dados **não** estão no repositório. Os scripts baixam sozinhos do Zenodo — são cerca
de 5 MB comprimidos.

---

## Os módulos

<table>
<tr><td width="60"><b>M0</b></td><td>Ambiente, pacotes, versões</td></tr>
<tr><td><b>M1</b></td><td>A API do GDC. Que grupos existem neste dado?</td></tr>
<tr><td><b>M2</b></td><td>Da contagem à matriz — e a checagem de alinhamento que salva a análise</td></tr>
<tr><td><b>M3</b></td><td>QC antes do teste: biblioteca, filtro de genes, VST, PCA</td></tr>
<tr><td><b>M4</b></td><td>DESeq2 em quatro etapas. Onde entra o limiar de fold change</td></tr>
<tr><td><b>M5</b></td><td>Visualização honesta — e o heatmap que engana</td></tr>
<tr><td><b>M6</b></td><td>Da lista ao significado. O universo de fundo correto</td></tr>
<tr><td><b>M7</b></td><td>Reprodutibilidade: o registro que acompanha o resultado</td></tr>
<tr><td><b>M8</b></td><td>Controle positivo por sexo, validado em XIST</td></tr>
<tr><td><b>M9</b></td><td><b>Concordância R × Python</b> — o módulo que justifica ter feito duas vezes</td></tr>
</table>

---

## O que esperar do resultado

Rodando com os parâmetros padrão, você deve chegar nisto:

| | |
|---|---|
| Matriz bruta | 60.660 genes × 60 amostras |
| Após o filtro de expressão | 19.453 genes |
| Segmentos V(D)J removidos | 170 |
| **Genes testados** | **19.283** |
| DEGs no R | 4 — CXCL9, CD4, IFI27, HLA-DQA2 |
| DEGs no Python | 1 — CXCL9 |
| Correlação dos log2FC | Pearson 0,9996 |
| Sobreposição das listas | Jaccard 0,250 |
| Controle positivo (XIST) | log2FC −8,15 · padj < 0,001 ✅ |

**Quatro genes.** Não é pouco por falha do pipeline — é o que existe entre dois momentos
da mesma doença, no mesmo paciente. Compare com um contraste tumor × tecido normal,
onde se acham 10 mil. A diferença é biológica, não técnica.

> [!TIP]
> **Ausência de mudança não é ausência de expressão.**
> CD38, SLAMF7 e TNFRSF17 — os alvos do daratumumabe, do elotuzumabe e do CAR-T —
> têm baseMean na casa das dezenas de milhares e log2FC praticamente zero.
> Continuam lá, expressos, na recidiva.

---

## O achado do Módulo 9

Os dois pipelines chegam a **fold changes idênticos** (diferença mediana de 0,00025) e a
**listas de DEGs diferentes**. A explicação não é ruído numérico:

```
razão padj_Python / padj_R  =  1,9983
coeficiente de variação      =  0,0%
```

Uma razão constante em todos os genes comparáveis. As duas bibliotecas calculam de forma
diferente o p-valor do teste com limiar — sob a hipótese nula composta `|LFC| ≤ limiar`
há mais de uma forma defensável de somar as caudas, e elas escolheram formas diferentes.

Qual das duas está certa? O dado não diz. Só o código-fonte.

**A lição:** um gene com `padj = 0,027` numa implementação e `0,055` na outra é o mesmo
gene, com o mesmo efeito. O que mudou foi o lado do corte em que ele caiu. Leve a sério
os que sobrevivem às duas análises com folga — e reporte como fronteira os que ficam na
fronteira.

---

## Quando algo dá errado

| Sintoma | O que fazer |
|---|---|
| `there is no package called 'DESeq2'` | Rodou o `00_instalar_pacotes.R`? Se ele criou uma biblioteca pessoal, **feche e reabra o RStudio** |
| `figure margins too large` | O painel *Plots* está pequeno. Arraste a divisória ou clique em *Zoom*. O PNG já foi salvo em `figuras/` |
| Figuras aparecem em Documentos | Você não abriu pelo `.Rproj` |
| `invalid multibyte string` ao ler o Zenodo | Rede instável. Rode a seção de novo — ela tem quatro tentativas |
| Colab travado ou com erro estranho | *Ambiente de execução → Reiniciar e executar tudo* |
| Módulo 9 diz "amostras diferentes" | Falta o `amostras_R.csv`. Suba pelo painel 📁 e rode desde o M1 |
| Módulo 9 diz "genes só no R/Python" | Os dois lados usaram `REMOVER_IG` diferente. Confira o painel de controle das duas trilhas |

---

## Os dados

Publicados no Zenodo, com DOI, e congelados: todo mundo analisa exatamente a mesma
matriz, hoje e daqui a dois anos.

**DOI:** [10.5281/zenodo.22399365](https://doi.org/10.5281/zenodo.22399365)

| Arquivo | Conteúdo |
|---|---|
| `mm_counts.csv.gz` | contagens brutas, Ensembl ID com versão |
| `mm_metadata.csv` | file_id, case_id, barcode, sample_type, condition |
| `gene_names.csv.gz` | mapa gene_id → símbolo e tipo |
| `mm_counts_symbols.csv` | contagens por símbolo, para o iDEP |
| `REGISTRO.txt` | procedência completa, somas MD5 e DOI |

O `REGISTRO.txt` traz tudo que é preciso para reconstruir a matriz a partir do GDC:
release, filtros, semente e a lista nominal dos 60 arquivos de origem. É o modelo do que
deveria acompanhar qualquer dado publicado.

---

## Se você quiser refazer o congelamento

Os scripts do instrutor não estão aqui — eles baixam do GDC, montam a matriz e publicam
no Zenodo. Peça se precisar.

O que vale saber: o GDC reprocessa os dados entre *releases*, e os rótulos de
`sample_type` mudam. Antes de qualquer edição nova do workshop, confirme os rótulos
literais na API.

---

## Limitações deste desenho

Estão no notebook também, mas merecem estar aqui:

1. **Não há plasmócito normal.** O contraste é entre dois momentos da mesma doença.
2. **O tratamento é um confundidor não controlado.** Entre o diagnóstico e a recidiva o
   paciente recebeu terapia; parte do que aparece como progressão é resposta à droga.
3. **Pureza variável.** A fração de plasmócitos na seleção CD138+ difere entre pacientes.
   Diferença de expressão pode ser diferença de composição celular.
4. **Sem covariáveis clínicas** — idade, sexo, subtipo citogenético, estádio ISS.
5. **Sem coorte de validação.**
6. **Subamostragem didática:** 30 por grupo, não a coorte inteira.

Este é material de ensino. Demonstra o método; não é um estudo com rigor metodológico.

---

## Créditos e referências

**DESeq2** — Love MI, Huber W, Anders S. *Moderated estimation of fold change and
dispersion for RNA-seq data with DESeq2.* Genome Biology 2014;15(12):550.

**apeglm** — Zhu A, Ibrahim JG, Love MI. *Heavy-tailed prior distributions for sequence
count data.* Bioinformatics 2019;35(12):2084–2092.

**PyDESeq2** — Muzellec B, Teleńczuk M, Cabeli V, Andreux M. *PyDESeq2: a python package
for bulk RNA-seq differential expression analysis.* Bioinformatics 2023;39(9):btad547.

**Dados** — MMRF CoMMpass Study, via NCI Genomic Data Commons. Os dados de acesso aberto
foram gerados pela Multiple Myeloma Research Foundation.

---

## Licença

Material didático sob [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/deed.pt-br).
Use, adapte e redistribua, citando a fonte.

<div align="center">
<br>
<sub>

**Tiago Fernando Chaves** · Universidade Federal de Santa Catarina

*A análise que não pode ser refeita não é um resultado — é uma opinião.*

</sub>
</div>

# Etapa 3 — Curadoria mitogenômica inspirada no MitoCurator

Este arquivo descreve a etapa `curate.py`, adicionada ao pipeline para fazer uma curadoria semelhante às funções centrais do MitoCurator.

A ideia é: depois da montagem e da escolha do melhor contig/FASTA candidato, você roda esta etapa para verificar se o candidato parece um mitogenoma coerente.

## O que esta etapa faz

A etapa de curadoria acrescenta:

1. leitura de FASTA ou GenBank;
2. estatísticas da sequência candidata;
3. busca de genes mitocondriais por BLASTN/TBLASTN contra os CDS de referência;
4. presença/ausência de genes esperados;
5. detecção de genes com múltiplos hits, que podem indicar duplicação, repeats ou redundância;
6. teste simples de circularidade por sobreposição terminal;
7. rotação de sequência circular para iniciar em um gene definido, por exemplo `COX1`;
8. quando a entrada é GenBank, diagnóstico de CDS, tRNA, rRNA e regiões intergênicas;
9. quando a entrada é GenBank, detecção de CDS com tamanho fora de múltiplo de 3, stops internos e ausência de stop terminal;
10. quando a entrada é GenBank, busca de regiões intergênicas AT-rich;
11. geração de tabelas TSV e relatório Markdown.

## Importante

Esta etapa **não monta** o genoma e **não substitui** anotadores como MITOS, MFannot ou MitoFinder.

Ela serve para curar e diagnosticar um candidato final ou semfinal.

Fluxo recomendado:

```text
reads brutas
  ↓
mito_pipeline.py ou mito_extract_reassemble.py
  ↓
contig/FASTA candidato
  ↓
anotação com MITOS/MFannot/MitoFinder, se possível
  ↓
curate.py
  ↓
relatório final de curadoria
```

## Instalação

Atualize o ambiente:

```bash
conda activate mito_assembly
mamba install -y -c conda-forge -c bioconda biopython blast
```

Se você recriar o ambiente do zero:

```bash
mamba env create -f environment.yml
conda activate mito_assembly
```

## Uso com FASTA candidato

Este é o uso que você provavelmente consegue rodar agora, usando o melhor contig candidato da montagem.

Exemplo usando o FASTA de candidatos gerado pelo primeiro pipeline:

```bash
python curate.py \
  --candidate mito_illumina/03_mito_candidates/mitochondrial_candidates.fasta \
  --input-format fasta \
  --nt-ref Agabis_H97_nt.fasta \
  --prot-ref Agabis_H97_prot.fasta \
  --expected-from-ref \
  --topology circular \
  --rotate-to COX1 \
  --genetic-code 4 \
  --outdir mito_illumina_curation \
  --threads 15
```

Para muitos fungos mitocondriais, o código genético usado frequentemente é o NCBI transl_table 4. Para invertebrados, use 5. Para vertebrados, use 2. Se a anotação final vier com outro `transl_table`, use o código indicado pela própria anotação.

## Uso com apenas o melhor contig

Como o seu relatório indicou que o melhor candidato original foi `NODE_131_length_46048_cov_410.638120`, uma boa prática é criar um FASTA só com esse contig e curar apenas ele.

Com `seqkit`:

```bash
seqkit grep -n -p 'NODE_131_length_46048_cov_410.638120' \
  mito_illumina/03_mito_candidates/mitochondrial_candidates.fasta \
  > NODE_131_mitogenome_candidate.fasta
```

Depois rode:

```bash
python curate.py \
  --candidate NODE_131_mitogenome_candidate.fasta \
  --input-format fasta \
  --nt-ref Agabis_H97_nt.fasta \
  --prot-ref Agabis_H97_prot.fasta \
  --expected-from-ref \
  --topology circular \
  --rotate-to COX1 \
  --genetic-code 4 \
  --outdir NODE_131_curation \
  --threads 15
```

Se seus headers não usam `COX1`, tente outros nomes comuns:

```bash
--rotate-to COI
```

ou:

```bash
--rotate-to COB
```

O script normaliza alguns sinônimos, como `COI` para `COX1` e `COB` para `CYTB`.

## Uso com GenBank anotado

Quando você tiver um arquivo `.gb` ou `.gbk`, por exemplo gerado por MITOS/MFannot/MitoFinder, rode:

```bash
python curate.py \
  --candidate candidato_anotado.gbk \
  --input-format genbank \
  --nt-ref Agabis_H97_nt.fasta \
  --prot-ref Agabis_H97_prot.fasta \
  --expected-from-ref \
  --expected-from-genbank \
  --topology circular \
  --rotate-to COX1 \
  --genetic-code 4 \
  --outdir candidato_anotado_curation \
  --threads 15
```

Com GenBank, o pipeline consegue fazer mais coisas do que com FASTA puro:

- contar CDS, rRNA, tRNA e outras features;
- traduzir CDS com o código genético informado;
- detectar CDS fora de frame;
- detectar stops internos;
- detectar CDS sem stop terminal;
- calcular regiões intergênicas;
- detectar candidatas AT-rich.

## Arquivos gerados

Dentro do diretório de saída, os principais arquivos são:

```text
CURATION_REPORT.md
commands.log
circularity_summary.tsv
gene_presence_absence.tsv
missing_genes.tsv
duplicate_or_repeat_gene_hits.tsv
blast_gene_search/gene_hits.filtered.tsv
blast_gene_search/gene_hits.all.tsv
gene_qc.tsv
problematic_features.tsv
intergenic_regions.tsv
at_rich_candidates.tsv
rotated_start_<GENE>.fasta
```

## Como interpretar os principais arquivos

### `CURATION_REPORT.md`

É o relatório principal. Ele resume:

- tamanho e N50 da sequência candidata;
- evidência de circularidade por overlap terminal;
- genes presentes/ausentes;
- genes com múltiplos hits;
- problemas em CDS, quando a entrada é GenBank;
- regiões AT-rich, quando a entrada é GenBank;
- arquivo FASTA rotacionado, se a rotação funcionou.

### `gene_presence_absence.tsv`

Tabela com os genes esperados e a decisão para cada um:

- `PRESENT`: gene encontrado;
- `MISSING`: gene esperado não encontrado;
- `MULTIPLE_HITS_CHECK_DUPLICATION_OR_REPEATS`: há múltiplos hits; revisar se é duplicação real, repeat, fragmentação ou redundância.

### `circularity_summary.tsv`

Mostra o melhor overlap entre o início e o fim de cada sequência.

Um overlap terminal forte precisa de identidade alta. Um overlap de 500 bp com 30–40% de identidade, por exemplo, **não prova circularidade**.

### `rotated_start_<GENE>.fasta`

FASTA rotacionado para iniciar no gene solicitado, caso o gene tenha sido encontrado.

Atenção: rotação **não prova circularidade**. Ela só muda o ponto inicial da sequência circular para padronizar a apresentação.

### `gene_qc.tsv`

Só é informativo quando a entrada é GenBank.

Mostra, por feature:

- tipo da feature;
- gene;
- coordenadas;
- strand;
- tamanho;
- se CDS é múltiplo de 3;
- número de stops internos;
- presença de stop terminal;
- status e sugestão de revisão.

### `at_rich_candidates.tsv`

Só aparece com conteúdo quando a entrada é GenBank e há regiões intergênicas grandes e ricas em A+T.

Por padrão, o filtro é:

```text
comprimento >= 500 bp
A+T >= 75%
```

Você pode mudar isso:

```bash
--at-rich-min-len 300 --at-rich-min-at 70
```

## Recomendações para o seu caso

Pelos relatórios anteriores, o melhor candidato ainda parece ser o contig original de aproximadamente 46 kb. Então eu começaria a curadoria por ele, não pela remontagem enriquecida fragmentada.

A ordem prática seria:

```bash
# 1. Extrair o melhor candidato
seqkit grep -n -p 'NODE_131_length_46048_cov_410.638120' \
  mito_illumina/03_mito_candidates/mitochondrial_candidates.fasta \
  > NODE_131_mitogenome_candidate.fasta

# 2. Curar com os CDS de Agabis como referência
python curate.py \
  --candidate NODE_131_mitogenome_candidate.fasta \
  --input-format fasta \
  --nt-ref Agabis_H97_nt.fasta \
  --prot-ref Agabis_H97_prot.fasta \
  --expected-from-ref \
  --topology circular \
  --rotate-to COX1 \
  --genetic-code 4 \
  --outdir NODE_131_curation \
  --threads 15

# 3. Ler o relatório
less NODE_131_curation/CURATION_REPORT.md
```

Depois, quando você tiver uma anotação em GenBank, rode novamente usando o `.gbk`, porque a curadoria fica mais completa.

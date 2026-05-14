# Pipeline Python para montagem de genoma mitocondrial

Este repositório contém um pipeline em Python para montar e avaliar um possível genoma mitocondrial usando **Illumina paired-end** ou **PacBio**. O tipo de dado é escolhido por argumento no comando (`--mode illumina` ou `--mode pacbio`).

O pipeline foi pensado para o seu caso com estes arquivos:

```text
Agabis_H97_nt.fasta              # CDS mitocondriais em nucleotídeos
Agabis_H97_prot.fasta            # proteínas codificadas pelos CDS mitocondriais
A13_20_illumina_1_trimmed.fastq.gz
A13_20_illumina_2_trimmed.fastq.gz
A13_20_pb.fastq.gz
```

## Ideia geral

O pipeline **não monta o genoma mitocondrial apenas a partir dos CDS de Agabis**. A montagem é feita a partir das reads. Os arquivos `Agabis_H97_nt.fasta` e `Agabis_H97_prot.fasta` são usados como **iscas/referências funcionais** para localizar, entre os contigs montados, quais contigs provavelmente são mitocondriais.

Fluxo geral:

```text
reads Illumina ou PacBio
        |
        v
montagem de novo
        |
        v
contigs/scaffolds montados
        |
        v
BLASTN: CDS nucleotídicos de Agabis vs assembly
TBLASTN: proteínas de Agabis vs assembly
        |
        v
seleção dos contigs candidatos mitocondriais
        |
        v
mapeamento das reads de volta nos candidatos
        |
        v
cobertura, profundidade, breadth, evidência de circularidade e relatório final
```

## Arquivos incluídos

```text
mito_pipeline.py      # pipeline principal
README.md             # este documento
environment.yml       # ambiente Conda/Mamba sugerido
```

## Dependências

O pipeline é escrito em Python, mas chama ferramentas externas de bioinformática. As dependências principais são:

- Python 3.10 ou superior
- SPAdes, para montagem com Illumina
- Flye, para montagem com PacBio
- BLAST+, para `blastn`, `tblastn` e `makeblastdb`
- minimap2, para mapear reads de volta nos candidatos
- samtools, para indexar BAM e calcular cobertura/profundidade

### Instalação recomendada com mamba

```bash
mamba env create -f environment.yml
mamba activate mito_assembly
```

Ou, manualmente:

```bash
mamba create -n mito_assembly -c conda-forge -c bioconda \
  python=3.11 spades flye minimap2 samtools blast

mamba activate mito_assembly
```

Teste se tudo está no `PATH`:

```bash
which python
which spades.py
which flye
which blastn
which tblastn
which makeblastdb
which minimap2
which samtools
```

## Como executar

Entre no diretório onde estão os arquivos FASTA/FASTQ e o script `mito_pipeline.py`.

### 1. Rodando com Illumina

```bash
python mito_pipeline.py \
  --mode illumina \
  --nt-ref Agabis_H97_nt.fasta \
  --prot-ref Agabis_H97_prot.fasta \
  --r1 A13_20_illumina_1_trimmed.fastq.gz \
  --r2 A13_20_illumina_2_trimmed.fastq.gz \
  --outdir mito_illumina \
  --threads 20
```

O modo Illumina usa o SPAdes com `--careful` e depois identifica os contigs mitocondriais por BLASTN/TBLASTN.

### 2. Rodando com PacBio raw/CLR

```bash
python mito_pipeline.py \
  --mode pacbio \
  --nt-ref Agabis_H97_nt.fasta \
  --prot-ref Agabis_H97_prot.fasta \
  --pacbio A13_20_pb.fastq.gz \
  --pacbio-type raw \
  --flye-genome-size 80k \
  --outdir mito_pacbio \
  --threads 20
```

Use `--pacbio-type raw` se as reads forem PacBio CLR/raw. Se forem PacBio HiFi, use:

```bash
--pacbio-type hifi
```

Exemplo HiFi:

```bash
python mito_pipeline.py \
  --mode pacbio \
  --nt-ref Agabis_H97_nt.fasta \
  --prot-ref Agabis_H97_prot.fasta \
  --pacbio A13_20_pb.fastq.gz \
  --pacbio-type hifi \
  --flye-genome-size 80k \
  --outdir mito_pacbio_hifi \
  --threads 20
```

### 3. Rodando apenas para conferir comandos

Antes de gastar tempo computacional, você pode fazer um teste seco:

```bash
python mito_pipeline.py \
  --mode illumina \
  --nt-ref Agabis_H97_nt.fasta \
  --prot-ref Agabis_H97_prot.fasta \
  --r1 A13_20_illumina_1_trimmed.fastq.gz \
  --r2 A13_20_illumina_2_trimmed.fastq.gz \
  --outdir teste_dry_run \
  --threads 20 \
  --dry-run
```

Isso registra os comandos em `commands.log`, mas não executa a montagem.

## Parâmetros importantes

### `--mode`

Escolhe a estratégia de montagem:

```text
illumina  -> usa SPAdes
pacbio    -> usa Flye
```

### `--flye-genome-size`

Tamanho esperado aproximado do mitogenoma. Para fungos, pode variar bastante. Um ponto inicial razoável é:

```bash
--flye-genome-size 80k
```

Se a montagem ficar ruim, teste outros valores, por exemplo:

```bash
--flye-genome-size 40k
--flye-genome-size 120k
--flye-genome-size 200k
```

### `--flye-meta`

Pode ser útil quando o DNA mitocondrial é uma fração pequena dentro do DNA total e a cobertura é desigual:

```bash
--flye-meta
```

Exemplo:

```bash
python mito_pipeline.py \
  --mode pacbio \
  --nt-ref Agabis_H97_nt.fasta \
  --prot-ref Agabis_H97_prot.fasta \
  --pacbio A13_20_pb.fastq.gz \
  --pacbio-type raw \
  --flye-genome-size 80k \
  --flye-meta \
  --outdir mito_pacbio_meta \
  --threads 20
```

### Filtros de BLAST

Os filtros padrão são conservadores o suficiente para encontrar contigs relacionados, mas podem ser ajustados:

```bash
--evalue 1e-5
--min-bitscore 50
--min-hit-qcov 0.30
--min-nt-pident 70
--min-prot-pident 30
```

Se a espécie for distante de Agabis, você pode reduzir a exigência do BLASTN e confiar mais no TBLASTN:

```bash
--min-nt-pident 60 --min-prot-pident 25
```

Se aparecerem muitos falsos positivos, aumente os filtros:

```bash
--min-nt-pident 80 --min-prot-pident 40 --min-hit-qcov 0.50
```

### Estatísticas de reads

Por padrão, o pipeline lê o FASTQ inteiro para calcular estatísticas. Se quiser acelerar apenas para teste:

```bash
--stats-max-reads 100000
```

Use `0` para estatística completa:

```bash
--stats-max-reads 0
```

## Estrutura de saída

Para `--outdir mito_illumina`, por exemplo:

```text
mito_illumina/
├── commands.log
├── REPORT.md
├── 01_assembly_spades/
│   ├── contigs.fasta
│   └── scaffolds.fasta
├── 02_blast_mito_hits/
│   ├── assembly_db.*
│   ├── blastn_nt_vs_assembly.tsv
│   └── tblastn_prot_vs_assembly.tsv
├── 03_mito_candidates/
│   ├── mitochondrial_candidates.fasta
│   └── mitochondrial_candidates_summary.tsv
└── 04_mapping_coverage/
    ├── reads_vs_mito_candidates.sorted.bam
    ├── reads_vs_mito_candidates.sorted.bam.bai
    ├── samtools_flagstat.txt
    ├── samtools_stats.txt
    ├── samtools_coverage.tsv
    ├── depth_per_base.tsv
    └── coverage_summary.tsv
```

No modo PacBio, a pasta de montagem será:

```text
01_assembly_flye/
```

## Como interpretar os principais resultados

### 1. `REPORT.md`

É o arquivo principal para apresentar ao professor. Ele contém:

- estatísticas das reads;
- estatísticas dos CDS de referência;
- estatísticas da montagem;
- lista dos contigs candidatos mitocondriais;
- hits de BLASTN/TBLASTN;
- cobertura média;
- profundidade mínima e máxima;
- breadth de cobertura em 1x, 5x, 10x e 20x;
- avaliação automática simples de circularidade.

### 2. `mitochondrial_candidates.fasta`

Este é o FASTA dos contigs que parecem mitocondriais.

Se houver **um único contig grande**, com muitos genes mitocondriais, cobertura alta e uniforme, e evidência de circularidade, você provavelmente está próximo de um mitogenoma completo.

Se houver **vários contigs**, a montagem está fragmentada ou existem regiões repetitivas/ambíguas. Nesse caso, o resultado ainda é útil, mas não deve ser chamado diretamente de genoma mitocondrial completo sem curadoria adicional.

### 3. `mitochondrial_candidates_summary.tsv`

Tabela com os candidatos. Campos importantes:

```text
contig                      nome do contig
contig_len                  tamanho do contig
unique_nt_queries           quantos CDS nucleotídicos bateram no contig
unique_prot_queries         quantas proteínas bateram no contig
total_bitscore              soma dos bitscores dos hits
best_nt_pident              melhor identidade no BLASTN
best_prot_pident            melhor identidade no TBLASTN
circular_overlap_len        tamanho da sobreposição terminal detectada
circular_by_terminal_overlap True/False para evidência simples de circularidade
circular_by_flye            True/False quando Flye informa circularidade
```

### 4. `coverage_summary.tsv`

Tabela de cobertura por contig candidato. Campos importantes:

```text
mean_depth       cobertura média
median_depth     cobertura mediana
min_depth        menor profundidade observada
max_depth        maior profundidade observada
breadth_1x       porcentagem do contig coberta por pelo menos 1 read
breadth_10x      porcentagem do contig coberta por pelo menos 10 reads
breadth_20x      porcentagem do contig coberta por pelo menos 20 reads
```

Uma montagem candidata forte deve ter, idealmente:

- breadth 1x próximo de 100%;
- cobertura média consistente;
- ausência de quedas longas para 0x;
- muitos genes mitocondriais detectados;
- contig único ou poucos contigs fáceis de conectar;
- evidência de circularidade, quando possível.

## Illumina ou PacBio: como decidir?

### Illumina

Vantagens:

- alta acurácia por base;
- boa para corrigir SNPs/indels;
- geralmente produz cobertura confiável.

Limitações:

- reads curtas podem fragmentar regiões repetitivas;
- mitogenomas com repeats ou regiões AT-rich podem sair em vários contigs;
- circularização pode ser difícil.

### PacBio

Vantagens:

- reads longas podem atravessar repeats;
- maior chance de obter um contig mitocondrial único;
- melhor para avaliar circularidade e estrutura.

Limitações:

- se for PacBio raw/CLR, o erro por base é maior;
- pode precisar de cobertura suficiente;
- se houver pouco DNA mitocondrial nas reads totais, a montagem pode ser instável.

### Estratégia prática recomendada

Rode os dois modos e compare:

```bash
# Illumina
python mito_pipeline.py \
  --mode illumina \
  --nt-ref Agabis_H97_nt.fasta \
  --prot-ref Agabis_H97_prot.fasta \
  --r1 A13_20_illumina_1_trimmed.fastq.gz \
  --r2 A13_20_illumina_2_trimmed.fastq.gz \
  --outdir mito_illumina \
  --threads 20

# PacBio
python mito_pipeline.py \
  --mode pacbio \
  --nt-ref Agabis_H97_nt.fasta \
  --prot-ref Agabis_H97_prot.fasta \
  --pacbio A13_20_pb.fastq.gz \
  --pacbio-type raw \
  --flye-genome-size 80k \
  --outdir mito_pacbio \
  --threads 20
```

Depois compare:

```text
mito_illumina/REPORT.md
mito_pacbio/REPORT.md
```

Escolha a melhor montagem com base em:

1. número de contigs mitocondriais candidatos;
2. tamanho do melhor candidato;
3. quantidade de genes/CDS detectados;
4. cobertura média e breadth;
5. uniformidade da cobertura;
6. evidência de circularidade;
7. coerência biológica após anotação.

## Próximos passos depois do pipeline

Este pipeline entrega candidatos mitocondriais e métricas de suporte. Para entregar um genoma mitocondrial final, recomenda-se ainda:

1. Abrir o BAM em IGV e verificar cobertura e regiões problemáticas.
2. Visualizar grafo de montagem quando possível, por exemplo no Bandage.
3. Fazer dotplot do melhor candidato para detectar repeats/sobreposição terminal.
4. Anotar o mitogenoma com ferramenta específica, como MITOS ou MFannot.
5. Conferir genes esperados: `cox1`, `cox2`, `cox3`, `cob`, `nad`, `atp`, rRNAs e tRNAs.
6. Comparar tamanho, conteúdo gênico e ordem gênica com mitogenomas próximos.

## Problemas comuns

### Nenhum contig mitocondrial foi encontrado

Possíveis causas:

- baixa cobertura mitocondrial;
- referência de Agabis muito distante;
- montagem fragmentada;
- filtros de BLAST rígidos demais.

Tente:

```bash
--min-nt-pident 60 --min-prot-pident 25 --min-hit-qcov 0.20
```

### Muitos contigs candidatos

Possíveis causas:

- montagem fragmentada;
- repeats;
- sequências mitocondriais transferidas para o núcleo;
- hits curtos pouco informativos.

Tente filtros mais rígidos:

```bash
--min-hit-qcov 0.50 --min-bitscore 100
```

### PacBio gerou montagem ruim

Teste outro tamanho esperado:

```bash
--flye-genome-size 40k
--flye-genome-size 120k
```

Ou teste o modo metagenômico:

```bash
--flye-meta
```

### Illumina gerou muitos contigs

Isso pode ser normal para regiões repetitivas. Compare com PacBio. Se PacBio gerar um contig único, ele pode ser melhor para estrutura; Illumina pode ser usada depois para validar/corrigir bases.

## Observação final

O arquivo mais importante para começar a discussão com o professor é:

```text
REPORT.md
```

Ele resume o que foi montado, quais contigs parecem mitocondriais, qual a cobertura e se há evidência de circularidade.

---

# Etapa adicional: extrair reads mitocondriais e remontar

Depois que o `REPORT.md` indicar contigs candidatos mitocondriais, você pode rodar uma segunda etapa para recuperar reads mitocondriais das reads originais e remontar apenas esse subconjunto enriquecido.

Script:

```text
mito_extract_reassemble.py
```

README específico:

```text
README_REASSEMBLY.md
```

Exemplo para Illumina, usando o FASTA de candidatos da primeira etapa como isca:

```bash
python mito_extract_reassemble.py \
  --mode illumina \
  --bait-fasta mito_illumina/03_mito_candidates/mitochondrial_candidates.fasta \
  --nt-ref Agabis_H97_nt.fasta \
  --prot-ref Agabis_H97_prot.fasta \
  --r1 A13_20_illumina_1_trimmed.fastq.gz \
  --r2 A13_20_illumina_2_trimmed.fastq.gz \
  --outdir mito_illumina_reassembly \
  --threads 20 \
  --evaluate-with-blast
```

Exemplo mais sensível, incluindo também os CDS nucleotídicos como isca:

```bash
python mito_extract_reassemble.py \
  --mode illumina \
  --bait-fasta mito_illumina/03_mito_candidates/mitochondrial_candidates.fasta \
  --nt-ref Agabis_H97_nt.fasta \
  --prot-ref Agabis_H97_prot.fasta \
  --include-nt-ref-as-bait \
  --r1 A13_20_illumina_1_trimmed.fastq.gz \
  --r2 A13_20_illumina_2_trimmed.fastq.gz \
  --outdir mito_illumina_reassembly_sensitive \
  --threads 20 \
  --evaluate-with-blast
```

No modo Illumina, essa etapa usa **BWA-MEM** para mapear contra a isca e **SPAdes** para remontar as reads extraídas. No modo PacBio, usa **minimap2** e **Flye**.

---

# Etapa adicional: curadoria mitogenômica

A versão atual do pacote também inclui o script `curate.py`, que adiciona funções inspiradas no MitoCurator para curar um mitogenoma candidato.

Ele deve ser usado depois que você tiver um FASTA candidato ou, idealmente, um GenBank anotado.

Exemplo rápido com o melhor contig candidato:

```bash
seqkit grep -n -p 'NODE_131_length_46048_cov_410.638120' \
  mito_illumina/03_mito_candidates/mitochondrial_candidates.fasta \
  > NODE_131_mitogenome_candidate.fasta

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

Leia o manual específico em:

```bash
less README_CURATION.md
```

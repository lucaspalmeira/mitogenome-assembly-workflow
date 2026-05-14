# Fluxo de Trabalho de Montagem do Mitogenoma

Pipeline para identificar scaffolds candidatos do genoma mitocondrial a partir de uma montagem de genoma híbrido, recrutando reads mitocondriais e realizando uma remontagem híbrida direcionada usando reads curtos da Illumina e reads longos da PacBio.

### Fluxo de trabalho geral

O fluxo de trabalho começa com uma montagem de genoma gerada com o SPAdes, representada aqui por `scaffolds_1000bp.fasta`, e um conjunto de sequências de CDS mitocondriais em formato FASTA de nucleotídeos e proteínas.

O primeiro script, `find_mitogenome_tblastn.py`, usa proteínas mitocondriais como consultas em uma busca `tBLASTn` contra a montagem do genoma. Esta etapa identifica scaffolds que contêm genes mitocondriais putativos. As ferramentas de linha de comando do NCBI BLAST+ suportam formatos de saída tabulares personalizados, que são usados ​​aqui para gerar tabelas detalhadas de hits, coordenadas, identidade, cobertura, valor E e bitscore.

O segundo script, `mito_recruit_reassemble_pipeline.py`, utiliza os scaffolds candidatos da primeira etapa para recrutar reads dos dados de sequenciamento originais. Reads curtos são mapeados com `bowtie2` ou `bwa`, enquanto reads PacBio são mapeados com `minimap2`. Para reads PacBio CLR, o minimap2 fornece a predefinição `map-pb` para mapear reads longos contínuos PacBio mais antigos a um genoma de referência. Os reads recrutados são então remontados com o SPAdes, utilizando reads curtos pareados e reads longos PacBio como dados suplementares. O SPAdes suporta montagem híbrida através da opção `--pacbio` para reads PacBio CLR.

### Diagrama de fluxo de trabalho

```sereia
fluxograma TD 

    A[Arquivos de entrada] -> B[Montagem do genoma híbrido] 
    A -> C [arquivos CDS FASTA mitocondriais] 

    B --> B1[scaffolds_1000bp.fasta] 
    C --> C1[Agabis_H97_prot.fasta] 
    C --> C2[Agabis_H97_nt.fasta] 

    B1 --> D[Etapa 1: pesquisa tBLASTn] 
    C1 --> D 
    C2 --> D 

    D -> E[find_mitogenoma_tblastn.py] 

    E -> F1[tblastn_raw.tsv] 
    E -> F2[tblastn_hits_all.csv] 
    E -> F3[tblastn_hits_filtered.csv] 
    E--> F4[best_hit_per_query.csv]

    E --> F5[scaffold_summary.csv]

    E --> F6[candidate_mitogene_scaffolds.fasta]

    F5 --> G[Etapa 2: Selecionar scaffolds candidatos]

    B1 --> G

    G --> H[mito_recruit_reassemble_pipeline.py]

    I1[FASTQ original Illumina R1] --> H

    I2[FASTQ original Illumina R2] --> H

    I3[Leituras PacBio originais] --> H

    H --> J1[Mapear leituras curtas contra scaffolds candidatos]

    H --> J2[Mapear leituras PacBio contra scaffolds candidatos]

    J1 --> K1[short_reads_vs_candidates.sorted.bam]

    J2 --> K2[pacbio_reads_vs_candidates.sorted.bam]

    K1 --> L1[Extrair IDs de leituras curtas mapeadas]

    K2 --> L2[Extrair IDs de leituras PacBio mapeadas]

    L1 --> M1[recruited_R1.fastq.gz]

    L1 --> M2[recruited_R2.fastq.gz]

    L2 --> M3[recruited_pacbio_reads.fastq.gz ou .fasta.gz]

    M1 --> N[Etapa 3: Remontagem híbrida direcionada]

    M2 --> N
    M3 --> N

    N --> O[Remontagem SPAdes com leituras recrutadas]

    O --> P1[spades_reassembly/scaffolds.fasta]

    O --> P2[spades_reassembly/contigs.fasta]

    P1 --> Q[Etapa 4: Validar candidatos remontados]
    C1 --> Q
    Q --> R[Executar find_mitogenom_tblastn.py novamente na nova montagem]

    R --> S[Scaffolds finais do mitogenoma candidato]
```

---

Ordem de execução

Este fluxo de trabalho deve ser executado na seguinte ordem:

1. Prepare os arquivos de entrada.

2. Execute find_mitogenom_tblastn.py.

3. Inspecione scaffold_summary.csv.

4. Execute mito_recruit_reassemble_pipeline.py.

5. Inspecione a remontagem das leituras recrutadas.

6. Execute find_mitogenom_tblastn.py novamente na nova montagem.

7. Avalie os scaffolds candidatos finais do mitogenoma.

---

| File                                               | Description                                                                                                                                                 |
| -------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `scaffolds_1000bp.fasta`                           | Montagem do genoma gerada com SPAdes. Este arquivo é usado como referência para o primeiro `tBLASTn`. search.                                                   |
| `Agabis_H97_prot.fasta`                            | Sequências de proteínas de CDS mitocondrial. Estas são usadas como consultas no `tBLASTn`.                                                                           |
| `Agabis_H97_nt.fasta`                              | Sequências de nucleotídeos de CDS mitocondrial. Estas não são usadas diretamente pelo `tBLASTn`, mas são extraídas posteriormente para comparação com as sequências de proteínas correspondentes. |
| `reads_R1.fastq.gz`                                | Leituras Illumina originais utilizadas na montagem híbrida.                                                                                                |
| `reads_R2.fastq.gz`                                | Sequências de leitura reversa originais da Illumina utilizadas na montagem híbrida.                                                                                                |
| `pacbio_reads.fastq.gz` or `pacbio_reads.fasta.gz` | Sequências originais de leitura PacBio utilizadas na montagem híbrida.                                                                                                          |

---

### Etapa 1: Identificar possíveis estruturas mitocondriais com o tBLASTn O primeiro passo é pesquisar as proteínas mitocondriais no conjunto do genoma.

```bash
python find_mitogenome_tblastn.py \
  --assembly scaffolds_1000bp.fasta \
  --cds-prot Agabis_H97_prot.fasta \
  --cds-nt Agabis_H97_nt.fasta \
  --threads 12 \
  --outdir tblastn_mitogenome_results
``

Esta etapa cria:

| Arquivo de saída | Descrição |
| ------------------------------------------------------------------------- | ------------------------------------------------------------------------------------- |
| `tblastn_mitogenome_results/tables/tblastn_raw.tsv` | Saída bruta do `tBLASTn`. |
| `tblastn_mitogenome_results/tables/tblastn_hits_all.csv` | Todos os resultados do BLAST em formato CSV. |
| `tblastn_mitogenome_results/tables/tblastn_hits_filtered.csv` | Resultados do BLAST filtrados após a aplicação de limiares de identidade, cobertura e pontuação de bits. |
| `tblastn_mitogenome_results/tables/best_hit_per_query.csv` | Melhor resultado para cada consulta de proteína mitocondrial. |
| `tblastn_mitogenom_results/tables/scaffold_summary.csv` | Resumo das evidências mitocondriais por scaffold. Este é o principal resultado desta etapa. |
| `tblastn_mitogenom_results/fasta/candidate_mitogenom_scaffolds.fasta` | Arquivo FASTA contendo os scaffolds candidatos identificados pelo `tBLASTn`. |
| `tblastn_mitogenom_results/fasta/tblastn_hit_regions_plus_padding.fasta` | Arquivo FASTA contendo as regiões ao redor de cada resultado do BLAST. |
| `tblastn_mitogenom_results/fasta/matched_cds_proteins.fasta` | Sequências de CDS de proteínas que corresponderam à montagem. |
| `tblastn_mitogenom_results/fasta/matched_cds_nucleotides.fasta` | Sequências de CDS de nucleotídeos correspondentes às proteínas correspondentes.

---


| Etapa                                           | Ferramenta                                                    |
| ----------------------------------------------- | ------------------------------------------------------------- |
| Escolher scaffolds candidatos ao mitogenoma     | resultado do `tBLASTn`, principalmente `scaffold_summary.csv` |
| Extrair os scaffolds candidatos                 | Python                                                        |
| Mapear short reads contra esses scaffolds       | `bowtie2` ou `bwa`                                            |
| Mapear long reads PacBio contra esses scaffolds | `minimap2 -ax map-pb`                                         |
| Extrair reads que mapearam                      | `samtools` + Python                                           |
| Remontar somente reads recrutadas               | `spades.py --pacbio`                                          |
| Tentar juntar contigs candidatos, opcional      | `CAP3`                                                        |

---

# Instalar/carregar dependências

No conda:

```bash
conda install -c bioconda bowtie2 bwa samtools minimap2 spades
```

O `CAP3` é opcional. O pipeline roda sem ele.

---

### Como rodar com PacBio

Você precisa dos arquivos originais usados na montagem híbrida:

```text
short reads R1
short reads R2
PacBio reads
scaffolds_1000bp.fasta
scaffold_summary.csv
```

Exemplo:

```bash
python mito_recruit_reassemble_pipeline.py \
  --assembly scaffolds_1000bp.fasta \
  --scaffold-summary tblastn_mitogenome_results/tables/scaffold_summary.csv \
  --short-r1 reads_R1.fastq.gz \
  --short-r2 reads_R2.fastq.gz \
  --pacbio-reads pacbio_reads.fastq.gz \
  --threads 20 \
  --memory 80 \
  --top-n 4 \
  --outdir mito_recruitment_pipeline_pacbio
```

Se suas PacBio reads estiverem em FASTA, também funciona:

```bash
python mito_recruit_reassemble_pipeline.py \
  --assembly scaffolds_1000bp.fasta \
  --scaffold-summary tblastn_mitogenome_results/tables/scaffold_summary.csv \
  --short-r1 reads_R1.fastq.gz \
  --short-r2 reads_R2.fastq.gz \
  --pacbio-reads pacbio_reads.fasta.gz \
  --threads 20 \
  --memory 80 \
  --top-n 4 \
  --outdir mito_recruitment_pipeline_pacbio
```

O script detecta automaticamente se o arquivo PacBio é `FASTQ` ou `FASTA`.


---

### O que o pipeline faz biologicamente

A lógica é:

```text
CDS mitocondriais
        ↓
tBLASTn contra scaffolds_1000bp.fasta
        ↓
identificação dos scaffolds com hits mitocondriais
        ↓
seleção dos melhores scaffolds candidatos
        ↓
mapeamento das reads originais contra esses scaffolds
        ↓
extração das reads que mapearam
        ↓
nova montagem usando apenas reads recrutadas
        ↓
candidatos mais limpos ao mitogenoma
```

Isso corresponde ao que foi dito no grupo: selecionar os contigs/scaffolds com hits mitocondriais, mapear as leituras contra eles, extrair essas leituras e remontar.

---

### Arquivos gerados

| Arquivo/pasta                                            | O que informa                                                                                        |
| -------------------------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| `selected_scaffolds.csv`                                 | Lista dos scaffolds selecionados a partir do `scaffold_summary.csv`. Por padrão, pega os 4 melhores. |
| `selected_candidate_scaffolds.fasta`                     | FASTA dos scaffolds candidatos usados como referência para recrutar reads.                           |
| `mapping_short/short_reads_vs_candidates.sorted.bam`     | Mapeamento das short reads contra os scaffolds candidatos.                                           |
| `mapping_short/short_reads_vs_candidates.flagstat.txt`   | Estatísticas do mapeamento das short reads.                                                          |
| `recruited_reads/mapped_short_read_ids.txt`              | IDs das short reads que mapearam nos scaffolds candidatos.                                           |
| `recruited_reads/recruited_R1.fastq.gz`                  | Reads R1 recrutadas para remontagem.                                                                 |
| `recruited_reads/recruited_R2.fastq.gz`                  | Reads R2 recrutadas para remontagem.                                                                 |
| `recruited_reads/recruited_short_reads_summary.csv`      | Quantidade de pares totais e pares recrutados.                                                       |
| `mapping_pacbio/pacbio_reads_vs_candidates.sorted.bam`   | Mapeamento das PacBio reads contra os scaffolds candidatos.                                          |
| `mapping_pacbio/pacbio_reads_vs_candidates.flagstat.txt` | Estatísticas do mapeamento das PacBio reads.                                                         |
| `recruited_reads/mapped_pacbio_read_ids.txt`             | IDs das PacBio reads que mapearam nos scaffolds candidatos.                                          |
| `recruited_reads/recruited_pacbio_reads.fastq.gz`        | PacBio reads recrutadas, se o arquivo original era FASTQ.                                            |
| `recruited_reads/recruited_pacbio_reads.fasta.gz`        | PacBio reads recrutadas, se o arquivo original era FASTA.                                            |
| `recruited_reads/recruited_pacbio_reads_summary.csv`     | Quantidade de PacBio reads analisadas e recrutadas.                                                  |
| `coverage/coverage_by_scaffold.csv`                      | Cobertura média e porcentagem coberta dos scaffolds candidatos.                                      |
| `coverage/depth_per_base.tsv`                            | Cobertura base a base dos scaffolds candidatos.                                                      |
| `spades_reassembly/scaffolds.fasta`                      | Nova montagem feita com as reads recrutadas. Principal resultado final.                              |
| `spades_reassembly/contigs.fasta`                        | Contigs da nova montagem com reads recrutadas.                                                       |
| `cap3/`                                                  | Saída opcional do CAP3, caso você use `--run-cap3`.                                                  |

---

### Como selecionar só 4 scaffolds

O padrão já é:

```bash
--top-n 4
```

Isso pega os 4 melhores scaffolds de acordo com:

```text
maior número de genes/CDS mitocondriais diferentes
maior cobertura somada dos hits
menor e-value
maior bitscore
```

Para usar só os 2 melhores:

```bash
--top-n 2
```

Para usar os 10 melhores:

```bash
--top-n 10
```

---

### Como selecionar scaffolds manualmente

Primeiro abra:

```bash
column -s, -t tblastn_mitogenome_results/tables/scaffold_summary.csv | less -S
```

Depois escolha os nomes dos scaffolds.

Exemplo:

```bash
python mito_recruit_reassemble_pipeline.py \
  --assembly scaffolds_1000bp.fasta \
  --scaffold-summary tblastn_mitogenome_results/tables/scaffold_summary.csv \
  --short-r1 reads_R1.fastq.gz \
  --short-r2 reads_R2.fastq.gz \
  --pacbio-reads pacbio_reads.fastq.gz \
  --manual-scaffolds NODE_12_length_45000_cov_80,NODE_33_length_12000_cov_65,NODE_41_length_8000_cov_70,NODE_90_length_5000_cov_60 \
  --threads 20 \
  --memory 80 \
  --outdir mito_recruitment_pipeline_manual_pacbio
```

---

### Rodar só recrutamento, sem remontar

Para testar primeiro se o mapeamento e a extração funcionam:

```bash
python mito_recruit_reassemble_pipeline.py \
  --assembly scaffolds_1000bp.fasta \
  --scaffold-summary tblastn_mitogenome_results/tables/scaffold_summary.csv \
  --short-r1 reads_R1.fastq.gz \
  --short-r2 reads_R2.fastq.gz \
  --pacbio-reads pacbio_reads.fastq.gz \
  --threads 20 \
  --top-n 4 \
  --skip-assembly \
  --outdir mito_recruitment_pipeline_pacbio_test
```

Depois olhe:

```bash
cat mito_recruitment_pipeline_pacbio_test/mapping_short/short_reads_vs_candidates.flagstat.txt
cat mito_recruitment_pipeline_pacbio_test/mapping_pacbio/pacbio_reads_vs_candidates.flagstat.txt
```

E veja quantas reads foram recrutadas:

```bash
cat mito_recruitment_pipeline_pacbio_test/recruited_reads/recruited_short_reads_summary.csv
cat mito_recruitment_pipeline_pacbio_test/recruited_reads/recruited_pacbio_reads_summary.csv
```

---

### Rodar CAP3 opcionalmente

Usar `CAP3` para juntar contigs que têm genes mitocondriais.

No pipeline, isso é opcional:

```bash
python mito_recruit_reassemble_pipeline.py \
  --assembly scaffolds_1000bp.fasta \
  --scaffold-summary tblastn_mitogenome_results/tables/scaffold_summary.csv \
  --short-r1 reads_R1.fastq.gz \
  --short-r2 reads_R2.fastq.gz \
  --pacbio-reads pacbio_reads.fastq.gz \
  --threads 20 \
  --memory 80 \
  --top-n 4 \
  --run-cap3 \
  --outdir mito_recruitment_pipeline_pacbio_cap3
```

O `CAP3` poderá ser uma etapa complementar:

```text
selecionar scaffolds mitocondriais
mapear reads
extrair reads mapeadas
remontar com SPAdes usando --pacbio
```

---

### Depois da remontagem

O resultado principal será:

```bash
mito_recruitment_pipeline_pacbio/spades_reassembly/scaffolds.fasta
```

Depois rode novamente o script de `tBLASTn` contra essa nova montagem:

```bash
python find_mitogenome_tblastn.py \
  --assembly mito_recruitment_pipeline_pacbio/spades_reassembly/scaffolds.fasta \
  --cds-prot Agabis_H97_prot.fasta \
  --cds-nt Agabis_H97_nt.fasta \
  --threads 12 \
  --outdir tblastn_after_pacbio_recruitment
```

O que você espera ver no novo `scaffold_summary.csv`:

```text
menos scaffolds candidatos
mais genes mitocondriais no mesmo scaffold
maior continuidade da montagem mitocondrial
melhor evidência de qual scaffold representa o mitogenoma
```

O arquivo mais importante depois disso será:

```bash
tblastn_after_pacbio_recruitment/tables/scaffold_summary.csv
```


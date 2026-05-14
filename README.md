

| Etapa                                           | Ferramenta                                                    |
| ----------------------------------------------- | ------------------------------------------------------------- |
| Escolher scaffolds candidatos ao mitogenoma     | resultado do `tBLASTn`, principalmente `scaffold_summary.csv` |
| Extrair os scaffolds candidatos                 | Python                                                        |
| Mapear short reads contra esses scaffolds       | `bowtie2` ou `bwa`                                            |
| Mapear long reads PacBio contra esses scaffolds | `minimap2 -ax map-pb`                                         |
| Extrair reads que mapearam                      | `samtools` + Python                                           |
| Remontar somente reads recrutadas               | `spades.py --pacbio`                                          |
| Tentar juntar contigs candidatos, opcional      | `CAP3`                                                        |


Baixe e execute o script <a href=''>mito_recruit_reassemble_pipeline.py</>


```bash
chmod +x mito_recruit_reassemble_pipeline.py
```

---

# Instalar/carregar dependências

No conda:

```bash
conda install -c bioconda bowtie2 bwa samtools minimap2 spades
```

O `CAP3` é opcional. O pipeline roda sem ele.

---

# Como rodar com PacBio

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

# O que o pipeline faz biologicamente

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

# Arquivos gerados

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

# Como selecionar só 4 scaffolds

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

# Como selecionar scaffolds manualmente

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

# Rodar só recrutamento, sem remontar

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

# Rodar CAP3 opcionalmente

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

# Depois da remontagem

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


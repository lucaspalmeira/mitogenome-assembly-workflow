# Etapa 2: enriquecimento de reads mitocondriais e remontagem

Este README explica o script `mito_extract_reassemble.py`, que deve ser usado **depois** da primeira etapa (`mito_pipeline.py`).

A ideia é simples:

```text
contigs mitocondriais candidatos da primeira montagem
        |
        v
isca mitocondrial
        |
        v
mapeamento das reads originais contra a isca
        |
        v
extração das reads que mapearam na mitocôndria
        |
        v
remontagem apenas das reads enriquecidas
        |
        v
avaliação de cobertura, BLAST/TBLASTN e circularidade
```

No seu caso, o principal arquivo de isca é:

```text
mito_illumina/03_mito_candidates/mitochondrial_candidates.fasta
```

Esse arquivo contém os contigs candidatos detectados na primeira montagem.

---

## Por que fazer esta segunda etapa?

A primeira etapa montou as reads totais e depois procurou os contigs mitocondriais. Isso funcionou, mas ainda pode deixar o mitogenoma fragmentado ou sem circularização clara.

A segunda etapa tenta enriquecer as reads mitocondriais antes da montagem. Isso reduz a complexidade do conjunto de reads e pode ajudar o assembler a montar uma molécula mitocondrial mais contínua.

---

## O pipeline “sabe” que o genoma é circular?

Biologicamente, sim, você pode informar e interpretar o resultado como uma molécula circular esperada. Computacionalmente, os assemblers normalmente escrevem uma molécula circular como uma sequência linear, com um ponto arbitrário de quebra.

Portanto, o objetivo não é obrigar o FASTA a parecer circular, mas demonstrar evidências de fechamento:

1. um contig principal com tamanho compatível;
2. genes mitocondriais esperados;
3. cobertura uniforme;
4. reads atravessando a junção entre final e início;
5. overlap terminal ou indicação de circularidade pelo assembler;
6. grafo de montagem compatível com círculo, quando disponível.

---

## Decisão de alinhador

### Illumina

Use **BWA-MEM**.

Justificativa: suas reads Illumina têm aproximadamente 150 bp paired-end. BWA-MEM é apropriado para reads Illumina nessa faixa e preserva bem a lógica paired-end. O script usa BWA-MEM para mapear as reads contra a isca e depois extrai o par completo quando pelo menos um mate mapeia.

### PacBio

Use **minimap2**.

Justificativa: para reads longas PacBio, minimap2 é a escolha mais adequada e possui presets específicos para PacBio raw (`map-pb`) e PacBio HiFi (`map-hifi`).

---

## Instalação

Atualize o ambiente com BWA:

```bash
mamba env create -f environment.yml
mamba activate mito_assembly
```

Ou, se o ambiente já existe:

```bash
mamba activate mito_assembly
mamba install -c bioconda bwa
```

Confira:

```bash
which bwa
which spades.py
which flye
which minimap2
which samtools
which blastn
which tblastn
```

---

## Uso recomendado para o seu resultado Illumina

Como seu `REPORT.md` mostrou um candidato muito forte (`NODE_131_length_46048_cov_410.638120`), o primeiro teste deve usar os contigs candidatos da primeira etapa como isca.

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

### Versão mais sensível

Se você quiser capturar reads de genes que talvez não tenham entrado nos contigs candidatos, inclua também os CDS de `Agabis_H97_nt.fasta` como isca:

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

Essa opção é mais sensível, mas pode puxar reads de regiões conservadas semelhantes à mitocôndria, então a interpretação da cobertura e do BLAST final fica ainda mais importante.

---

## Uso com PacBio

Se você quiser fazer a mesma estratégia usando as reads PacBio:

```bash
python mito_extract_reassemble.py \
  --mode pacbio \
  --bait-fasta mito_illumina/03_mito_candidates/mitochondrial_candidates.fasta \
  --nt-ref Agabis_H97_nt.fasta \
  --prot-ref Agabis_H97_prot.fasta \
  --pacbio A13_20_pb.fastq.gz \
  --pacbio-type raw \
  --flye-genome-size 80k \
  --outdir mito_pacbio_reassembly \
  --threads 20 \
  --evaluate-with-blast
```

Se suas reads forem HiFi:

```bash
--pacbio-type hifi
```

---

## Estrutura de saída

Exemplo para `--outdir mito_illumina_reassembly`:

```text
mito_illumina_reassembly/
├── commands.log
├── REPORT_REASSEMBLY.md
├── 00_bait_reference/
│   └── bait_reference.fasta
├── 01_map_reads_to_bait/
│   ├── illumina_reads_vs_bait.sorted.bam
│   ├── illumina_reads_vs_bait.sorted.bam.bai
│   └── samtools_flagstat.txt
├── 02_extracted_mito_reads/
│   ├── mapped_read_names.txt
│   ├── mito_enriched_R1.fastq.gz
│   ├── mito_enriched_R2.fastq.gz
│   └── extraction_summary.tsv
├── 03_reassembly_spades_mito_reads/
│   ├── contigs.fasta
│   └── scaffolds.fasta
└── 04_evaluate_reassembly/
    ├── circularity_summary.tsv
    ├── blast/
    │   ├── blastn_nt_vs_reassembly.tsv
    │   ├── tblastn_prot_vs_reassembly.tsv
    │   └── reassembly_mito_candidates_summary.tsv
    └── coverage_original_reads_vs_reassembly/
        ├── original_reads_vs_reassembly.sorted.bam
        ├── original_reads_vs_reassembly.sorted.bam.bai
        ├── depth_per_base.tsv
        ├── samtools_coverage.tsv
        └── coverage_summary.tsv
```

---

## Arquivo principal para olhar

```text
REPORT_REASSEMBLY.md
```

Ele resume:

- quantas reads foram extraídas;
- qual porcentagem das reads originais foi usada;
- estatísticas da remontagem;
- candidatos detectados por BLAST/TBLASTN;
- cobertura das reads originais contra a remontagem;
- evidência automática de circularidade.

---

## Como interpretar

### Caso ideal

Você quer ver algo como:

```text
1 contig principal de ~46 kb ou próximo disso
muitos genes mitocondriais detectados
cobertura média alta
breadth_1x = 100%
breadth_20x = 100% ou próximo
circular_by_terminal_overlap = True
```

Nesse caso, você tem forte evidência de mitogenoma completo ou quase completo.

### Se ainda sair fragmentado

Se a remontagem ainda gerar vários contigs, isso pode indicar:

- repeats longos que Illumina não atravessa;
- regiões com viés de cobertura;
- mistura de variantes mitocondriais;
- contigs alternativos;
- NUMTs ou regiões semelhantes à mitocôndria no genoma nuclear;
- necessidade de usar PacBio para fechar a estrutura.

Nesse cenário, rode a versão PacBio ou faça montagem híbrida/curadoria com grafo.

---

## Sobre circularização final

Mesmo se o contig é circular, o FASTA final ainda será uma sequência linear. Para entrega, normalmente você pode:

1. escolher um ponto de início biológico, por exemplo `cox1`, `cox2`, `cob` ou outro gene de referência;
2. rotacionar a sequência para começar nesse gene;
3. remover redundância terminal se houver overlap duplicado;
4. validar a junção com reads mapeadas.

Esse script detecta overlap terminal simples, mas não substitui inspeção manual em IGV, dotplot ou Bandage.

---

## Comando de teste seco

Antes de rodar tudo:

```bash
python mito_extract_reassemble.py \
  --mode illumina \
  --bait-fasta mito_illumina/03_mito_candidates/mitochondrial_candidates.fasta \
  --nt-ref Agabis_H97_nt.fasta \
  --prot-ref Agabis_H97_prot.fasta \
  --r1 A13_20_illumina_1_trimmed.fastq.gz \
  --r2 A13_20_illumina_2_trimmed.fastq.gz \
  --outdir teste_reassembly_dryrun \
  --threads 20 \
  --evaluate-with-blast \
  --dry-run
```

Isso cria `commands.log` sem executar a análise pesada.

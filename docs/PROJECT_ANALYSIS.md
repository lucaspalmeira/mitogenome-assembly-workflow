# Análise do Projeto e Estratégia de Integração

Este projeto contém dois fluxos relacionados para recuperar e avaliar candidatos a mitogenoma. Eles partem de premissas diferentes, mas são complementares:

- o workflow da raiz começa com uma montagem genômica já existente, identifica scaffolds mitocondriais por tBLASTn, recruta reads contra esses scaffolds e faz uma remontagem direcionada;
- `toolkit/` monta reads de novo, avalia contigs candidatos com BLASTN/TBLASTN, calcula estatísticas/cobertura/circularidade e adiciona uma etapa de enriquecimento/remontagem e curadoria.

A integração implementada preserva os scripts existentes e adiciona uma camada orientada a objetos em `workflow/`, com uma CLI unificada em `run_workflow.py`.

## Pipeline da Raiz

### `find_mitogenome_tblastn.py`

Função principal: detectar scaffolds candidatos ao mitogenoma usando proteínas mitocondriais como query em `tblastn`.

Fluxo interno:

1. lê a assembly FASTA;
2. cria um banco BLAST nucleotídico com `makeblastdb`;
3. roda `tblastn` contra a assembly;
4. parseia a tabela BLAST;
5. filtra hits por identidade, cobertura do HSP e bitscore;
6. agrupa hits por scaffold;
7. calcula evidência por scaffold: número de hits, queries únicas, intervalo coberto, melhor e-value e maior bitscore;
8. extrai FASTA dos scaffolds candidatos;
9. extrai regiões de hits com padding;
10. extrai CDS/proteínas de referência que tiveram hit.

Papel no workflow: etapa de descoberta funcional. Ela cria o arquivo ponte `tables/scaffold_summary.csv`, usado pela etapa de recrutamento.

### `mito_recruit_reassemble_pipeline.py`

Função principal: selecionar scaffolds candidatos, mapear reads originais contra esses scaffolds, extrair reads recrutadas e remontar esse subconjunto com SPAdes.

Fluxo interno:

1. carrega `scaffold_summary.csv`;
2. ranqueia scaffolds por queries únicas, bases cobertas, e-value, bitscore e número de hits;
3. permite seleção manual de scaffolds;
4. extrai os scaffolds selecionados da assembly original;
5. mapeia Illumina com Bowtie2 ou BWA;
6. mapeia PacBio com minimap2;
7. extrai IDs das reads mapeadas;
8. extrai pares Illumina completos quando pelo menos um mate mapeou;
9. extrai reads PacBio em FASTQ ou FASTA;
10. calcula cobertura das short reads sobre os scaffolds candidatos;
11. roda SPAdes com `--careful` e, quando houver PacBio, `--pacbio`;
12. opcionalmente roda CAP3.

Ajustes adicionados:

- `--pacbio-type raw|hifi|corrected`, usando `map-hifi` para HiFi e `map-pb` para raw/corrected;
- `--min-mapq`, aplicado ao `samtools view` na extração dos IDs de reads recrutadas.

## Toolkit Paralelo em `toolkit/`

### `mito_pipeline.py`

Fluxo de primeira montagem e avaliação:

1. monta Illumina com SPAdes ou PacBio com Flye;
2. calcula estatísticas de FASTQ/FASTA;
3. roda BLASTN dos CDS nucleotídicos e TBLASTN das proteínas contra a assembly;
4. seleciona candidatos mitocondriais;
5. avalia circularidade por overlap terminal e flags do Flye;
6. remapeia reads contra candidatos;
7. calcula cobertura média, mediana, mínima, máxima e breadth em 1x/5x/10x/20x;
8. gera `REPORT.md`.

Funcionalidades úteis para a raiz:

- relatório Markdown;
- BLASTN + TBLASTN como evidência complementar;
- métricas de cobertura mais completas;
- circularidade por overlap terminal;
- suporte explícito a PacBio HiFi/raw/corrected.

### `mito_extract_reassemble.py`

Fluxo de enriquecimento:

1. combina FASTAs de isca;
2. opcionalmente adiciona CDS nucleotídicos como isca;
3. mapeia reads originais contra as iscas;
4. extrai reads mitocondriais candidatas com MAPQ mínimo;
5. remonta Illumina com SPAdes ou PacBio com Flye;
6. avalia a remontagem com circularidade, BLAST e cobertura;
7. gera `REPORT_REASSEMBLY.md`.

Funcionalidades úteis para a raiz:

- conceito de isca flexível;
- MAPQ mínimo para extração;
- avaliação automática da remontagem;
- relatório final da etapa enriquecida.

### `curate.py`

Fluxo de curadoria:

1. aceita FASTA ou GenBank candidato;
2. calcula estatísticas da sequência;
3. testa circularidade por overlap terminal;
4. busca genes por BLASTN/TBLASTN;
5. resume presença/ausência de genes esperados;
6. detecta múltiplos hits que podem sugerir duplicações/repeats;
7. rotaciona sequência circular para um gene definido;
8. em GenBank, avalia CDS, tRNAs, rRNAs, regiões intergênicas e regiões AT-rich;
9. gera `CURATION_REPORT.md`.

Funcionalidades úteis para a raiz:

- etapa final de diagnóstico biológico;
- rotação para padronizar o candidato circular;
- checagem de genes ausentes, duplicados ou problemáticos.

## Lacunas Identificadas

O workflow da raiz era funcional, mas exigia passos manuais após a remontagem:

- a validação pós-SPAdes precisava ser chamada manualmente com `find_mitogenome_tblastn.py`;
- não havia relatório único amarrando descoberta, recrutamento, remontagem e validação;
- a curadoria final do pipeline paralelo não estava conectada ao workflow da raiz;
- PacBio HiFi usava o mesmo preset de PacBio raw;
- a extração de reads mapeadas não tinha filtro MAPQ configurável.

## Revisão da Lógica Bioinformática

A lógica geral do workflow está coerente para recuperação de mitogenoma a partir de uma assembly nuclear/organellar mista:

1. a descoberta inicial por `tBLASTn` é adequada quando a distância evolutiva pode ser grande, pois proteínas conservadas toleram mais divergência do que CDS nucleotídicos;
2. o `scaffold_summary.csv` é a ponte correta entre evidência funcional e recrutamento de reads;
3. mapear reads originais contra scaffolds candidatos e extrair o par Illumina completo quando um mate mapeia é biologicamente razoável, porque preserva informação paired-end para a remontagem;
4. a validação pós-remontagem com novo `tBLASTn` não é repetição inútil: o alvo mudou de assembly original para assembly enriquecida;
5. BLASTN e TBLASTN no `toolkit/` são complementares, não redundantes, porque capturam evidência nucleotídica mais específica e evidência proteica mais sensível;
6. a circularidade por overlap terminal deve ser tratada como triagem, não como prova final sem cobertura na junção, grafo de montagem ou inspeção manual;
7. `--min-mapq` melhora a seletividade do recrutamento quando há NUMTs, repeats ou scaffolds conservados demais;
8. `--pacbio-type hifi` agora evita usar preset de PacBio raw em reads HiFi.

Duplicações removidas nesta revisão:

- parsers e escritores FASTA/FASTQ dos scripts da raiz foram movidos para `workflow/bioio.py`;
- escrita FASTA, teste de overlap terminal e parsing de circularidade do Flye foram movidos para `toolkit/common.py`;
- os diretórios de execução agora usam nomes curtos e não repetem o objetivo do projeto.

As rotinas de BLAST, seleção de candidatos e cobertura ainda parecem semelhantes entre scripts, mas têm contratos diferentes de tabela, filtros e relatório. Elas não foram fundidas nesta etapa para não misturar evidências biologicamente distintas nem quebrar compatibilidade dos arquivos gerados.

## Integração Implementada

A nova camada fica em:

```text
workflow/
├── __init__.py
├── bioio.py
└── orchestrator.py
toolkit/
├── common.py
├── curate.py
├── mito_extract_reassemble.py
└── mito_pipeline.py
run_workflow.py
```

Classes principais:

- `WorkflowConfig`: parâmetros globais;
- `WorkflowOutputs`: arquivos importantes produzidos;
- `CommandRunner`: executor com `commands.log` reprodutível;
- `TblastnDiscoveryStage`: roda a descoberta inicial ou usa um `scaffold_summary.csv` existente;
- `RecruitReassembleStage`: roda o recrutamento/remontagem da raiz;
- `ReassemblyValidationStage`: valida automaticamente a assembly remontada por tBLASTn;
- `CurationStage`: opcionalmente chama `toolkit/curate.py`;
- `AssemblyWorkflow`: orquestra as etapas;
- `WorkflowReport`: gera `WORKFLOW_REPORT.md`.
- `workflow/bioio.py`: centraliza leitura/escrita FASTA/FASTQ e normalização de IDs de reads.
- `toolkit/common.py`: centraliza escrita FASTA, parsing de circularidade do Flye e teste de overlap terminal.

## Comando Unificado Recomendado

```bash
python run_workflow.py \
  --assembly scaffolds_1000bp.fasta \
  --cds-prot Agabis_H97_prot.fasta \
  --cds-nt Agabis_H97_nt.fasta \
  --short-r1 reads_R1.fastq.gz \
  --short-r2 reads_R2.fastq.gz \
  --pacbio-reads pacbio_reads.fastq.gz \
  --pacbio-type raw \
  --threads 20 \
  --memory 80 \
  --top-n 4 \
  --min-mapq 20 \
  --outdir analysis_run
```

Com curadoria final:

```bash
python run_workflow.py \
  --assembly scaffolds_1000bp.fasta \
  --cds-prot Agabis_H97_prot.fasta \
  --cds-nt Agabis_H97_nt.fasta \
  --short-r1 reads_R1.fastq.gz \
  --short-r2 reads_R2.fastq.gz \
  --pacbio-reads pacbio_reads.fastq.gz \
  --threads 20 \
  --memory 80 \
  --top-n 4 \
  --min-mapq 20 \
  --curate-final \
  --expected-from-ref \
  --rotate-to COX1 \
  --outdir curated_run
```

## Estrutura de Saída

```text
analysis_run/
├── workflow_config.json
├── commands.log
├── WORKFLOW_REPORT.md
├── 01_discovery/
├── 02_recruitment/
├── 03_validation/
└── 04_curate/                # apenas com --curate-final
```

## Decisão de Arquitetura

A integração foi feita como camada de orquestração, não como reescrita total, por três motivos:

1. os scripts existentes já codificam decisões biológicas e formatos de saída importantes;
2. o risco de regressão seria alto se toda a lógica de BLAST, mapeamento e extração fosse reimplementada de uma vez;
3. a camada OOP já cria pontos claros para evoluir o projeto sem quebrar os comandos antigos.

Próximos candidatos naturais para modularização interna:

- consolidar funções de BLAST e cobertura;
- padronizar relatórios Markdown;
- adicionar testes unitários para seleção de scaffolds, normalização de IDs e parsing de BLAST.

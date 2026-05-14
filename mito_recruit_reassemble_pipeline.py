#!/usr/bin/env python3

import argparse
import csv
import gzip
import shutil
import subprocess
from pathlib import Path
from itertools import zip_longest


def smart_open(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def check_executable(program, required=True):
    found = shutil.which(program)
    if found is None and required:
        raise RuntimeError(
            f"ERRO: '{program}' não foi encontrado no PATH. "
            f"Instale ou carregue esse programa antes de rodar o pipeline."
        )
    return found


def run_cmd(cmd, log_file=None, dry_run=False):
    print("\n[CMD]")
    print(cmd)

    if dry_run:
        return

    if log_file is not None:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)

        with open(log_file, "w") as log:
            log.write(cmd + "\n\n")
            subprocess.run(
                cmd,
                shell=True,
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,
                executable="/bin/bash",
            )
    else:
        subprocess.run(
            cmd,
            shell=True,
            check=True,
            executable="/bin/bash",
        )


def read_fasta(fasta_path):
    records = {}
    header = None
    seq_chunks = []

    with open(fasta_path, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    seq_id = header.split()[0]
                    records[seq_id] = {
                        "header": header,
                        "seq": "".join(seq_chunks).upper(),
                    }

                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)

        if header is not None:
            seq_id = header.split()[0]
            records[seq_id] = {
                "header": header,
                "seq": "".join(seq_chunks).upper(),
            }

    return records


def write_fasta(records, output_path, width=80):
    with open(output_path, "w") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")


def detect_sequence_format(path):
    with smart_open(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith("@"):
                return "fastq"
            if line.startswith(">"):
                return "fasta"

    raise ValueError(f"Não foi possível detectar o formato de {path}")


def normalize_read_id(header):
    read_id = header.strip()

    if read_id.startswith("@") or read_id.startswith(">"):
        read_id = read_id[1:]

    read_id = read_id.split()[0]

    if read_id.endswith("/1") or read_id.endswith("/2"):
        read_id = read_id[:-2]

    return read_id


def fastq_iter(fastq_path):
    with smart_open(fastq_path, "rt") as handle:
        while True:
            h = handle.readline()
            if not h:
                break

            s = handle.readline()
            p = handle.readline()
            q = handle.readline()

            if not q:
                raise ValueError(f"FASTQ incompleto ou corrompido: {fastq_path}")

            yield h, s, p, q, normalize_read_id(h)


def fasta_iter(fasta_path):
    with smart_open(fasta_path, "rt") as handle:
        header = None
        seq_chunks = []

        for line in handle:
            line = line.rstrip("\n")

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_chunks)
                    yield header, seq, normalize_read_id(header)

                header = line
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())

        if header is not None:
            seq = "".join(seq_chunks)
            yield header, seq, normalize_read_id(header)


def write_fastq_record(handle, record):
    h, s, p, q, _ = record
    handle.write(h)
    handle.write(s)
    handle.write(p)
    handle.write(q)


def write_fasta_record(handle, header, seq, width=80):
    handle.write(header + "\n")
    for i in range(0, len(seq), width):
        handle.write(seq[i:i + width] + "\n")


def load_read_ids(ids_path):
    ids = set()
    with open(ids_path, "r") as handle:
        for line in handle:
            line = line.strip()
            if line:
                ids.add(normalize_read_id(line))
    return ids


def to_int(value, default=0):
    try:
        return int(float(value))
    except Exception:
        return default


def to_float(value, default=0.0):
    try:
        return float(value)
    except Exception:
        return default


def load_scaffold_summary(summary_csv):
    rows = []

    with open(summary_csv, "r", newline="") as handle:
        reader = csv.DictReader(handle)

        required = {
            "sseqid",
            "n_hits",
            "n_unique_queries",
            "scaffold_length_bp",
            "covered_bp_merged_hits",
            "best_evalue",
            "max_bitscore",
        }

        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"O arquivo {summary_csv} não tem as colunas esperadas: {missing}"
            )

        for row in reader:
            row["_n_hits"] = to_int(row.get("n_hits"))
            row["_n_unique_queries"] = to_int(row.get("n_unique_queries"))
            row["_scaffold_length_bp"] = to_int(row.get("scaffold_length_bp"))
            row["_covered_bp_merged_hits"] = to_int(row.get("covered_bp_merged_hits"))
            row["_best_evalue"] = to_float(row.get("best_evalue"), default=999999.0)
            row["_max_bitscore"] = to_float(row.get("max_bitscore"))
            rows.append(row)

    rows.sort(
        key=lambda r: (
            -r["_n_unique_queries"],
            -r["_covered_bp_merged_hits"],
            r["_best_evalue"],
            -r["_max_bitscore"],
            -r["_n_hits"],
        )
    )

    return rows


def select_scaffolds(rows, top_n, min_unique_queries, manual_scaffolds):
    if manual_scaffolds:
        wanted = [x.strip() for x in manual_scaffolds.split(",") if x.strip()]
        by_id = {row["sseqid"]: row for row in rows}

        selected = []
        for scaffold in wanted:
            if scaffold in by_id:
                selected.append(by_id[scaffold])
            else:
                selected.append({
                    "sseqid": scaffold,
                    "n_hits": "NA",
                    "n_unique_queries": "NA",
                    "scaffold_length_bp": "NA",
                    "covered_bp_merged_hits": "NA",
                    "best_evalue": "NA",
                    "max_bitscore": "NA",
                    "queries": "manual_selection",
                    "_n_hits": 0,
                    "_n_unique_queries": 0,
                    "_scaffold_length_bp": 0,
                    "_covered_bp_merged_hits": 0,
                    "_best_evalue": 0,
                    "_max_bitscore": 0,
                })

        return selected

    filtered = [
        row for row in rows
        if row["_n_unique_queries"] >= min_unique_queries
    ]

    if not filtered:
        print(
            "[AVISO] Nenhum scaffold passou em --min-unique-queries. "
            "Usando os melhores scaffolds mesmo assim."
        )
        filtered = rows

    return filtered[:top_n]


def write_selected_scaffolds_csv(selected, output_csv):
    fieldnames = [
        "rank",
        "sseqid",
        "n_hits",
        "n_unique_queries",
        "scaffold_length_bp",
        "covered_bp_merged_hits",
        "best_evalue",
        "max_bitscore",
        "queries",
    ]

    with open(output_csv, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()

        for i, row in enumerate(selected, start=1):
            writer.writerow({
                "rank": i,
                "sseqid": row.get("sseqid"),
                "n_hits": row.get("n_hits"),
                "n_unique_queries": row.get("n_unique_queries"),
                "scaffold_length_bp": row.get("scaffold_length_bp"),
                "covered_bp_merged_hits": row.get("covered_bp_merged_hits"),
                "best_evalue": row.get("best_evalue"),
                "max_bitscore": row.get("max_bitscore"),
                "queries": row.get("queries", ""),
            })


def extract_selected_scaffolds(assembly_fasta, selected, output_fasta):
    assembly = read_fasta(assembly_fasta)
    fasta_records = []

    for row in selected:
        scaffold_id = row["sseqid"]

        if scaffold_id not in assembly:
            print(f"[AVISO] Scaffold não encontrado no assembly: {scaffold_id}")
            continue

        rec = assembly[scaffold_id]

        header = (
            f"{scaffold_id} "
            f"n_unique_queries={row.get('n_unique_queries', 'NA')} "
            f"n_hits={row.get('n_hits', 'NA')} "
            f"covered_bp_merged_hits={row.get('covered_bp_merged_hits', 'NA')} "
            f"best_evalue={row.get('best_evalue', 'NA')} "
            f"max_bitscore={row.get('max_bitscore', 'NA')}"
        )

        fasta_records.append((header, rec["seq"]))

    if not fasta_records:
        raise RuntimeError(
            "Nenhum scaffold foi extraído. Verifique se os IDs do "
            "scaffold_summary.csv batem com os IDs do FASTA da montagem."
        )

    write_fasta(fasta_records, output_fasta)
    return [header.split()[0] for header, seq in fasta_records]


def extract_paired_fastq_by_ids(r1, r2, ids_path, out_r1, out_r2):
    wanted = load_read_ids(ids_path)

    total_pairs = 0
    kept_pairs = 0
    mismatched_ids = 0

    with smart_open(out_r1, "wt") as out1, smart_open(out_r2, "wt") as out2:
        for rec1, rec2 in zip_longest(fastq_iter(r1), fastq_iter(r2)):
            if rec1 is None or rec2 is None:
                raise ValueError(
                    "Os arquivos R1 e R2 têm quantidades diferentes de reads."
                )

            total_pairs += 1

            id1 = rec1[4]
            id2 = rec2[4]

            if id1 != id2:
                mismatched_ids += 1

            if id1 in wanted or id2 in wanted:
                write_fastq_record(out1, rec1)
                write_fastq_record(out2, rec2)
                kept_pairs += 1

    return {
        "total_pairs": total_pairs,
        "kept_pairs": kept_pairs,
        "mismatched_ids": mismatched_ids,
    }


def extract_single_reads_by_ids(reads_path, ids_path, output_prefix):
    wanted = load_read_ids(ids_path)
    seq_format = detect_sequence_format(reads_path)

    total_reads = 0
    kept_reads = 0

    if seq_format == "fastq":
        output_reads = Path(str(output_prefix) + ".fastq.gz")

        with smart_open(output_reads, "wt") as out:
            for rec in fastq_iter(reads_path):
                total_reads += 1

                if rec[4] in wanted:
                    write_fastq_record(out, rec)
                    kept_reads += 1

    elif seq_format == "fasta":
        output_reads = Path(str(output_prefix) + ".fasta.gz")

        with smart_open(output_reads, "wt") as out:
            for header, seq, read_id in fasta_iter(reads_path):
                total_reads += 1

                if read_id in wanted:
                    write_fasta_record(out, header, seq)
                    kept_reads += 1

    else:
        raise ValueError(f"Formato não suportado: {seq_format}")

    return {
        "sequence_format": seq_format,
        "total_reads": total_reads,
        "kept_reads": kept_reads,
        "output_reads": str(output_reads),
    }


def write_dict_csv(row, output_csv):
    with open(output_csv, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)


def map_short_reads_bowtie2(ref_fasta, r1, r2, outdir, threads, dry_run=False):
    check_executable("bowtie2-build")
    check_executable("bowtie2")
    check_executable("samtools")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    index_prefix = outdir / "candidate_scaffolds.bt2"
    bam = outdir / "short_reads_vs_candidates.sorted.bam"
    flagstat = outdir / "short_reads_vs_candidates.flagstat.txt"
    log = outdir / "bowtie2_mapping.log"

    run_cmd(
        f"bowtie2-build {ref_fasta} {index_prefix}",
        log_file=outdir / "bowtie2_build.log",
        dry_run=dry_run,
    )

    cmd = (
        f"bowtie2 --very-sensitive-local "
        f"-x {index_prefix} "
        f"-1 {r1} "
        f"-2 {r2} "
        f"-p {threads} "
        f"2> {log} "
        f"| samtools sort -@ {threads} -o {bam} -"
    )

    run_cmd(cmd, dry_run=dry_run)
    run_cmd(f"samtools index {bam}", dry_run=dry_run)
    run_cmd(f"samtools flagstat {bam} > {flagstat}", dry_run=dry_run)

    return bam, flagstat


def map_short_reads_bwa(ref_fasta, r1, r2, outdir, threads, dry_run=False):
    check_executable("bwa")
    check_executable("samtools")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    bam = outdir / "short_reads_vs_candidates.sorted.bam"
    flagstat = outdir / "short_reads_vs_candidates.flagstat.txt"
    log = outdir / "bwa_mapping.log"

    run_cmd(f"bwa index {ref_fasta}", log_file=outdir / "bwa_index.log", dry_run=dry_run)

    cmd = (
        f"bwa mem -t {threads} {ref_fasta} {r1} {r2} "
        f"2> {log} "
        f"| samtools sort -@ {threads} -o {bam} -"
    )

    run_cmd(cmd, dry_run=dry_run)
    run_cmd(f"samtools index {bam}", dry_run=dry_run)
    run_cmd(f"samtools flagstat {bam} > {flagstat}", dry_run=dry_run)

    return bam, flagstat


def map_pacbio_reads_minimap2(ref_fasta, pacbio_reads, outdir, threads, dry_run=False):
    check_executable("minimap2")
    check_executable("samtools")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    bam = outdir / "pacbio_reads_vs_candidates.sorted.bam"
    flagstat = outdir / "pacbio_reads_vs_candidates.flagstat.txt"
    log = outdir / "minimap2_pacbio_mapping.log"

    cmd = (
        f"minimap2 -t {threads} -ax map-pb {ref_fasta} {pacbio_reads} "
        f"2> {log} "
        f"| samtools sort -@ {threads} -o {bam} -"
    )

    run_cmd(cmd, dry_run=dry_run)
    run_cmd(f"samtools index {bam}", dry_run=dry_run)
    run_cmd(f"samtools flagstat {bam} > {flagstat}", dry_run=dry_run)

    return bam, flagstat


def extract_mapped_read_ids_from_bam(bam, output_ids, dry_run=False):
    check_executable("samtools")

    cmd = (
        f"samtools view -F 4 {bam} "
        f"| cut -f 1 "
        f"| LC_ALL=C sort -u "
        f"> {output_ids}"
    )

    run_cmd(cmd, dry_run=dry_run)


def calculate_coverage_by_scaffold(bam, ref_fasta, output_depth, output_csv, dry_run=False):
    check_executable("samtools")

    if not dry_run:
        ref_records = read_fasta(ref_fasta)
    else:
        ref_records = {}

    run_cmd(f"samtools depth -a {bam} > {output_depth}", dry_run=dry_run)

    if dry_run:
        return

    stats = {}

    for scaffold_id, rec in ref_records.items():
        length = len(rec["seq"])
        stats[scaffold_id] = {
            "scaffold": scaffold_id,
            "length_bp": length,
            "sum_depth": 0,
            "covered_bases_depth_ge1": 0,
            "max_depth": 0,
        }

    with open(output_depth, "r") as handle:
        for line in handle:
            scaffold, pos, depth = line.rstrip("\n").split("\t")
            depth = int(depth)

            if scaffold not in stats:
                continue

            stats[scaffold]["sum_depth"] += depth

            if depth >= 1:
                stats[scaffold]["covered_bases_depth_ge1"] += 1

            if depth > stats[scaffold]["max_depth"]:
                stats[scaffold]["max_depth"] = depth

    fieldnames = [
        "scaffold",
        "length_bp",
        "mean_depth",
        "covered_bases_depth_ge1",
        "breadth_coverage_percent",
        "max_depth",
    ]

    with open(output_csv, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()

        for scaffold_id, row in stats.items():
            length = row["length_bp"]
            mean_depth = row["sum_depth"] / length if length else 0
            breadth = (
                row["covered_bases_depth_ge1"] / length * 100
                if length else 0
            )

            writer.writerow({
                "scaffold": scaffold_id,
                "length_bp": length,
                "mean_depth": round(mean_depth, 4),
                "covered_bases_depth_ge1": row["covered_bases_depth_ge1"],
                "breadth_coverage_percent": round(breadth, 4),
                "max_depth": row["max_depth"],
            })


def run_spades_pacbio(recruited_r1, recruited_r2, recruited_pacbio,
                      outdir, threads, memory, dry_run=False):
    check_executable("spades.py")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = (
        f"spades.py "
        f"-1 {recruited_r1} "
        f"-2 {recruited_r2} "
        f"-o {outdir} "
        f"-t {threads} "
        f"-m {memory} "
        f"--careful"
    )

    if recruited_pacbio is not None:
        cmd += f" --pacbio {recruited_pacbio}"

    run_cmd(cmd, log_file=outdir / "spades_pipeline.log", dry_run=dry_run)


def run_cap3(selected_scaffolds_fasta, cap3_outdir, dry_run=False):
    cap3 = check_executable("cap3", required=False)

    if cap3 is None:
        print(
            "\n[AVISO] CAP3 não encontrado no PATH. "
            "Pulando etapa opcional de CAP3."
        )
        return None

    cap3_outdir = Path(cap3_outdir)
    cap3_outdir.mkdir(parents=True, exist_ok=True)

    input_copy = cap3_outdir / "selected_scaffolds_for_cap3.fasta"
    shutil.copy2(selected_scaffolds_fasta, input_copy)

    log = cap3_outdir / "cap3.log"

    cmd = f"cd {cap3_outdir} && cap3 {input_copy.name} > {log.name}"
    run_cmd(cmd, dry_run=dry_run)

    return log


def print_final_report(outdir):
    outdir = Path(outdir)

    print("\n" + "=" * 80)
    print("PIPELINE FINALIZADO")
    print("=" * 80)

    files = [
        ("Scaffolds selecionados", outdir / "selected_scaffolds.csv"),
        ("FASTA dos scaffolds candidatos", outdir / "selected_candidate_scaffolds.fasta"),
        ("BAM short reads", outdir / "mapping_short" / "short_reads_vs_candidates.sorted.bam"),
        ("Flagstat short reads", outdir / "mapping_short" / "short_reads_vs_candidates.flagstat.txt"),
        ("IDs das short reads recrutadas", outdir / "recruited_reads" / "mapped_short_read_ids.txt"),
        ("Short reads R1 recrutadas", outdir / "recruited_reads" / "recruited_R1.fastq.gz"),
        ("Short reads R2 recrutadas", outdir / "recruited_reads" / "recruited_R2.fastq.gz"),
        ("BAM PacBio", outdir / "mapping_pacbio" / "pacbio_reads_vs_candidates.sorted.bam"),
        ("Flagstat PacBio", outdir / "mapping_pacbio" / "pacbio_reads_vs_candidates.flagstat.txt"),
        ("IDs das PacBio reads recrutadas", outdir / "recruited_reads" / "mapped_pacbio_read_ids.txt"),
        ("PacBio reads recrutadas FASTQ", outdir / "recruited_reads" / "recruited_pacbio_reads.fastq.gz"),
        ("PacBio reads recrutadas FASTA", outdir / "recruited_reads" / "recruited_pacbio_reads.fasta.gz"),
        ("Cobertura por scaffold", outdir / "coverage" / "coverage_by_scaffold.csv"),
        ("Depth por base", outdir / "coverage" / "depth_per_base.tsv"),
        ("Montagem SPAdes scaffolds", outdir / "spades_reassembly" / "scaffolds.fasta"),
        ("Montagem SPAdes contigs", outdir / "spades_reassembly" / "contigs.fasta"),
        ("Log CAP3", outdir / "cap3" / "cap3.log"),
    ]

    for desc, path in files:
        if path.exists:
            try:
                if path.exists():
                    print(f"{desc}: {path}")
            except Exception:
                pass

    print("=" * 80)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Pipeline para recrutar reads mitocondriais a partir dos scaffolds "
            "com hits de tBLASTn e remontar usando short reads + PacBio."
        )
    )

    parser.add_argument(
        "--assembly",
        required=True,
        help="FASTA da montagem original. Exemplo: scaffolds_1000bp.fasta",
    )

    parser.add_argument(
        "--scaffold-summary",
        required=True,
        help="CSV gerado antes: tblastn_mitogenome_results/tables/scaffold_summary.csv",
    )

    parser.add_argument(
        "--short-r1",
        required=True,
        help="FASTQ R1 das short reads originais.",
    )

    parser.add_argument(
        "--short-r2",
        required=True,
        help="FASTQ R2 das short reads originais.",
    )

    parser.add_argument(
        "--pacbio-reads",
        default=None,
        help="Arquivo FASTQ ou FASTA das reads PacBio originais. Pode ser .gz.",
    )

    parser.add_argument(
        "--outdir",
        default="mitogenome_recruitment_pipeline",
        help="Diretório de saída.",
    )

    parser.add_argument(
        "--top-n",
        type=int,
        default=4,
        help="Número de scaffolds candidatos a selecionar. Padrão: 4.",
    )

    parser.add_argument(
        "--min-unique-queries",
        type=int,
        default=1,
        help="Mínimo de genes/CDS mitocondriais diferentes batendo no scaffold.",
    )

    parser.add_argument(
        "--manual-scaffolds",
        default=None,
        help=(
            "Lista manual de scaffolds separados por vírgula. "
            "Exemplo: --manual-scaffolds scaffold_1,scaffold_8"
        ),
    )

    parser.add_argument(
        "--short-aligner",
        choices=["bowtie2", "bwa"],
        default="bowtie2",
        help="Alinhador para short reads.",
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=12,
        help="Número de threads.",
    )

    parser.add_argument(
        "--memory",
        type=int,
        default=64,
        help="Memória máxima para SPAdes, em GB.",
    )

    parser.add_argument(
        "--skip-assembly",
        action="store_true",
        help="Executa recrutamento das reads, mas não roda SPAdes.",
    )

    parser.add_argument(
        "--run-cap3",
        action="store_true",
        help="Tenta rodar CAP3 nos scaffolds candidatos selecionados.",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Mostra comandos, mas não executa.",
    )

    args = parser.parse_args()

    assembly = Path(args.assembly).resolve()
    summary = Path(args.scaffold_summary).resolve()
    short_r1 = Path(args.short_r1).resolve()
    short_r2 = Path(args.short_r2).resolve()
    pacbio_reads = Path(args.pacbio_reads).resolve() if args.pacbio_reads else None

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    mapping_short_dir = outdir / "mapping_short"
    mapping_pacbio_dir = outdir / "mapping_pacbio"
    recruited_dir = outdir / "recruited_reads"
    coverage_dir = outdir / "coverage"
    spades_dir = outdir / "spades_reassembly"
    cap3_dir = outdir / "cap3"

    recruited_dir.mkdir(parents=True, exist_ok=True)
    coverage_dir.mkdir(parents=True, exist_ok=True)

    print("\n[1/8] Lendo scaffold_summary.csv e selecionando candidatos...")
    rows = load_scaffold_summary(summary)

    selected = select_scaffolds(
        rows=rows,
        top_n=args.top_n,
        min_unique_queries=args.min_unique_queries,
        manual_scaffolds=args.manual_scaffolds,
    )

    selected_csv = outdir / "selected_scaffolds.csv"
    selected_fasta = outdir / "selected_candidate_scaffolds.fasta"

    write_selected_scaffolds_csv(selected, selected_csv)

    print("\nScaffolds selecionados:")
    for i, row in enumerate(selected, start=1):
        print(
            f"{i}. {row['sseqid']} | "
            f"n_unique_queries={row.get('n_unique_queries')} | "
            f"n_hits={row.get('n_hits')} | "
            f"covered_bp={row.get('covered_bp_merged_hits')}"
        )

    print("\n[2/8] Extraindo FASTA dos scaffolds candidatos...")
    extracted_ids = extract_selected_scaffolds(
        assembly_fasta=assembly,
        selected=selected,
        output_fasta=selected_fasta,
    )

    print(f"Scaffolds extraídos: {len(extracted_ids)}")
    print(f"FASTA: {selected_fasta}")

    print("\n[3/8] Mapeando short reads contra scaffolds candidatos...")

    if args.short_aligner == "bowtie2":
        short_bam, short_flagstat = map_short_reads_bowtie2(
            ref_fasta=selected_fasta,
            r1=short_r1,
            r2=short_r2,
            outdir=mapping_short_dir,
            threads=args.threads,
            dry_run=args.dry_run,
        )
    else:
        short_bam, short_flagstat = map_short_reads_bwa(
            ref_fasta=selected_fasta,
            r1=short_r1,
            r2=short_r2,
            outdir=mapping_short_dir,
            threads=args.threads,
            dry_run=args.dry_run,
        )

    print("\n[4/8] Extraindo IDs das short reads mapeadas...")

    short_ids = recruited_dir / "mapped_short_read_ids.txt"

    extract_mapped_read_ids_from_bam(
        bam=short_bam,
        output_ids=short_ids,
        dry_run=args.dry_run,
    )

    recruited_r1 = recruited_dir / "recruited_R1.fastq.gz"
    recruited_r2 = recruited_dir / "recruited_R2.fastq.gz"

    if not args.dry_run:
        print("\n[5/8] Extraindo pares FASTQ recrutados...")

        paired_stats = extract_paired_fastq_by_ids(
            r1=short_r1,
            r2=short_r2,
            ids_path=short_ids,
            out_r1=recruited_r1,
            out_r2=recruited_r2,
        )

        paired_stats_csv = recruited_dir / "recruited_short_reads_summary.csv"
        write_dict_csv(paired_stats, paired_stats_csv)

        print(f"Pares totais analisados: {paired_stats['total_pairs']}")
        print(f"Pares recrutados:        {paired_stats['kept_pairs']}")
        print(f"IDs R1/R2 divergentes:   {paired_stats['mismatched_ids']}")
    else:
        print("\n[5/8] Dry-run: pulando extração real dos FASTQ.")

    recruited_pacbio = None

    if pacbio_reads is not None:
        print("\n[6/8] Mapeando PacBio reads contra scaffolds candidatos...")

        pacbio_bam, pacbio_flagstat = map_pacbio_reads_minimap2(
            ref_fasta=selected_fasta,
            pacbio_reads=pacbio_reads,
            outdir=mapping_pacbio_dir,
            threads=args.threads,
            dry_run=args.dry_run,
        )

        pacbio_ids = recruited_dir / "mapped_pacbio_read_ids.txt"

        extract_mapped_read_ids_from_bam(
            bam=pacbio_bam,
            output_ids=pacbio_ids,
            dry_run=args.dry_run,
        )

        if not args.dry_run:
            print("\nExtraindo PacBio reads recrutadas...")

            pacbio_prefix = recruited_dir / "recruited_pacbio_reads"

            pacbio_stats = extract_single_reads_by_ids(
                reads_path=pacbio_reads,
                ids_path=pacbio_ids,
                output_prefix=pacbio_prefix,
            )

            pacbio_stats_csv = recruited_dir / "recruited_pacbio_reads_summary.csv"
            write_dict_csv(pacbio_stats, pacbio_stats_csv)

            recruited_pacbio = Path(pacbio_stats["output_reads"])

            print(f"Formato PacBio detectado:       {pacbio_stats['sequence_format']}")
            print(f"PacBio reads totais analisadas: {pacbio_stats['total_reads']}")
            print(f"PacBio reads recrutadas:        {pacbio_stats['kept_reads']}")
            print(f"Arquivo PacBio recrutado:       {recruited_pacbio}")
    else:
        print("\n[6/8] Nenhuma PacBio read fornecida. Pulando etapa de PacBio.")

    print("\n[7/8] Calculando cobertura das short reads nos scaffolds candidatos...")

    calculate_coverage_by_scaffold(
        bam=short_bam,
        ref_fasta=selected_fasta,
        output_depth=coverage_dir / "depth_per_base.tsv",
        output_csv=coverage_dir / "coverage_by_scaffold.csv",
        dry_run=args.dry_run,
    )

    if args.run_cap3:
        print("\n[CAP3] Rodando etapa opcional de CAP3...")

        run_cap3(
            selected_scaffolds_fasta=selected_fasta,
            cap3_outdir=cap3_dir,
            dry_run=args.dry_run,
        )

    if args.skip_assembly:
        print("\n[8/8] --skip-assembly usado. Pulando SPAdes.")
    else:
        print("\n[8/8] Rodando SPAdes com short reads recrutadas + PacBio recrutadas...")

        run_spades_pacbio(
            recruited_r1=recruited_r1,
            recruited_r2=recruited_r2,
            recruited_pacbio=recruited_pacbio,
            outdir=spades_dir,
            threads=args.threads,
            memory=args.memory,
            dry_run=args.dry_run,
        )

    print_final_report(outdir)


if __name__ == "__main__":
    main()

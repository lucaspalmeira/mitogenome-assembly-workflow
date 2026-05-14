#!/usr/bin/env python3
"""
mito_pipeline.py

Pipeline transparente para montagem e avaliação de genoma mitocondrial a partir de
Illumina paired-end ou PacBio, usando CDS mitocondriais de referência como isca
para identificar contigs mitocondriais após a montagem.

Autor: gerado para Lucas Palmeira

Principais saídas:
  - montagem de novo
  - BLASTN/TBLASTN dos CDS mitocondriais contra contigs montados
  - FASTA com contigs candidatos mitocondriais
  - BAM/coverage/depth do mapeamento das reads contra os candidatos
  - relatório Markdown com comandos, estatísticas e interpretação

Dependências externas esperadas no PATH:
  - Illumina: spades.py
  - PacBio: flye
  - Ambos: makeblastdb, blastn, tblastn, minimap2, samtools
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import shutil
import statistics
import subprocess
import sys
import textwrap
import time
from collections import defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:
    from .common import parse_flye_circular_flags, terminal_overlap, write_fasta
except ImportError:
    from common import parse_flye_circular_flags, terminal_overlap, write_fasta


# -------------------------
# Utilitários gerais
# -------------------------

@dataclass
class FastqStats:
    file: str
    reads_counted: int
    total_bases: int
    min_len: int
    max_len: int
    mean_len: float
    median_len: float
    n50: int
    exact: bool


@dataclass
class FastaStats:
    file: str
    sequences: int
    total_bases: int
    min_len: int
    max_len: int
    mean_len: float
    median_len: float
    n50: int


@dataclass
class CandidateSummary:
    contig: str
    contig_len: int
    blastn_hits: int
    tblastn_hits: int
    unique_nt_queries: int
    unique_prot_queries: int
    blastn_bitscore: float
    tblastn_bitscore: float
    total_bitscore: float
    best_nt_pident: float
    best_prot_pident: float
    circular_overlap_len: int
    circular_overlap_identity: float
    circular_by_terminal_overlap: bool
    circular_by_flye: Optional[bool]


@dataclass
class CoverageSummary:
    contig: str
    length_from_depth: int
    mean_depth: float
    median_depth: float
    min_depth: int
    max_depth: int
    breadth_1x: float
    breadth_5x: float
    breadth_10x: float
    breadth_20x: float


class PipelineError(RuntimeError):
    pass


def now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def open_text_auto(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def which_or_fail(program: str) -> str:
    hit = shutil.which(program)
    if not hit:
        raise PipelineError(
            f"Dependência não encontrada no PATH: {program}. "
            "Instale via conda/mamba ou carregue o módulo correspondente."
        )
    return hit


def check_dependencies(mode: str) -> None:
    common = ["makeblastdb", "blastn", "tblastn", "minimap2", "samtools"]
    mode_specific = []
    if mode == "illumina":
        mode_specific = ["spades.py"]
    elif mode == "pacbio":
        mode_specific = ["flye"]
    for program in common + mode_specific:
        which_or_fail(program)


def run_cmd(
    cmd: Sequence[str],
    log_file: Path,
    cwd: Optional[Path] = None,
    stdout_file: Optional[Path] = None,
    stderr_file: Optional[Path] = None,
    dry_run: bool = False,
) -> None:
    """Executa comando e registra a chamada no log."""
    cmd_str = " ".join(str(x) for x in cmd)
    with log_file.open("a", encoding="utf-8") as log:
        log.write(f"\n[{now()}] $ {cmd_str}\n")
        if cwd:
            log.write(f"# cwd: {cwd}\n")

    if dry_run:
        return

    stdout_handle = None
    stderr_handle = None
    try:
        if stdout_file:
            stdout_file.parent.mkdir(parents=True, exist_ok=True)
            stdout_handle = stdout_file.open("w", encoding="utf-8")
        if stderr_file:
            stderr_file.parent.mkdir(parents=True, exist_ok=True)
            stderr_handle = stderr_file.open("w", encoding="utf-8")

        result = subprocess.run(
            list(map(str, cmd)),
            cwd=str(cwd) if cwd else None,
            stdout=stdout_handle if stdout_handle else subprocess.PIPE,
            stderr=stderr_handle if stderr_handle else subprocess.PIPE,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            extra = ""
            if not stderr_file and result.stderr:
                extra = f"\nSTDERR:\n{result.stderr[-4000:]}"
            raise PipelineError(
                f"Comando falhou com código {result.returncode}: {cmd_str}{extra}"
            )
        if not stdout_file and result.stdout:
            with log_file.open("a", encoding="utf-8") as log:
                log.write(result.stdout[-4000:] + "\n")
        if not stderr_file and result.stderr:
            with log_file.open("a", encoding="utf-8") as log:
                log.write(result.stderr[-4000:] + "\n")
    finally:
        if stdout_handle:
            stdout_handle.close()
        if stderr_handle:
            stderr_handle.close()


def n50(lengths: Sequence[int]) -> int:
    if not lengths:
        return 0
    total = sum(lengths)
    half = total / 2
    acc = 0
    for length in sorted(lengths, reverse=True):
        acc += length
        if acc >= half:
            return int(length)
    return 0


# -------------------------
# Leitura de FASTA/FASTQ
# -------------------------


def iter_fasta(path: Path) -> Iterable[Tuple[str, str]]:
    header = None
    chunks: List[str] = []
    with open_text_auto(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks)


def read_fasta_dict(path: Path) -> Dict[str, str]:
    return {h: s for h, s in iter_fasta(path)}


def fasta_stats(path: Path) -> FastaStats:
    lengths = [len(seq) for _, seq in iter_fasta(path)]
    if not lengths:
        return FastaStats(str(path), 0, 0, 0, 0, 0.0, 0.0, 0)
    return FastaStats(
        file=str(path),
        sequences=len(lengths),
        total_bases=sum(lengths),
        min_len=min(lengths),
        max_len=max(lengths),
        mean_len=float(statistics.mean(lengths)),
        median_len=float(statistics.median(lengths)),
        n50=n50(lengths),
    )


def fastq_lengths(path: Path, max_reads: int = 0) -> Tuple[List[int], bool]:
    """Retorna comprimentos. max_reads=0 significa arquivo completo."""
    lengths: List[int] = []
    with open_text_auto(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                break
            lengths.append(len(seq.strip()))
            if max_reads and len(lengths) >= max_reads:
                return lengths, False
    return lengths, True


def fastq_stats(path: Path, max_reads: int = 0) -> FastqStats:
    lengths, exact = fastq_lengths(path, max_reads=max_reads)
    if not lengths:
        return FastqStats(str(path), 0, 0, 0, 0, 0.0, 0.0, 0, exact)
    return FastqStats(
        file=str(path),
        reads_counted=len(lengths),
        total_bases=sum(lengths),
        min_len=min(lengths),
        max_len=max(lengths),
        mean_len=float(statistics.mean(lengths)),
        median_len=float(statistics.median(lengths)),
        n50=n50(lengths),
        exact=exact,
    )


# -------------------------
# Montagem
# -------------------------


def run_illumina_assembly(args, log_file: Path) -> Path:
    spades_dir = args.outdir / "01_assembly_spades"
    ensure_dir(spades_dir)
    cmd = [
        "spades.py",
        "--careful",
        "-1",
        str(args.r1),
        "-2",
        str(args.r2),
        "-o",
        str(spades_dir),
        "--threads",
        str(args.threads),
    ]
    if args.memory_gb:
        cmd.extend(["--memory", str(args.memory_gb)])
    if args.spades_k:
        cmd.extend(["-k", args.spades_k])
    run_cmd(cmd, log_file, dry_run=args.dry_run)

    scaffolds = spades_dir / "scaffolds.fasta"
    contigs = spades_dir / "contigs.fasta"
    if args.dry_run:
        return scaffolds
    if args.spades_output == "contigs":
        assembly = contigs
    else:
        assembly = scaffolds if scaffolds.exists() else contigs
    if not assembly.exists():
        raise PipelineError(
            f"SPAdes terminou, mas não encontrei {scaffolds} nem {contigs}."
        )
    return assembly


def run_pacbio_assembly(args, log_file: Path) -> Path:
    flye_dir = args.outdir / "01_assembly_flye"
    ensure_dir(flye_dir)
    pacbio_flag = {
        "raw": "--pacbio-raw",
        "hifi": "--pacbio-hifi",
        "corrected": "--pacbio-corr",
    }[args.pacbio_type]
    cmd = [
        "flye",
        pacbio_flag,
        str(args.pacbio),
        "--out-dir",
        str(flye_dir),
        "--threads",
        str(args.threads),
        "--genome-size",
        args.flye_genome_size,
    ]
    if args.flye_meta:
        cmd.append("--meta")
    run_cmd(cmd, log_file, dry_run=args.dry_run)

    assembly = flye_dir / "assembly.fasta"
    if args.dry_run:
        return assembly
    if not assembly.exists():
        raise PipelineError(f"Flye terminou, mas não encontrei {assembly}.")
    return assembly


# -------------------------
# BLAST e seleção de candidatos
# -------------------------


def make_blast_db(assembly: Path, db_prefix: Path, log_file: Path, dry_run: bool) -> None:
    ensure_dir(db_prefix.parent)
    run_cmd(
        ["makeblastdb", "-in", str(assembly), "-dbtype", "nucl", "-out", str(db_prefix)],
        log_file,
        dry_run=dry_run,
    )


def run_blasts(args, assembly: Path, log_file: Path) -> Tuple[Path, Path]:
    blast_dir = args.outdir / "02_blast_mito_hits"
    ensure_dir(blast_dir)
    db_prefix = blast_dir / "assembly_db"
    make_blast_db(assembly, db_prefix, log_file, args.dry_run)

    fields = (
        "qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"
    )
    blastn_out = blast_dir / "blastn_nt_vs_assembly.tsv"
    tblastn_out = blast_dir / "tblastn_prot_vs_assembly.tsv"

    run_cmd(
        [
            "blastn",
            "-query",
            str(args.nt_ref),
            "-db",
            str(db_prefix),
            "-outfmt",
            f"6 {fields}",
            "-evalue",
            str(args.evalue),
            "-max_target_seqs",
            str(args.max_target_seqs),
            "-num_threads",
            str(args.threads),
            "-out",
            str(blastn_out),
        ],
        log_file,
        dry_run=args.dry_run,
    )

    run_cmd(
        [
            "tblastn",
            "-query",
            str(args.prot_ref),
            "-db",
            str(db_prefix),
            "-outfmt",
            f"6 {fields}",
            "-evalue",
            str(args.evalue),
            "-max_target_seqs",
            str(args.max_target_seqs),
            "-num_threads",
            str(args.threads),
            "-out",
            str(tblastn_out),
        ],
        log_file,
        dry_run=args.dry_run,
    )
    return blastn_out, tblastn_out


def parse_blast_table(path: Path, source: str) -> List[dict]:
    hits: List[dict] = []
    if not path.exists():
        return hits
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 12:
                continue
            qseqid, sseqid = parts[0], parts[1]
            pident = float(parts[2])
            length = int(float(parts[3]))
            qlen = int(float(parts[4]))
            slen = int(float(parts[5]))
            evalue = float(parts[10])
            bitscore = float(parts[11])
            qcov = length / qlen if qlen else 0.0
            hits.append(
                {
                    "source": source,
                    "qseqid": qseqid,
                    "sseqid": sseqid,
                    "pident": pident,
                    "length": length,
                    "qlen": qlen,
                    "slen": slen,
                    "qcov": qcov,
                    "evalue": evalue,
                    "bitscore": bitscore,
                    "raw": parts,
                }
            )
    return hits


def filter_hits(args, hits: List[dict]) -> List[dict]:
    filtered = []
    for hit in hits:
        if hit["evalue"] > args.evalue:
            continue
        if hit["bitscore"] < args.min_bitscore:
            continue
        if hit["qcov"] < args.min_hit_qcov:
            continue
        if hit["source"] == "blastn" and hit["pident"] < args.min_nt_pident:
            continue
        if hit["source"] == "tblastn" and hit["pident"] < args.min_prot_pident:
            continue
        filtered.append(hit)
    return filtered


def select_candidates(args, assembly: Path, blastn_out: Path, tblastn_out: Path) -> Tuple[Path, Path, List[CandidateSummary]]:
    all_hits = parse_blast_table(blastn_out, "blastn") + parse_blast_table(tblastn_out, "tblastn")
    passing = filter_hits(args, all_hits)

    seqs = read_fasta_dict(assembly)
    per_contig = defaultdict(list)
    for hit in passing:
        per_contig[hit["sseqid"]].append(hit)

    # Caso os nomes da assembly tenham sido truncados por BLAST ou FASTA parser,
    # usa o primeiro token e cria um índice auxiliar.
    seq_by_first_token = {k.split()[0]: v for k, v in seqs.items()}

    flye_flags: Dict[str, Optional[bool]] = {}
    if args.mode == "pacbio":
        flye_flags = parse_flye_circular_flags(args.outdir / "01_assembly_flye")

    summaries: List[CandidateSummary] = []
    candidate_records: Dict[str, str] = {}

    for contig, hits in per_contig.items():
        seq = seqs.get(contig) or seq_by_first_token.get(contig)
        if not seq:
            continue
        if len(seq) < args.min_candidate_len:
            continue
        nt_hits = [h for h in hits if h["source"] == "blastn"]
        prot_hits = [h for h in hits if h["source"] == "tblastn"]
        nt_queries = {h["qseqid"] for h in nt_hits}
        prot_queries = {h["qseqid"] for h in prot_hits}
        nt_bits = sum(h["bitscore"] for h in nt_hits)
        prot_bits = sum(h["bitscore"] for h in prot_hits)
        best_nt_pid = max((h["pident"] for h in nt_hits), default=0.0)
        best_prot_pid = max((h["pident"] for h in prot_hits), default=0.0)
        ov_len, ov_id, ov_bool = terminal_overlap(
            seq, args.min_circular_overlap, args.max_circular_overlap, args.min_circular_identity
        )
        circ_flye = flye_flags.get(contig, flye_flags.get(contig.split()[0]))
        summaries.append(
            CandidateSummary(
                contig=contig,
                contig_len=len(seq),
                blastn_hits=len(nt_hits),
                tblastn_hits=len(prot_hits),
                unique_nt_queries=len(nt_queries),
                unique_prot_queries=len(prot_queries),
                blastn_bitscore=nt_bits,
                tblastn_bitscore=prot_bits,
                total_bitscore=nt_bits + prot_bits,
                best_nt_pident=best_nt_pid,
                best_prot_pident=best_prot_pid,
                circular_overlap_len=ov_len,
                circular_overlap_identity=ov_id,
                circular_by_terminal_overlap=ov_bool,
                circular_by_flye=circ_flye,
            )
        )
        candidate_records[contig] = seq

    summaries.sort(key=lambda x: (x.total_bitscore, x.contig_len), reverse=True)
    ordered_records = {s.contig: candidate_records[s.contig] for s in summaries}

    candidates_fasta = args.outdir / "03_mito_candidates" / "mitochondrial_candidates.fasta"
    candidate_summary_tsv = args.outdir / "03_mito_candidates" / "mitochondrial_candidates_summary.tsv"
    ensure_dir(candidates_fasta.parent)
    write_fasta(ordered_records, candidates_fasta)
    write_candidate_summary_tsv(summaries, candidate_summary_tsv)
    return candidates_fasta, candidate_summary_tsv, summaries


def write_candidate_summary_tsv(summaries: List[CandidateSummary], path: Path) -> None:
    ensure_dir(path.parent)
    fields = list(asdict(summaries[0]).keys()) if summaries else [
        "contig",
        "contig_len",
        "blastn_hits",
        "tblastn_hits",
        "unique_nt_queries",
        "unique_prot_queries",
        "blastn_bitscore",
        "tblastn_bitscore",
        "total_bitscore",
        "best_nt_pident",
        "best_prot_pident",
        "circular_overlap_len",
        "circular_overlap_identity",
        "circular_by_terminal_overlap",
        "circular_by_flye",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for s in summaries:
            writer.writerow(asdict(s))


# -------------------------
# Cobertura e mapeamento
# -------------------------


def map_reads_and_coverage(args, candidates_fasta: Path, log_file: Path) -> Tuple[Optional[Path], Optional[Path], Optional[Path], List[CoverageSummary]]:
    cov_dir = args.outdir / "04_mapping_coverage"
    ensure_dir(cov_dir)

    if args.dry_run:
        bam = cov_dir / "reads_vs_mito_candidates.sorted.bam"
        depth_tsv = cov_dir / "depth_per_base.tsv"
        coverage_tsv = cov_dir / "samtools_coverage.tsv"
        return bam, depth_tsv, coverage_tsv, []

    if not candidates_fasta.exists() or candidates_fasta.stat().st_size == 0:
        return None, None, None, []

    # Indexa FASTA para samtools depth -a funcionar de forma previsível.
    run_cmd(["samtools", "faidx", str(candidates_fasta)], log_file, dry_run=args.dry_run)

    bam = cov_dir / "reads_vs_mito_candidates.sorted.bam"
    minimap_preset = "sr" if args.mode == "illumina" else ("map-hifi" if args.pacbio_type == "hifi" else "map-pb")

    if args.mode == "illumina":
        map_cmd = [
            "minimap2",
            "-ax",
            minimap_preset,
            "-t",
            str(args.threads),
            str(candidates_fasta),
            str(args.r1),
            str(args.r2),
        ]
    else:
        map_cmd = [
            "minimap2",
            "-ax",
            minimap_preset,
            "-t",
            str(args.threads),
            str(candidates_fasta),
            str(args.pacbio),
        ]

    # Usa shell apenas aqui para stream minimap2 -> samtools sort sem salvar SAM gigante.
    shell_cmd = " ".join(map(str, map_cmd)) + f" | samtools sort -@ {args.threads} -o {bam} -"
    with log_file.open("a", encoding="utf-8") as log:
        log.write(f"\n[{now()}] $ {shell_cmd}\n")
    result = subprocess.run(shell_cmd, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise PipelineError(f"Mapeamento/sort falhou:\n{result.stderr[-4000:]}")
    with (cov_dir / "minimap2_samtools_sort.stderr.log").open("w", encoding="utf-8") as handle:
        handle.write(result.stderr)

    run_cmd(["samtools", "index", str(bam)], log_file, dry_run=args.dry_run)

    flagstat = cov_dir / "samtools_flagstat.txt"
    run_cmd(["samtools", "flagstat", str(bam)], log_file, stdout_file=flagstat, dry_run=args.dry_run)

    stats = cov_dir / "samtools_stats.txt"
    run_cmd(["samtools", "stats", str(bam)], log_file, stdout_file=stats, dry_run=args.dry_run)

    coverage_tsv = cov_dir / "samtools_coverage.tsv"
    run_cmd(["samtools", "coverage", str(bam)], log_file, stdout_file=coverage_tsv, dry_run=args.dry_run)

    depth_tsv = cov_dir / "depth_per_base.tsv"
    run_cmd(["samtools", "depth", "-a", str(bam)], log_file, stdout_file=depth_tsv, dry_run=args.dry_run)

    summaries = parse_depth(depth_tsv)
    cov_summary_tsv = cov_dir / "coverage_summary.tsv"
    write_coverage_summary_tsv(summaries, cov_summary_tsv)
    return bam, depth_tsv, coverage_tsv, summaries


def parse_depth(depth_tsv: Path) -> List[CoverageSummary]:
    depths: Dict[str, List[int]] = defaultdict(list)
    if not depth_tsv.exists():
        return []
    with depth_tsv.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            contig = parts[0]
            depth = int(float(parts[2]))
            depths[contig].append(depth)
    summaries: List[CoverageSummary] = []
    for contig, values in depths.items():
        length = len(values)
        if length == 0:
            continue
        summaries.append(
            CoverageSummary(
                contig=contig,
                length_from_depth=length,
                mean_depth=float(statistics.mean(values)),
                median_depth=float(statistics.median(values)),
                min_depth=min(values),
                max_depth=max(values),
                breadth_1x=100.0 * sum(v >= 1 for v in values) / length,
                breadth_5x=100.0 * sum(v >= 5 for v in values) / length,
                breadth_10x=100.0 * sum(v >= 10 for v in values) / length,
                breadth_20x=100.0 * sum(v >= 20 for v in values) / length,
            )
        )
    summaries.sort(key=lambda x: x.mean_depth, reverse=True)
    return summaries


def write_coverage_summary_tsv(summaries: List[CoverageSummary], path: Path) -> None:
    ensure_dir(path.parent)
    fields = list(asdict(summaries[0]).keys()) if summaries else [
        "contig",
        "length_from_depth",
        "mean_depth",
        "median_depth",
        "min_depth",
        "max_depth",
        "breadth_1x",
        "breadth_5x",
        "breadth_10x",
        "breadth_20x",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for s in summaries:
            writer.writerow(asdict(s))


# -------------------------
# Relatório
# -------------------------


def md_table(rows: List[Dict[str, object]], columns: List[str]) -> str:
    if not rows:
        return "_Nenhum registro._\n"
    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join(["---"] * len(columns)) + " |"
    body = []
    for row in rows:
        vals = []
        for col in columns:
            value = row.get(col, "")
            if isinstance(value, float):
                value = f"{value:.3f}"
            vals.append(str(value))
        body.append("| " + " | ".join(vals) + " |")
    return "\n".join([header, sep] + body) + "\n"


def write_report(
    args,
    read_stats: List[FastqStats],
    ref_nt_stats: FastaStats,
    ref_prot_stats: FastaStats,
    assembly: Path,
    assembly_stats_obj: Optional[FastaStats],
    candidate_summaries: List[CandidateSummary],
    coverage_summaries: List[CoverageSummary],
    paths: Dict[str, Optional[Path]],
) -> Path:
    report = args.outdir / "REPORT.md"
    cmd_log = args.outdir / "commands.log"

    read_rows = []
    for st in read_stats:
        read_rows.append(
            {
                "file": Path(st.file).name,
                "reads_counted": st.reads_counted,
                "total_bases": st.total_bases,
                "min_len": st.min_len,
                "mean_len": st.mean_len,
                "median_len": st.median_len,
                "max_len": st.max_len,
                "N50": st.n50,
                "exact": st.exact,
            }
        )

    ref_rows = [
        {
            "file": Path(ref_nt_stats.file).name,
            "type": "CDS nucleotide",
            "sequences": ref_nt_stats.sequences,
            "total_bases/residues": ref_nt_stats.total_bases,
            "min": ref_nt_stats.min_len,
            "mean": ref_nt_stats.mean_len,
            "median": ref_nt_stats.median_len,
            "max": ref_nt_stats.max_len,
            "N50": ref_nt_stats.n50,
        },
        {
            "file": Path(ref_prot_stats.file).name,
            "type": "protein",
            "sequences": ref_prot_stats.sequences,
            "total_bases/residues": ref_prot_stats.total_bases,
            "min": ref_prot_stats.min_len,
            "mean": ref_prot_stats.mean_len,
            "median": ref_prot_stats.median_len,
            "max": ref_prot_stats.max_len,
            "N50": ref_prot_stats.n50,
        },
    ]

    candidate_rows = []
    for s in candidate_summaries[:30]:
        candidate_rows.append(
            {
                "contig": s.contig,
                "len": s.contig_len,
                "nt_genes": s.unique_nt_queries,
                "prot_genes": s.unique_prot_queries,
                "bitscore_total": s.total_bitscore,
                "best_nt_%id": s.best_nt_pident,
                "best_prot_%id": s.best_prot_pident,
                "terminal_overlap": s.circular_overlap_len,
                "overlap_%id": 100.0 * s.circular_overlap_identity,
                "circular_overlap": s.circular_by_terminal_overlap,
                "circular_flye": s.circular_by_flye,
            }
        )

    cov_rows = []
    for c in coverage_summaries[:30]:
        cov_rows.append(
            {
                "contig": c.contig,
                "len": c.length_from_depth,
                "mean_depth": c.mean_depth,
                "median_depth": c.median_depth,
                "min_depth": c.min_depth,
                "max_depth": c.max_depth,
                "breadth_1x_%": c.breadth_1x,
                "breadth_5x_%": c.breadth_5x,
                "breadth_10x_%": c.breadth_10x,
                "breadth_20x_%": c.breadth_20x,
            }
        )

    interpretation = []
    if not candidate_summaries:
        interpretation.append(
            "- Nenhum contig passou os filtros de BLAST/TBLASTN. Isso pode indicar baixa cobertura mitocondrial, montagem fragmentada, referência muito distante, ou filtros muito rígidos."
        )
    else:
        top = candidate_summaries[0]
        interpretation.append(
            f"- Melhor candidato: `{top.contig}` com {top.contig_len:,} bp e bitscore total {top.total_bitscore:.1f}.".replace(",", ".")
        )
        if len(candidate_summaries) == 1:
            interpretation.append("- Apenas um contig mitocondrial candidato foi encontrado; isso é compatível com uma montagem mitocondrial mais contínua.")
        else:
            interpretation.append(
                f"- Foram encontrados {len(candidate_summaries)} contigs candidatos; isso pode representar fragmentação, repeats, contigs alternativos ou regiões nucleares semelhantes a mitocôndria."
            )
        if top.circular_by_terminal_overlap or top.circular_by_flye is True:
            interpretation.append("- Há evidência automática de circularidade para o melhor candidato, mas ainda é recomendado validar manualmente com inspeção do BAM, Bandage quando houver grafo de montagem, e/ou dotplot.")
        else:
            interpretation.append("- Não houve evidência automática forte de circularidade no melhor candidato; isso não descarta circularidade biológica, apenas indica que o contig final não mostrou sobreposição terminal clara nos parâmetros usados.")

    if coverage_summaries:
        top_cov = coverage_summaries[0]
        interpretation.append(
            f"- Maior cobertura média observada: `{top_cov.contig}` com {top_cov.mean_depth:.2f}x e breadth >=1x de {top_cov.breadth_1x:.2f}%."
        )
        low = [c for c in coverage_summaries if c.breadth_1x < 95.0]
        if low:
            interpretation.append("- Pelo menos um contig candidato tem regiões sem cobertura completa; confira `depth_per_base.tsv` antes de chamar o genoma como final.")
    else:
        interpretation.append("- Cobertura não foi calculada ou não houve candidatos para mapear reads.")

    asm_block = "_Montagem ainda não avaliada, possivelmente por `--dry-run`._"
    if assembly_stats_obj:
        asm_block = md_table(
            [
                {
                    "assembly": str(assembly),
                    "contigs": assembly_stats_obj.sequences,
                    "total_bp": assembly_stats_obj.total_bases,
                    "min": assembly_stats_obj.min_len,
                    "mean": assembly_stats_obj.mean_len,
                    "median": assembly_stats_obj.median_len,
                    "max": assembly_stats_obj.max_len,
                    "N50": assembly_stats_obj.n50,
                }
            ],
            ["assembly", "contigs", "total_bp", "min", "mean", "median", "max", "N50"],
        )

    text = f"""
# Relatório do pipeline de montagem mitocondrial

**Data:** {now()}  
**Modo:** `{args.mode}`  
**Diretório de saída:** `{args.outdir}`

## 1. Objetivo

Montar contigs a partir das reads, identificar contigs mitocondriais usando os CDS de referência (`nt` e `prot`) e avaliar a sustentação dos candidatos por mapeamento das reads e cobertura por base.

## 2. Entradas

### Reads
{md_table(read_rows, ["file", "reads_counted", "total_bases", "min_len", "mean_len", "median_len", "max_len", "N50", "exact"])}

### Referências mitocondriais usadas como isca
{md_table(ref_rows, ["file", "type", "sequences", "total_bases/residues", "min", "mean", "median", "max", "N50"])}

## 3. Estatísticas da montagem

{asm_block}

## 4. Contigs mitocondriais candidatos

Os candidatos foram selecionados por similaridade contra os CDS mitocondriais de `Agabis_H97_nt.fasta` via BLASTN e contra `Agabis_H97_prot.fasta` via TBLASTN.

{md_table(candidate_rows, ["contig", "len", "nt_genes", "prot_genes", "bitscore_total", "best_nt_%id", "best_prot_%id", "terminal_overlap", "overlap_%id", "circular_overlap", "circular_flye"])}

Arquivo completo: `{paths.get('candidate_summary')}`  
FASTA dos candidatos: `{paths.get('candidates_fasta')}`

## 5. Cobertura dos candidatos

As reads foram mapeadas de volta contra os contigs mitocondriais candidatos. A cobertura por base foi calculada com `samtools depth -a`.

{md_table(cov_rows, ["contig", "len", "mean_depth", "median_depth", "min_depth", "max_depth", "breadth_1x_%", "breadth_5x_%", "breadth_10x_%", "breadth_20x_%"])}

BAM: `{paths.get('bam')}`  
Cobertura por base: `{paths.get('depth_tsv')}`  
Resumo `samtools coverage`: `{paths.get('coverage_tsv')}`

## 6. Interpretação automática

{chr(10).join(interpretation)}

## 7. Arquivos importantes

- Log dos comandos: `{cmd_log}`
- BLASTN nt vs assembly: `{paths.get('blastn')}`
- TBLASTN prot vs assembly: `{paths.get('tblastn')}`
- FASTA dos contigs candidatos: `{paths.get('candidates_fasta')}`
- Resumo dos candidatos: `{paths.get('candidate_summary')}`
- BAM de mapeamento: `{paths.get('bam')}`
- Profundidade por base: `{paths.get('depth_tsv')}`

## 8. Observações importantes

1. O FASTA final `mitochondrial_candidates.fasta` contém os contigs que parecem mitocondriais. Se houver mais de um contig, ainda não é uma molécula mitocondrial única e circularizada.
2. A evidência de circularidade calculada aqui é simples: busca sobreposição entre início e fim do contig. Valide o resultado com inspeção manual, dotplot, Bandage/assembly graph quando possível e cobertura uniforme.
3. CDS mitocondriais de uma espécie próxima ajudam a encontrar o contig mitocondrial, mas não garantem que todos os genes ou regiões intergênicas estejam completos.
4. Para entrega final, recomenda-se anotar o candidato com ferramentas específicas para mitogenoma, como MITOS/MFannot, e comparar ordem gênica e completude.
"""
    write_text(report, textwrap.dedent(text).strip() + "\n")
    return report


# -------------------------
# CLI
# -------------------------


def parse_args(argv: Optional[Sequence[str]] = None):
    parser = argparse.ArgumentParser(
        prog="mito_pipeline.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Pipeline Python para montagem e avaliação de genoma mitocondrial com Illumina ou PacBio.",
    )

    parser.add_argument("--mode", choices=["illumina", "pacbio"], required=True, help="Tipo de dado usado para montar.")
    parser.add_argument("--nt-ref", type=Path, required=True, help="FASTA de CDS mitocondriais nucleotídicos, ex.: Agabis_H97_nt.fasta.")
    parser.add_argument("--prot-ref", type=Path, required=True, help="FASTA de proteínas mitocondriais, ex.: Agabis_H97_prot.fasta.")
    parser.add_argument("--outdir", type=Path, required=True, help="Diretório de saída.")
    parser.add_argument("--threads", type=int, default=8, help="Número de threads.")
    parser.add_argument("--dry-run", action="store_true", help="Não executa comandos externos; apenas escreve comandos no log.")

    # Illumina
    parser.add_argument("--r1", type=Path, help="FASTQ R1 Illumina paired-end.")
    parser.add_argument("--r2", type=Path, help="FASTQ R2 Illumina paired-end.")
    parser.add_argument("--memory-gb", type=int, default=0, help="Memória máxima para SPAdes em GB. 0 deixa SPAdes decidir.")
    parser.add_argument("--spades-k", default="", help="Lista de k-mers para SPAdes, ex.: 21,33,55,77,99,127. Vazio usa padrão do SPAdes.")
    parser.add_argument("--spades-output", choices=["scaffolds", "contigs"], default="scaffolds", help="Arquivo de saída do SPAdes a usar como assembly.")

    # PacBio
    parser.add_argument("--pacbio", type=Path, help="FASTQ/FASTA PacBio.")
    parser.add_argument("--pacbio-type", choices=["raw", "hifi", "corrected"], default="raw", help="Tipo de reads PacBio para escolher o preset do Flye/minimap2.")
    parser.add_argument("--flye-genome-size", default="80k", help="Tamanho esperado aproximado do mitogenoma para Flye, ex.: 40k, 80k, 120k.")
    parser.add_argument("--flye-meta", action="store_true", help="Usa modo metagenômico do Flye; útil quando o DNA mitocondrial é mistura pequena dentro de DNA total.")

    # Estatísticas
    parser.add_argument("--stats-max-reads", type=int, default=0, help="Máximo de reads para estatísticas FASTQ. 0 = arquivo completo.")

    # Filtros BLAST
    parser.add_argument("--evalue", type=float, default=1e-5, help="E-value máximo para BLAST/TBLASTN.")
    parser.add_argument("--max-target-seqs", type=int, default=20, help="Máximo de alvos por query no BLAST.")
    parser.add_argument("--min-bitscore", type=float, default=50.0, help="Bitscore mínimo por hit.")
    parser.add_argument("--min-hit-qcov", type=float, default=0.30, help="Cobertura mínima da query por hit, de 0 a 1.")
    parser.add_argument("--min-nt-pident", type=float, default=70.0, help="Identidade mínima para hits BLASTN.")
    parser.add_argument("--min-prot-pident", type=float, default=30.0, help="Identidade mínima para hits TBLASTN.")
    parser.add_argument("--min-candidate-len", type=int, default=1000, help="Comprimento mínimo do contig candidato.")

    # Circularidade
    parser.add_argument("--min-circular-overlap", type=int, default=500, help="Sobreposição terminal mínima para teste simples de circularidade.")
    parser.add_argument("--max-circular-overlap", type=int, default=5000, help="Sobreposição terminal máxima para teste simples de circularidade.")
    parser.add_argument("--min-circular-identity", type=float, default=0.98, help="Identidade mínima da sobreposição terminal, de 0 a 1.")

    args = parser.parse_args(argv)

    if args.mode == "illumina":
        if not args.r1 or not args.r2:
            parser.error("--mode illumina requer --r1 e --r2.")
    if args.mode == "pacbio":
        if not args.pacbio:
            parser.error("--mode pacbio requer --pacbio.")
    for p in [args.nt_ref, args.prot_ref, args.r1, args.r2, args.pacbio]:
        if p is not None and not p.exists():
            parser.error(f"Arquivo não encontrado: {p}")
    if args.threads < 1:
        parser.error("--threads deve ser >= 1.")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    ensure_dir(args.outdir)
    log_file = args.outdir / "commands.log"
    write_text(log_file, f"# Log de comandos - {now()}\n")

    try:
        check_dependencies(args.mode)
        with log_file.open("a", encoding="utf-8") as log:
            log.write("\n# Argumentos\n")
            log.write(json.dumps({k: str(v) for k, v in vars(args).items()}, indent=2) + "\n")

        print(f"[{now()}] Calculando estatísticas das entradas...")
        read_stats: List[FastqStats] = []
        if args.mode == "illumina":
            read_stats.append(fastq_stats(args.r1, args.stats_max_reads))
            read_stats.append(fastq_stats(args.r2, args.stats_max_reads))
        else:
            read_stats.append(fastq_stats(args.pacbio, args.stats_max_reads))
        ref_nt_stats = fasta_stats(args.nt_ref)
        ref_prot_stats = fasta_stats(args.prot_ref)

        print(f"[{now()}] Rodando montagem {args.mode}...")
        if args.mode == "illumina":
            assembly = run_illumina_assembly(args, log_file)
        else:
            assembly = run_pacbio_assembly(args, log_file)

        assembly_stats_obj = None
        if not args.dry_run and assembly.exists():
            assembly_stats_obj = fasta_stats(assembly)

        print(f"[{now()}] Buscando contigs mitocondriais com BLASTN/TBLASTN...")
        blastn_out, tblastn_out = run_blasts(args, assembly, log_file)

        print(f"[{now()}] Selecionando candidatos mitocondriais...")
        if args.dry_run:
            candidates_fasta = args.outdir / "03_mito_candidates" / "mitochondrial_candidates.fasta"
            candidate_summary_tsv = args.outdir / "03_mito_candidates" / "mitochondrial_candidates_summary.tsv"
            candidate_summaries: List[CandidateSummary] = []
        else:
            candidates_fasta, candidate_summary_tsv, candidate_summaries = select_candidates(
                args, assembly, blastn_out, tblastn_out
            )

        print(f"[{now()}] Mapeando reads de volta aos candidatos e calculando cobertura...")
        bam, depth_tsv, coverage_tsv, coverage_summaries = map_reads_and_coverage(args, candidates_fasta, log_file)

        paths = {
            "blastn": blastn_out,
            "tblastn": tblastn_out,
            "candidates_fasta": candidates_fasta,
            "candidate_summary": candidate_summary_tsv,
            "bam": bam,
            "depth_tsv": depth_tsv,
            "coverage_tsv": coverage_tsv,
        }
        report = write_report(
            args,
            read_stats,
            ref_nt_stats,
            ref_prot_stats,
            assembly,
            assembly_stats_obj,
            candidate_summaries,
            coverage_summaries,
            paths,
        )

        print("\nPipeline finalizado.")
        print(f"Relatório: {report}")
        print(f"Log de comandos: {log_file}")
        print(f"FASTA candidatos: {candidates_fasta}")
        return 0

    except PipelineError as exc:
        sys.stderr.write(f"\nERRO: {exc}\n")
        sys.stderr.write(f"Consulte também o log: {log_file}\n")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

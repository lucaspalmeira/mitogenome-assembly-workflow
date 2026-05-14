#!/usr/bin/env python3
"""
mito_extract_reassemble.py

Segunda etapa do pipeline de montagem mitocondrial.

Objetivo:
  1. Usar contigs mitocondriais candidatos ou CDS mitocondriais como isca.
  2. Mapear reads brutas/originais contra essas iscas.
  3. Extrair as reads que mapearam nas iscas.
  4. Remontar apenas o subconjunto enriquecido de reads mitocondriais.
  5. Avaliar cobertura, circularidade simples e hits contra CDS de referência.

Decisão de alinhador:
  - Illumina paired-end: BWA-MEM
  - PacBio: minimap2

Autor: gerado para Lucas Palmeira
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import shutil
import statistics
import subprocess
import sys
import textwrap
import time
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

try:
    from .common import parse_flye_circular_flags, terminal_overlap, write_fasta
except ImportError:
    from common import parse_flye_circular_flags, terminal_overlap, write_fasta


# =============================================================================
# Utilitários gerais
# =============================================================================


class PipelineError(RuntimeError):
    pass


@dataclass
class FastqStats:
    file: str
    reads_counted: int
    total_bases: int
    min_len: int
    mean_len: float
    median_len: float
    max_len: int
    n50: int
    exact: bool


@dataclass
class FastaStats:
    file: str
    sequences: int
    total_bases: int
    min_len: int
    mean_len: float
    median_len: float
    max_len: int
    n50: int


@dataclass
class ExtractionStats:
    mode: str
    bait_fasta: str
    mapped_read_names: int
    input_reads_r1_or_single: int
    input_reads_r2: int
    extracted_reads_r1_or_single: int
    extracted_reads_r2: int
    extracted_fraction_r1_or_single: float
    extracted_fraction_r2: float
    extracted_bases_r1_or_single: int
    extracted_bases_r2: int


@dataclass
class CircularityCall:
    contig: str
    length: int
    best_terminal_overlap: int
    best_terminal_identity_pct: float
    circular_by_terminal_overlap: bool
    circular_by_assembler: Optional[bool]


@dataclass
class BlastCandidate:
    contig: str
    length: int
    nt_genes: int
    prot_genes: int
    blastn_hits: int
    tblastn_hits: int
    total_bitscore: float
    best_nt_pident: float
    best_prot_pident: float


@dataclass
class CoverageSummary:
    contig: str
    length_from_depth: int
    mean_depth: float
    median_depth: float
    min_depth: int
    max_depth: int
    breadth_1x_pct: float
    breadth_5x_pct: float
    breadth_10x_pct: float
    breadth_20x_pct: float


def now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def open_text_auto(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def open_write_text_auto(path: Path):
    ensure_dir(path.parent)
    if str(path).endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8")
    return path.open("wt", encoding="utf-8")


def write_text(path: Path, text: str) -> None:
    ensure_dir(path.parent)
    path.write_text(text, encoding="utf-8")


def which_or_fail(program: str) -> None:
    if shutil.which(program) is None:
        raise PipelineError(
            f"Dependência não encontrada no PATH: {program}. "
            "Instale via mamba/conda ou carregue o módulo correspondente."
        )


def run_cmd(
    cmd: Sequence[str],
    log_file: Path,
    cwd: Optional[Path] = None,
    stdout_file: Optional[Path] = None,
    stderr_file: Optional[Path] = None,
    dry_run: bool = False,
) -> None:
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
            ensure_dir(stdout_file.parent)
            stdout_handle = stdout_file.open("w", encoding="utf-8")
        if stderr_file:
            ensure_dir(stderr_file.parent)
            stderr_handle = stderr_file.open("w", encoding="utf-8")

        result = subprocess.run(
            [str(x) for x in cmd],
            cwd=str(cwd) if cwd else None,
            stdout=stdout_handle if stdout_handle else subprocess.PIPE,
            stderr=stderr_handle if stderr_handle else subprocess.PIPE,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            msg = f"Comando falhou com código {result.returncode}: {cmd_str}"
            if not stderr_file and result.stderr:
                msg += "\nSTDERR:\n" + result.stderr[-5000:]
            raise PipelineError(msg)
        if not stdout_file and result.stdout:
            with log_file.open("a", encoding="utf-8") as log:
                log.write(result.stdout[-5000:] + "\n")
        if not stderr_file and result.stderr:
            with log_file.open("a", encoding="utf-8") as log:
                log.write(result.stderr[-5000:] + "\n")
    finally:
        if stdout_handle:
            stdout_handle.close()
        if stderr_handle:
            stderr_handle.close()


def run_shell(shell_cmd: str, log_file: Path, stderr_file: Optional[Path] = None, dry_run: bool = False) -> None:
    with log_file.open("a", encoding="utf-8") as log:
        log.write(f"\n[{now()}] $ {shell_cmd}\n")
    if dry_run:
        return
    result = subprocess.run(shell_cmd, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if stderr_file:
        write_text(stderr_file, result.stderr)
    if result.returncode != 0:
        raise PipelineError(f"Comando falhou com código {result.returncode}:\n{shell_cmd}\n\nSTDERR:\n{result.stderr[-5000:]}")


def n50(lengths: Sequence[int]) -> int:
    if not lengths:
        return 0
    total = sum(lengths)
    cutoff = total / 2
    acc = 0
    for length in sorted(lengths, reverse=True):
        acc += length
        if acc >= cutoff:
            return int(length)
    return 0


# =============================================================================
# FASTA / FASTQ
# =============================================================================


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
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks)


def fasta_dict(path: Path) -> Dict[str, str]:
    return {h.split()[0]: seq for h, seq in iter_fasta(path)}


def fasta_stats(path: Path) -> FastaStats:
    lengths = [len(seq) for _, seq in iter_fasta(path)]
    if not lengths:
        return FastaStats(str(path), 0, 0, 0, 0.0, 0.0, 0, 0)
    return FastaStats(
        file=str(path),
        sequences=len(lengths),
        total_bases=sum(lengths),
        min_len=min(lengths),
        mean_len=float(statistics.mean(lengths)),
        median_len=float(statistics.median(lengths)),
        max_len=max(lengths),
        n50=n50(lengths),
    )


def fastq_name(header: str) -> str:
    """Normaliza o nome da read para comparar FASTQ e BAM/SAM."""
    h = header.strip()
    if h.startswith("@"):
        h = h[1:]
    h = h.split()[0]
    if h.endswith("/1") or h.endswith("/2"):
        h = h[:-2]
    return h


def read_fastq_record(handle) -> Optional[Tuple[str, str, str, str]]:
    h = handle.readline()
    if not h:
        return None
    s = handle.readline()
    p = handle.readline()
    q = handle.readline()
    if not q:
        raise PipelineError("FASTQ truncado ou mal formatado.")
    return h, s, p, q


def fastq_stats(path: Path, max_reads: int = 0) -> FastqStats:
    lengths: List[int] = []
    exact = True
    with open_text_auto(path) as handle:
        while True:
            rec = read_fastq_record(handle)
            if rec is None:
                break
            _, seq, _, _ = rec
            lengths.append(len(seq.strip()))
            if max_reads and len(lengths) >= max_reads:
                exact = False
                break
    if not lengths:
        return FastqStats(str(path), 0, 0, 0, 0.0, 0.0, 0, 0, exact)
    return FastqStats(
        file=str(path),
        reads_counted=len(lengths),
        total_bases=sum(lengths),
        min_len=min(lengths),
        mean_len=float(statistics.mean(lengths)),
        median_len=float(statistics.median(lengths)),
        max_len=max(lengths),
        n50=n50(lengths),
        exact=exact,
    )


def combine_baits(args) -> Path:
    bait_dir = args.outdir / "00_bait_reference"
    ensure_dir(bait_dir)
    bait_out = bait_dir / "bait_reference.fasta"
    records: List[Tuple[str, str]] = []

    for idx, path in enumerate(args.bait_fasta or [], start=1):
        for header, seq in iter_fasta(path):
            records.append((f"bait{idx}|{header}", seq))

    if args.nt_ref and args.include_nt_ref_as_bait:
        for header, seq in iter_fasta(args.nt_ref):
            records.append((f"nt_ref|{header}", seq))

    if not records:
        raise PipelineError("Nenhuma sequência de isca foi fornecida. Use --bait-fasta e/ou --nt-ref --include-nt-ref-as-bait.")

    write_fasta(records, bait_out)
    return bait_out


# =============================================================================
# Mapeamento contra isca e extração de reads
# =============================================================================


def check_dependencies(args) -> None:
    common = ["samtools"]
    if args.mode == "illumina":
        common += ["bwa", "spades.py"]
    else:
        common += ["minimap2", "flye"]
    if args.evaluate_with_blast:
        common += ["makeblastdb", "blastn", "tblastn"]
    # Para mapear reads no assembly final em modo PacBio usa minimap2; para Illumina usa BWA.
    for program in common:
        which_or_fail(program)


def map_illumina_to_bait(args, bait_fasta: Path, log_file: Path) -> Path:
    map_dir = args.outdir / "01_map_reads_to_bait"
    ensure_dir(map_dir)
    bam = map_dir / "illumina_reads_vs_bait.sorted.bam"

    run_cmd(["bwa", "index", str(bait_fasta)], log_file, dry_run=args.dry_run)
    shell_cmd = (
        f"bwa mem -t {args.threads} {bait_fasta} {args.r1} {args.r2} "
        f"| samtools sort -@ {args.threads} -o {bam} -"
    )
    run_shell(shell_cmd, log_file, stderr_file=map_dir / "bwa_mem_samtools_sort.stderr.log", dry_run=args.dry_run)
    run_cmd(["samtools", "index", str(bam)], log_file, dry_run=args.dry_run)
    run_cmd(["samtools", "flagstat", str(bam)], log_file, stdout_file=map_dir / "samtools_flagstat.txt", dry_run=args.dry_run)
    return bam


def map_pacbio_to_bait(args, bait_fasta: Path, log_file: Path) -> Path:
    map_dir = args.outdir / "01_map_reads_to_bait"
    ensure_dir(map_dir)
    bam = map_dir / "pacbio_reads_vs_bait.sorted.bam"
    preset = "map-hifi" if args.pacbio_type == "hifi" else "map-pb"
    shell_cmd = (
        f"minimap2 -ax {preset} -t {args.threads} {bait_fasta} {args.pacbio} "
        f"| samtools sort -@ {args.threads} -o {bam} -"
    )
    run_shell(shell_cmd, log_file, stderr_file=map_dir / "minimap2_samtools_sort.stderr.log", dry_run=args.dry_run)
    run_cmd(["samtools", "index", str(bam)], log_file, dry_run=args.dry_run)
    run_cmd(["samtools", "flagstat", str(bam)], log_file, stdout_file=map_dir / "samtools_flagstat.txt", dry_run=args.dry_run)
    return bam


def mapped_read_names_from_bam(bam: Path, min_mapq: int, log_file: Path, out_txt: Path) -> Set[str]:
    """Coleta nomes das reads alinhadas com MAPQ mínimo."""
    cmd = ["samtools", "view", "-F", "4", "-q", str(min_mapq), str(bam)]
    with log_file.open("a", encoding="utf-8") as log:
        log.write(f"\n[{now()}] $ {' '.join(cmd)} > {out_txt}\n")
    proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if proc.returncode != 0:
        raise PipelineError(f"Falha ao listar reads mapeadas:\n{proc.stderr[-5000:]}")
    names = set()
    for line in proc.stdout.splitlines():
        if not line:
            continue
        qname = line.split("\t", 1)[0]
        names.add(fastq_name(qname))
    ensure_dir(out_txt.parent)
    with out_txt.open("w", encoding="utf-8") as handle:
        for name in sorted(names):
            handle.write(name + "\n")
    return names


def extract_illumina_pairs(args, names: Set[str]) -> Tuple[Path, Path, ExtractionStats]:
    out_dir = args.outdir / "02_extracted_mito_reads"
    ensure_dir(out_dir)
    out_r1 = out_dir / "mito_enriched_R1.fastq.gz"
    out_r2 = out_dir / "mito_enriched_R2.fastq.gz"

    total_pairs = 0
    kept_pairs = 0
    bases_r1 = 0
    bases_r2 = 0

    with open_text_auto(args.r1) as h1, open_text_auto(args.r2) as h2, open_write_text_auto(out_r1) as o1, open_write_text_auto(out_r2) as o2:
        while True:
            rec1 = read_fastq_record(h1)
            rec2 = read_fastq_record(h2)
            if rec1 is None and rec2 is None:
                break
            if rec1 is None or rec2 is None:
                raise PipelineError("R1 e R2 têm número diferente de reads.")
            n1 = fastq_name(rec1[0])
            n2 = fastq_name(rec2[0])
            if n1 != n2:
                raise PipelineError(f"R1/R2 fora de sincronia: {n1} != {n2}")
            total_pairs += 1
            if n1 in names:
                kept_pairs += 1
                o1.write("".join(rec1))
                o2.write("".join(rec2))
                bases_r1 += len(rec1[1].strip())
                bases_r2 += len(rec2[1].strip())

    stats = ExtractionStats(
        mode="illumina",
        bait_fasta=str(args._bait_fasta_for_stats),
        mapped_read_names=len(names),
        input_reads_r1_or_single=total_pairs,
        input_reads_r2=total_pairs,
        extracted_reads_r1_or_single=kept_pairs,
        extracted_reads_r2=kept_pairs,
        extracted_fraction_r1_or_single=100.0 * kept_pairs / total_pairs if total_pairs else 0.0,
        extracted_fraction_r2=100.0 * kept_pairs / total_pairs if total_pairs else 0.0,
        extracted_bases_r1_or_single=bases_r1,
        extracted_bases_r2=bases_r2,
    )
    return out_r1, out_r2, stats


def extract_pacbio_reads(args, names: Set[str]) -> Tuple[Path, ExtractionStats]:
    out_dir = args.outdir / "02_extracted_mito_reads"
    ensure_dir(out_dir)
    out_pb = out_dir / "mito_enriched_pacbio.fastq.gz"
    total = 0
    kept = 0
    bases = 0
    with open_text_auto(args.pacbio) as hin, open_write_text_auto(out_pb) as hout:
        while True:
            rec = read_fastq_record(hin)
            if rec is None:
                break
            total += 1
            name = fastq_name(rec[0])
            if name in names:
                kept += 1
                hout.write("".join(rec))
                bases += len(rec[1].strip())
    stats = ExtractionStats(
        mode="pacbio",
        bait_fasta=str(args._bait_fasta_for_stats),
        mapped_read_names=len(names),
        input_reads_r1_or_single=total,
        input_reads_r2=0,
        extracted_reads_r1_or_single=kept,
        extracted_reads_r2=0,
        extracted_fraction_r1_or_single=100.0 * kept / total if total else 0.0,
        extracted_fraction_r2=0.0,
        extracted_bases_r1_or_single=bases,
        extracted_bases_r2=0,
    )
    return out_pb, stats


def write_extraction_stats(stats: ExtractionStats, path: Path) -> None:
    ensure_dir(path.parent)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(stats).keys()), delimiter="\t")
        writer.writeheader()
        writer.writerow(asdict(stats))


# =============================================================================
# Remontagem
# =============================================================================


def run_spades_reassembly(args, r1: Path, r2: Path, log_file: Path) -> Path:
    asm_dir = args.outdir / "03_reassembly_spades_mito_reads"
    ensure_dir(asm_dir)
    cmd = [
        "spades.py",
        "--careful",
        "-1",
        str(r1),
        "-2",
        str(r2),
        "-o",
        str(asm_dir),
        "--threads",
        str(args.threads),
    ]
    if args.memory_gb:
        cmd += ["--memory", str(args.memory_gb)]
    if args.spades_k:
        cmd += ["-k", args.spades_k]
    if args.use_bait_as_trusted_contigs:
        cmd += ["--trusted-contigs", str(args._bait_fasta_for_stats)]
    run_cmd(cmd, log_file, dry_run=args.dry_run)
    assembly = asm_dir / ("contigs.fasta" if args.spades_output == "contigs" else "scaffolds.fasta")
    if not args.dry_run and not assembly.exists():
        fallback = asm_dir / "contigs.fasta"
        if fallback.exists():
            assembly = fallback
        else:
            raise PipelineError(f"SPAdes terminou, mas não encontrei {assembly} nem {fallback}.")
    return assembly


def run_flye_reassembly(args, pacbio_reads: Path, log_file: Path) -> Path:
    asm_dir = args.outdir / "03_reassembly_flye_mito_reads"
    ensure_dir(asm_dir)
    flag = {"raw": "--pacbio-raw", "hifi": "--pacbio-hifi", "corrected": "--pacbio-corr"}[args.pacbio_type]
    cmd = [
        "flye",
        flag,
        str(pacbio_reads),
        "--out-dir",
        str(asm_dir),
        "--threads",
        str(args.threads),
        "--genome-size",
        args.flye_genome_size,
    ]
    if args.flye_meta:
        cmd.append("--meta")
    run_cmd(cmd, log_file, dry_run=args.dry_run)
    assembly = asm_dir / "assembly.fasta"
    if not args.dry_run and not assembly.exists():
        raise PipelineError(f"Flye terminou, mas não encontrei {assembly}.")
    return assembly


# =============================================================================
# Avaliação: circularidade, BLAST, cobertura
# =============================================================================


def evaluate_circularity(args, assembly: Path) -> List[CircularityCall]:
    seqs = fasta_dict(assembly)
    flye_flags: Dict[str, Optional[bool]] = {}
    if args.mode == "pacbio":
        flye_flags = parse_flye_circular_flags(args.outdir / "03_reassembly_flye_mito_reads")
    calls: List[CircularityCall] = []
    for contig, seq in seqs.items():
        ov_len, ov_id, is_circ = terminal_overlap(seq, args.min_circular_overlap, args.max_circular_overlap, args.min_circular_identity)
        calls.append(
            CircularityCall(
                contig=contig,
                length=len(seq),
                best_terminal_overlap=ov_len,
                best_terminal_identity_pct=100.0 * ov_id,
                circular_by_terminal_overlap=is_circ,
                circular_by_assembler=flye_flags.get(contig),
            )
        )
    calls.sort(key=lambda x: (x.circular_by_terminal_overlap or x.circular_by_assembler is True, x.length), reverse=True)
    return calls


def write_circularity_tsv(calls: List[CircularityCall], path: Path) -> None:
    ensure_dir(path.parent)
    fields = list(asdict(calls[0]).keys()) if calls else [
        "contig", "length", "best_terminal_overlap", "best_terminal_identity_pct", "circular_by_terminal_overlap", "circular_by_assembler"
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for call in calls:
            writer.writerow(asdict(call))


def run_blast_evaluation(args, assembly: Path, log_file: Path) -> Tuple[Optional[Path], Optional[Path], List[BlastCandidate]]:
    if not args.evaluate_with_blast:
        return None, None, []
    if not args.nt_ref or not args.prot_ref:
        return None, None, []

    blast_dir = args.outdir / "04_evaluate_reassembly" / "blast"
    ensure_dir(blast_dir)
    db_prefix = blast_dir / "reassembly_db"
    run_cmd(["makeblastdb", "-in", str(assembly), "-dbtype", "nucl", "-out", str(db_prefix)], log_file, dry_run=args.dry_run)
    fields = "qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"
    blastn_out = blast_dir / "blastn_nt_vs_reassembly.tsv"
    tblastn_out = blast_dir / "tblastn_prot_vs_reassembly.tsv"
    run_cmd(
        [
            "blastn", "-query", str(args.nt_ref), "-db", str(db_prefix), "-outfmt", f"6 {fields}",
            "-evalue", str(args.evalue), "-max_target_seqs", str(args.max_target_seqs),
            "-num_threads", str(args.threads), "-out", str(blastn_out),
        ],
        log_file,
        dry_run=args.dry_run,
    )
    run_cmd(
        [
            "tblastn", "-query", str(args.prot_ref), "-db", str(db_prefix), "-outfmt", f"6 {fields}",
            "-evalue", str(args.evalue), "-max_target_seqs", str(args.max_target_seqs),
            "-num_threads", str(args.threads), "-out", str(tblastn_out),
        ],
        log_file,
        dry_run=args.dry_run,
    )
    if args.dry_run:
        return blastn_out, tblastn_out, []
    candidates = summarize_blast_candidates(args, assembly, blastn_out, tblastn_out)
    write_blast_candidates_tsv(candidates, blast_dir / "reassembly_mito_candidates_summary.tsv")
    return blastn_out, tblastn_out, candidates


def parse_blast_table(path: Path, source: str) -> List[dict]:
    hits = []
    if not path or not path.exists():
        return hits
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 12:
                continue
            qlen = int(float(parts[4])) if parts[4] else 0
            aln_len = int(float(parts[3])) if parts[3] else 0
            hits.append(
                {
                    "source": source,
                    "qseqid": parts[0],
                    "sseqid": parts[1],
                    "pident": float(parts[2]),
                    "length": aln_len,
                    "qlen": qlen,
                    "qcov": aln_len / qlen if qlen else 0.0,
                    "evalue": float(parts[10]),
                    "bitscore": float(parts[11]),
                }
            )
    return hits


def summarize_blast_candidates(args, assembly: Path, blastn_out: Path, tblastn_out: Path) -> List[BlastCandidate]:
    seqs = fasta_dict(assembly)
    all_hits = parse_blast_table(blastn_out, "blastn") + parse_blast_table(tblastn_out, "tblastn")
    by_contig: Dict[str, List[dict]] = defaultdict(list)
    for hit in all_hits:
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
        by_contig[hit["sseqid"]].append(hit)

    candidates: List[BlastCandidate] = []
    for contig, hits in by_contig.items():
        seq = seqs.get(contig)
        if not seq:
            continue
        nt_hits = [h for h in hits if h["source"] == "blastn"]
        prot_hits = [h for h in hits if h["source"] == "tblastn"]
        candidates.append(
            BlastCandidate(
                contig=contig,
                length=len(seq),
                nt_genes=len({h["qseqid"] for h in nt_hits}),
                prot_genes=len({h["qseqid"] for h in prot_hits}),
                blastn_hits=len(nt_hits),
                tblastn_hits=len(prot_hits),
                total_bitscore=sum(h["bitscore"] for h in hits),
                best_nt_pident=max([h["pident"] for h in nt_hits] or [0.0]),
                best_prot_pident=max([h["pident"] for h in prot_hits] or [0.0]),
            )
        )
    candidates.sort(key=lambda x: (x.total_bitscore, x.length), reverse=True)
    return candidates


def write_blast_candidates_tsv(candidates: List[BlastCandidate], path: Path) -> None:
    ensure_dir(path.parent)
    fields = list(asdict(candidates[0]).keys()) if candidates else [
        "contig", "length", "nt_genes", "prot_genes", "blastn_hits", "tblastn_hits", "total_bitscore", "best_nt_pident", "best_prot_pident"
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for c in candidates:
            writer.writerow(asdict(c))


def map_original_reads_to_reassembly(args, assembly: Path, log_file: Path) -> Tuple[Path, Path, Path, List[CoverageSummary]]:
    cov_dir = args.outdir / "04_evaluate_reassembly" / "coverage_original_reads_vs_reassembly"
    ensure_dir(cov_dir)
    bam = cov_dir / "original_reads_vs_reassembly.sorted.bam"

    run_cmd(["samtools", "faidx", str(assembly)], log_file, dry_run=args.dry_run)

    if args.mode == "illumina":
        run_cmd(["bwa", "index", str(assembly)], log_file, dry_run=args.dry_run)
        shell_cmd = (
            f"bwa mem -t {args.threads} {assembly} {args.r1} {args.r2} "
            f"| samtools sort -@ {args.threads} -o {bam} -"
        )
    else:
        preset = "map-hifi" if args.pacbio_type == "hifi" else "map-pb"
        shell_cmd = (
            f"minimap2 -ax {preset} -t {args.threads} {assembly} {args.pacbio} "
            f"| samtools sort -@ {args.threads} -o {bam} -"
        )
    run_shell(shell_cmd, log_file, stderr_file=cov_dir / "mapping.stderr.log", dry_run=args.dry_run)
    run_cmd(["samtools", "index", str(bam)], log_file, dry_run=args.dry_run)
    run_cmd(["samtools", "flagstat", str(bam)], log_file, stdout_file=cov_dir / "samtools_flagstat.txt", dry_run=args.dry_run)
    coverage_tsv = cov_dir / "samtools_coverage.tsv"
    depth_tsv = cov_dir / "depth_per_base.tsv"
    run_cmd(["samtools", "coverage", str(bam)], log_file, stdout_file=coverage_tsv, dry_run=args.dry_run)
    run_cmd(["samtools", "depth", "-a", str(bam)], log_file, stdout_file=depth_tsv, dry_run=args.dry_run)
    summaries = parse_depth(depth_tsv) if not args.dry_run else []
    write_coverage_summary_tsv(summaries, cov_dir / "coverage_summary.tsv")
    return bam, coverage_tsv, depth_tsv, summaries


def parse_depth(depth_tsv: Path) -> List[CoverageSummary]:
    values_by_contig: Dict[str, List[int]] = defaultdict(list)
    if not depth_tsv.exists():
        return []
    with depth_tsv.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            values_by_contig[parts[0]].append(int(float(parts[2])))
    summaries = []
    for contig, values in values_by_contig.items():
        if not values:
            continue
        length = len(values)
        summaries.append(
            CoverageSummary(
                contig=contig,
                length_from_depth=length,
                mean_depth=float(statistics.mean(values)),
                median_depth=float(statistics.median(values)),
                min_depth=min(values),
                max_depth=max(values),
                breadth_1x_pct=100.0 * sum(v >= 1 for v in values) / length,
                breadth_5x_pct=100.0 * sum(v >= 5 for v in values) / length,
                breadth_10x_pct=100.0 * sum(v >= 10 for v in values) / length,
                breadth_20x_pct=100.0 * sum(v >= 20 for v in values) / length,
            )
        )
    summaries.sort(key=lambda x: (x.mean_depth, x.length_from_depth), reverse=True)
    return summaries


def write_coverage_summary_tsv(summaries: List[CoverageSummary], path: Path) -> None:
    ensure_dir(path.parent)
    fields = list(asdict(summaries[0]).keys()) if summaries else [
        "contig", "length_from_depth", "mean_depth", "median_depth", "min_depth", "max_depth", "breadth_1x_pct", "breadth_5x_pct", "breadth_10x_pct", "breadth_20x_pct"
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for s in summaries:
            writer.writerow(asdict(s))


# =============================================================================
# Relatório
# =============================================================================


def md_table(rows: List[Dict[str, object]], columns: List[str]) -> str:
    if not rows:
        return "_Nenhum registro._\n"
    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join(["---"] * len(columns)) + " |"
    body = []
    for row in rows:
        vals = []
        for col in columns:
            val = row.get(col, "")
            if isinstance(val, float):
                val = f"{val:.3f}"
            vals.append(str(val))
        body.append("| " + " | ".join(vals) + " |")
    return "\n".join([header, sep] + body) + "\n"


def write_report(
    args,
    bait_fasta: Path,
    bait_stats: FastaStats,
    input_read_stats: List[FastqStats],
    extraction_stats: ExtractionStats,
    assembly: Path,
    assembly_stats_obj: Optional[FastaStats],
    circularity_calls: List[CircularityCall],
    blast_candidates: List[BlastCandidate],
    coverage_summaries: List[CoverageSummary],
    paths: Dict[str, Optional[Path]],
) -> Path:
    report = args.outdir / "REPORT_REASSEMBLY.md"

    read_rows = [asdict(x) for x in input_read_stats]
    for row in read_rows:
        row["file"] = Path(row["file"]).name

    extract_row = asdict(extraction_stats)
    bait_row = asdict(bait_stats)
    bait_row["file"] = Path(bait_row["file"]).name

    asm_block = "_Montagem ainda não avaliada._\n"
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

    circ_rows = [asdict(c) for c in circularity_calls[:30]]
    blast_rows = [asdict(c) for c in blast_candidates[:30]]
    cov_rows = [asdict(c) for c in coverage_summaries[:30]]

    interpretation: List[str] = []
    interpretation.append(
        f"- Foram extraídas {extraction_stats.extracted_reads_r1_or_single} reads R1/single "
        f"({extraction_stats.extracted_fraction_r1_or_single:.3f}% do total) a partir do mapeamento contra a isca."
    )
    if args.mode == "illumina":
        interpretation.append(
            "- Como o dado é Illumina paired-end, o script extraiu o par completo quando pelo menos um mate mapeou contra a isca mitocondrial."
        )
    if assembly_stats_obj:
        interpretation.append(
            f"- A remontagem enriquecida gerou {assembly_stats_obj.sequences} contigs/scaffolds, totalizando {assembly_stats_obj.total_bases} bp, com N50 de {assembly_stats_obj.n50} bp."
        )
    if blast_candidates:
        top = blast_candidates[0]
        interpretation.append(
            f"- Melhor candidato por BLAST na remontagem: `{top.contig}` com {top.length} bp, {top.nt_genes} CDS nucleotídicos, {top.prot_genes} proteínas e bitscore total {top.total_bitscore:.1f}."
        )
    elif args.evaluate_with_blast:
        interpretation.append("- Nenhum candidato passou os filtros de BLAST/TBLASTN na remontagem; avalie reduzir filtros ou conferir o FASTA de isca.")
    if circularity_calls:
        circ_yes = [c for c in circularity_calls if c.circular_by_terminal_overlap or c.circular_by_assembler is True]
        if circ_yes:
            best = circ_yes[0]
            interpretation.append(
                f"- Há evidência automática de circularidade em `{best.contig}`: overlap terminal de {best.best_terminal_overlap} bp com {best.best_terminal_identity_pct:.2f}% de identidade."
            )
        else:
            best = circularity_calls[0]
            interpretation.append(
                f"- Não foi detectada circularidade terminal forte. O melhor overlap observado foi em `{best.contig}`: {best.best_terminal_overlap} bp com {best.best_terminal_identity_pct:.2f}% de identidade."
            )
    if coverage_summaries:
        top_cov = coverage_summaries[0]
        interpretation.append(
            f"- Maior cobertura média com reads originais contra a remontagem: `{top_cov.contig}` com {top_cov.mean_depth:.2f}x e breadth >=1x de {top_cov.breadth_1x_pct:.2f}%."
        )

    text = f"""
# Relatório da etapa de enriquecimento e remontagem mitocondrial

**Data:** {now()}  
**Modo:** `{args.mode}`  
**Diretório de saída:** `{args.outdir}`

## 1. O que esta etapa faz

Esta etapa usa uma ou mais sequências mitocondriais candidatas como **isca** para recuperar, a partir das reads originais, as reads que provavelmente vêm da mitocôndria. Em seguida, ela remonta apenas esse subconjunto enriquecido.

No modo Illumina, o mapeamento contra a isca é feito com **BWA-MEM**. No modo PacBio, é feito com **minimap2**.

## 2. Isca usada

Arquivo combinado de isca: `{bait_fasta}`

{md_table([bait_row], ["file", "sequences", "total_bases", "min_len", "mean_len", "median_len", "max_len", "n50"])}

## 3. Reads originais

{md_table(read_rows, ["file", "reads_counted", "total_bases", "min_len", "mean_len", "median_len", "max_len", "n50", "exact"])}

## 4. Extração das reads mitocondriais

{md_table([extract_row], ["mode", "mapped_read_names", "input_reads_r1_or_single", "input_reads_r2", "extracted_reads_r1_or_single", "extracted_reads_r2", "extracted_fraction_r1_or_single", "extracted_fraction_r2", "extracted_bases_r1_or_single", "extracted_bases_r2"])}

Arquivos extraídos:

- R1/single: `{paths.get('extracted_r1_or_single')}`
- R2: `{paths.get('extracted_r2')}`
- Lista de nomes mapeados: `{paths.get('mapped_read_names')}`

## 5. Estatísticas da remontagem enriquecida

{asm_block}

Assembly enriquecida: `{assembly}`

## 6. Circularidade da remontagem

{md_table(circ_rows, ["contig", "length", "best_terminal_overlap", "best_terminal_identity_pct", "circular_by_terminal_overlap", "circular_by_assembler"])}

Tabela completa: `{paths.get('circularity_tsv')}`

## 7. BLAST/TBLASTN contra CDS mitocondriais

{md_table(blast_rows, ["contig", "length", "nt_genes", "prot_genes", "blastn_hits", "tblastn_hits", "total_bitscore", "best_nt_pident", "best_prot_pident"])}

BLASTN: `{paths.get('blastn')}`  
TBLASTN: `{paths.get('tblastn')}`

## 8. Cobertura com reads originais contra a remontagem

{md_table(cov_rows, ["contig", "length_from_depth", "mean_depth", "median_depth", "min_depth", "max_depth", "breadth_1x_pct", "breadth_5x_pct", "breadth_10x_pct", "breadth_20x_pct"])}

BAM: `{paths.get('coverage_bam')}`  
Cobertura por base: `{paths.get('depth_tsv')}`  
Resumo `samtools coverage`: `{paths.get('coverage_tsv')}`

## 9. Interpretação automática

{chr(10).join(interpretation)}

## 10. Como decidir se virou mitogenoma final

Um candidato forte deve combinar:

1. um contig/scaffold principal próximo do tamanho esperado do mitogenoma;
2. muitos genes mitocondriais detectados por BLASTN/TBLASTN;
3. cobertura média alta e breadth próximo de 100%;
4. ausência de quedas longas de cobertura para 0x;
5. evidência de circularidade por overlap terminal, pelo grafo de montagem ou por reads atravessando a junção;
6. anotação final coerente com genes mitocondriais esperados.

Se o genoma é biologicamente circular, mas o FASTA aparece linear, isso é normal: a maioria dos assemblers escreve uma sequência circular como uma sequência linear com um ponto arbitrário de quebra. O importante é demonstrar fechamento da junção e cobertura contínua.
"""
    write_text(report, textwrap.dedent(text).strip() + "\n")
    return report


# =============================================================================
# CLI e main
# =============================================================================


def parse_args(argv: Optional[Sequence[str]] = None):
    parser = argparse.ArgumentParser(
        prog="mito_extract_reassemble.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Extrai reads mitocondriais por mapeamento contra iscas e remonta o subconjunto enriquecido.",
    )
    parser.add_argument("--mode", choices=["illumina", "pacbio"], required=True, help="Tipo de reads a processar.")
    parser.add_argument("--bait-fasta", type=Path, action="append", help="FASTA de isca. Pode repetir. Ex.: mitochondrial_candidates.fasta.")
    parser.add_argument("--nt-ref", type=Path, help="CDS mitocondriais nucleotídicos. Usado como isca opcional e para BLAST.")
    parser.add_argument("--prot-ref", type=Path, help="Proteínas mitocondriais. Usado para TBLASTN na avaliação.")
    parser.add_argument("--include-nt-ref-as-bait", action="store_true", help="Inclui os CDS nucleotídicos no FASTA de isca além dos contigs candidatos.")
    parser.add_argument("--outdir", type=Path, required=True, help="Diretório de saída.")
    parser.add_argument("--threads", type=int, default=8, help="Threads.")
    parser.add_argument("--dry-run", action="store_true", help="Apenas registra comandos, sem executar.")

    # Illumina
    parser.add_argument("--r1", type=Path, help="FASTQ R1 Illumina.")
    parser.add_argument("--r2", type=Path, help="FASTQ R2 Illumina.")
    parser.add_argument("--memory-gb", type=int, default=0, help="Memória para SPAdes em GB. 0 = padrão do SPAdes.")
    parser.add_argument("--spades-k", default="", help="K-mers para SPAdes, ex.: 21,33,55,77,99,127.")
    parser.add_argument("--spades-output", choices=["scaffolds", "contigs"], default="scaffolds", help="Arquivo SPAdes a usar como assembly final.")
    parser.add_argument("--use-bait-as-trusted-contigs", action="store_true", help="Usa a isca como trusted-contigs no SPAdes. Use com cautela, pois pode enviesar pela referência.")

    # PacBio
    parser.add_argument("--pacbio", type=Path, help="FASTQ PacBio.")
    parser.add_argument("--pacbio-type", choices=["raw", "hifi", "corrected"], default="raw", help="Tipo de PacBio para minimap2/Flye.")
    parser.add_argument("--flye-genome-size", default="80k", help="Tamanho esperado do mitogenoma para Flye.")
    parser.add_argument("--flye-meta", action="store_true", help="Usa modo --meta do Flye.")

    # Extração
    parser.add_argument("--min-mapq", type=int, default=20, help="MAPQ mínimo para considerar uma read como mitocondrial.")
    parser.add_argument("--stats-max-reads", type=int, default=0, help="Máximo de reads para estatísticas. 0 = arquivo inteiro.")

    # Avaliação BLAST
    parser.add_argument("--evaluate-with-blast", action="store_true", help="Roda BLASTN/TBLASTN contra a remontagem usando --nt-ref e --prot-ref.")
    parser.add_argument("--evalue", type=float, default=1e-5, help="E-value máximo.")
    parser.add_argument("--max-target-seqs", type=int, default=20, help="Máximo de alvos por query.")
    parser.add_argument("--min-bitscore", type=float, default=50.0, help="Bitscore mínimo por hit.")
    parser.add_argument("--min-hit-qcov", type=float, default=0.30, help="Cobertura mínima da query por hit.")
    parser.add_argument("--min-nt-pident", type=float, default=70.0, help="Identidade mínima para BLASTN.")
    parser.add_argument("--min-prot-pident", type=float, default=30.0, help="Identidade mínima para TBLASTN.")

    # Circularidade
    parser.add_argument("--min-circular-overlap", type=int, default=500, help="Overlap terminal mínimo para chamada simples de circularidade.")
    parser.add_argument("--max-circular-overlap", type=int, default=5000, help="Overlap terminal máximo para busca de circularidade.")
    parser.add_argument("--min-circular-identity", type=float, default=0.98, help="Identidade mínima do overlap terminal, de 0 a 1.")

    args = parser.parse_args(argv)

    if args.mode == "illumina":
        if not args.r1 or not args.r2:
            parser.error("--mode illumina requer --r1 e --r2.")
    else:
        if not args.pacbio:
            parser.error("--mode pacbio requer --pacbio.")

    if not args.bait_fasta and not (args.nt_ref and args.include_nt_ref_as_bait):
        parser.error("Forneça --bait-fasta e/ou use --nt-ref com --include-nt-ref-as-bait.")

    if args.evaluate_with_blast and (not args.nt_ref or not args.prot_ref):
        parser.error("--evaluate-with-blast requer --nt-ref e --prot-ref.")

    paths = []
    if args.bait_fasta:
        paths.extend(args.bait_fasta)
    for p in [args.nt_ref, args.prot_ref, args.r1, args.r2, args.pacbio]:
        if p:
            paths.append(p)
    for p in paths:
        if not p.exists():
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
        check_dependencies(args)
        with log_file.open("a", encoding="utf-8") as log:
            log.write("\n# Argumentos\n")
            log.write(json.dumps({k: str(v) for k, v in vars(args).items()}, indent=2, ensure_ascii=False) + "\n")

        print(f"[{now()}] Combinando iscas mitocondriais...")
        bait_fasta = combine_baits(args)
        args._bait_fasta_for_stats = bait_fasta
        bait_stats = fasta_stats(bait_fasta)

        print(f"[{now()}] Calculando estatísticas das reads originais...")
        input_read_stats: List[FastqStats] = []
        if args.mode == "illumina":
            input_read_stats.append(fastq_stats(args.r1, args.stats_max_reads))
            input_read_stats.append(fastq_stats(args.r2, args.stats_max_reads))
        else:
            input_read_stats.append(fastq_stats(args.pacbio, args.stats_max_reads))

        print(f"[{now()}] Mapeando reads originais contra as iscas...")
        if args.mode == "illumina":
            bait_bam = map_illumina_to_bait(args, bait_fasta, log_file)
        else:
            bait_bam = map_pacbio_to_bait(args, bait_fasta, log_file)

        print(f"[{now()}] Coletando nomes das reads mitocondriais candidatas...")
        mapped_names_path = args.outdir / "02_extracted_mito_reads" / "mapped_read_names.txt"
        if args.dry_run:
            mapped_names: Set[str] = set()
        else:
            mapped_names = mapped_read_names_from_bam(bait_bam, args.min_mapq, log_file, mapped_names_path)

        print(f"[{now()}] Extraindo reads enriquecidas...")
        extracted_r2: Optional[Path] = None
        if args.dry_run:
            if args.mode == "illumina":
                extracted_r1_or_single = args.outdir / "02_extracted_mito_reads" / "mito_enriched_R1.fastq.gz"
                extracted_r2 = args.outdir / "02_extracted_mito_reads" / "mito_enriched_R2.fastq.gz"
                extraction_stats = ExtractionStats("illumina", str(bait_fasta), 0, 0, 0, 0, 0, 0.0, 0.0, 0, 0)
            else:
                extracted_r1_or_single = args.outdir / "02_extracted_mito_reads" / "mito_enriched_pacbio.fastq.gz"
                extraction_stats = ExtractionStats("pacbio", str(bait_fasta), 0, 0, 0, 0, 0, 0.0, 0.0, 0, 0)
        elif args.mode == "illumina":
            extracted_r1_or_single, extracted_r2, extraction_stats = extract_illumina_pairs(args, mapped_names)
        else:
            extracted_r1_or_single, extraction_stats = extract_pacbio_reads(args, mapped_names)
        write_extraction_stats(extraction_stats, args.outdir / "02_extracted_mito_reads" / "extraction_summary.tsv")

        if not args.dry_run and extraction_stats.extracted_reads_r1_or_single == 0:
            raise PipelineError(
                "Nenhuma read foi extraída. Tente reduzir --min-mapq, incluir --include-nt-ref-as-bait, ou conferir se o FASTA de isca está correto."
            )

        print(f"[{now()}] Remontando reads enriquecidas...")
        if args.mode == "illumina":
            if extracted_r2 is None:
                raise PipelineError("R2 extraído não definido para modo Illumina.")
            assembly = run_spades_reassembly(args, extracted_r1_or_single, extracted_r2, log_file)
        else:
            assembly = run_flye_reassembly(args, extracted_r1_or_single, log_file)

        assembly_stats_obj = None
        circularity_calls: List[CircularityCall] = []
        blastn_out = None
        tblastn_out = None
        blast_candidates: List[BlastCandidate] = []
        cov_bam = None
        coverage_tsv = None
        depth_tsv = None
        coverage_summaries: List[CoverageSummary] = []

        if not args.dry_run and assembly.exists():
            print(f"[{now()}] Avaliando assembly enriquecida...")
            assembly_stats_obj = fasta_stats(assembly)
            circularity_calls = evaluate_circularity(args, assembly)
            circularity_tsv = args.outdir / "04_evaluate_reassembly" / "circularity_summary.tsv"
            write_circularity_tsv(circularity_calls, circularity_tsv)
            blastn_out, tblastn_out, blast_candidates = run_blast_evaluation(args, assembly, log_file)
            cov_bam, coverage_tsv, depth_tsv, coverage_summaries = map_original_reads_to_reassembly(args, assembly, log_file)
        else:
            circularity_tsv = args.outdir / "04_evaluate_reassembly" / "circularity_summary.tsv"

        paths = {
            "bait_bam": bait_bam,
            "mapped_read_names": mapped_names_path,
            "extracted_r1_or_single": extracted_r1_or_single,
            "extracted_r2": extracted_r2,
            "circularity_tsv": circularity_tsv,
            "blastn": blastn_out,
            "tblastn": tblastn_out,
            "coverage_bam": cov_bam,
            "coverage_tsv": coverage_tsv,
            "depth_tsv": depth_tsv,
        }
        report = write_report(
            args,
            bait_fasta,
            bait_stats,
            input_read_stats,
            extraction_stats,
            assembly,
            assembly_stats_obj,
            circularity_calls,
            blast_candidates,
            coverage_summaries,
            paths,
        )

        print("\nPipeline de enriquecimento/remontagem finalizado.")
        print(f"Relatório: {report}")
        print(f"Reads extraídas: {extracted_r1_or_single}" + (f" e {extracted_r2}" if extracted_r2 else ""))
        print(f"Assembly enriquecida: {assembly}")
        return 0

    except PipelineError as exc:
        sys.stderr.write(f"\nERRO: {exc}\n")
        sys.stderr.write(f"Consulte também o log: {log_file}\n")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

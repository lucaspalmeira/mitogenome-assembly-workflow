#!/usr/bin/env python3
"""
curate.py

Etapa de curadoria mitogenômica inspirada nas funções centrais do MitoCurator.

Ela NÃO substitui uma anotação manual ou ferramentas dedicadas como MITOS, MFannot
ou MitoFinder, mas acrescenta ao pipeline:

1. leitura de FASTA ou GenBank;
2. diagnóstico gene-a-gene quando há GenBank anotado;
3. detecção de CDS fora de frame e stops internos usando código genético definido;
4. busca de genes mitocondriais por BLASTN/TBLASTN contra CDS de referência;
5. checagem de genes esperados, ausentes e duplicados;
6. regiões intergênicas e candidatas AT-rich quando há GenBank;
7. rotação de sequência circular para iniciar em um gene definido pelo usuário;
8. relatório Markdown + tabelas TSV.

Autor: gerado para Lucas Palmeira
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import shutil
import statistics
import subprocess
import sys
import textwrap
import time
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:
    from .common import write_fasta
except ImportError:
    from common import write_fasta

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except Exception:  # pragma: no cover
    SeqIO = None
    Seq = None


# =============================================================================
# Estruturas de dados
# =============================================================================


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
class BlastGeneHit:
    source: str
    query_id: str
    gene: str
    contig: str
    pident: float
    aln_len: int
    query_len: int
    query_cov_pct: float
    qstart: int
    qend: int
    sstart: int
    send: int
    strand: str
    evalue: str
    bitscore: float
    status: str


@dataclass
class ExpectedGeneCall:
    gene: str
    nt_hits: int
    prot_hits: int
    best_pident: float
    best_qcov_pct: float
    best_bitscore: float
    best_contig: str
    status: str


@dataclass
class FeatureQC:
    seqid: str
    feature_type: str
    gene: str
    gene_normalized: str
    product: str
    start: int
    end: int
    strand: str
    length_nt: int
    length_aa: str
    multiple_of_three: str
    internal_stop_count: str
    internal_stop_positions: str
    terminal_stop: str
    status: str
    decision_hint: str
    comment: str


@dataclass
class IntergenicRegion:
    seqid: str
    start: int
    end: int
    length: int
    AT_percent: float
    upstream_feature: str
    downstream_feature: str
    circular_gap: str


@dataclass
class CircularityCall:
    contig: str
    length: int
    best_terminal_overlap: int
    best_terminal_identity_pct: float
    circular_by_terminal_overlap: bool


class PipelineError(RuntimeError):
    pass


# =============================================================================
# Utilitários gerais
# =============================================================================


def now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_text(path: Path, text: str) -> None:
    ensure_dir(path.parent)
    path.write_text(text, encoding="utf-8")


def open_text_auto(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def which_or_fail(program: str) -> None:
    if shutil.which(program) is None:
        raise PipelineError(
            f"Dependência não encontrada no PATH: {program}. "
            "Instale via mamba/conda ou carregue o módulo correspondente."
        )


def which(program: str) -> Optional[str]:
    return shutil.which(program)


def run_cmd(
    cmd: Sequence[str],
    log_file: Path,
    cwd: Optional[Path] = None,
    stdout_file: Optional[Path] = None,
    stderr_file: Optional[Path] = None,
    dry_run: bool = False,
) -> None:
    cmd_str = " ".join(str(x) for x in cmd)
    ensure_dir(log_file.parent)
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


def write_tsv(path: Path, rows: List[dict], fieldnames: Optional[List[str]] = None) -> None:
    ensure_dir(path.parent)
    if not fieldnames:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def markdown_table(rows: List[dict], max_rows: int = 30) -> str:
    if not rows:
        return "_Sem registros._\n"
    fieldnames = list(rows[0].keys())
    out = []
    out.append("| " + " | ".join(fieldnames) + " |")
    out.append("| " + " | ".join(["---"] * len(fieldnames)) + " |")
    for row in rows[:max_rows]:
        out.append("| " + " | ".join(str(row.get(k, "")) for k in fieldnames) + " |")
    if len(rows) > max_rows:
        out.append(f"\n_Mostrando {max_rows} de {len(rows)} linhas._")
    return "\n".join(out) + "\n"


# =============================================================================
# FASTA / GenBank
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
    return {h.split()[0]: s for h, s in iter_fasta(path)}


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


def detect_input_format(path: Path, explicit: str = "auto") -> str:
    if explicit != "auto":
        return explicit
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith((".gb", ".gbk", ".genbank")):
        return "genbank"
    with open_text_auto(path) as handle:
        first = handle.read(200)
    if first.lstrip().startswith(">"):
        return "fasta"
    if first.startswith("LOCUS") or "FEATURES" in first:
        return "genbank"
    raise PipelineError(f"Não consegui inferir o formato do arquivo: {path}. Use --input-format fasta ou genbank.")


def require_biopython() -> None:
    if SeqIO is None or Seq is None:
        raise PipelineError("Biopython não está instalado. Instale com: mamba install -c conda-forge biopython")


def genbank_to_fasta(genbank: Path, out_fasta: Path) -> Path:
    require_biopython()
    records = []
    for rec in SeqIO.parse(str(genbank), "genbank"):
        records.append((rec.id, str(rec.seq)))
    if not records:
        raise PipelineError(f"Nenhum registro GenBank lido em {genbank}")
    write_fasta(records, out_fasta)
    return out_fasta


# =============================================================================
# Normalização de nomes de genes
# =============================================================================


GENE_SYNONYMS = {
    "COI": "COX1",
    "COXI": "COX1",
    "CO1": "COX1",
    "COX I": "COX1",
    "COX1": "COX1",
    "COII": "COX2",
    "CO2": "COX2",
    "COXII": "COX2",
    "COX2": "COX2",
    "COIII": "COX3",
    "CO3": "COX3",
    "COXIII": "COX3",
    "COX3": "COX3",
    "COB": "CYTB",
    "CYTB": "CYTB",
    "CYTOCHROMEB": "CYTB",
    "CYTOCHROMEB": "CYTB",
    "ATPASE6": "ATP6",
    "ATPASE8": "ATP8",
    "ATPASE9": "ATP9",
    "ATP6": "ATP6",
    "ATP8": "ATP8",
    "ATP9": "ATP9",
    "NAD1": "ND1",
    "NAD2": "ND2",
    "NAD3": "ND3",
    "NAD4": "ND4",
    "NAD4L": "ND4L",
    "NAD5": "ND5",
    "NAD6": "ND6",
    "ND1": "ND1",
    "ND2": "ND2",
    "ND3": "ND3",
    "ND4": "ND4",
    "ND4L": "ND4L",
    "ND5": "ND5",
    "ND6": "ND6",
    "RNL": "RNL",
    "RRNL": "RNL",
    "LSURRNA": "RNL",
    "RRNLARGE": "RNL",
    "RNS": "RNS",
    "RRNS": "RNS",
    "SSURRNA": "RNS",
    "RRNSMALL": "RNS",
    "RPS3": "RPS3",
    "VAR1": "RPS3",
}

KNOWN_GENE_PATTERNS = [
    "cox1", "cox2", "cox3", "coi", "coii", "coiii",
    "cob", "cytb",
    "atp6", "atp8", "atp9", "atpase6", "atpase8", "atpase9",
    "nad1", "nad2", "nad3", "nad4l", "nad4", "nad5", "nad6",
    "nd1", "nd2", "nd3", "nd4l", "nd4", "nd5", "nd6",
    "rnl", "rns", "rrnl", "rrns", "lsu", "ssu", "rps3", "var1",
]


def normalize_gene_name(value: str) -> str:
    if value is None:
        return "."
    text = str(value).strip()
    if not text or text == ".":
        return "."

    # Remove qualificadores comuns.
    text = re.sub(r"^(gene|product|locus_tag|protein_id|Name)[:=]", "", text, flags=re.I).strip()

    # Tenta capturar genes conhecidos dentro de descrições longas.
    lowered = text.lower()
    for pat in sorted(KNOWN_GENE_PATTERNS, key=len, reverse=True):
        if re.search(rf"(^|[^a-z0-9]){re.escape(pat)}([^a-z0-9]|$)", lowered):
            key = re.sub(r"[^A-Z0-9]", "", pat.upper())
            return GENE_SYNONYMS.get(key, key)

    # tRNAs: trnK, tRNA-Lys, trnL1 etc.
    m = re.search(r"\btrn([A-Za-z][0-9]?)\b", text, flags=re.I)
    if m:
        return "TRN" + m.group(1).upper()
    m = re.search(r"tRNA[-_ ]?([A-Za-z]{3}|[A-Za-z])", text, flags=re.I)
    if m:
        aa = m.group(1).upper()
        aa3_to_1 = {
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
            "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
            "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        }
        aa1 = aa3_to_1.get(aa, aa[0])
        return "TRN" + aa1

    key = re.sub(r"[^A-Za-z0-9]", "", text).upper()
    return GENE_SYNONYMS.get(key, key if key else ".")


def gene_name_from_header(header: str) -> str:
    candidates = []
    h = header.strip()
    for pattern in [
        r"\bgene[=:]([^\s;|,]+)",
        r"\bproduct[=:]([^;|]+)",
        r"\[gene=([^\]]+)\]",
        r"\[product=([^\]]+)\]",
    ]:
        m = re.search(pattern, h, flags=re.I)
        if m:
            candidates.append(m.group(1))
    # Também tenta o primeiro token e a descrição completa.
    candidates.append(h.split()[0])
    candidates.append(h)
    for c in candidates:
        norm = normalize_gene_name(c)
        if norm != ".":
            return norm
    return normalize_gene_name(h)


def build_query_gene_map(paths: Sequence[Path]) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for path in paths:
        if not path:
            continue
        for header, _seq in iter_fasta(path):
            qid = header.split()[0]
            mapping[qid] = gene_name_from_header(header)
    return mapping


def expected_genes_from_reference(paths: Sequence[Path]) -> List[str]:
    genes = []
    for path in paths:
        if not path:
            continue
        for header, _ in iter_fasta(path):
            genes.append(gene_name_from_header(header))
    return sorted(set(g for g in genes if g and g != "."))


def expected_genes_from_file(path: Path) -> List[str]:
    genes = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # aceita TSV/CSV/lista simples
            token = re.split(r"[\t,; ]+", line)[0]
            genes.append(normalize_gene_name(token))
    return sorted(set(g for g in genes if g and g != "."))


def expected_genes_from_preset(preset: str) -> List[str]:
    preset = preset.lower()
    cds13 = ["ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"]
    fungal_common = ["ATP6", "ATP8", "ATP9", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "RNL", "RNS", "RPS3"]
    if preset in {"animal", "animal_mito", "insect", "insect_mito", "metazoa"}:
        return sorted(cds13 + ["RNL", "RNS"])
    if preset in {"fungal", "fungal_mito", "fungi"}:
        return sorted(fungal_common)
    if preset in {"cds13", "standard"}:
        return sorted(cds13)
    if preset in {"none", "off"}:
        return []
    raise PipelineError(f"Preset de genes esperados não reconhecido: {preset}")


# =============================================================================
# BLAST/TBLASTN e presença gênica
# =============================================================================


def make_blast_db(target_fasta: Path, out_prefix: Path, log_file: Path, dry_run: bool = False) -> Path:
    ensure_dir(out_prefix.parent)
    run_cmd(
        ["makeblastdb", "-in", str(target_fasta), "-dbtype", "nucl", "-out", str(out_prefix)],
        log_file,
        stderr_file=out_prefix.parent / "makeblastdb.stderr.log",
        dry_run=dry_run,
    )
    return out_prefix


def run_blast_gene_search(args, target_fasta: Path, log_file: Path) -> Tuple[List[BlastGeneHit], Path, Path]:
    blast_dir = args.outdir / "blast_gene_search"
    ensure_dir(blast_dir)
    db_prefix = blast_dir / "target_db"
    make_blast_db(target_fasta, db_prefix, log_file, dry_run=args.dry_run)

    all_hits: List[BlastGeneHit] = []
    blastn_out = blast_dir / "blastn_nt_vs_candidate.tsv"
    tblastn_out = blast_dir / "tblastn_prot_vs_candidate.tsv"

    outfmt = "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"

    if args.nt_ref:
        run_cmd(
            [
                "blastn", "-query", str(args.nt_ref), "-db", str(db_prefix),
                "-out", str(blastn_out), "-outfmt", outfmt,
                "-evalue", str(args.evalue), "-num_threads", str(args.threads),
                "-max_target_seqs", str(args.max_target_seqs),
            ],
            log_file,
            stderr_file=blast_dir / "blastn.stderr.log",
            dry_run=args.dry_run,
        )
        qmap = build_query_gene_map([args.nt_ref])
        all_hits.extend(parse_blast_table(blastn_out, "blastn", qmap, args.min_identity, args.min_query_cov))
    else:
        write_text(blastn_out, "")

    if args.prot_ref:
        run_cmd(
            [
                "tblastn", "-query", str(args.prot_ref), "-db", str(db_prefix),
                "-out", str(tblastn_out), "-outfmt", outfmt,
                "-evalue", str(args.evalue), "-num_threads", str(args.threads),
                "-max_target_seqs", str(args.max_target_seqs),
            ],
            log_file,
            stderr_file=blast_dir / "tblastn.stderr.log",
            dry_run=args.dry_run,
        )
        qmap = build_query_gene_map([args.prot_ref])
        all_hits.extend(parse_blast_table(tblastn_out, "tblastn", qmap, args.min_identity, args.min_query_cov))
    else:
        write_text(tblastn_out, "")

    return all_hits, blastn_out, tblastn_out


def parse_blast_table(path: Path, source: str, qmap: Dict[str, str], min_identity: float, min_qcov: float) -> List[BlastGeneHit]:
    hits = []
    if not path.exists() or path.stat().st_size == 0:
        return hits
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            qid, sid = parts[0], parts[1]
            pident = float(parts[2])
            aln_len = int(float(parts[3]))
            qlen = int(float(parts[4])) if parts[4] not in {"0", ""} else 0
            qstart = int(float(parts[6]))
            qend = int(float(parts[7]))
            sstart = int(float(parts[8]))
            send = int(float(parts[9]))
            evalue = parts[10]
            bitscore = float(parts[11])
            qcov = 100.0 * aln_len / qlen if qlen else 0.0
            status = "PASS" if pident >= min_identity and qcov >= min_qcov else "LOW_CONFIDENCE"
            gene = qmap.get(qid, normalize_gene_name(qid))
            hits.append(
                BlastGeneHit(
                    source=source,
                    query_id=qid,
                    gene=gene,
                    contig=sid,
                    pident=pident,
                    aln_len=aln_len,
                    query_len=qlen,
                    query_cov_pct=qcov,
                    qstart=qstart,
                    qend=qend,
                    sstart=sstart,
                    send=send,
                    strand="+" if sstart <= send else "-",
                    evalue=evalue,
                    bitscore=bitscore,
                    status=status,
                )
            )
    return hits


def summarize_expected_gene_calls(hits: List[BlastGeneHit], expected_genes: Sequence[str]) -> List[ExpectedGeneCall]:
    by_gene: Dict[str, List[BlastGeneHit]] = defaultdict(list)
    for h in hits:
        if h.status == "PASS":
            by_gene[h.gene].append(h)

    genes = sorted(set(expected_genes) | set(by_gene.keys()))
    rows: List[ExpectedGeneCall] = []
    for gene in genes:
        ghits = by_gene.get(gene, [])
        nt_hits = sum(1 for h in ghits if h.source == "blastn")
        prot_hits = sum(1 for h in ghits if h.source == "tblastn")
        if not ghits:
            rows.append(ExpectedGeneCall(gene, 0, 0, 0.0, 0.0, 0.0, ".", "MISSING"))
            continue
        best = sorted(ghits, key=lambda h: (h.bitscore, h.query_cov_pct, h.pident), reverse=True)[0]
        # Duplicado se houver hits PASS em mais de uma região/contig para o mesmo gene.
        loci = {(h.contig, min(h.sstart, h.send), max(h.sstart, h.send)) for h in ghits}
        status = "PRESENT"
        if len(loci) > 1:
            status = "MULTIPLE_HITS_CHECK_DUPLICATION_OR_REPEATS"
        rows.append(
            ExpectedGeneCall(
                gene=gene,
                nt_hits=nt_hits,
                prot_hits=prot_hits,
                best_pident=round(best.pident, 3),
                best_qcov_pct=round(best.query_cov_pct, 3),
                best_bitscore=round(best.bitscore, 3),
                best_contig=best.contig,
                status=status,
            )
        )
    return rows


# =============================================================================
# Diagnóstico de GenBank
# =============================================================================


def feature_bounds(feat) -> Tuple[int, int]:
    return int(feat.location.start) + 1, int(feat.location.end)


def strand_symbol(feat) -> str:
    return "+" if feat.location.strand == 1 else "-" if feat.location.strand == -1 else "."


def feature_name(feat) -> str:
    for key in ["gene", "locus_tag", "product", "note"]:
        if key in feat.qualifiers and feat.qualifiers[key]:
            return str(feat.qualifiers[key][0])
    return "."


def summarize_genbank_features(genbank: Path, genetic_code: int) -> Tuple[List[FeatureQC], List[IntergenicRegion], int, Dict[str, int]]:
    require_biopython()
    records = list(SeqIO.parse(str(genbank), "genbank"))
    if not records:
        raise PipelineError(f"Nenhum registro GenBank lido em {genbank}")

    all_feature_rows: List[FeatureQC] = []
    all_intergenic_rows: List[IntergenicRegion] = []
    feature_counts: Counter[str] = Counter()
    total_len = 0

    for rec in records:
        total_len += len(rec.seq)
        feature_counts.update([feat.type for feat in rec.features if feat.type != "source"])
        for feat in rec.features:
            if feat.type == "source":
                continue
            start, end = feature_bounds(feat)
            gene = feat.qualifiers.get("gene", [feature_name(feat)])[0]
            product = feat.qualifiers.get("product", ["."])[0]
            norm = normalize_gene_name(gene if gene != "." else product)
            seq = feat.extract(rec.seq)
            length_nt = len(seq)
            row = FeatureQC(
                seqid=rec.id,
                feature_type=feat.type,
                gene=str(gene),
                gene_normalized=norm,
                product=str(product),
                start=start,
                end=end,
                strand=strand_symbol(feat),
                length_nt=length_nt,
                length_aa=".",
                multiple_of_three=".",
                internal_stop_count=".",
                internal_stop_positions=".",
                terminal_stop=".",
                status="OK",
                decision_hint="OK",
                comment=".",
            )
            if feat.type == "CDS":
                nt = str(seq).upper()
                aa = str(Seq(nt).translate(table=genetic_code, to_stop=False))
                stops = []
                for i, a in enumerate(aa, start=1):
                    if a == "*":
                        status = "terminal" if i == len(aa) else "internal"
                        codon = nt[(i - 1) * 3 : i * 3]
                        stops.append((i, codon, status))
                internal = [s for s in stops if s[2] == "internal"]
                terminal = any(s[2] == "terminal" for s in stops)
                row.length_aa = str(len(aa))
                row.multiple_of_three = "yes" if length_nt % 3 == 0 else "no"
                row.internal_stop_count = str(len(internal))
                row.internal_stop_positions = ",".join(str(s[0]) for s in internal) if internal else "."
                row.terminal_stop = "yes" if terminal else "no"
                if length_nt % 3 != 0 and internal:
                    row.status = "PROBLEM"
                    row.decision_hint = "CHECK_FRAMESHIFT_AND_INTERNAL_STOP"
                elif length_nt % 3 != 0:
                    row.status = "PROBLEM"
                    row.decision_hint = "CHECK_LENGTH_NOT_MULTIPLE_OF_THREE"
                elif internal:
                    row.status = "PROBLEM"
                    row.decision_hint = "CHECK_INTERNAL_STOP"
                elif not terminal:
                    row.status = "WARNING"
                    row.decision_hint = "NO_TERMINAL_STOP_OR_PARTIAL_CDS"
            all_feature_rows.append(row)

        all_intergenic_rows.extend(find_intergenic_regions(rec))

    return all_feature_rows, all_intergenic_rows, total_len, dict(feature_counts)


def find_intergenic_regions(record) -> List[IntergenicRegion]:
    intervals = []
    for feat in record.features:
        if feat.type == "source":
            continue
        start, end = feature_bounds(feat)
        if start <= end:
            intervals.append((start, end, feat.type, feature_name(feat)))
    intervals.sort(key=lambda x: (x[0], x[1]))
    rows: List[IntergenicRegion] = []
    seq_len = len(record.seq)

    for i in range(len(intervals) - 1):
        a = intervals[i]
        b = intervals[i + 1]
        gap_start = a[1] + 1
        gap_end = b[0] - 1
        if gap_end >= gap_start:
            gap_seq = record.seq[gap_start - 1 : gap_end]
            rows.append(make_intergenic_row(record.id, gap_start, gap_end, str(gap_seq), a, b, "no"))

    # Gap circular entre última e primeira feature.
    if len(intervals) > 1:
        first = intervals[0]
        last = intervals[-1]
        if last[1] < seq_len or first[0] > 1:
            seq = str(record.seq[last[1] :]) + str(record.seq[: first[0] - 1])
            length = len(seq)
            if length > 0:
                rows.append(
                    IntergenicRegion(
                        seqid=record.id,
                        start=last[1] + 1,
                        end=first[0] - 1,
                        length=length,
                        AT_percent=round(at_percent(seq), 2),
                        upstream_feature=f"{last[2]}:{last[3]}",
                        downstream_feature=f"{first[2]}:{first[3]}",
                        circular_gap="yes",
                    )
                )
    return rows


def make_intergenic_row(seqid: str, start: int, end: int, seq: str, upstream, downstream, circular_gap: str) -> IntergenicRegion:
    return IntergenicRegion(
        seqid=seqid,
        start=start,
        end=end,
        length=len(seq),
        AT_percent=round(at_percent(seq), 2),
        upstream_feature=f"{upstream[2]}:{upstream[3]}",
        downstream_feature=f"{downstream[2]}:{downstream[3]}",
        circular_gap=circular_gap,
    )


def at_percent(seq: str) -> float:
    seq = seq.upper()
    valid = sum(1 for x in seq if x in "ACGT")
    if valid == 0:
        return 0.0
    return 100.0 * (seq.count("A") + seq.count("T")) / valid


# =============================================================================
# Circularidade e rotação
# =============================================================================


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def best_terminal_overlap(seq: str, min_overlap: int = 100, max_overlap: int = 5000) -> Tuple[int, float, bool]:
    seq = seq.upper()
    max_o = min(max_overlap, len(seq) // 2)
    best_len = 0
    best_id = 0.0
    for o in range(min_overlap, max_o + 1):
        left = seq[:o]
        right = seq[-o:]
        matches = sum(1 for a, b in zip(left, right) if a == b)
        identity = 100.0 * matches / o if o else 0.0
        if identity > best_id or (identity == best_id and o > best_len):
            best_len = o
            best_id = identity
    circular = best_len >= min_overlap and best_id >= 95.0
    return best_len, best_id, circular


def circularity_for_fasta(path: Path, min_overlap: int, max_overlap: int) -> List[CircularityCall]:
    calls = []
    for header, seq in iter_fasta(path):
        contig = header.split()[0]
        olen, ident, circ = best_terminal_overlap(seq, min_overlap=min_overlap, max_overlap=max_overlap)
        calls.append(CircularityCall(contig, len(seq), olen, round(ident, 3), circ))
    return calls


def rotate_sequence(seq: str, start_1based: int) -> str:
    if start_1based <= 1:
        return seq
    idx = start_1based - 1
    return seq[idx:] + seq[:idx]


def rotate_fasta_by_coordinate(input_fasta: Path, output_fasta: Path, contig: str, start_1based: int, rotate_to: str) -> Optional[Path]:
    records = []
    did = False
    for header, seq in iter_fasta(input_fasta):
        cid = header.split()[0]
        if cid == contig:
            rseq = rotate_sequence(seq, start_1based)
            records.append((f"{cid}_rotated_start_{rotate_to}_at_{start_1based}", rseq))
            did = True
        else:
            records.append((header, seq))
    if not did:
        return None
    write_fasta(records, output_fasta)
    return output_fasta


def choose_rotation_from_blast(hits: List[BlastGeneHit], rotate_to: str) -> Optional[BlastGeneHit]:
    target = normalize_gene_name(rotate_to)
    candidates = [h for h in hits if h.gene == target and h.status == "PASS"]
    if not candidates:
        return None
    return sorted(candidates, key=lambda h: (h.bitscore, h.query_cov_pct, h.pident), reverse=True)[0]


def choose_rotation_from_genbank(genbank: Path, rotate_to: str) -> Optional[Tuple[str, int]]:
    require_biopython()
    target = normalize_gene_name(rotate_to)
    for rec in SeqIO.parse(str(genbank), "genbank"):
        for feat in rec.features:
            if feat.type == "source":
                continue
            gene = feat.qualifiers.get("gene", [feature_name(feat)])[0]
            product = feat.qualifiers.get("product", ["."])[0]
            norm = normalize_gene_name(gene if gene != "." else product)
            if norm == target:
                start = int(feat.location.start) + 1
                return rec.id, start
    return None


# =============================================================================
# Relatório
# =============================================================================


def write_curation_report(
    args,
    input_format: str,
    target_fasta: Path,
    stats: FastaStats,
    expected_rows: List[ExpectedGeneCall],
    blast_hits: List[BlastGeneHit],
    feature_rows: List[FeatureQC],
    intergenic_rows: List[IntergenicRegion],
    at_rich_rows: List[IntergenicRegion],
    circularity_rows: List[CircularityCall],
    rotated_fasta: Optional[Path],
    feature_counts: Dict[str, int],
    blastn_out: Optional[Path],
    tblastn_out: Optional[Path],
) -> Path:
    report = args.outdir / "CURATION_REPORT.md"
    missing = [r for r in expected_rows if r.status == "MISSING"]
    duplicated = [r for r in expected_rows if r.status == "MULTIPLE_HITS_CHECK_DUPLICATION_OR_REPEATS"]
    present = [r for r in expected_rows if r.status != "MISSING"]
    problems = [r for r in feature_rows if r.status == "PROBLEM"]
    warnings = [r for r in feature_rows if r.status == "WARNING"]

    best_circ = None
    if circularity_rows:
        best_circ = sorted(circularity_rows, key=lambda r: (r.circular_by_terminal_overlap, r.best_terminal_identity_pct, r.best_terminal_overlap), reverse=True)[0]

    text = f"""# Relatório de curadoria mitogenômica

**Data:** {now()}  
**Entrada:** `{args.candidate}`  
**Formato da entrada:** `{input_format}`  
**Diretório de saída:** `{args.outdir}`  
**Topologia informada:** `{args.topology}`  
**Código genético usado para CDS:** `{args.genetic_code}`

## 1. Objetivo

Esta etapa acrescenta ao pipeline funções de curadoria semelhantes às do MitoCurator: leitura de FASTA/GenBank, diagnóstico gene-a-gene quando há anotação, detecção de stops internos em CDS, checagem de genes esperados por BLASTN/TBLASTN, avaliação simples de circularidade terminal, busca de regiões intergênicas/AT-rich em GenBank e rotação de sequência circular para iniciar em um gene escolhido.

## 2. Estatísticas da sequência candidata

{markdown_table([asdict(stats)])}

## 3. Circularidade por sobreposição terminal

{markdown_table([asdict(r) for r in circularity_rows])}
"""

    if best_circ:
        if best_circ.circular_by_terminal_overlap:
            text += f"\n**Interpretação:** há evidência de sobreposição terminal forte em `{best_circ.contig}`: {best_circ.best_terminal_overlap} bp com {best_circ.best_terminal_identity_pct}% de identidade. Isso sugere uma sequência circular redundante nas pontas, mas ainda deve ser validado por reads atravessando a junção/grafo.\n"
        else:
            text += f"\n**Interpretação:** não foi detectada sobreposição terminal forte. O melhor caso foi `{best_circ.contig}`, com {best_circ.best_terminal_overlap} bp e {best_circ.best_terminal_identity_pct}% de identidade. Isso não descarta circularidade biológica; apenas indica que o FASTA não traz uma redundância terminal clara.\n"

    text += "\n## 4. Presença de genes esperados por BLASTN/TBLASTN\n\n"
    if expected_rows:
        text += f"- Genes avaliados: **{len(expected_rows)}**\n"
        text += f"- Presentes ou com hit: **{len(present)}**\n"
        text += f"- Ausentes: **{len(missing)}**\n"
        text += f"- Com múltiplos hits/repeats/possível duplicação: **{len(duplicated)}**\n\n"
        text += markdown_table([asdict(r) for r in expected_rows], max_rows=args.report_max_rows)
    else:
        text += "Nenhuma lista de genes esperados foi definida. Use `--expected-from-ref`, `--expected-genes-file` ou `--expected-preset`.\n"

    if blastn_out or tblastn_out:
        text += "\nArquivos BLAST:\n\n"
        if blastn_out:
            text += f"- BLASTN: `{blastn_out}`\n"
        if tblastn_out:
            text += f"- TBLASTN: `{tblastn_out}`\n"

    text += "\n## 5. Diagnóstico de anotação GenBank\n\n"
    if feature_rows:
        text += "### Contagem de features\n\n"
        fc_rows = [{"feature_type": k, "count": v} for k, v in sorted(feature_counts.items())]
        text += markdown_table(fc_rows)
        text += f"\n### Problemas em CDS\n\n- CDS/features problemáticas: **{len(problems)}**\n- Avisos: **{len(warnings)}**\n\n"
        if problems:
            text += markdown_table([asdict(r) for r in problems], max_rows=args.report_max_rows)
        else:
            text += "Nenhum CDS problemático detectado pelas regras atuais.\n"
    else:
        text += "A entrada não era GenBank anotado, então esta etapa não avaliou features, stops internos por anotação, tRNAs/rRNAs anotados nem regiões intergênicas. Para essa parte, forneça um `.gb/.gbk` anotado por MITOS, MFannot, MitoFinder ou curadoria manual.\n"

    text += "\n## 6. Regiões intergênicas e candidatas AT-rich\n\n"
    if intergenic_rows:
        text += f"- Regiões intergênicas detectadas: **{len(intergenic_rows)}**\n"
        text += f"- Candidatas AT-rich: **{len(at_rich_rows)}** usando mínimo de {args.at_rich_min_len} bp e {args.at_rich_min_at}% A+T.\n\n"
        if at_rich_rows:
            text += markdown_table([asdict(r) for r in at_rich_rows], max_rows=args.report_max_rows)
        else:
            text += "Nenhuma região intergênica passou os filtros de AT-rich.\n"
    else:
        text += "Sem regiões intergênicas calculadas. Isso é esperado para FASTA sem anotação.\n"

    text += "\n## 7. Rotação da sequência circular\n\n"
    if rotated_fasta:
        text += f"Sequência rotacionada gerada em: `{rotated_fasta}`\n"
    elif args.rotate_to:
        text += f"Foi solicitada rotação para `{args.rotate_to}`, mas o gene não foi encontrado com confiança suficiente para definir a coordenada de início.\n"
    else:
        text += "Rotação não solicitada. Use `--rotate-to COX1` ou outro gene para gerar um FASTA rotacionado.\n"

    text += f"""
## 8. Arquivos gerados

- Relatório: `{report}`
- Hits BLAST filtrados: `{args.outdir / 'blast_gene_search' / 'gene_hits.filtered.tsv'}`
- Presença/ausência de genes: `{args.outdir / 'gene_presence_absence.tsv'}`
- Circularidade terminal: `{args.outdir / 'circularity_summary.tsv'}`
- Diagnóstico de features GenBank: `{args.outdir / 'gene_qc.tsv'}`
- Features problemáticas: `{args.outdir / 'problematic_features.tsv'}`
- Regiões intergênicas: `{args.outdir / 'intergenic_regions.tsv'}`
- Candidatas AT-rich: `{args.outdir / 'at_rich_candidates.tsv'}`
- Log dos comandos: `{args.outdir / 'commands.log'}`

## 9. Como interpretar

Um mitogenoma final forte deve combinar: contig único próximo ao tamanho esperado, muitos genes mitocondriais presentes, poucos ou nenhum CDS com stop interno, cobertura contínua, ausência de quebras longas com 0x, anotação coerente e evidência independente de circularidade por reads atravessando a junção, grafo ou assembler. A rotação para um gene inicial, como COX1/COX1-like, padroniza a representação do genoma circular, mas não prova circularidade por si só.
"""
    write_text(report, text)
    return report


# =============================================================================
# Pipeline principal
# =============================================================================


def check_tools(args, input_format: str) -> None:
    # BLAST só é obrigatório se houver referência para busca gênica.
    if args.nt_ref or args.prot_ref:
        for p in ["makeblastdb"]:
            which_or_fail(p)
        if args.nt_ref:
            which_or_fail("blastn")
        if args.prot_ref:
            which_or_fail("tblastn")
    if input_format == "genbank":
        require_biopython()


def run_curation(args) -> Path:
    ensure_dir(args.outdir)
    log_file = args.outdir / "commands.log"
    write_text(log_file, f"[{now()}] Início da curadoria\n")

    input_format = detect_input_format(args.candidate, args.input_format)
    check_tools(args, input_format)

    work_dir = args.outdir / "00_input"
    ensure_dir(work_dir)
    if input_format == "genbank":
        target_fasta = genbank_to_fasta(args.candidate, work_dir / "candidate_from_genbank.fasta")
    else:
        target_fasta = args.candidate

    stats = fasta_stats(target_fasta)

    # Circularidade por sobreposição terminal.
    circularity = circularity_for_fasta(target_fasta, args.min_terminal_overlap, args.max_terminal_overlap)
    write_tsv(args.outdir / "circularity_summary.tsv", [asdict(r) for r in circularity])

    # Diagnóstico GenBank.
    feature_rows: List[FeatureQC] = []
    intergenic_rows: List[IntergenicRegion] = []
    feature_counts: Dict[str, int] = {}
    if input_format == "genbank":
        feature_rows, intergenic_rows, _total_len, feature_counts = summarize_genbank_features(args.candidate, args.genetic_code)
        write_tsv(args.outdir / "gene_qc.tsv", [asdict(r) for r in feature_rows])
        problem_rows = [r for r in feature_rows if r.status != "OK"]
        write_tsv(args.outdir / "problematic_features.tsv", [asdict(r) for r in problem_rows], fieldnames=list(asdict(feature_rows[0]).keys()) if feature_rows else [])
        write_tsv(args.outdir / "intergenic_regions.tsv", [asdict(r) for r in intergenic_rows])
    else:
        write_tsv(args.outdir / "gene_qc.tsv", [])
        write_tsv(args.outdir / "problematic_features.tsv", [])
        write_tsv(args.outdir / "intergenic_regions.tsv", [])

    at_rich_rows = [
        r for r in intergenic_rows
        if r.length >= args.at_rich_min_len and r.AT_percent >= args.at_rich_min_at
    ]
    at_rich_rows.sort(key=lambda r: (r.length, r.AT_percent), reverse=True)
    write_tsv(args.outdir / "at_rich_candidates.tsv", [asdict(r) for r in at_rich_rows])

    # BLAST gene search.
    hits: List[BlastGeneHit] = []
    blastn_out: Optional[Path] = None
    tblastn_out: Optional[Path] = None
    if args.nt_ref or args.prot_ref:
        hits, blastn_out, tblastn_out = run_blast_gene_search(args, target_fasta, log_file)
        pass_hits = [asdict(h) for h in hits if h.status == "PASS"]
        write_tsv(args.outdir / "blast_gene_search" / "gene_hits.filtered.tsv", pass_hits)
        write_tsv(args.outdir / "blast_gene_search" / "gene_hits.all.tsv", [asdict(h) for h in hits])
    else:
        ensure_dir(args.outdir / "blast_gene_search")
        write_tsv(args.outdir / "blast_gene_search" / "gene_hits.filtered.tsv", [])
        write_tsv(args.outdir / "blast_gene_search" / "gene_hits.all.tsv", [])

    # Genes esperados.
    expected: List[str] = []
    if args.expected_genes_file:
        expected.extend(expected_genes_from_file(args.expected_genes_file))
    if args.expected_preset and args.expected_preset.lower() not in {"none", "off"}:
        expected.extend(expected_genes_from_preset(args.expected_preset))
    if args.expected_from_ref:
        expected.extend(expected_genes_from_reference([p for p in [args.nt_ref, args.prot_ref] if p]))
    if input_format == "genbank" and args.expected_from_genbank:
        expected.extend([r.gene_normalized for r in feature_rows if r.gene_normalized != "."])
    expected = sorted(set(g for g in expected if g and g != "."))
    expected_rows = summarize_expected_gene_calls(hits, expected) if expected or hits else []
    write_tsv(args.outdir / "gene_presence_absence.tsv", [asdict(r) for r in expected_rows])
    missing = [asdict(r) for r in expected_rows if r.status == "MISSING"]
    duplicated = [asdict(r) for r in expected_rows if r.status == "MULTIPLE_HITS_CHECK_DUPLICATION_OR_REPEATS"]
    write_tsv(args.outdir / "missing_genes.tsv", missing)
    write_tsv(args.outdir / "duplicate_or_repeat_gene_hits.tsv", duplicated)

    # Rotação.
    rotated_fasta: Optional[Path] = None
    if args.rotate_to and args.topology == "circular":
        rotation_hit: Optional[Tuple[str, int]] = None
        if input_format == "genbank":
            rotation_hit = choose_rotation_from_genbank(args.candidate, args.rotate_to)
        if rotation_hit is None and hits:
            h = choose_rotation_from_blast(hits, args.rotate_to)
            if h is not None:
                rotation_hit = (h.contig, min(h.sstart, h.send))
        if rotation_hit is not None:
            contig, start = rotation_hit
            rotated_fasta = rotate_fasta_by_coordinate(
                target_fasta,
                args.outdir / f"rotated_start_{normalize_gene_name(args.rotate_to)}.fasta",
                contig,
                start,
                normalize_gene_name(args.rotate_to),
            )
            if rotated_fasta:
                with log_file.open("a", encoding="utf-8") as log:
                    log.write(f"\n[{now()}] Rotação: contig={contig}; start={start}; gene={args.rotate_to}\n")

    report = write_curation_report(
        args=args,
        input_format=input_format,
        target_fasta=target_fasta,
        stats=stats,
        expected_rows=expected_rows,
        blast_hits=hits,
        feature_rows=feature_rows,
        intergenic_rows=intergenic_rows,
        at_rich_rows=at_rich_rows,
        circularity_rows=circularity,
        rotated_fasta=rotated_fasta,
        feature_counts=feature_counts,
        blastn_out=blastn_out,
        tblastn_out=tblastn_out,
    )
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Curadoria mitogenômica: diagnóstico de genes, CDS, circularidade, rotação e relatórios.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--candidate", type=Path, required=True, help="FASTA ou GenBank do mitogenoma candidato.")
    parser.add_argument("--input-format", choices=["auto", "fasta", "genbank"], default="auto", help="Formato da entrada.")
    parser.add_argument("--outdir", type=Path, required=True, help="Diretório de saída da curadoria.")

    parser.add_argument("--nt-ref", type=Path, help="CDS nucleotídicos mitocondriais de referência, em FASTA.")
    parser.add_argument("--prot-ref", type=Path, help="Proteínas mitocondriais de referência, em FASTA.")
    parser.add_argument("--expected-from-ref", action="store_true", help="Usar nomes de genes dos arquivos --nt-ref/--prot-ref como conjunto esperado.")
    parser.add_argument("--expected-from-genbank", action="store_true", help="Usar genes anotados no GenBank como conjunto esperado adicional.")
    parser.add_argument("--expected-genes-file", type=Path, help="Arquivo com genes esperados, um por linha ou TSV simples.")
    parser.add_argument("--expected-preset", default="none", help="Preset: none, cds13, animal_mito/insect_mito, fungal_mito.")

    parser.add_argument("--topology", choices=["circular", "linear", "unknown"], default="circular", help="Topologia esperada da molécula.")
    parser.add_argument("--rotate-to", default="", help="Gene para iniciar a sequência circular, por exemplo COX1, COX2, ATP6, COB/CYTB.")
    parser.add_argument("--genetic-code", type=int, default=4, help="Código genético NCBI para traduzir CDS. Para muitos fungos mitocondriais, use 4; invertebrados, 5; vertebrados, 2.")

    parser.add_argument("--min-identity", type=float, default=70.0, help="Identidade mínima para aceitar hits BLAST como presença gênica.")
    parser.add_argument("--min-query-cov", type=float, default=50.0, help="Cobertura mínima da query para aceitar hits BLAST como presença gênica.")
    parser.add_argument("--evalue", default="1e-10", help="E-value para BLASTN/TBLASTN.")
    parser.add_argument("--max-target-seqs", type=int, default=20, help="Máximo de alvos por query no BLAST.")

    parser.add_argument("--min-terminal-overlap", type=int, default=100, help="Sobreposição terminal mínima para testar circularidade simples.")
    parser.add_argument("--max-terminal-overlap", type=int, default=5000, help="Sobreposição terminal máxima para testar circularidade simples.")
    parser.add_argument("--at-rich-min-len", type=int, default=500, help="Comprimento mínimo para chamar região intergênica AT-rich.")
    parser.add_argument("--at-rich-min-at", type=float, default=75.0, help="Percentual mínimo de A+T para região AT-rich.")

    parser.add_argument("--threads", type=int, default=8, help="Threads para BLAST.")
    parser.add_argument("--report-max-rows", type=int, default=40, help="Número máximo de linhas exibidas em tabelas do relatório.")
    parser.add_argument("--dry-run", action="store_true", help="Registra comandos, mas não executa ferramentas externas.")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        if not args.candidate.exists():
            raise PipelineError(f"Arquivo candidato não encontrado: {args.candidate}")
        if args.nt_ref and not args.nt_ref.exists():
            raise PipelineError(f"Arquivo --nt-ref não encontrado: {args.nt_ref}")
        if args.prot_ref and not args.prot_ref.exists():
            raise PipelineError(f"Arquivo --prot-ref não encontrado: {args.prot_ref}")
        if args.expected_genes_file and not args.expected_genes_file.exists():
            raise PipelineError(f"Arquivo --expected-genes-file não encontrado: {args.expected_genes_file}")
        report = run_curation(args)
        print(f"\nCuradoria concluída. Relatório: {report}")
        return 0
    except PipelineError as exc:
        print(f"\nERRO: {exc}", file=sys.stderr)
        try:
            if 'args' in locals():
                print(f"Consulte também o log: {args.outdir / 'commands.log'}", file=sys.stderr)
        except Exception:
            pass
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

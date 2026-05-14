#!/usr/bin/env python3

import argparse
import csv
import shutil
import subprocess
from pathlib import Path
from collections import defaultdict


BLAST_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
    "qcovhsp",
    "qcovs",
    "sframe",
]


def check_executable(program):
    path = shutil.which(program)
    if path is None:
        raise RuntimeError(
            f"ERROR: '{program}' não foi encontrado no PATH. "
            f"Instale o BLAST+ ou carregue o módulo correto."
        )
    return path


def run_command(cmd):
    print("\n[CMD]", " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True)


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


def reverse_complement(seq):
    table = str.maketrans(
        "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
        "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
    )
    return seq.translate(table)[::-1].upper()


def make_blast_db(assembly, db_prefix, force=False):
    db_prefix = Path(db_prefix)
    existing_db_files = list(db_prefix.parent.glob(db_prefix.name + ".*"))

    if existing_db_files and not force:
        print(f"\n[INFO] Banco BLAST já existe: {db_prefix}")
        print("[INFO] Use --force-db para recriar o banco.")
        return

    cmd = [
        "makeblastdb",
        "-in", str(assembly),
        "-dbtype", "nucl",
        "-out", str(db_prefix),
    ]

    run_command(cmd)


def run_tblastn(query_proteins, db_prefix, output_tsv, threads, evalue,
                max_target_seqs, seg, db_gencode):
    outfmt = "6 " + " ".join(BLAST_COLUMNS)

    cmd = [
        "tblastn",
        "-query", str(query_proteins),
        "-db", str(db_prefix),
        "-out", str(output_tsv),
        "-outfmt", outfmt,
        "-evalue", str(evalue),
        "-num_threads", str(threads),
        "-max_target_seqs", str(max_target_seqs),
        "-seg", str(seg),
        "-db_gencode", str(db_gencode),
    ]

    run_command(cmd)


def parse_blast_tsv(blast_tsv):
    hits = []

    if not Path(blast_tsv).exists() or Path(blast_tsv).stat().st_size == 0:
        return hits

    with open(blast_tsv, "r") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) != len(BLAST_COLUMNS):
                raise ValueError(
                    f"Linha BLAST com número inesperado de colunas: {line}"
                )

            hit = dict(zip(BLAST_COLUMNS, parts))

            for col in ["pident", "evalue", "bitscore", "qcovhsp", "qcovs"]:
                hit[col] = float(hit[col])

            for col in [
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "qlen",
                "slen",
                "sframe",
            ]:
                hit[col] = int(float(hit[col]))

            hit_start = min(hit["sstart"], hit["send"])
            hit_end = max(hit["sstart"], hit["send"])
            hit["hit_start"] = hit_start
            hit["hit_end"] = hit_end
            hit["hit_span_bp"] = hit_end - hit_start + 1

            if hit["sframe"] < 0 or hit["sstart"] > hit["send"]:
                hit["strand"] = "-"
            else:
                hit["strand"] = "+"

            hits.append(hit)

    return hits


def write_csv(rows, output_csv, fieldnames):
    with open(output_csv, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def filter_hits(hits, min_pident, min_qcov_hsp, min_bitscore):
    filtered = []

    for hit in hits:
        if hit["pident"] < min_pident:
            continue
        if hit["qcovhsp"] < min_qcov_hsp:
            continue
        if hit["bitscore"] < min_bitscore:
            continue

        filtered.append(hit)

    return filtered


def merge_intervals(intervals):
    if not intervals:
        return []

    intervals = sorted(intervals)
    merged = [list(intervals[0])]

    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]

        if start <= last_end + 1:
            merged[-1][1] = max(last_end, end)
        else:
            merged.append([start, end])

    return [(start, end) for start, end in merged]


def summarize_by_scaffold(hits):
    grouped = defaultdict(list)

    for hit in hits:
        grouped[hit["sseqid"]].append(hit)

    summary = []

    for scaffold, scaffold_hits in grouped.items():
        intervals = [
            (hit["hit_start"], hit["hit_end"])
            for hit in scaffold_hits
        ]
        merged = merge_intervals(intervals)

        covered_bp = sum(end - start + 1 for start, end in merged)
        min_coord = min(start for start, end in intervals)
        max_coord = max(end for start, end in intervals)
        span_bp = max_coord - min_coord + 1

        query_ids = sorted(set(hit["qseqid"] for hit in scaffold_hits))

        row = {
            "sseqid": scaffold,
            "n_hits": len(scaffold_hits),
            "n_unique_queries": len(query_ids),
            "scaffold_length_bp": scaffold_hits[0]["slen"],
            "min_hit_coord": min_coord,
            "max_hit_coord": max_coord,
            "span_bp": span_bp,
            "covered_bp_merged_hits": covered_bp,
            "best_evalue": min(hit["evalue"] for hit in scaffold_hits),
            "max_bitscore": max(hit["bitscore"] for hit in scaffold_hits),
            "mean_pident": round(
                sum(hit["pident"] for hit in scaffold_hits) / len(scaffold_hits),
                3,
            ),
            "mean_qcovhsp": round(
                sum(hit["qcovhsp"] for hit in scaffold_hits) / len(scaffold_hits),
                3,
            ),
            "queries": ";".join(query_ids),
        }

        summary.append(row)

    summary.sort(
        key=lambda r: (
            -r["n_unique_queries"],
            -r["covered_bp_merged_hits"],
            r["best_evalue"],
            -r["max_bitscore"],
        )
    )

    return summary


def best_hit_per_query(hits):
    grouped = defaultdict(list)

    for hit in hits:
        grouped[hit["qseqid"]].append(hit)

    best_rows = []

    for query, query_hits in grouped.items():
        query_hits.sort(
            key=lambda h: (
                h["evalue"],
                -h["bitscore"],
                -h["qcovhsp"],
                -h["pident"],
            )
        )
        best_rows.append(query_hits[0])

    best_rows.sort(
        key=lambda h: (
            h["qseqid"],
            h["evalue"],
            -h["bitscore"],
        )
    )

    return best_rows


def print_summary(total_hits, filtered_hits, scaffold_summary, args):
    print("\n" + "=" * 80)
    print("RESUMO DO tBLASTn")
    print("=" * 80)
    print(f"Assembly:              {args.assembly}")
    print(f"Proteínas query:       {args.cds_prot}")
    print(f"Total de hits BLAST:   {total_hits}")
    print(f"Hits após filtros:     {filtered_hits}")
    print(f"Filtros aplicados:     pident >= {args.min_pident}; "
          f"qcovhsp >= {args.min_qcov_hsp}; bitscore >= {args.min_bitscore}")
    print(f"Código genético alvo:  {args.db_gencode}")
    print("=" * 80)

    if not scaffold_summary:
        print("\nNenhum scaffold passou nos filtros.")
        return

    print("\nTOP SCAFFOLDS CANDIDATOS")
    print("-" * 80)

    header = (
        f"{'rank':<5}"
        f"{'scaffold':<28}"
        f"{'queries':>8}"
        f"{'hits':>8}"
        f"{'len_bp':>12}"
        f"{'covered_bp':>14}"
        f"{'best_eval':>14}"
        f"{'bitscore':>12}"
    )
    print(header)
    print("-" * 80)

    for i, row in enumerate(scaffold_summary[:args.top_scaffolds], start=1):
        print(
            f"{i:<5}"
            f"{row['sseqid']:<28.28}"
            f"{row['n_unique_queries']:>8}"
            f"{row['n_hits']:>8}"
            f"{row['scaffold_length_bp']:>12}"
            f"{row['covered_bp_merged_hits']:>14}"
            f"{row['best_evalue']:>14.2e}"
            f"{row['max_bitscore']:>12.1f}"
        )

    print("-" * 80)


def extract_candidate_scaffolds(assembly_records, scaffold_summary,
                                output_fasta, min_queries_per_scaffold,
                                top_scaffolds):
    fasta_records = []

    selected = [
        row for row in scaffold_summary
        if row["n_unique_queries"] >= min_queries_per_scaffold
    ]

    if not selected and scaffold_summary:
        print(
            "\n[AVISO] Nenhum scaffold atingiu --min-queries-per-scaffold. "
            "Extraindo o melhor scaffold mesmo assim."
        )
        selected = scaffold_summary[:1]

    selected = selected[:top_scaffolds]

    for row in selected:
        scaffold_id = row["sseqid"]

        if scaffold_id not in assembly_records:
            print(f"[AVISO] Scaffold não encontrado no FASTA: {scaffold_id}")
            continue

        seq = assembly_records[scaffold_id]["seq"]

        header = (
            f"{scaffold_id} "
            f"n_unique_queries={row['n_unique_queries']} "
            f"n_hits={row['n_hits']} "
            f"covered_bp={row['covered_bp_merged_hits']} "
            f"best_evalue={row['best_evalue']} "
            f"max_bitscore={row['max_bitscore']}"
        )

        fasta_records.append((header, seq))

    write_fasta(fasta_records, output_fasta)


def extract_hit_regions(assembly_records, hits, output_fasta, padding):
    fasta_records = []

    for i, hit in enumerate(hits, start=1):
        scaffold_id = hit["sseqid"]

        if scaffold_id not in assembly_records:
            print(f"[AVISO] Scaffold não encontrado no FASTA: {scaffold_id}")
            continue

        scaffold_seq = assembly_records[scaffold_id]["seq"]
        scaffold_len = len(scaffold_seq)

        start = max(1, hit["hit_start"] - padding)
        end = min(scaffold_len, hit["hit_end"] + padding)

        subseq = scaffold_seq[start - 1:end]

        if hit["strand"] == "-":
            subseq = reverse_complement(subseq)

        header = (
            f"hit_{i}|query={hit['qseqid']}|scaffold={scaffold_id}|"
            f"region={start}-{end}|hit={hit['hit_start']}-{hit['hit_end']}|"
            f"strand={hit['strand']}|pident={hit['pident']}|"
            f"qcovhsp={hit['qcovhsp']}|evalue={hit['evalue']}|"
            f"bitscore={hit['bitscore']}"
        )

        fasta_records.append((header, subseq))

    write_fasta(fasta_records, output_fasta)


def write_matched_queries(query_fasta, matched_ids, output_fasta):
    if query_fasta is None:
        return

    records = read_fasta(query_fasta)
    fasta_records = []

    for seq_id, rec in records.items():
        if seq_id in matched_ids:
            fasta_records.append((rec["header"], rec["seq"]))

    write_fasta(fasta_records, output_fasta)

    if len(fasta_records) == 0:
        print(
            f"[AVISO] Nenhuma sequência de {query_fasta} foi extraída. "
            f"Possível diferença entre IDs do FASTA proteico e do FASTA nucleotídico."
        )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Roda tBLASTn de CDS/proteínas mitocondriais contra uma montagem "
            "e extrai scaffolds/regiões candidatos ao mitogenoma."
        )
    )

    parser.add_argument(
        "--assembly",
        required=True,
        help="FASTA da montagem. Exemplo: scaffolds_1000bp.fasta",
    )
    parser.add_argument(
        "--cds-prot",
        required=True,
        help="FASTA das proteínas dos CDS mitocondriais. Exemplo: Agabis_H97_prot.fasta",
    )
    parser.add_argument(
        "--cds-nt",
        default=None,
        help="FASTA nucleotídico dos CDS. Opcional. Exemplo: Agabis_H97_nt.fasta",
    )
    parser.add_argument(
        "--outdir",
        default="tblastn_mitogenome_results",
        help="Diretório de saída.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Número de threads para o BLAST.",
    )
    parser.add_argument(
        "--evalue",
        type=float,
        default=1e-5,
        help="E-value máximo usado no tBLASTn.",
    )
    parser.add_argument(
        "--min-pident",
        type=float,
        default=25.0,
        help="Identidade mínima percentual para manter hit na tabela filtrada.",
    )
    parser.add_argument(
        "--min-qcov-hsp",
        type=float,
        default=30.0,
        help="Cobertura mínima do HSP sobre a query.",
    )
    parser.add_argument(
        "--min-bitscore",
        type=float,
        default=50.0,
        help="Bitscore mínimo para manter hit na tabela filtrada.",
    )
    parser.add_argument(
        "--max-target-seqs",
        type=int,
        default=5000,
        help="Número máximo de alvos reportados por query no BLAST.",
    )
    parser.add_argument(
        "--seg",
        default="no",
        choices=["yes", "no"],
        help="Usar filtro SEG para baixa complexidade nas proteínas.",
    )
    parser.add_argument(
        "--db-gencode",
        type=int,
        default=4,
        help=(
            "Código genético usado para traduzir o banco nucleotídico no tBLASTn. "
            "Para muitos mitogenomas fúngicos, use 4. Para código padrão, use 1."
        ),
    )
    parser.add_argument(
        "--padding",
        type=int,
        default=1000,
        help="Bases extras antes/depois de cada hit ao extrair regiões FASTA.",
    )
    parser.add_argument(
        "--min-queries-per-scaffold",
        type=int,
        default=2,
        help=(
            "Número mínimo de CDS/proteínas diferentes batendo no scaffold "
            "para extrair como candidato."
        ),
    )
    parser.add_argument(
        "--top-scaffolds",
        type=int,
        default=10,
        help="Número máximo de scaffolds candidatos a extrair.",
    )
    parser.add_argument(
        "--force-db",
        action="store_true",
        help="Recria o banco BLAST mesmo se ele já existir.",
    )

    args = parser.parse_args()

    check_executable("makeblastdb")
    check_executable("tblastn")

    assembly = Path(args.assembly).resolve()
    cds_prot = Path(args.cds_prot).resolve()
    cds_nt = Path(args.cds_nt).resolve() if args.cds_nt else None

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    db_dir = outdir / "blastdb"
    db_dir.mkdir(parents=True, exist_ok=True)

    tables_dir = outdir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    fasta_dir = outdir / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)

    db_prefix = db_dir / assembly.stem

    raw_blast_tsv = tables_dir / "tblastn_raw.tsv"
    blast_hits_csv = tables_dir / "tblastn_hits_all.csv"
    filtered_hits_csv = tables_dir / "tblastn_hits_filtered.csv"
    best_hits_csv = tables_dir / "best_hit_per_query.csv"
    scaffold_summary_csv = tables_dir / "scaffold_summary.csv"

    candidate_scaffolds_fasta = fasta_dir / "candidate_mitogenome_scaffolds.fasta"
    hit_regions_fasta = fasta_dir / "tblastn_hit_regions_plus_padding.fasta"
    matched_prot_fasta = fasta_dir / "matched_cds_proteins.fasta"
    matched_nt_fasta = fasta_dir / "matched_cds_nucleotides.fasta"

    print("\n[INFO] Lendo assembly FASTA...")
    assembly_records = read_fasta(assembly)
    print(f"[INFO] Scaffolds carregados: {len(assembly_records)}")

    print("\n[INFO] Criando banco BLAST da montagem...")
    make_blast_db(assembly, db_prefix, force=args.force_db)

    print("\n[INFO] Rodando tBLASTn...")
    run_tblastn(
        query_proteins=cds_prot,
        db_prefix=db_prefix,
        output_tsv=raw_blast_tsv,
        threads=args.threads,
        evalue=args.evalue,
        max_target_seqs=args.max_target_seqs,
        seg=args.seg,
        db_gencode=args.db_gencode,
    )

    print("\n[INFO] Processando resultados do BLAST...")
    hits = parse_blast_tsv(raw_blast_tsv)

    all_fieldnames = BLAST_COLUMNS + [
        "hit_start",
        "hit_end",
        "hit_span_bp",
        "strand",
    ]

    write_csv(hits, blast_hits_csv, all_fieldnames)

    filtered_hits = filter_hits(
        hits,
        min_pident=args.min_pident,
        min_qcov_hsp=args.min_qcov_hsp,
        min_bitscore=args.min_bitscore,
    )

    write_csv(filtered_hits, filtered_hits_csv, all_fieldnames)

    best_hits = best_hit_per_query(filtered_hits)
    write_csv(best_hits, best_hits_csv, all_fieldnames)

    scaffold_summary = summarize_by_scaffold(filtered_hits)

    scaffold_summary_fieldnames = [
        "sseqid",
        "n_hits",
        "n_unique_queries",
        "scaffold_length_bp",
        "min_hit_coord",
        "max_hit_coord",
        "span_bp",
        "covered_bp_merged_hits",
        "best_evalue",
        "max_bitscore",
        "mean_pident",
        "mean_qcovhsp",
        "queries",
    ]

    write_csv(scaffold_summary, scaffold_summary_csv, scaffold_summary_fieldnames)

    print_summary(
        total_hits=len(hits),
        filtered_hits=len(filtered_hits),
        scaffold_summary=scaffold_summary,
        args=args,
    )

    print("\n[INFO] Extraindo FASTA dos scaffolds candidatos...")
    extract_candidate_scaffolds(
        assembly_records=assembly_records,
        scaffold_summary=scaffold_summary,
        output_fasta=candidate_scaffolds_fasta,
        min_queries_per_scaffold=args.min_queries_per_scaffold,
        top_scaffolds=args.top_scaffolds,
    )

    print("\n[INFO] Extraindo FASTA das regiões com hits...")
    extract_hit_regions(
        assembly_records=assembly_records,
        hits=filtered_hits,
        output_fasta=hit_regions_fasta,
        padding=args.padding,
    )

    matched_query_ids = sorted(set(hit["qseqid"] for hit in filtered_hits))

    print("\n[INFO] Extraindo proteínas query que tiveram hits...")
    write_matched_queries(
        query_fasta=cds_prot,
        matched_ids=matched_query_ids,
        output_fasta=matched_prot_fasta,
    )

    if cds_nt is not None:
        print("\n[INFO] Extraindo CDS nucleotídicos correspondentes às queries com hits...")
        write_matched_queries(
            query_fasta=cds_nt,
            matched_ids=matched_query_ids,
            output_fasta=matched_nt_fasta,
        )

    print("\n" + "=" * 80)
    print("ARQUIVOS GERADOS")
    print("=" * 80)
    print(f"BLAST bruto TSV:              {raw_blast_tsv}")
    print(f"Todos os hits CSV:            {blast_hits_csv}")
    print(f"Hits filtrados CSV:           {filtered_hits_csv}")
    print(f"Melhor hit por query CSV:     {best_hits_csv}")
    print(f"Resumo por scaffold CSV:      {scaffold_summary_csv}")
    print(f"Scaffolds candidatos FASTA:   {candidate_scaffolds_fasta}")
    print(f"Regiões dos hits FASTA:       {hit_regions_fasta}")
    print(f"Proteínas com hits FASTA:     {matched_prot_fasta}")

    if cds_nt is not None:
        print(f"CDS nucleotídicos com hits:   {matched_nt_fasta}")

    print("=" * 80)
    print("\n[OK] Análise finalizada.")


if __name__ == "__main__":
    main()

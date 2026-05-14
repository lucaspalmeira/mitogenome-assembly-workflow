#!/usr/bin/env python3
"""Unified CLI for the object-oriented analysis workflow."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

from workflow import AssemblyWorkflow, WorkflowConfig
from workflow.orchestrator import WorkflowError


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_workflow.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Runs the complete analysis workflow: tBLASTn discovery, read "
            "recruitment, targeted reassembly, post-reassembly validation, and "
            "optional curation."
        ),
    )

    parser.add_argument("--assembly", type=Path, required=True, help="Original assembly FASTA.")
    parser.add_argument("--cds-prot", type=Path, required=True, help="Mitochondrial CDS proteins FASTA.")
    parser.add_argument("--cds-nt", type=Path, help="Mitochondrial CDS nucleotide FASTA.")
    parser.add_argument("--short-r1", type=Path, required=True, help="Original Illumina R1 FASTQ.")
    parser.add_argument("--short-r2", type=Path, required=True, help="Original Illumina R2 FASTQ.")
    parser.add_argument("--pacbio-reads", type=Path, help="Original PacBio reads, FASTQ or FASTA, optionally gzipped.")
    parser.add_argument("--scaffold-summary", type=Path, help="Existing scaffold_summary.csv to skip initial tBLASTn.")
    parser.add_argument("--outdir", type=Path, default=Path("analysis_run"), help="Workflow output directory.")
    parser.add_argument("--threads", type=int, default=12, help="Threads for external tools.")
    parser.add_argument("--memory", type=int, default=64, help="SPAdes memory limit in GB.")
    parser.add_argument("--dry-run", action="store_true", help="Write the workflow commands without executing them.")

    recruitment = parser.add_argument_group("Recruitment and reassembly")
    recruitment.add_argument("--top-n", type=int, default=4, help="Number of candidate scaffolds to select.")
    recruitment.add_argument("--min-unique-queries", type=int, default=1, help="Minimum distinct mitochondrial queries per scaffold.")
    recruitment.add_argument("--manual-scaffolds", default="", help="Comma-separated scaffold IDs to use instead of automatic ranking.")
    recruitment.add_argument("--short-aligner", choices=["bowtie2", "bwa"], default="bowtie2", help="Short-read aligner.")
    recruitment.add_argument("--min-mapq", type=int, default=0, help="Minimum MAPQ for read recruitment.")
    recruitment.add_argument("--pacbio-type", choices=["raw", "hifi", "corrected"], default="raw", help="PacBio type for minimap2 preset.")
    recruitment.add_argument("--skip-assembly", action="store_true", help="Recruit reads but skip SPAdes reassembly.")
    recruitment.add_argument("--run-cap3", action="store_true", help="Run optional CAP3 on selected candidate scaffolds.")

    validation = parser.add_argument_group("BLAST/tBLASTn filters")
    validation.add_argument("--skip-validation", action="store_true", help="Skip tBLASTn validation after reassembly.")
    validation.add_argument("--evalue", type=float, default=1e-5, help="E-value threshold for tBLASTn.")
    validation.add_argument("--min-pident", type=float, default=25.0, help="Minimum protein identity for tBLASTn hits.")
    validation.add_argument("--min-qcov-hsp", type=float, default=30.0, help="Minimum HSP query coverage percentage.")
    validation.add_argument("--min-bitscore", type=float, default=50.0, help="Minimum BLAST bitscore.")
    validation.add_argument("--max-target-seqs", type=int, default=5000, help="Maximum BLAST targets per query.")
    validation.add_argument("--seg", choices=["yes", "no"], default="no", help="Use SEG low-complexity filtering in tBLASTn.")
    validation.add_argument("--db-gencode", type=int, default=4, help="Genetic code used by tBLASTn to translate the target database.")
    validation.add_argument("--padding", type=int, default=1000, help="Padding around tBLASTn hit regions.")
    validation.add_argument("--min-queries-per-scaffold", type=int, default=2, help="Minimum queries to extract a candidate scaffold.")
    validation.add_argument("--top-scaffolds", type=int, default=10, help="Maximum candidate scaffolds to extract.")
    validation.add_argument("--force-db", action="store_true", help="Recreate BLAST databases even when files exist.")

    curation = parser.add_argument_group("Optional curation")
    curation.add_argument("--curate-final", action="store_true", help="Run toolkit/curate.py on final candidates.")
    curation.add_argument("--curation-topology", choices=["circular", "linear", "unknown"], default="circular", help="Expected final molecule topology.")
    curation.add_argument("--rotate-to", default="", help="Gene used to rotate the circular sequence, e.g. COX1.")
    curation.add_argument("--genetic-code", type=int, default=4, help="NCBI genetic code for curation.")
    curation.add_argument("--expected-from-ref", action="store_true", help="Use reference FASTA gene names as expected genes during curation.")

    return parser


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = build_parser()
    args = parser.parse_args(argv)

    required_paths = [args.assembly, args.cds_prot, args.short_r1, args.short_r2]
    optional_paths = [args.cds_nt, args.pacbio_reads, args.scaffold_summary]
    for path in required_paths + [p for p in optional_paths if p is not None]:
        if not path.exists():
            parser.error(f"File not found: {path}")
    if args.threads < 1:
        parser.error("--threads must be >= 1.")
    if args.memory < 1:
        parser.error("--memory must be >= 1.")
    if args.min_mapq < 0:
        parser.error("--min-mapq must be >= 0.")
    return args


def config_from_args(args: argparse.Namespace) -> WorkflowConfig:
    return WorkflowConfig(
        assembly=args.assembly.resolve(),
        cds_prot=args.cds_prot.resolve(),
        cds_nt=args.cds_nt.resolve() if args.cds_nt else None,
        short_r1=args.short_r1.resolve(),
        short_r2=args.short_r2.resolve(),
        pacbio_reads=args.pacbio_reads.resolve() if args.pacbio_reads else None,
        scaffold_summary=args.scaffold_summary.resolve() if args.scaffold_summary else None,
        outdir=args.outdir.resolve(),
        threads=args.threads,
        memory=args.memory,
        top_n=args.top_n,
        min_unique_queries=args.min_unique_queries,
        manual_scaffolds=args.manual_scaffolds,
        short_aligner=args.short_aligner,
        skip_assembly=args.skip_assembly,
        run_cap3=args.run_cap3,
        validate_reassembly=not args.skip_validation,
        curate_final=args.curate_final,
        curation_topology=args.curation_topology,
        rotate_to=args.rotate_to,
        genetic_code=args.genetic_code,
        expected_from_ref=args.expected_from_ref,
        dry_run=args.dry_run,
        evalue=args.evalue,
        min_pident=args.min_pident,
        min_qcov_hsp=args.min_qcov_hsp,
        min_bitscore=args.min_bitscore,
        max_target_seqs=args.max_target_seqs,
        seg=args.seg,
        db_gencode=args.db_gencode,
        padding=args.padding,
        min_queries_per_scaffold=args.min_queries_per_scaffold,
        top_scaffolds=args.top_scaffolds,
        force_db=args.force_db,
        min_mapq=args.min_mapq,
        pacbio_type=args.pacbio_type,
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    config = config_from_args(args)

    try:
        outputs = AssemblyWorkflow(config).run()
    except WorkflowError as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 1

    print("\nWorkflow finished.")
    print(f"Report: {outputs.workflow_report}")
    print(f"Commands log: {outputs.commands_log}")
    if outputs.final_candidates_fasta:
        print(f"Final candidates: {outputs.final_candidates_fasta}")
    elif not config.skip_assembly and outputs.reassembly_scaffolds:
        print(f"Reassembly scaffolds: {outputs.reassembly_scaffolds}")
    elif outputs.selected_scaffolds_fasta:
        print(f"Selected scaffolds: {outputs.selected_scaffolds_fasta}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""High-level orchestration for the assembly analysis.

This module intentionally keeps the original scripts as the execution units and
adds a small object-oriented layer around them. That makes the full workflow
easier to run end-to-end without duplicating the bioinformatics logic already
present in the repository.
"""

from __future__ import annotations

import json
import shlex
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import List, Optional, Sequence


class WorkflowError(RuntimeError):
    """Raised when a workflow stage cannot continue."""


@dataclass
class WorkflowConfig:
    """Configuration shared by all workflow stages."""

    assembly: Path
    cds_prot: Path
    short_r1: Path
    short_r2: Path
    outdir: Path
    cds_nt: Optional[Path] = None
    pacbio_reads: Optional[Path] = None
    scaffold_summary: Optional[Path] = None
    threads: int = 12
    memory: int = 64
    top_n: int = 4
    min_unique_queries: int = 1
    manual_scaffolds: str = ""
    short_aligner: str = "bowtie2"
    skip_assembly: bool = False
    run_cap3: bool = False
    validate_reassembly: bool = True
    curate_final: bool = False
    curation_topology: str = "circular"
    rotate_to: str = ""
    genetic_code: int = 4
    expected_from_ref: bool = False
    dry_run: bool = False
    evalue: float = 1e-5
    min_pident: float = 25.0
    min_qcov_hsp: float = 30.0
    min_bitscore: float = 50.0
    max_target_seqs: int = 5000
    seg: str = "no"
    db_gencode: int = 4
    padding: int = 1000
    min_queries_per_scaffold: int = 2
    top_scaffolds: int = 10
    force_db: bool = False
    min_mapq: int = 0
    pacbio_type: str = "raw"

    @property
    def project_root(self) -> Path:
        return Path(__file__).resolve().parents[1]


@dataclass
class WorkflowOutputs:
    """Important files produced or expected by the workflow."""

    initial_tblastn_dir: Optional[Path] = None
    scaffold_summary: Optional[Path] = None
    recruitment_dir: Optional[Path] = None
    selected_scaffolds_csv: Optional[Path] = None
    selected_scaffolds_fasta: Optional[Path] = None
    reassembly_scaffolds: Optional[Path] = None
    validation_tblastn_dir: Optional[Path] = None
    validation_summary: Optional[Path] = None
    final_candidates_fasta: Optional[Path] = None
    curation_dir: Optional[Path] = None
    curation_report: Optional[Path] = None
    workflow_report: Optional[Path] = None
    commands_log: Optional[Path] = None


@dataclass
class StageRecord:
    name: str
    status: str
    outputs: List[str] = field(default_factory=list)
    note: str = ""


class CommandRunner:
    """Runs external commands and records a reproducible command log."""

    def __init__(self, log_file: Path, dry_run: bool = False) -> None:
        self.log_file = log_file
        self.dry_run = dry_run
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        self.log_file.write_text(f"# Workflow command log - {now()}\n", encoding="utf-8")

    def note(self, text: str) -> None:
        with self.log_file.open("a", encoding="utf-8") as handle:
            handle.write(f"\n# {text}\n")

    def run(self, cmd: Sequence[object], cwd: Optional[Path] = None) -> None:
        cmd_list = [str(part) for part in cmd]
        cmd_text = " ".join(shlex.quote(part) for part in cmd_list)
        with self.log_file.open("a", encoding="utf-8") as handle:
            handle.write(f"\n[{now()}] $ {cmd_text}\n")
            if cwd:
                handle.write(f"# cwd: {cwd}\n")

        if self.dry_run:
            return

        result = subprocess.run(
            cmd_list,
            cwd=str(cwd) if cwd else None,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        with self.log_file.open("a", encoding="utf-8") as handle:
            if result.stdout:
                handle.write(result.stdout)
                if not result.stdout.endswith("\n"):
                    handle.write("\n")

        if result.returncode != 0:
            tail = (result.stdout or "")[-4000:]
            raise WorkflowError(
                f"Stage command failed with exit code {result.returncode}: {cmd_text}\n{tail}"
            )


class WorkflowStage:
    """Base class for composable workflow stages."""

    name = "stage"

    def run(self, config: WorkflowConfig, runner: CommandRunner, outputs: WorkflowOutputs) -> StageRecord:
        raise NotImplementedError

    @staticmethod
    def _existing_or_expected(paths: Sequence[Optional[Path]]) -> List[str]:
        values = []
        for path in paths:
            if path is not None:
                values.append(str(path))
        return values


class TblastnDiscoveryStage(WorkflowStage):
    name = "01_discovery"

    def run(self, config: WorkflowConfig, runner: CommandRunner, outputs: WorkflowOutputs) -> StageRecord:
        if config.scaffold_summary:
            outputs.scaffold_summary = config.scaffold_summary.resolve()
            return StageRecord(
                name=self.name,
                status="skipped",
                outputs=[str(outputs.scaffold_summary)],
                note="Using scaffold summary supplied by --scaffold-summary.",
            )

        outdir = config.outdir / "01_discovery"
        script = config.project_root / "find_mitogenome_tblastn.py"
        cmd: List[object] = [
            sys.executable,
            script,
            "--assembly",
            config.assembly,
            "--cds-prot",
            config.cds_prot,
            "--outdir",
            outdir,
            "--threads",
            config.threads,
            "--evalue",
            config.evalue,
            "--min-pident",
            config.min_pident,
            "--min-qcov-hsp",
            config.min_qcov_hsp,
            "--min-bitscore",
            config.min_bitscore,
            "--max-target-seqs",
            config.max_target_seqs,
            "--seg",
            config.seg,
            "--db-gencode",
            config.db_gencode,
            "--padding",
            config.padding,
            "--min-queries-per-scaffold",
            config.min_queries_per_scaffold,
            "--top-scaffolds",
            config.top_scaffolds,
        ]
        if config.cds_nt:
            cmd.extend(["--cds-nt", config.cds_nt])
        if config.force_db:
            cmd.append("--force-db")

        runner.run(cmd, cwd=config.project_root)
        outputs.initial_tblastn_dir = outdir
        outputs.scaffold_summary = outdir / "tables" / "scaffold_summary.csv"
        return StageRecord(
            name=self.name,
            status="planned" if config.dry_run else "completed",
            outputs=self._existing_or_expected([outputs.scaffold_summary]),
        )


class RecruitReassembleStage(WorkflowStage):
    name = "02_recruitment"

    def run(self, config: WorkflowConfig, runner: CommandRunner, outputs: WorkflowOutputs) -> StageRecord:
        if outputs.scaffold_summary is None:
            raise WorkflowError("No scaffold summary is available for recruitment.")

        outdir = config.outdir / "02_recruitment"
        script = config.project_root / "mito_recruit_reassemble_pipeline.py"
        cmd: List[object] = [
            sys.executable,
            script,
            "--assembly",
            config.assembly,
            "--scaffold-summary",
            outputs.scaffold_summary,
            "--short-r1",
            config.short_r1,
            "--short-r2",
            config.short_r2,
            "--outdir",
            outdir,
            "--threads",
            config.threads,
            "--memory",
            config.memory,
            "--top-n",
            config.top_n,
            "--min-unique-queries",
            config.min_unique_queries,
            "--short-aligner",
            config.short_aligner,
            "--min-mapq",
            config.min_mapq,
            "--pacbio-type",
            config.pacbio_type,
        ]
        if config.pacbio_reads:
            cmd.extend(["--pacbio-reads", config.pacbio_reads])
        if config.manual_scaffolds:
            cmd.extend(["--manual-scaffolds", config.manual_scaffolds])
        if config.skip_assembly:
            cmd.append("--skip-assembly")
        if config.run_cap3:
            cmd.append("--run-cap3")

        runner.run(cmd, cwd=config.project_root)
        outputs.recruitment_dir = outdir
        outputs.selected_scaffolds_csv = outdir / "selected_scaffolds.csv"
        outputs.selected_scaffolds_fasta = outdir / "selected_candidate_scaffolds.fasta"
        outputs.reassembly_scaffolds = None if config.skip_assembly else outdir / "spades_reassembly" / "scaffolds.fasta"

        return StageRecord(
            name=self.name,
            status="planned" if config.dry_run else "completed",
            outputs=self._existing_or_expected(
                [
                    outputs.selected_scaffolds_csv,
                    outputs.selected_scaffolds_fasta,
                    outputs.reassembly_scaffolds,
                ]
            ),
        )


class ReassemblyValidationStage(WorkflowStage):
    name = "03_validation"

    def run(self, config: WorkflowConfig, runner: CommandRunner, outputs: WorkflowOutputs) -> StageRecord:
        if not config.validate_reassembly:
            return StageRecord(name=self.name, status="skipped", note="Validation disabled.")
        if config.skip_assembly:
            return StageRecord(name=self.name, status="skipped", note="Reassembly was skipped.")
        if outputs.reassembly_scaffolds is None:
            raise WorkflowError("No reassembly FASTA path is available for validation.")
        if not config.dry_run and not outputs.reassembly_scaffolds.exists():
            raise WorkflowError(f"Expected reassembly FASTA was not found: {outputs.reassembly_scaffolds}")

        outdir = config.outdir / "03_validation"
        script = config.project_root / "find_mitogenome_tblastn.py"
        cmd: List[object] = [
            sys.executable,
            script,
            "--assembly",
            outputs.reassembly_scaffolds,
            "--cds-prot",
            config.cds_prot,
            "--outdir",
            outdir,
            "--threads",
            config.threads,
            "--evalue",
            config.evalue,
            "--min-pident",
            config.min_pident,
            "--min-qcov-hsp",
            config.min_qcov_hsp,
            "--min-bitscore",
            config.min_bitscore,
            "--max-target-seqs",
            config.max_target_seqs,
            "--seg",
            config.seg,
            "--db-gencode",
            config.db_gencode,
            "--padding",
            config.padding,
            "--min-queries-per-scaffold",
            config.min_queries_per_scaffold,
            "--top-scaffolds",
            config.top_scaffolds,
        ]
        if config.cds_nt:
            cmd.extend(["--cds-nt", config.cds_nt])
        if config.force_db:
            cmd.append("--force-db")

        runner.run(cmd, cwd=config.project_root)
        outputs.validation_tblastn_dir = outdir
        outputs.validation_summary = outdir / "tables" / "scaffold_summary.csv"
        outputs.final_candidates_fasta = outdir / "fasta" / "candidate_mitogenome_scaffolds.fasta"
        return StageRecord(
            name=self.name,
            status="planned" if config.dry_run else "completed",
            outputs=self._existing_or_expected([outputs.validation_summary, outputs.final_candidates_fasta]),
        )


class CurationStage(WorkflowStage):
    name = "04_curate"

    def run(self, config: WorkflowConfig, runner: CommandRunner, outputs: WorkflowOutputs) -> StageRecord:
        if not config.curate_final:
            return StageRecord(name=self.name, status="skipped", note="Curation disabled.")

        candidate_options = [outputs.final_candidates_fasta]
        if not config.skip_assembly:
            candidate_options.append(outputs.reassembly_scaffolds)
        candidate_options.append(outputs.selected_scaffolds_fasta)

        candidate = None
        for option in candidate_options:
            if option is None:
                continue
            if config.dry_run or option.exists():
                candidate = option
                break
        if candidate is None:
            candidate = next((option for option in candidate_options if option is not None), None)
        if candidate is None:
            raise WorkflowError("No candidate FASTA is available for curation.")
        if not config.dry_run and not candidate.exists():
            raise WorkflowError(f"Expected candidate FASTA was not found: {candidate}")

        script = config.project_root / "toolkit" / "curate.py"
        if not script.exists():
            raise WorkflowError(f"Curation script not found: {script}")

        outdir = config.outdir / "04_curate"
        cmd: List[object] = [
            sys.executable,
            script,
            "--candidate",
            candidate,
            "--input-format",
            "fasta",
            "--outdir",
            outdir,
            "--topology",
            config.curation_topology,
            "--genetic-code",
            config.genetic_code,
            "--threads",
            config.threads,
        ]
        if config.cds_nt:
            cmd.extend(["--nt-ref", config.cds_nt])
        if config.cds_prot:
            cmd.extend(["--prot-ref", config.cds_prot])
        if config.expected_from_ref:
            cmd.append("--expected-from-ref")
        if config.rotate_to:
            cmd.extend(["--rotate-to", config.rotate_to])
        if config.dry_run:
            cmd.append("--dry-run")

        runner.run(cmd, cwd=config.project_root)
        outputs.curation_dir = outdir
        outputs.curation_report = outdir / "CURATION_REPORT.md"
        return StageRecord(
            name=self.name,
            status="planned" if config.dry_run else "completed",
            outputs=self._existing_or_expected([outputs.curation_report]),
        )


class WorkflowReport:
    """Writes a compact report tying the stage outputs together."""

    def __init__(self, config: WorkflowConfig, outputs: WorkflowOutputs, records: Sequence[StageRecord]) -> None:
        self.config = config
        self.outputs = outputs
        self.records = list(records)

    def write(self) -> Path:
        report = self.config.outdir / "WORKFLOW_REPORT.md"
        self.outputs.workflow_report = report
        self.outputs.commands_log = self.config.outdir / "commands.log"
        report.parent.mkdir(parents=True, exist_ok=True)
        lines = [
            "# Assembly Workflow Report",
            "",
            f"Generated: {now()}",
            "",
            "## Configuration",
            "",
            "```json",
            json.dumps(_jsonable(asdict(self.config)), indent=2, ensure_ascii=False),
            "```",
            "",
            "## Stage Summary",
            "",
            "| Stage | Status | Main outputs | Note |",
            "| --- | --- | --- | --- |",
        ]
        for record in self.records:
            outputs = "<br>".join(record.outputs) if record.outputs else ""
            lines.append(f"| {record.name} | {record.status} | {outputs} | {record.note} |")

        lines.extend(["", "## Key Files", ""])
        key_files = [
            ("Commands log", self.outputs.commands_log),
            ("Initial scaffold summary", self.outputs.scaffold_summary),
            ("Selected scaffolds", self.outputs.selected_scaffolds_fasta),
            ("Reassembly scaffolds", self.outputs.reassembly_scaffolds),
            ("Reassembly validation summary", self.outputs.validation_summary),
            ("Final candidates FASTA", self.outputs.final_candidates_fasta),
            ("Curation report", self.outputs.curation_report),
        ]
        for label, path in key_files:
            if path is not None:
                lines.append(f"- {label}: `{path}`")

        lines.extend(
            [
                "",
                "## Interpretation Checklist",
                "",
                "- Prefer candidates with many distinct mitochondrial genes/proteins and high merged hit coverage.",
                "- Check that recruited-read coverage is broad and not concentrated in only short regions.",
                "- After reassembly, expect mitochondrial hits to collapse into fewer scaffolds than in the original assembly.",
                "- Treat circularity and rotation as evidence to review, not as proof without coverage/junction inspection.",
            ]
        )
        report.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return report


class AssemblyWorkflow:
    """Runs the complete modular workflow."""

    def __init__(self, config: WorkflowConfig) -> None:
        self.config = config
        self.outputs = WorkflowOutputs()
        self.records: List[StageRecord] = []

    def run(self) -> WorkflowOutputs:
        self.config.outdir.mkdir(parents=True, exist_ok=True)
        runner = CommandRunner(self.config.outdir / "commands.log", dry_run=self.config.dry_run)
        self.outputs.commands_log = runner.log_file
        self._write_config()

        stages: Sequence[WorkflowStage] = [
            TblastnDiscoveryStage(),
            RecruitReassembleStage(),
            ReassemblyValidationStage(),
            CurationStage(),
        ]

        for stage in stages:
            runner.note(f"Stage: {stage.name}")
            record = stage.run(self.config, runner, self.outputs)
            self.records.append(record)

        report = WorkflowReport(self.config, self.outputs, self.records).write()
        self.outputs.workflow_report = report
        return self.outputs

    def _write_config(self) -> None:
        config_path = self.config.outdir / "workflow_config.json"
        config_path.write_text(
            json.dumps(_jsonable(asdict(self.config)), indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )


def now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def _jsonable(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _jsonable(val) for key, val in value.items()}
    if isinstance(value, list):
        return [_jsonable(item) for item in value]
    return value

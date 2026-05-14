"""Small shared helpers for the standalone toolkit scripts."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Tuple


def write_fasta(records, path: Path, width: int = 80) -> None:
    """Write FASTA from either a dict or an iterable of ``(header, seq)`` pairs."""
    path.parent.mkdir(parents=True, exist_ok=True)
    items = records.items() if isinstance(records, dict) else records
    with path.open("w", encoding="utf-8") as handle:
        for header, seq in items:
            handle.write(f">{header}\n")
            for i in range(0, len(seq), width):
                handle.write(seq[i : i + width] + "\n")


def terminal_overlap(seq: str, min_overlap: int, max_overlap: int, min_identity: float) -> Tuple[int, float, bool]:
    """Compare contig prefix/suffix as a lightweight circularity screen."""
    seq = seq.upper()
    max_ov = min(max_overlap, len(seq) // 2)
    if max_ov < min_overlap:
        return 0, 0.0, False

    best_len = 0
    best_id = 0.0
    for overlap in range(max_ov, min_overlap - 1, -1):
        prefix = seq[:overlap]
        suffix = seq[-overlap:]
        matches = sum(1 for a, b in zip(prefix, suffix) if a == b)
        identity = matches / overlap if overlap else 0.0
        if identity > best_id:
            best_id = identity
            best_len = overlap
        if identity >= min_identity:
            return overlap, identity, True
    return best_len, best_id, False


def parse_flye_circular_flags(flye_dir: Path) -> Dict[str, Optional[bool]]:
    """Parse Flye assembly_info.txt circular flags when available."""
    info = flye_dir / "assembly_info.txt"
    flags: Dict[str, Optional[bool]] = {}
    if not info.exists():
        return flags

    with info.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            lower = [part.lower() for part in parts]
            if "seq_name" in lower or "#seq_name" in lower:
                continue
            if len(parts) < 4:
                continue

            value = parts[3].strip().lower()
            if value in {"y", "yes", "true", "circular"}:
                flags[parts[0]] = True
            elif value in {"n", "no", "false", "linear"}:
                flags[parts[0]] = False
            else:
                flags[parts[0]] = None

    return flags

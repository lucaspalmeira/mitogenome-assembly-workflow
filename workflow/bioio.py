"""Shared FASTA/FASTQ helpers used by the workflow scripts."""

from __future__ import annotations

import gzip
from typing import Iterable, Iterator, Sequence, TextIO, Tuple


FastqRecord = Tuple[str, str, str, str, str]


def smart_open(path, mode: str = "rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_fasta_records(fasta_path) -> dict:
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


def write_fasta(records: Iterable[Sequence[str]], output_path, width: int = 80) -> None:
    with open(output_path, "w") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")


def reverse_complement(seq: str) -> str:
    table = str.maketrans(
        "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
        "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
    )
    return seq.translate(table)[::-1].upper()


def detect_sequence_format(path) -> str:
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


def normalize_read_id(header: str) -> str:
    read_id = header.strip()

    if read_id.startswith("@") or read_id.startswith(">"):
        read_id = read_id[1:]

    read_id = read_id.split()[0]

    if read_id.endswith("/1") or read_id.endswith("/2"):
        read_id = read_id[:-2]

    return read_id


def fastq_iter(fastq_path) -> Iterator[FastqRecord]:
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


def fasta_iter(fasta_path) -> Iterator[Tuple[str, str, str]]:
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


def write_fastq_record(handle: TextIO, record: FastqRecord) -> None:
    h, s, p, q, _ = record
    handle.write(h)
    handle.write(s)
    handle.write(p)
    handle.write(q)


def write_fasta_record(handle: TextIO, header: str, seq: str, width: int = 80) -> None:
    handle.write(header + "\n")
    for i in range(0, len(seq), width):
        handle.write(seq[i:i + width] + "\n")


def load_read_ids(ids_path) -> set:
    ids = set()
    with open(ids_path, "r") as handle:
        for line in handle:
            line = line.strip()
            if line:
                ids.add(normalize_read_id(line))
    return ids

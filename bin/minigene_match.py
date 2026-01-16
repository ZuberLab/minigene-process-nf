#!/usr/bin/env python3
"""
Count minigene barcode–sequence matches with fuzzy matching.

FASTQ read layout (read orientation):
    [barcode][STOP][minigene][START]

Assumptions:
- Barcode and minigene in the FASTQ are the *reverse complement* of the
  corresponding sequences in the design table.
- The design table (e.g. DART2OS_pat001.txt) has:
    - 'sequence'          : barcode (forward)
    - 'sequence_minigene' : minigene sequence (forward)

This script:
  1. Extracts barcode and minigene segments from each read.
  2. Reverse-complements both segments.
  3. Fuzzy-matches barcode (<= 2 mismatches).
  4. For the matched barcode, compares minigene (<= 9 mismatches).
  5. Writes:
     - <prefix>.per_minigene_counts.tsv
     - <prefix>.summary.tsv

Output files
============

This script produces two tab-separated output files per run:

1) <prefix>.per_minigene_counts.tsv
-----------------------------------

One row per designed minigene (i.e. per row in the design table with a valid
barcode + minigene).

Columns:

- barcode
    Barcode sequence from the design file (column "sequence"), in forward
    orientation. This is the canonical barcode to which all fuzzy-matched
    barcode reads are assigned.

- minigene_seq
    Minigene sequence from the design file (column "sequence_minigene"), in
    forward orientation. This is what the reverse-complemented minigene
    segment extracted from the read is compared to (with up to
    max_mg_mismatches mismatches allowed).

- design_id
    Identifier for the minigene. If the design file contains a column "id",
    that value is used. Otherwise, if a column "name" exists, it is used.
    If neither column exists, the barcode itself is used as the design_id.
    This field is for reporting only and does not affect matching.

- reads_barcode_assigned
    Number of reads whose barcode segment (after reverse complement) was
    successfully assigned to this canonical barcode under the barcode
    fuzzy-matching rules:

      * The barcode segment has the expected length.
      * EITHER it matches this barcode exactly (0 mismatches),
        OR it matches this barcode with Hamming distance
        <= max_barc_mismatches (default 2),
        AND this match is unique (not ambiguous with any other barcode).

    Once a read's barcode is assigned to this barcode, it contributes +1
    here, regardless of whether the corresponding minigene sequence later
    passes or fails its mismatch threshold.

- reads_valid_pair
    Number of reads that:

      1) Have their barcode assigned to this canonical barcode (as above), AND
      2) Have a minigene segment (after reverse complement) whose Hamming
         distance to the expected minigene sequence for this barcode is
         <= max_mg_mismatches (default 9).

    These reads support a "valid" barcode–minigene pair under the specified
    mismatch thresholds. This is typically the main usable count per
    minigene.

- reads_seq_fail
    Number of reads that:

      1) Have their barcode assigned to this canonical barcode, BUT
      2) The minigene segment (after reverse complement) differs from the
         expected minigene sequence for this barcode by MORE than
         max_mg_mismatches.

    These reads indicate cases where the barcode looks correct (within the
    barcode mismatch threshold) but the associated minigene sequence is too
    divergent from the design.

Note: For each row (i.e. each barcode/minigene), we expect:
    reads_barcode_assigned = reads_valid_pair + reads_seq_fail
with the current script logic.


2) <prefix>.summary.tsv
------------------------

A small table with global counts and QC statistics for the entire FASTQ.

Columns:

- metric
    Name of the metric.

- value
    Integer count for that metric.

Metrics:

- total_reads
    Total number of reads processed from the input FASTQ file (i.e. number of
    complete 4-line FASTQ records read).

- too_short
    Number of reads whose total length is shorter than the expected layout:
        barcode_len + STOP_len + minigene_len + START_len
    For these reads, the script cannot reliably slice the read into barcode,
    STOP, minigene, and START segments, so they are discarded from further
    analysis.

- bc_exact
    Number of reads where:

      * The extracted barcode segment (after reverse complement) has the
        expected length; and
      * It matches a design barcode exactly (0 mismatches).

    These reads are considered barcode-assigned via exact matching and are
    used for downstream minigene comparison.

- bc_fuzzy
    Number of reads where:

      * The barcode segment (after reverse complement) has the expected
        length; and
      * It does NOT match any design barcode exactly, but
      * There is a unique design barcode within Hamming distance
        <= max_bc_mismatches (default 2).

    These reads are rescued by fuzzy matching and treated as if they had the
    canonical barcode they matched.

- bc_ambiguous
    Number of reads where:

      * The barcode segment (after reverse complement) has the expected
        length; and
      * The barcode sequence falls within the allowed mismatch radius of
        more than one design barcode.

    In this case, the correct assignment cannot be decided unambiguously,
    so these reads are not assigned to any barcode and are excluded from
    minigene comparison.

- bc_no_match
    Number of reads where:

      * The barcode segment (after reverse complement) has the expected
        length; and
      * It is not an exact match to any design barcode; and
      * It is not within the allowed mismatch radius of any barcode (i.e.
        no valid fuzzy match is found).

    These barcodes are considered off-design or too divergent and are
    discarded.

- bc_length_mismatch
    Number of reads where the extracted barcode segment length does not match
    the barcode length from the design table. These reads are treated as
    having invalid barcodes and are not used for minigene comparison.

- pair_valid
    Number of reads where:

      * The barcode is successfully assigned (either bc_exact or bc_fuzzy),
        AND
      * The minigene segment (after reverse complement) has Hamming
        distance <= max_mg_mismatches (default 9) to the expected minigene
        sequence for that assigned barcode.

    This is the total number of reads that support a valid barcode–minigene
    pair under the configured mismatch thresholds.

- mg_too_many_mismatches
    Number of reads where:

      * The barcode is successfully assigned to a design barcode, BUT
      * The minigene segment (after reverse complement) differs from the
        expected minigene sequence for that barcode by more than
        max_mg_mismatches.

    These correspond to all reads counted as reads_seq_fail in the
    per-minigene file, aggregated across all barcodes/minigenes.
"""

import argparse
import csv
import gzip
import sys
from collections import Counter

DNA_ALPHABET = "ACGT"


def open_maybe_gz(path, mode="rt"):
    """Open plain or gzipped file; '-' means stdin/stdout."""
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def revcomp(seq: str) -> str:
    """Reverse-complement a DNA sequence."""
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def hamming(a: str, b: str) -> int:
    """Hamming distance (length differences contribute as mismatches)."""
    if len(a) != len(b):
        m = min(len(a), len(b))
        return sum(1 for i in range(m) if a[i] != b[i]) + abs(len(a) - len(b))
    return sum(1 for x, y in zip(a, b) if x != y)


def generate_barcode_variants(seq: str, max_mismatches: int):
    """
    Generate all sequences within <= max_mismatches of seq.
    Only safe for short barcodes (e.g. 8 nt).
    """
    variants = {seq}
    n = len(seq)

    if max_mismatches >= 1:
        for i in range(n):
            for base in DNA_ALPHABET:
                if base != seq[i]:
                    variants.add(seq[:i] + base + seq[i + 1 :])

    if max_mismatches >= 2:
        for i in range(n):
            for j in range(i + 1, n):
                for b1 in DNA_ALPHABET:
                    if b1 == seq[i]:
                        continue
                    for b2 in DNA_ALPHABET:
                        if b2 == seq[j]:
                            continue
                        s_list = list(seq)
                        s_list[i] = b1
                        s_list[j] = b2
                        variants.add("".join(s_list))

    return variants


def load_design(path):
    """
    Load design table with columns:
      - sequence           (barcode, forward)
      - sequence_minigene  (minigene, forward)
    Optionally an 'id' or 'name' column for reporting.
    """
    with open(path, "r") as f:
        sample = f.readline()
        if not sample:
            raise ValueError("Design file appears to be empty.")
        f.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,;")
        reader = csv.DictReader(f, dialect=dialect)

        design = []
        for row in reader:
            if not row:
                continue
            try:
                bc = row["sequence"].strip().upper()
                mg = row["sequence_minigene"].strip().upper()
            except KeyError as e:
                raise KeyError(
                    "Design file must have columns 'sequence' (barcode) and "
                    "'sequence_minigene' (minigene); missing: " + str(e)
                )
            if not bc or not mg:
                continue
            design.append(
                {
                    "barcode": bc,
                    "minigene": mg,
                    "id": row.get("id", row.get("name", bc)),
                }
            )

    if not design:
        raise ValueError("Design file contains no usable rows.")

    # Check consistency of lengths
    bc_len = len(design[0]["barcode"])
    mg_len = len(design[0]["minigene"])
    for rec in design:
        if len(rec["barcode"]) != bc_len:
            raise ValueError("Inconsistent barcode lengths in design file.")
        if len(rec["minigene"]) != mg_len:
            raise ValueError("Inconsistent minigene lengths in design file.")

    return design, bc_len, mg_len


def build_barcode_maps(design, max_bc_mismatches: int):
    """
    Build:
      - barcode_to_design: exact barcode -> design record
      - variant_map: barcode variant (<= mism) -> canonical barcode or None if ambiguous
      - ambiguous_variants: set of variants shared by multiple barcodes
    """
    barcode_to_design = {}
    for rec in design:
        bc = rec["barcode"]
        if bc in barcode_to_design:
            raise ValueError(f"Duplicate barcode in design: {bc}")
        barcode_to_design[bc] = rec

    if max_bc_mismatches < 0:
        max_bc_mismatches = 0

    variant_map = {}
    ambiguous_variants = set()

    for rec in design:
        bc = rec["barcode"]
        variants = generate_barcode_variants(bc, max_bc_mismatches)
        for v in variants:
            if v in variant_map:
                if variant_map[v] != bc:
                    # conflict: mark as ambiguous
                    variant_map[v] = None
                    ambiguous_variants.add(v)
            else:
                variant_map[v] = bc

    return barcode_to_design, variant_map, ambiguous_variants


def classify_barcode(
    seq: str,
    barcode_to_design,
    variant_map,
    ambiguous_variants,
    bc_len: int,
    max_bc_mismatches: int,
):
    """
    Return (canonical_barcode, mismatches, status) where status is one of:
      {"exact", "fuzzy", "ambiguous", "no_match", "length_mismatch"}
    """
    seq = seq.upper()
    if len(seq) != bc_len:
        return None, None, "length_mismatch"

    # exact
    if seq in barcode_to_design:
        return seq, 0, "exact"

    # fuzzy
    canonical = variant_map.get(seq)
    if canonical is None:
        if seq in ambiguous_variants:
            return None, None, "ambiguous"
        return None, None, "no_match"

    mism = hamming(seq, canonical)
    if mism <= max_bc_mismatches:
        return canonical, mism, "fuzzy"
    else:
        return None, None, "no_match"


def process_fastq(
    fastq_path: str,
    design_path: str,
    out_prefix: str,
    stop_len: int = 3,
    start_len: int = 3,
    max_bc_mismatches: int = 2,
    max_mg_mismatches: int = 9,
):
    # Load design + barcode fuzzy index
    design, bc_len, mg_len = load_design(design_path)
    barcode_to_design, variant_map, ambiguous_variants = build_barcode_maps(
        design, max_bc_mismatches
    )

    # Counters
    pair_counts = Counter()        # canonical_barcode -> reads where bc+mg both OK
    bc_assigned_counts = Counter() # canonical_barcode -> reads with barcode assigned
    seq_fail_counts = Counter()    # canonical_barcode -> reads with bc OK, mg too many mism
    status_counts = Counter()
    total_reads = 0

    min_read_len = bc_len + stop_len + mg_len + start_len

    with open_maybe_gz(fastq_path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().strip()
            plus = fh.readline()
            qual = fh.readline()
            if not qual:
                break

            total_reads += 1

            if len(seq) < min_read_len:
                status_counts["too_short"] += 1
                continue

            # Extract segments in read orientation
            bc_seq = seq[:bc_len]
            # stop_seq = seq[bc_len : bc_len + stop_len]  # unused but could be checked
            mg_seq = seq[bc_len + stop_len : bc_len + stop_len + mg_len]
            # start_seq = seq[bc_len + stop_len + mg_len : bc_len + stop_len + mg_len + start_len]

            # Reverse complement barcode and minigene to match design orientation
            bc_seq = revcomp(bc_seq)
            mg_seq = revcomp(mg_seq)

            # Classify barcode (fuzzy <= max_bc_mismatches)
            canonical_bc, bc_mism, bc_status = classify_barcode(
                bc_seq,
                barcode_to_design,
                variant_map,
                ambiguous_variants,
                bc_len,
                max_bc_mismatches,
            )
            status_counts[f"bc_{bc_status}"] += 1

            # Without a barcode assignment, do not attempt minigene matching
            if canonical_bc is None:
                continue

            bc_assigned_counts[canonical_bc] += 1

            expected_mg = barcode_to_design[canonical_bc]["minigene"]
            mg_seq = mg_seq.upper()
            mg_mism = hamming(mg_seq, expected_mg)

            if mg_mism <= max_mg_mismatches:
                pair_counts[canonical_bc] += 1
                status_counts["pair_valid"] += 1
            else:
                seq_fail_counts[canonical_bc] += 1
                status_counts["mg_too_many_mismatches"] += 1

    # Per-minigene output
    per_path = f"{out_prefix}.per_minigene_counts.tsv"
    with open(per_path, "w") as out:
        out.write(
            "\t".join(
                [
                    "barcode",
                    "minigene_seq",
                    "design_id",
                    "reads_barcode_assigned",
                    "reads_valid_pair",
                    "reads_seq_fail",
                ]
            )
            + "\n"
        )
        for rec in design:
            bc = rec["barcode"]
            mg = rec["minigene"]
            did = rec["id"]
            out.write(
                "\t".join(
                    map(
                        str,
                        [
                            bc,
                            mg,
                            did,
                            bc_assigned_counts.get(bc, 0),
                            pair_counts.get(bc, 0),
                            seq_fail_counts.get(bc, 0),
                        ],
                    )
                )
                + "\n"
            )

    # Global summary
    summary_path = f"{out_prefix}.summary.tsv"
    with open(summary_path, "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"total_reads\t{total_reads}\n")
        for k, v in sorted(status_counts.items()):
            out.write(f"{k}\t{v}\n")

    return per_path, summary_path


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Count minigene barcode–sequence matches with fuzzy matching.\n"
            "FASTQ layout (read): [barcode][STOP][minigene][START]\n"
            "Barcode and minigene from the read are reverse-complemented\n"
            "to match the design table orientation."
        )
    )
    ap.add_argument("--fastq", required=True, help="Input FASTQ (optionally .gz)")
    ap.add_argument(
        "--design",
        required=True,
        help=(
            "Design table (e.g. DART2OS_pat001.txt) with columns "
            "'sequence' (barcode) and 'sequence_minigene' (minigene)"
        ),
    )
    ap.add_argument(
        "--out-prefix",
        required=True,
        help="Prefix for output files (TSV)",
    )
    ap.add_argument(
        "--stop-len",
        type=int,
        default=3,
        help="Length of STOP codon segment (default: 3)",
    )
    ap.add_argument(
        "--start-len",
        type=int,
        default=3,
        help="Length of START codon segment (default: 3)",
    )
    ap.add_argument(
        "--max-bc-mismatches",
        type=int,
        default=2,
        help="Max mismatches allowed in barcode (default: 2)",
    )
    ap.add_argument(
        "--max-mg-mismatches",
        type=int,
        default=9,
        help="Max mismatches allowed in minigene sequence (default: 9)",
    )

    args = ap.parse_args()

    process_fastq(
        fastq_path=args.fastq,
        design_path=args.design,
        out_prefix=args.out_prefix,
        stop_len=args.stop_len,
        start_len=args.start_len,
        max_bc_mismatches=args.max_bc_mismatches,
        max_mg_mismatches=args.max_mg_mismatches,
    )


if __name__ == "__main__":
    main()

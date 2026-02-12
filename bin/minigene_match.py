#!/usr/bin/env python3
"""
Minigene matching with:
- barcode fuzzy match (<= max_bc_mismatches, default 2)
- minigene mismatch sweep (default: 0,1,2,4,6,9,12)
- STOP codon check (default expected: TAA)
- START codon check (default expected: ATG)
- empty plasmid detection with <= empty_plasmid_max_mismatches (default 9)
  performed BEFORE any barcode/minigene assignment

Input FASTQ layout (READ ORIENTATION):
    [barcode][STOP][minigene][START]

Important:
- The read segments are assumed to be reverse-complements of the design orientation.
- Therefore, barcode/STOP/minigene/START are reverse-complemented (segment-wise) before matching.

Design file (e.g. DART2OS_pat001.txt) must contain:
    - sequence            : barcode (forward)
    - sequence_minigene   : minigene sequence (forward)

Outputs:
1) <prefix>.per_minigene_counts.tsv
   - One row per design entry plus an extra EMPTY_PLASMID row.
   - Columns include barcode assignment counts, START/STOP codon match counts+%, and
     reads_valid_pair_mmX / reads_seq_fail_mmX for each mismatch level.

2) <prefix>.summary.tsv
   - Columns: metric, value, pct_total_reads
   - Includes empty plasmid metric (empty_plasmid_mm9 by default),
     bc_* status counts, codon match counts, and mismatch sweep metrics.
"""

import argparse
import csv
import gzip
import sys
from collections import Counter

DNA_ALPHABET = "ACGT"

DEFAULT_EMPTY_PLASMID_SEQ = (
    "GGGATGTCGAAGAAAATCCTGGCCCAATGCGAGACGCTATCAGTCAGCAGCAGGTGAGAGTATACAAAGTGGCTTTATAATTAGCTTGGTCCATACGTAACGTCTCC"
).upper()


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


def parse_int_list(csv_like: str):
    """Parse a comma-separated list of non-negative integers."""
    parts = [p.strip() for p in csv_like.split(",") if p.strip() != ""]
    vals = []
    for p in parts:
        try:
            vals.append(int(p))
        except ValueError:
            raise ValueError(f"Invalid integer in --mg-mismatch-levels: '{p}'")
    vals = sorted(set(vals))
    if not vals:
        raise ValueError("--mg-mismatch-levels must contain at least one integer")
    if any(v < 0 for v in vals):
        raise ValueError("--mg-mismatch-levels cannot contain negative values")
    return vals


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

    if seq in barcode_to_design:
        return seq, 0, "exact"

    canonical = variant_map.get(seq)
    if canonical is None:
        if seq in ambiguous_variants:
            return None, None, "ambiguous"
        return None, None, "no_match"

    mism = hamming(seq, canonical)
    if mism <= max_bc_mismatches:
        return canonical, mism, "fuzzy"

    return None, None, "no_match"


def safe_pct(numer: int, denom: int) -> float:
    return 0.0 if denom <= 0 else 100.0 * float(numer) / float(denom)


def process_fastq(
    fastq_path: str,
    design_path: str,
    out_prefix: str,
    stop_len: int,
    start_len: int,
    max_bc_mismatches: int,
    mg_mismatch_levels,
    expected_stop: str,
    expected_start: str,
    empty_plasmid_seq: str,
    empty_plasmid_max_mismatches: int,
):
    design, bc_len, mg_len = load_design(design_path)
    barcode_to_design, variant_map, ambiguous_variants = build_barcode_maps(
        design, max_bc_mismatches
    )

    mg_levels = sorted(set(mg_mismatch_levels))

    expected_stop = expected_stop.upper()
    expected_start = expected_start.upper()
    empty_plasmid_seq = empty_plasmid_seq.upper()

    if len(expected_stop) != stop_len:
        raise ValueError(
            f"--expected-stop length ({len(expected_stop)}) must equal --stop-len ({stop_len})"
        )
    if len(expected_start) != start_len:
        raise ValueError(
            f"--expected-start length ({len(expected_start)}) must equal --start-len ({start_len})"
        )

    expected_len = bc_len + stop_len + mg_len + start_len
    if len(empty_plasmid_seq) != expected_len:
        raise ValueError(
            f"Empty plasmid sequence length is {len(empty_plasmid_seq)} but expected {expected_len} "
            f"(barcode_len {bc_len} + stop_len {stop_len} + minigene_len {mg_len} + start_len {start_len})."
        )

    # Empty plasmid counters
    empty_plasmid_count = 0
    empty_stop_ok = 0
    empty_start_ok = 0
    empty_both_ok = 0

    # Per-barcode counters
    bc_assigned_counts = Counter()

    # Codon QC per barcode (among barcode-assigned reads)
    stop_ok_counts = Counter()
    start_ok_counts = Counter()
    both_ok_counts = Counter()

    # Minigene mismatch sweeps
    pair_counts = {mm: Counter() for mm in mg_levels}
    seq_fail_counts = {mm: Counter() for mm in mg_levels}

    # Global counters
    status_counts = Counter()
    pair_valid_global = {mm: 0 for mm in mg_levels}
    mg_fail_global = {mm: 0 for mm in mg_levels}
    total_reads = 0

    # Codon QC global across non-empty reads (excluding too_short and empty plasmid reads)
    global_stop_ok = 0
    global_start_ok = 0
    global_both_ok = 0

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

            if len(seq) < expected_len:
                status_counts["too_short"] += 1
                continue

            # Extract segments in READ orientation
            bc_read = seq[:bc_len]
            stop_read = seq[bc_len : bc_len + stop_len]
            mg_read = seq[bc_len + stop_len : bc_len + stop_len + mg_len]
            start_read = seq[bc_len + stop_len + mg_len : expected_len]

            # Reverse complement segment-wise to DESIGN orientation
            bc_fwd = revcomp(bc_read).upper()
            stop_fwd = revcomp(stop_read).upper()
            mg_fwd = revcomp(mg_read).upper()
            start_fwd = revcomp(start_read).upper()

            full_fwd = start_fwd + mg_fwd + stop_fwd + bc_fwd

            # Empty plasmid check (<= 9 mismatches by default), BEFORE any barcode/minigene assignment
            empty_mism = hamming(full_fwd, empty_plasmid_seq)
            if empty_mism <= empty_plasmid_max_mismatches:
                empty_plasmid_count += 1

                # Codon stats for EMPTY_PLASMID row only
                stop_ok = (stop_fwd == expected_stop)
                start_ok = (start_fwd == expected_start)
                both_ok = stop_ok and start_ok
                if stop_ok:
                    empty_stop_ok += 1
                if start_ok:
                    empty_start_ok += 1
                if both_ok:
                    empty_both_ok += 1

                continue

            # Codon checks for non-empty reads (global & per-barcode if assigned)
            stop_ok = (stop_fwd == expected_stop)
            start_ok = (start_fwd == expected_start)
            both_ok = stop_ok and start_ok

            if stop_ok:
                global_stop_ok += 1
            if start_ok:
                global_start_ok += 1
            if both_ok:
                global_both_ok += 1

            # Barcode assignment
            canonical_bc, bc_mism, bc_status = classify_barcode(
                bc_fwd,
                barcode_to_design,
                variant_map,
                ambiguous_variants,
                bc_len,
                max_bc_mismatches,
            )
            status_counts[f"bc_{bc_status}"] += 1

            if canonical_bc is None:
                continue

            bc_assigned_counts[canonical_bc] += 1

            # Codon checks per barcode (among barcode-assigned reads)
            if stop_ok:
                stop_ok_counts[canonical_bc] += 1
            if start_ok:
                start_ok_counts[canonical_bc] += 1
            if both_ok:
                both_ok_counts[canonical_bc] += 1

            # Minigene mismatch distance
            expected_mg = barcode_to_design[canonical_bc]["minigene"]
            mg_mism = hamming(mg_fwd, expected_mg)

            for mm in mg_levels:
                if mg_mism <= mm:
                    pair_counts[mm][canonical_bc] += 1
                    pair_valid_global[mm] += 1
                else:
                    seq_fail_counts[mm][canonical_bc] += 1
                    mg_fail_global[mm] += 1

    # Write per-minigene output (design rows + EMPTY_PLASMID row)
    per_path = f"{out_prefix}.per_minigene_counts.tsv"
    with open(per_path, "w") as out:
        header = [
            "barcode",
            "minigene_seq",
            "design_id",
            "reads_barcode_assigned",
            "stop_codon_matches",
            "stop_codon_match_pct",
            "start_codon_matches",
            "start_codon_match_pct",
            "start_stop_both_matches",
            "start_stop_both_match_pct",
        ]
        for mm in mg_levels:
            header.append(f"reads_valid_pair_mm{mm}")
            header.append(f"reads_seq_fail_mm{mm}")
        out.write("\t".join(header) + "\n")

        # Normal design entries
        for rec in design:
            bc = rec["barcode"]
            mg = rec["minigene"]
            did = rec["id"]

            assigned = bc_assigned_counts.get(bc, 0)
            stop_m = stop_ok_counts.get(bc, 0)
            start_m = start_ok_counts.get(bc, 0)
            both_m = both_ok_counts.get(bc, 0)

            row = [
                bc,
                mg,
                did,
                str(assigned),
                str(stop_m),
                f"{safe_pct(stop_m, assigned):.6f}",
                str(start_m),
                f"{safe_pct(start_m, assigned):.6f}",
                str(both_m),
                f"{safe_pct(both_m, assigned):.6f}",
            ]
            for mm in mg_levels:
                row.append(str(pair_counts[mm].get(bc, 0)))
                row.append(str(seq_fail_counts[mm].get(bc, 0)))
            out.write("\t".join(row) + "\n")

        # Append EMPTY_PLASMID row
        empty_assigned = empty_plasmid_count
        empty_row = [
            "EMPTY_PLASMID",
            empty_plasmid_seq,  # full 107nt sequence in design orientation
            "EMPTY_PLASMID",
            str(empty_assigned),
            str(empty_stop_ok),
            f"{safe_pct(empty_stop_ok, empty_assigned):.6f}",
            str(empty_start_ok),
            f"{safe_pct(empty_start_ok, empty_assigned):.6f}",
            str(empty_both_ok),
            f"{safe_pct(empty_both_ok, empty_assigned):.6f}",
        ]
        # minigene mismatch sweep columns are not applicable for empty plasmid; keep as 0
        for mm in mg_levels:
            empty_row.append("0")
            empty_row.append("0")
        out.write("\t".join(empty_row) + "\n")

    # Write summary with % relative to total_reads
    summary_path = f"{out_prefix}.summary.tsv"
    with open(summary_path, "w") as out:
        out.write("metric\tvalue\tpct_total_reads\n")

        def write_metric(metric: str, value: int):
            pct = safe_pct(value, total_reads)
            out.write(f"{metric}\t{value}\t{pct:.6f}\n")

        out.write(f"total_reads\t{total_reads}\t{(100.0 if total_reads > 0 else 0.0):.6f}\n")

        write_metric("too_short", status_counts.get("too_short", 0))
        write_metric(f"empty_plasmid_mm{empty_plasmid_max_mismatches}", empty_plasmid_count)

        # barcode statuses (only for non-empty reads)
        for m in ["bc_exact", "bc_fuzzy", "bc_ambiguous", "bc_no_match", "bc_length_mismatch"]:
            write_metric(m, status_counts.get(m, 0))

        # codon QC across non-empty reads (excluding too_short and empty plasmid reads)
        write_metric("stop_codon_match", global_stop_ok)
        write_metric("start_codon_match", global_start_ok)
        write_metric("start_stop_both_match", global_both_ok)

        for mm in mg_levels:
            write_metric(f"pair_valid_mm{mm}", pair_valid_global[mm])
            write_metric(f"mg_too_many_mismatches_mm{mm}", mg_fail_global[mm])

    return per_path, summary_path


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Minigene matching with barcode fuzzy matching, multi-threshold minigene matching, "
            "START/STOP codon checks, and empty plasmid detection (<= mismatches)."
        )
    )
    ap.add_argument("--fastq", required=True, help="Input FASTQ (optionally .gz)")
    ap.add_argument(
        "--design",
        required=True,
        help="Design file with columns: 'sequence' (barcode) and 'sequence_minigene' (minigene)",
    )
    ap.add_argument("--out-prefix", required=True, help="Prefix for output files")

    ap.add_argument("--stop-len", type=int, default=3, help="STOP length (default: 3)")
    ap.add_argument("--start-len", type=int, default=3, help="START length (default: 3)")
    ap.add_argument(
        "--expected-stop",
        type=str,
        default="TAA",
        help="Expected STOP codon after reverse-complementing (default: TAA)",
    )
    ap.add_argument(
        "--expected-start",
        type=str,
        default="ATG",
        help="Expected START codon after reverse-complementing (default: ATG)",
    )

    ap.add_argument(
        "--empty-plasmid-seq",
        type=str,
        default=DEFAULT_EMPTY_PLASMID_SEQ,
        help="Expected empty plasmid full 107nt sequence in design orientation (default: built-in sequence)",
    )
    ap.add_argument(
        "--empty-plasmid-max-mismatches",
        type=int,
        default=9,
        help="Max mismatches allowed for empty plasmid detection (default: 9)",
    )

    ap.add_argument(
        "--max-bc-mismatches",
        type=int,
        default=2,
        help="Max mismatches allowed in barcode (default: 2)",
    )
    ap.add_argument(
        "--mg-mismatch-levels",
        type=str,
        default="0,1,2,4,6,9,12",
        help="Comma-separated minigene mismatch thresholds to evaluate (default: 0,1,2,4,6,9,12)",
    )

    args = ap.parse_args()
    mg_levels = parse_int_list(args.mg_mismatch_levels)

    process_fastq(
        fastq_path=args.fastq,
        design_path=args.design,
        out_prefix=args.out_prefix,
        stop_len=args.stop_len,
        start_len=args.start_len,
        max_bc_mismatches=args.max_bc_mismatches,
        mg_mismatch_levels=mg_levels,
        expected_stop=args.expected_stop,
        expected_start=args.expected_start,
        empty_plasmid_seq=args.empty_plasmid_seq,
        empty_plasmid_max_mismatches=args.empty_plasmid_max_mismatches,
    )


if __name__ == "__main__":
    main()

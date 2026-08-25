#!/usr/bin/env python3

import argparse
import re
import sys


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
MD_RE = re.compile(r"(\d+|\^[A-Z]+|[A-Z])")


ALIGNMENT_TAG_PREFIXES_TO_REMOVE = {
    "MD",
    "NM",
    "AS",
    "XS",
    "XN",
    "XM",
    "XO",
    "XG",
    "YT",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Mark SAM alignments as unmapped if they fail positional mismatch "
            "rules: first X read bases allow N mismatches, remaining read bases "
            "allow M mismatches."
        )
    )

    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Input SAM file. Use '-' for stdin.",
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="Output SAM file. Use '-' for stdout.",
    )
    parser.add_argument(
        "--prefix-length",
        "-x",
        type=int,
        required=True,
        help="Number of leading read bases to treat as constrained prefix.",
    )
    parser.add_argument(
        "--max-prefix-mismatches",
        type=int,
        default=1,
        help="Maximum mismatches allowed in the first X read bases. Default: 1.",
    )
    parser.add_argument(
        "--max-suffix-mismatches",
        type=int,
        default=6,
        help="Maximum mismatches allowed after the first X read bases. Default: 6.",
    )
    parser.add_argument(
        "--fail-missing-md",
        action="store_true",
        help=(
            "If an aligned record lacks an MD tag, mark it unmapped. "
            "Default: keep records without MD unchanged."
        ),
    )
    parser.add_argument(
        "--allow-indels",
        action="store_true",
        help=(
            "Allow alignments with insertions or deletions. "
            "Default: mark alignments containing I or D in the CIGAR as unmapped."
        ),
    )
    parser.add_argument(
        "--count-soft-clipped-in-read-position",
        action="store_true",
        help=(
            "Count soft-clipped bases as read positions when deciding where the "
            "first X bases are. Default: ignore soft-clipped bases."
        ),
    )
    parser.add_argument(
        "--keep-alignment-tags-on-unmapped",
        action="store_true",
        help=(
            "Keep optional alignment tags such as MD, NM, AS on records that are "
            "marked unmapped. Default: remove common alignment-specific tags."
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print summary statistics to stderr.",
    )

    return parser.parse_args()


def open_maybe_stdin(path):
    if path == "-":
        return sys.stdin
    return open(path, "r")


def open_maybe_stdout(path):
    if path == "-":
        return sys.stdout
    return open(path, "w")


def parse_cigar(cigar):
    if cigar == "*":
        return []

    parsed = []
    consumed = 0

    for length, op in CIGAR_RE.findall(cigar):
        length = int(length)
        parsed.append((length, op))
        consumed += len(str(length)) + 1

    if consumed != len(cigar):
        raise ValueError(f"Could not fully parse CIGAR string: {cigar}")

    return parsed


def has_indel(cigar_ops):
    return any(op in {"I", "D"} for _, op in cigar_ops)


def get_md_tag(fields):
    for field in fields[11:]:
        if field.startswith("MD:Z:"):
            return field[5:]
    return None


def build_ref_to_read_map(cigar_ops, count_soft_clipped=False):
    """
    Build a mapping from reference-consuming alignment positions to read positions.

    The returned list has one entry per reference-consuming aligned base.
    Each entry is either:
      - read position, 0-based, for M/= /X operations
      - None, for D/N operations

    CIGAR operation behavior:
      M, =, X consume read and reference
      I, S consume read only
      D, N consume reference only
      H, P consume neither
    """
    ref_to_read = []
    read_pos = 0

    for length, op in cigar_ops:
        if op in {"M", "=", "X"}:
            for _ in range(length):
                ref_to_read.append(read_pos)
                read_pos += 1

        elif op == "I":
            read_pos += length

        elif op in {"D", "N"}:
            for _ in range(length):
                ref_to_read.append(None)

        elif op == "S":
            if count_soft_clipped:
                read_pos += length
            else:
                pass

        elif op in {"H", "P"}:
            pass

        else:
            raise ValueError(f"Unsupported CIGAR operation: {op}")

    return ref_to_read


def mismatch_read_positions_from_md(md, ref_to_read):
    """
    Parse the MD tag and return read positions of mismatches.

    MD examples:
      10A5
      5^AC3T10
      0A0C10

    MD describes reference positions. Insertions are not represented in MD.
    """
    tokens = MD_RE.findall(md)

    consumed = "".join(tokens)
    if consumed != md:
        raise ValueError(f"Could not fully parse MD tag: {md}")

    ref_pos = 0
    mismatch_read_positions = []

    for token in tokens:
        if token.isdigit():
            ref_pos += int(token)

        elif token.startswith("^"):
            deleted_bases = token[1:]
            ref_pos += len(deleted_bases)

        else:
            if ref_pos >= len(ref_to_read):
                raise ValueError(
                    f"MD tag references position outside CIGAR-derived alignment. "
                    f"MD={md}, ref_pos={ref_pos}, ref_to_read_len={len(ref_to_read)}"
                )

            read_pos = ref_to_read[ref_pos]

            if read_pos is not None:
                mismatch_read_positions.append(read_pos)

            ref_pos += 1

    return mismatch_read_positions


def count_prefix_suffix_mismatches(mismatch_positions, prefix_length):
    prefix_mismatches = 0
    suffix_mismatches = 0

    for read_pos in mismatch_positions:
        if read_pos < prefix_length:
            prefix_mismatches += 1
        else:
            suffix_mismatches += 1

    return prefix_mismatches, suffix_mismatches


def is_unmapped(flag):
    return bool(flag & 0x4)


def remove_alignment_specific_tags(optional_fields):
    kept = []

    for field in optional_fields:
        tag = field.split(":", 1)[0]

        if tag in ALIGNMENT_TAG_PREFIXES_TO_REMOVE:
            continue

        kept.append(field)

    return kept


def mark_record_unmapped(fields, keep_alignment_tags=False):
    """
    Convert an aligned SAM record to an unmapped SAM record while preserving
    QNAME, FLAG metadata, SEQ, QUAL, and non-alignment optional tags.
    """
    flag = int(fields[1])

    # Add unmapped bit.
    flag |= 0x4

    # If this is paired-end data, keep mate information only if meaningful.
    # For single-end data, these are standard unmapped values.
    fields[1] = str(flag)
    fields[2] = "*"
    fields[3] = "0"
    fields[4] = "0"
    fields[5] = "*"

    # Mate fields. For your shown module this is likely single-end, so use unmapped defaults.
    fields[6] = "*"
    fields[7] = "0"
    fields[8] = "0"

    if not keep_alignment_tags and len(fields) > 11:
        fields[11:] = remove_alignment_specific_tags(fields[11:])

    return fields


def record_passes_mismatch_rules(fields, args):
    flag = int(fields[1])
    cigar = fields[5]

    if is_unmapped(flag):
        return True, "already_unmapped"

    md = get_md_tag(fields)

    if md is None:
        if args.fail_missing_md:
            return False, "missing_md"
        else:
            return True, "missing_md_kept"

    cigar_ops = parse_cigar(cigar)

    if not args.allow_indels and has_indel(cigar_ops):
        return False, "indel"

    ref_to_read = build_ref_to_read_map(
        cigar_ops,
        count_soft_clipped=args.count_soft_clipped_in_read_position,
    )

    mismatch_positions = mismatch_read_positions_from_md(md, ref_to_read)

    prefix_mismatches, suffix_mismatches = count_prefix_suffix_mismatches(
        mismatch_positions,
        args.prefix_length,
    )

    if prefix_mismatches > args.max_prefix_mismatches:
        return False, "prefix_mismatches"

    if suffix_mismatches > args.max_suffix_mismatches:
        return False, "suffix_mismatches"

    return True, "passed"


def process_sam(args):
    total_records = 0
    kept_aligned = 0
    already_unmapped = 0
    marked_unmapped = 0

    failed_missing_md = 0
    failed_indel = 0
    failed_prefix = 0
    failed_suffix = 0
    failed_parse_error = 0
    missing_md_kept = 0

    with open_maybe_stdin(args.input) as in_fh, open_maybe_stdout(args.output) as out_fh:
        for line in in_fh:
            if line.startswith("@"):
                out_fh.write(line)
                continue

            line = line.rstrip("\n")

            if not line:
                continue

            total_records += 1
            fields = line.split("\t")

            if len(fields) < 11:
                # Malformed record. Keep unchanged rather than silently dropping it.
                out_fh.write(line + "\n")
                failed_parse_error += 1
                continue

            try:
                passed, reason = record_passes_mismatch_rules(fields, args)

            except Exception:
                # If parsing fails for an apparently aligned read, mark it unmapped.
                fields = mark_record_unmapped(
                    fields,
                    keep_alignment_tags=args.keep_alignment_tags_on_unmapped,
                )
                out_fh.write("\t".join(fields) + "\n")
                marked_unmapped += 1
                failed_parse_error += 1
                continue

            if passed:
                out_fh.write("\t".join(fields) + "\n")

                if reason == "already_unmapped":
                    already_unmapped += 1
                elif reason == "missing_md_kept":
                    missing_md_kept += 1
                    kept_aligned += 1
                else:
                    kept_aligned += 1

            else:
                if reason == "missing_md":
                    failed_missing_md += 1
                elif reason == "indel":
                    failed_indel += 1
                elif reason == "prefix_mismatches":
                    failed_prefix += 1
                elif reason == "suffix_mismatches":
                    failed_suffix += 1

                fields = mark_record_unmapped(
                    fields,
                    keep_alignment_tags=args.keep_alignment_tags_on_unmapped,
                )
                out_fh.write("\t".join(fields) + "\n")
                marked_unmapped += 1

    if args.verbose:
        print(f"Total SAM records: {total_records}", file=sys.stderr)
        print(f"Kept aligned: {kept_aligned}", file=sys.stderr)
        print(f"Already unmapped: {already_unmapped}", file=sys.stderr)
        print(f"Marked unmapped: {marked_unmapped}", file=sys.stderr)
        print(f"Failed due to missing MD: {failed_missing_md}", file=sys.stderr)
        print(f"Missing MD kept unchanged: {missing_md_kept}", file=sys.stderr)
        print(f"Failed due to indel: {failed_indel}", file=sys.stderr)
        print(f"Failed due to prefix mismatches: {failed_prefix}", file=sys.stderr)
        print(f"Failed due to suffix mismatches: {failed_suffix}", file=sys.stderr)
        print(f"Failed due to parse error: {failed_parse_error}", file=sys.stderr)


def main():
    args = parse_args()

    if args.prefix_length < 0:
        raise ValueError("--prefix-length must be >= 0")

    if args.max_prefix_mismatches < 0:
        raise ValueError("--max-prefix-mismatches must be >= 0")

    if args.max_suffix_mismatches < 0:
        raise ValueError("--max-suffix-mismatches must be >= 0")

    process_sam(args)


if __name__ == "__main__":
    main()

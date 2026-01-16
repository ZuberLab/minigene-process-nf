#!/usr/bin/env python3
"""
Combine minigene_match per-sample summary files into one table.

Input:
  - One or more summary files as positional arguments.
    Each file is a TSV with columns:
        metric  value

Output:
  - A single TSV to stdout with:
        sample  <metric1>  <metric1_pct>  <metric2>  <metric2_pct> ...

    One row per sample. Percentages are computed as:
        metric_pct = 100 * metric_value / total_reads

    Metrics missing in a sample are reported as 0 and 0.0%.
"""

import sys
import os
import csv
import gzip
from collections import OrderedDict


def open_maybe_gz(path, mode="rt"):
    """Open plain or gzipped file."""
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def infer_sample_name(path: str) -> str:
    """
    Infer sample name from a summary file path.

    Examples:
      /path/to/sample1.summary.tsv   -> sample1
      /path/to/sampleA.summary.txt   -> sampleA
      /path/to/foo.bar.summary.tsv   -> foo.bar.summary (keeps dots)
      /path/to/sample1.txt           -> sample1
    """
    base = os.path.basename(path)
    # strip typical suffixes used in this pipeline
    for suffix in (".summary.tsv", ".summary.txt"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    # fall back: remove only the final extension
    name, _ext = os.path.splitext(base)
    return name


def read_summary_file(path: str):
    """
    Read a single summary file.
    Returns:
      sample_name (str),
      metrics (dict: metric -> int)
    """
    sample = infer_sample_name(path)
    metrics = {}
    with open_maybe_gz(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "metric" not in reader.fieldnames or "value" not in reader.fieldnames:
            raise ValueError(f"{path}: expected columns 'metric' and 'value'")
        for row in reader:
            m = row["metric"]
            v = row["value"]
            if m is None or m == "":
                continue
            try:
                val = int(v)
            except ValueError:
                # if it's not an int, you can either skip or raise
                raise ValueError(f"{path}: value for metric '{m}' is not an integer: {v}")
            metrics[m] = val
    return sample, metrics


def main():
    if len(sys.argv) < 2:
        sys.stderr.write(
            "Usage: combine_minigene_match_summary.py summary1.tsv [summary2.tsv ...]\n"
        )
        sys.exit(1)

    summary_paths = sys.argv[1:]

    # Read all summaries
    sample_to_metrics = OrderedDict()  # sample -> {metric: value}
    all_metrics = set()

    for path in summary_paths:
        sample, metrics = read_summary_file(path)
        if sample in sample_to_metrics:
            raise ValueError(f"Duplicate sample name inferred: {sample} (file: {path})")
        sample_to_metrics[sample] = metrics
        all_metrics.update(metrics.keys())

    if "total_reads" not in all_metrics:
        raise ValueError("No 'total_reads' metric found in any summary file.")

    # Define column order: total_reads first, then all other metrics sorted
    metrics_sorted = ["total_reads"] + sorted(
        m for m in all_metrics if m != "total_reads"
    )

    # Write header
    header = ["sample"]
    for m in metrics_sorted:
        header.append(m)
        header.append(f"{m}_pct")
    out = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    out.writerow(header)

    # Write one row per sample
    for sample, metrics in sample_to_metrics.items():
        total_reads = metrics.get("total_reads", 0)
        row = [sample]
        for m in metrics_sorted:
            count = metrics.get(m, 0)
            # avoid division by zero; if total_reads == 0, define % as 0.0
            if total_reads > 0:
                pct = 100.0 * count / total_reads
            else:
                pct = 0.0
            row.append(str(count))
            # format percentage with 3 decimal places
            row.append(f"{pct:.2f}")
        out.writerow(row)


if __name__ == "__main__":
    main()

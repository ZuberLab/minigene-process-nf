#!/usr/bin/env python3
"""
Combine minigene_match per-sample summary files into one table and make a saturation plot.

Accepted per-sample summary formats:
A) metric   value   pct_total_reads
B) metric   value  (pct computed from total_reads)

Output:
- Combined wide TSV to stdout:
    sample  <metric> <metric_pct> ...

- Saturation plot via --plot:
    x-axis: mismatch allowance (from pair_valid_mmX)
    y-axis: pair_valid_mmX_pct (% of total reads)
"""

import sys
import os
import csv
import gzip
import re
from collections import OrderedDict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PAIR_VALID_RE = re.compile(r"^pair_valid_mm(\d+)$")


def open_maybe_gz(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def infer_sample_name(path: str) -> str:
    base = os.path.basename(path)
    for suffix in (".summary.tsv", ".summary.txt"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    name, _ext = os.path.splitext(base)
    return name


def parse_number(s: str):
    s = s.strip()
    if s == "":
        return 0
    try:
        return int(s)
    except ValueError:
        return float(s)


def fmt_value(v):
    if isinstance(v, int):
        return str(v)
    if abs(v - round(v)) < 1e-12:
        return str(int(round(v)))
    return f"{v:.6f}"


def read_summary_file(path: str):
    sample = infer_sample_name(path)
    metrics_value = {}
    metrics_pct = {}

    with open_maybe_gz(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or "metric" not in reader.fieldnames or "value" not in reader.fieldnames:
            raise ValueError(f"{path}: expected at least columns 'metric' and 'value'")

        has_pct = "pct_total_reads" in reader.fieldnames

        for row in reader:
            m = row.get("metric", "").strip()
            if not m:
                continue
            v = parse_number(row.get("value", "0"))
            metrics_value[m] = v

            if has_pct:
                p = float(row.get("pct_total_reads", "0") or 0.0)
                metrics_pct[m] = p

    if "total_reads" not in metrics_value:
        raise ValueError(f"{path}: missing required metric 'total_reads'")

    total_reads = float(metrics_value["total_reads"]) if metrics_value["total_reads"] else 0.0
    if total_reads <= 0:
        for m in metrics_value.keys():
            metrics_pct.setdefault(m, 0.0)
    else:
        for m, v in metrics_value.items():
            if m not in metrics_pct:
                metrics_pct[m] = 100.0 * float(v) / total_reads

    return sample, metrics_value, metrics_pct


def metric_order(all_metrics):
    core = [
        "total_reads",
        "too_short",
        # new empty plasmid metric name used by minigene_match.py (default threshold 9):
        "empty_plasmid_mm9",
        "bc_exact",
        "bc_fuzzy",
        "bc_ambiguous",
        "bc_no_match",
        "bc_length_mismatch",
        "stop_codon_match",
        "start_codon_match",
        "start_stop_both_match",
    ]
    present_core = [m for m in core if m in all_metrics]

    mm_levels = sorted(
        {int(PAIR_VALID_RE.match(m).group(1)) for m in all_metrics if PAIR_VALID_RE.match(m)}
    )

    pair_valid = [f"pair_valid_mm{mm}" for mm in mm_levels if f"pair_valid_mm{mm}" in all_metrics]
    mg_fail = [
        f"mg_too_many_mismatches_mm{mm}"
        for mm in mm_levels
        if f"mg_too_many_mismatches_mm{mm}" in all_metrics
    ]

    used = set(present_core + pair_valid + mg_fail)
    rest = sorted(m for m in all_metrics if m not in used)

    return present_core + pair_valid + mg_fail + rest, mm_levels


def write_combined_table(sample_to_value, sample_to_pct, metrics_sorted):
    header = ["sample"]
    for m in metrics_sorted:
        header.append(m)
        header.append(f"{m}_pct")

    out = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    out.writerow(header)

    for sample in sample_to_value.keys():
        row = [sample]
        for m in metrics_sorted:
            v = sample_to_value[sample].get(m, 0)
            p = sample_to_pct[sample].get(m, 0.0)
            row.append(fmt_value(v))
            row.append(f"{p:.6f}")
        out.writerow(row)


def make_saturation_plot(sample_to_pct, mm_levels, plot_path):
    if not mm_levels:
        sys.stderr.write("WARNING: No pair_valid_mmX metrics found; skipping plot.\n")
        return

    x = mm_levels
    ys = []
    labels = []

    for sample, pct_map in sample_to_pct.items():
        y = []
        for mm in mm_levels:
            y.append(float(pct_map.get(f"pair_valid_mm{mm}", 0.0)))
        ys.append(y)
        labels.append(sample)

    mean_y = []
    for j in range(len(x)):
        vals = [ys[i][j] for i in range(len(ys))]
        mean_y.append(sum(vals) / len(vals) if vals else 0.0)

    plt.figure()
    for y, label in zip(ys, labels):
        plt.plot(x, y, marker="o", label=label)
    plt.plot(x, mean_y, marker="o", linestyle="--", linewidth=2, label="mean")

    plt.xlabel("Max allowed mismatches in minigene sequence")
    plt.ylabel("pair_valid_pct (% of total reads)")
    plt.title("Minigene matching saturation (pair_valid_pct vs mismatch allowance)")
    plt.xticks(x)

    if len(labels) <= 15:
        plt.legend()

    plt.tight_layout()
    plt.savefig(plot_path, dpi=200)
    plt.close()


def main():
    if len(sys.argv) < 2:
        sys.stderr.write(
            "Usage: combine_minigene_match_summary.py [--plot out.png] summary1.tsv [summary2.tsv ...]\n"
        )
        sys.exit(1)

    args = sys.argv[1:]
    plot_path = None

    if len(args) >= 2 and args[0] == "--plot":
        plot_path = args[1]
        args = args[2:]

    summary_paths = args
    if not summary_paths:
        sys.stderr.write("ERROR: No summary files provided.\n")
        sys.exit(1)

    sample_to_value = OrderedDict()
    sample_to_pct = OrderedDict()
    all_metrics = set()

    for path in summary_paths:
        sample, metrics_value, metrics_pct = read_summary_file(path)
        if sample in sample_to_value:
            raise ValueError(f"Duplicate sample name inferred: {sample} (file: {path})")
        sample_to_value[sample] = metrics_value
        sample_to_pct[sample] = metrics_pct
        all_metrics.update(metrics_value.keys())

    metrics_sorted, mm_levels = metric_order(all_metrics)

    write_combined_table(sample_to_value, sample_to_pct, metrics_sorted)

    if plot_path is not None:
        make_saturation_plot(sample_to_pct, mm_levels, plot_path)


if __name__ == "__main__":
    main()

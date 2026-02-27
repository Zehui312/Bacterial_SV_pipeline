#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pysam
from collections import defaultdict
import csv
import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


# ====== argument parsing ======
def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Analyze breakpoints in BAM file. Supports breakpoints on the same contig "
            "or on different contigs."
        )
    )

    parser.add_argument("--bam", required=True, help="Path to BAM file")

    # Backward compatible: old scripts used --contig for both A and B
    parser.add_argument(
        "--contig",
        required=False,
        help="Contig/chromosome name (legacy: used for both A and B unless --contig-a/--contig-b are provided)",
    )
    parser.add_argument("--contig-a", required=False, help="Contig/chromosome for breakpoint A")
    parser.add_argument("--contig-b", required=False, help="Contig/chromosome for breakpoint B")

    parser.add_argument("--breakpoint-a", type=int, required=True, help="Position of breakpoint A")
    parser.add_argument("--breakpoint-b", type=int, required=True, help="Position of breakpoint B")

    parser.add_argument("--slop", type=int, default=5, help="Slop distance for breakpoint detection (default: 5)")
    parser.add_argument(
        "--span-cutoff",
        type=int,
        default=0,
        help=(
            "Minimum bp on EACH side of the breakpoint to count as spanning "
            "(0 = no cutoff, default: 0)"
        ),
    )
    parser.add_argument("--sample-name", required=True, help="Sample name for output files")
    parser.add_argument(
        "--margin",
        type=int,
        default=2000,
        help="Fetch margin around each breakpoint when loading reads (default: 2000)",
    )

    args = parser.parse_args()

    # Resolve contigs
    contig_a = args.contig_a or args.contig
    contig_b = args.contig_b or args.contig

    if not contig_a or not contig_b:
        parser.error("You must provide either --contig (for both breakpoints) or both --contig-a and --contig-b.")

    args.contig_a = contig_a
    args.contig_b = contig_b
    return args


args = parse_args()
bam_path = args.bam
contig_a = args.contig_a
contig_b = args.contig_b
A = args.breakpoint_a
B = args.breakpoint_b
slop_num = args.slop
span_cutoff = args.span_cutoff
sample_name = args.sample_name
margin = args.margin

# Output names
if contig_a == contig_b:
    out_tag = f"{contig_a}-{A}_{B}"
else:
    out_tag = f"{contig_a}-{A}__{contig_b}-{B}"

csv_out = f"{sample_name}_{out_tag}.csv"
fig_out = f"{sample_name}_{out_tag}.png"
# ===========================


def load_alignments_near_breakpoints(bam_path, contig_a, A, contig_b, B, margin=2000):
    """
    Load alignments from BAM around each breakpoint region, potentially on different contigs.
    Group them by read (query_name).

    If contig_a == contig_b, fetch a single region spanning both breakpoints (plus margin).
    Otherwise, fetch two separate regions (one per contig) and merge reads by query_name.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    read_dict = defaultdict(list)

    def add_fetch(contig, pos):
        region_start = max(0, pos - margin)
        region_end = pos + margin
        for aln in bam.fetch(contig, region_start, region_end):
            if aln.is_unmapped:
                continue
            if aln.reference_end is None:
                continue
            read_dict[aln.query_name].append(aln)

    if contig_a == contig_b:
        region_start = max(0, min(A, B) - margin)
        region_end = max(A, B) + margin
        for aln in bam.fetch(contig_a, region_start, region_end):
            if aln.is_unmapped:
                continue
            if aln.reference_end is None:
                continue
            read_dict[aln.query_name].append(aln)
    else:
        add_fetch(contig_a, A)
        add_fetch(contig_b, B)

    bam.close()
    return read_dict


def classify_reads_for_breakpoint(read_dict, contig, P, slop, span_cutoff):
    """
    For a given breakpoint position P on `contig`:

    - spanning_reads:
        reads with an alignment that covers P (start <= P < end) ON THIS CONTIG and
        have at least `span_cutoff` bp on EACH side of P:
            left_len  = P - start
            right_len = end - P
            min(left_len, right_len) >= span_cutoff

      If span_cutoff == 0, then any read that covers P is counted as spanning.

    - containing_reads:
        split/chimeric reads where ANY alignment ON THIS CONTIG starts or ends within `slop`
        bp of the breakpoint AND the read has >=2 alignments total (can be on same contig
        or different contigs). This supports inter-contig breakpoints.
    """
    spanning_reads = set()
    containing_reads = set()

    for qname, alns in read_dict.items():
        # All valid alignments for this read
        alns_valid = [a for a in alns if (not a.is_unmapped) and (a.reference_end is not None)]
        if not alns_valid:
            continue

        # Alignments on the target contig
        alns_on_contig = [a for a in alns_valid if a.reference_name == contig]
        if not alns_on_contig:
            continue

        # ----- spanning reads -----
        for a in alns_on_contig:
            start = a.reference_start
            end = a.reference_end  # 0-based, end exclusive
            if start <= P < end:
                if span_cutoff > 0:
                    left_len = P - start
                    right_len = end - P
                    if min(left_len, right_len) < span_cutoff:
                        continue
                spanning_reads.add(qname)
                break

        # ----- chimeric / breakpoint-containing -----
        if len(alns_valid) < 2:
            continue

        for a in alns_on_contig:
            start = a.reference_start
            end = a.reference_end
            if abs(start - P) <= slop or abs(end - P) <= slop:
                containing_reads.add(qname)
                break

    return containing_reads, spanning_reads


def save_csv(summary, csv_out):
    """
    summary: list of dicts with keys:
      ['breakpoint_label', 'contig', 'position', 'reads_containing', 'reads_spanning']
    """
    fieldnames = ["breakpoint_label", "contig", "position", "reads_containing", "reads_spanning"]
    with open(csv_out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary:
            writer.writerow(row)
    print(f"[INFO] CSV saved to: {csv_out}")


def make_plot(summary, fig_out):
    labels = [f'{row["breakpoint_label"]}\n{row["contig"]}:{row["position"]}' for row in summary]
    spanning = [row["reads_spanning"] for row in summary]
    containing = [row["reads_containing"] for row in summary]

    x = range(len(labels))

    total_spanning = sum(spanning)
    total_chimeric = sum(containing)
    total_reads = total_spanning + total_chimeric
    chim_pct = (total_chimeric / total_reads) * 100 if total_reads > 0 else 0.0

    plt.figure(figsize=(6.0, 5))
    plt.bar(x, spanning, label="Spanning reads")
    plt.bar(x, containing, bottom=spanning, label="Chimeric reads")

    plt.xticks(x, labels, fontsize=12)
    plt.yticks(fontsize=14)
    plt.ylabel("Read Counts", fontsize=18)
    plt.xlabel("Breakpoint", fontsize=18)

    plt.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.15),
        ncol=2,
        frameon=False,
    )

    pct_text = f"SV percentage: {total_chimeric}/{total_reads} = {chim_pct:.2f}%"
    plt.text(
        0.5,
        1.08,
        pct_text,
        ha="center",
        va="bottom",
        transform=plt.gca().transAxes,
        fontsize=12,
    )

    plt.gca().yaxis.set_major_locator(MaxNLocator(5))

    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(2)

    plt.tight_layout()
    plt.savefig(fig_out, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[INFO] Figure saved to: {fig_out}")
    print(f"[INFO] Chimeric reads percentage: {chim_pct:.2f}%")


def main():
    read_dict = load_alignments_near_breakpoints(
        bam_path, contig_a, A, contig_b, B, margin=margin
    )

    contain_A, span_A = classify_reads_for_breakpoint(
        read_dict, contig_a, A, slop=slop_num, span_cutoff=span_cutoff
    )
    contain_B, span_B = classify_reads_for_breakpoint(
        read_dict, contig_b, B, slop=slop_num, span_cutoff=span_cutoff
    )

    print(f"BAM file: {bam_path}\n")

    print(f"Breakpoint A = {contig_a}:{A}")
    print(f"  Reads containing A breakpoint (chimeric around A): {len(contain_A)}")
    print(f"  Reads spanning A (continuous alignment covers A): {len(span_A)}\n")

    print(f"Breakpoint B = {contig_b}:{B}")
    print(f"  Reads containing B breakpoint (chimeric around B): {len(contain_B)}")
    print(f"  Reads spanning B (continuous alignment covers B): {len(span_B)}\n")

    summary = [
        {
            "breakpoint_label": "A",
            "contig": contig_a,
            "position": A,
            "reads_containing": len(contain_A),
            "reads_spanning": len(span_A),
        },
        {
            "breakpoint_label": "B",
            "contig": contig_b,
            "position": B,
            "reads_containing": len(contain_B),
            "reads_spanning": len(span_B),
        },
    ]

    save_csv(summary, csv_out)
    make_plot(summary, fig_out)


if __name__ == "__main__":
    main()
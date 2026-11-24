#!/usr/bin/env python3

import argparse
import csv
import json
from collections import Counter


def save_json(filename, data):
    with open(filename, "w") as file:
        json.dump(data, file, indent=2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("primer")
    parser.add_argument("db")
    args = parser.parse_args()

    # Read OBI-PCR TSV
    with open(args.input) as file:
        rows = list(csv.DictReader(file, delimiter="\t"))

    prefix = f"{args.primer}_{args.db}"
    name = f"{args.primer} ({args.db})"

    # Amplicon recovery
    number_of_amplicons = len(rows)

    amplicon_recovery = {
        "id": "obipcr_amplicon_recovery",
        "parent_id": "obipcr",
        "parent_name": "OBI-PCR",
        "section_name": "OBI-PCR Amplicon Recovery",
        "description": "Number of in-silico PCR amplicons recovered for each primer and database.",
        "plot_type": "bargraph",
        "pconfig": {
            "id": "obipcr_amplicon_recovery_plot",
            "title": "OBI-PCR: Amplicon Recovery",
            "ylab": "Recovered amplicons",
            "cpswitch": False
        },
        "data": {
            name: {
                "Amplicons": number_of_amplicons
            }
        }
    }

    save_json(
        f"{prefix}_obipcr_amplicon_recovery_mqc.json",
        amplicon_recovery
    )

    if len(rows) == 0:
        return

    # Get values from TSV
    lengths = []
    forward_errors = []
    reverse_errors = []
    forward_positions = []
    reverse_positions = []
    forward_binding_gc = []
    reverse_binding_gc = []
    amplicon_gc_values = []

    for row in rows:
        lengths.append(int(row["Amplicon_Length"]))
        forward_errors.append(int(row["Forward_Errors"]))
        reverse_errors.append(int(row["Reverse_Errors"]))
        forward_binding_gc.append(float(row["Forward_Match_GC"]))
        reverse_binding_gc.append(float(row["Reverse_Match_GC"]))
        amplicon_gc_values.append(float(row["Amplicon_GC"]))

        if row["Forward_Mismatch_Positions_Primer"] != "None":
            positions = row["Forward_Mismatch_Positions_Primer"].split(",")
            for position in positions:
                forward_positions.append(int(position))

        if row["Reverse_Mismatch_Positions_Primer"] != "None":
            positions = row["Reverse_Mismatch_Positions_Primer"].split(",")
            for position in positions:
                reverse_positions.append(int(position))

    # Perfect primer binding
    perfect = 0

    for i in range(len(rows)):
        if forward_errors[i] == 0 and reverse_errors[i] == 0:
            perfect += 1

    perfect_percentage = perfect / len(rows) * 100

    # Most common mismatch positions
    forward_counts = Counter(forward_positions)
    reverse_counts = Counter(reverse_positions)

    most_common_forward = forward_counts.most_common(1)
    most_common_reverse = reverse_counts.most_common(1)

    if most_common_forward:
        most_common_forward = most_common_forward[0][0]
    else:
        most_common_forward = "None"

    if most_common_reverse:
        most_common_reverse = most_common_reverse[0][0]
    else:
        most_common_reverse = "None"

    # Summary table
    summary = {
        "id": "obipcr_summary",
        "parent_id": "obipcr",
        "parent_name": "OBI-PCR",
        "section_name": "OBI-PCR Summary",
        "plot_type": "table",
        "data": {
            name: {
                "Database": args.db,
                "Amplicons": len(rows),
                "Mean length (bp)": round(sum(lengths) / len(lengths), 1),
                "Min length (bp)": min(lengths),
                "Max length (bp)": max(lengths),
                "Perfect binding (%)": round(perfect_percentage, 1),
                "Forward primer GC (%)": float(rows[0]["Forward_Primer_GC"]),
                "Reverse primer GC (%)": float(rows[0]["Reverse_Primer_GC"]),
                "Most common FWD mismatch": most_common_forward,
                "Most common REV mismatch": most_common_reverse
            }
        }
    }

    save_json(
        f"{prefix}_obipcr_summary_mqc.json",
        summary
    )

    # Primer and binding-site GC
    gc_data = {
        "id": "obipcr_gc",
        "parent_id": "obipcr",
        "parent_name": "OBI-PCR",
        "section_name": "OBI-PCR GC Content",
        "plot_type": "table",
        "data": {
            name: {
                "Forward primer": float(rows[0]["Forward_Primer_GC"]),
                "Reverse primer": float(rows[0]["Reverse_Primer_GC"]),
                "Forward binding site": round(
                    sum(forward_binding_gc) / len(forward_binding_gc), 2
                ),
                "Reverse binding site": round(
                    sum(reverse_binding_gc) / len(reverse_binding_gc), 2
                )
            }
        }
    }

    save_json(
        f"{prefix}_obipcr_gc_mqc.json",
        gc_data
    )

    # Amplicon GC distribution
    gc_counts = Counter()

    for gc in amplicon_gc_values:
        gc_counts[round(gc)] += 1

    amplicon_gc = {
        "id": "obipcr_amplicon_gc",
        "parent_id": "obipcr",
        "parent_name": "OBI-PCR",
        "section_name": "OBI-PCR Amplicon GC",
        "plot_type": "linegraph",
        "pconfig": {
            "xlab": "GC (%)",
            "ylab": "Amplicons"
        },
        "data": {
            name: [
                [gc, gc_counts[gc]]
                for gc in sorted(gc_counts)
            ]
        }
    }

    save_json(
        f"{prefix}_obipcr_amplicon_gc_mqc.json",
        amplicon_gc
    )

    # Amplicon length distribution
    length_counts = Counter(lengths)

    amplicon_lengths = {
        "id": "obipcr_amplicon_lengths",
        "parent_id": "obipcr",
        "parent_name": "OBI-PCR",
        "section_name": "OBI-PCR Amplicon Lengths",
        "plot_type": "linegraph",
        "pconfig": {
            "xlab": "Amplicon length (bp)",
            "ylab": "Amplicons"
        },
        "data": {
            name: [
                [length, length_counts[length]]
                for length in sorted(length_counts)
            ]
        }
    }

    save_json(
        f"{prefix}_obipcr_amplicon_lengths_mqc.json",
        amplicon_lengths
    )

    # Number of mismatches
    maximum_error = max(forward_errors + reverse_errors)

    errors = {
        "id": "obipcr_errors",
        "parent_id": "obipcr",
        "parent_name": "OBI-PCR",
        "section_name": "OBI-PCR Primer Mismatches",
        "plot_type": "bargraph",
        "pconfig": {
            "ylab": "Amplicons"
        },
        "data": {
            f"{name} Forward": {
                str(error): forward_errors.count(error)
                for error in range(maximum_error + 1)
            },
            f"{name} Reverse": {
                str(error): reverse_errors.count(error)
                for error in range(maximum_error + 1)
            }
        }
    }

    save_json(
        f"{prefix}_obipcr_errors_mqc.json",
        errors
    )

    # Mismatch positions
    forward_primer_length = len(rows[0]["Forward_Primer"])
    reverse_primer_length = len(rows[0]["Reverse_Primer"])

    positions = {
        "id": "obipcr_positions",
        "parent_id": "obipcr",
        "parent_name": "OBI-PCR",
        "section_name": "OBI-PCR Mismatch Positions",
        "plot_type": "linegraph",
        "pconfig": {
            "xlab": "Position from 5' end",
            "ylab": "Mismatch count"
        },
        "data": {
            f"{name} Forward": [
                [position, forward_counts[position]]
                for position in range(1, forward_primer_length + 1)
            ],
            f"{name} Reverse": [
                [position, reverse_counts[position]]
                for position in range(1, reverse_primer_length + 1)
            ]
        }
    }

    save_json(
        f"{prefix}_obipcr_positions_mqc.json",
        positions
    )


if __name__ == "__main__":
    main()
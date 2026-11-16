#!/usr/bin/env python3

import argparse
import csv
from collections import Counter, defaultdict
from itertools import chain

CONSENSUS_COLUMNS = [
    "cluster_id",
    "accession",
    "accession_taxid",
    "accession_name",
    "consensus_name",
    "consensus_taxid",
    "consensus_rank",
    "disambiguation",
]

REPORT_CATEGORIES = [
    "Species",
    "Genus",
    "Family",
    "Order",
    "Class",
    "Phylum",
    "Kingdom",
    "Unclassified",
    "Other",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create MultiQC summaries from cluster consensus results."
    )
    parser.add_argument("--input", required=True, help="Input cluster_consensus.tsv file.")
    parser.add_argument("--output-prefix", required=True, help="Prefix for output files.")
    parser.add_argument(
        "--sample",
        required=True,
        help="Sample name shown in MultiQC, e.g. ITS2 (NCBI).",
    )
    parser.add_argument(
        "--cluster-column",
        default="cluster_id",
        help="Column containing cluster IDs.",
    )
    parser.add_argument(
        "--rank-column",
        default="consensus_rank",
        help="Column containing consensus taxonomic ranks.",
    )
    return parser.parse_args()


def rank_category(rank):
    """Collapse detailed taxonomic ranks into four reporting groups."""
    rank = rank.lower().strip()

    species_ranks = {
        "strain",
        "isolate",
        "forma specialis",
        "forma",
        "varietas",
        "subspecies",
        "species",
    }
    genus_ranks = {
        "species subgroup",
        "species group",
        "series",
        "section",
        "subgenus",
        "genus",
    }
    family_ranks = {
        "subtribe",
        "tribe",
        "subfamily",
        "family",
    }
    order_ranks = {
        "superfamily",
        "parvorder",
        "infraorder",
        "suborder",
        "order",
    }
    class_ranks = {
        "superorder",
        "subcohort",
        "cohort",
        "infraclass",
        "subclass",
        "class",
    }
    phylum_ranks = {
        "superclass",
        "subphylum",
        "phylum",
    }
    kingdom_ranks = {
        "superphylum",
        "subkingdom",
        "kingdom",
        "superkingdom",
    }

    if rank in species_ranks:
        return "Species"
    if rank in genus_ranks:
        return "Genus"
    if rank in family_ranks:
        return "Family"
    if rank in order_ranks:
        return "Order"
    if rank in class_ranks:
        return "Class"
    if rank in phylum_ranks:
        return "Phylum"
    if rank in kingdom_ranks:
        return "Kingdom"
    if rank in {"no rank", "unranked"}:
        return "Unclassified"
    return "Other"


def read_clusters(input_file, cluster_column, rank_column):
    """Read either the pipeline's headerless TSV or a TSV with named columns."""
    cluster_sizes = Counter()
    cluster_ranks = defaultdict(set)

    with open(input_file, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        first_row = next(reader, None)

        if first_row is None:
            raise ValueError("Input file is empty.")

        if cluster_column in first_row and rank_column in first_row:
            cluster_index = first_row.index(cluster_column)
            rank_index = first_row.index(rank_column)
            rows = reader
            first_data_row = 2
        else:
            missing_columns = {
                cluster_column,
                rank_column,
            } - set(CONSENSUS_COLUMNS)
            if missing_columns:
                raise ValueError(
                    f"Missing required columns: {sorted(missing_columns)}\n"
                    f"Available columns: {CONSENSUS_COLUMNS}"
                )

            cluster_index = CONSENSUS_COLUMNS.index(cluster_column)
            rank_index = CONSENSUS_COLUMNS.index(rank_column)
            rows = chain([first_row], reader)
            first_data_row = 1

        required_index = max(cluster_index, rank_index)
        for row_number, row in enumerate(rows, start=first_data_row):
            if len(row) <= required_index:
                raise ValueError(
                    f"Row {row_number} has {len(row)} columns; "
                    f"expected at least {required_index + 1}."
                )

            cluster_id = row[cluster_index].strip()
            rank = row[rank_index].strip().lower() or "unranked"

            if not cluster_id:
                continue

            cluster_sizes[cluster_id] += 1
            cluster_ranks[cluster_id].add(rank)

    final_cluster_ranks = {}
    for cluster_id, ranks in cluster_ranks.items():
        if len(ranks) != 1:
            raise ValueError(
                f"Cluster {cluster_id} has multiple consensus ranks: {sorted(ranks)}"
            )
        final_cluster_ranks[cluster_id] = next(iter(ranks))

    return cluster_sizes, final_cluster_ranks


def calculate_taxonomic_resolution(cluster_ranks):
    rank_counts = Counter(rank_category(rank) for rank in cluster_ranks.values())
    total_clusters = len(cluster_ranks)

    if total_clusters == 0:
        raise ValueError("No clusters were found in the input file.")

    percentages = {
        category: rank_counts[category] / total_clusters * 100
        for category in REPORT_CATEGORIES
    }
    return rank_counts, percentages


def size_category(size):
    if size == 1:
        return "1"
    if size <= 5:
        return "2-5"
    if size <= 10:
        return "6-10"
    if size <= 50:
        return "11-50"
    return ">50"


def calculate_cluster_sizes(cluster_sizes):
    return Counter(size_category(size) for size in cluster_sizes.values())


def write_mqc_tsv(filename, config_lines, header, values):
    with open(filename, "w", newline="") as handle:
        handle.writelines(f"# {line}\n" for line in config_lines)

        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerow(values)


def write_resolution_report(prefix, sample, rank_counts):
    config = [
        'id: "cluster_consensus_resolution"',
        'parent_id: "cluster_consensus"',
        'parent_name: "Cluster Consensus"',
        'section_name: "Overall Taxonomic Resolution"',
        (
            'description: "Number and percentage of clusters resolving at each '
            'taxonomic level, plus unclassified and other ranks. Use the plot '
            'toggle to switch between percentages and cluster counts."'
        ),
        'plot_type: "bargraph"',
        'pconfig:',
        '  id: "cluster_consensus_resolution_plot"',
        '  title: "Cluster Consensus: Taxonomic Resolution"',
        '  ylab: "Clusters"',
        '  cpswitch: true',
        '  cpswitch_c_active: false',
        '  cpswitch_counts_label: "Number of clusters"',
        '  cpswitch_percent_label: "Percentage of clusters"',
    ]
    header = ["Primer / Database", *REPORT_CATEGORIES]
    values = [sample] + [
        rank_counts[category]
        for category in REPORT_CATEGORIES
    ]
    write_mqc_tsv(
        f"{prefix}.taxonomic_resolution_mqc.tsv", config, header, values
    )


def write_summary_report(prefix, sample, cluster_ranks, rank_counts):
    total_clusters = len(cluster_ranks)

    def count_with_percentage(category):
        count = rank_counts[category]
        percentage = count / total_clusters * 100
        return f"{count} ({percentage:.1f}%)"

    config = [
        'id: "cluster_consensus_summary"',
        'parent_id: "cluster_consensus"',
        'parent_name: "Cluster Consensus"',
        'section_name: "Overall Summary"',
        'description: "Summary of cluster-level taxonomic resolution."',
        'plot_type: "table"',
        'pconfig:',
        '  id: "cluster_consensus_summary_table"',
        '  title: "Cluster Consensus: Overall Summary"',
    ]
    header = [
        "Primer / Database",
        "Clusters",
        *[f"{category} clusters" for category in REPORT_CATEGORIES],
    ]
    values = [
        sample,
        total_clusters,
        *[
            count_with_percentage(category)
            for category in REPORT_CATEGORIES
        ],
    ]
    write_mqc_tsv(f"{prefix}.overall_summary_mqc.tsv", config, header, values)


def write_cluster_size_report(prefix, sample, size_counts):
    config = [
        'id: "cluster_consensus_size"',
        'parent_id: "cluster_consensus"',
        'parent_name: "Cluster Consensus"',
        'section_name: "Cluster Size"',
        'description: "Distribution of the number of sequences per cluster."',
        'plot_type: "bargraph"',
        'pconfig:',
        '  id: "cluster_consensus_size_plot"',
        '  title: "Cluster Consensus: Cluster Size"',
        '  ylab: "Number of clusters"',
    ]
    header = ["Primer / Database", "1 sequence", "2-5", "6-10", "11-50", ">50"]
    values = [
        sample,
        size_counts["1"],
        size_counts["2-5"],
        size_counts["6-10"],
        size_counts["11-50"],
        size_counts[">50"],
    ]
    write_mqc_tsv(f"{prefix}.cluster_size_mqc.tsv", config, header, values)


def main():
    args = parse_args()
    cluster_sizes, cluster_ranks = read_clusters(
        args.input, args.cluster_column, args.rank_column
    )
    rank_counts, percentages = calculate_taxonomic_resolution(cluster_ranks)
    size_counts = calculate_cluster_sizes(cluster_sizes)

    print(f"Sample: {args.sample}")
    print(f"Total clusters: {len(cluster_ranks)}")
    print(f"Taxonomic resolution counts: {dict(rank_counts)}")
    print(f"Taxonomic resolution percentages: {dict(percentages)}")
    print(f"Cluster size distribution: {dict(size_counts)}")

    write_resolution_report(args.output_prefix, args.sample, rank_counts)
    write_summary_report(args.output_prefix, args.sample, cluster_ranks, rank_counts)
    write_cluster_size_report(args.output_prefix, args.sample, size_counts)


if __name__ == "__main__":
    main()

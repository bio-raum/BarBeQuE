#!/usr/bin/env python3

import argparse
import os

import taxidTools


def consensus_fraction(value):
    value = float(value)
    if not 0.5 < value <= 1.0:
        raise argparse.ArgumentTypeError(
            "consensus fraction must be greater than 0.5 and at most 1.0"
        )
    return value


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Assign taxonomy to sequence clusters")
    parser.add_argument(
        "--input",
        required=True,
        help="Path to cluster_taxonomy.tsv from JOIN_ACCESSION_TAXONOMY",
    )
    parser.add_argument("--taxdump", required=True, help="Path to the taxdump folder")
    parser.add_argument("--output", required=True, help="Path to output table")
    parser.add_argument(
        "--min-consensus",
        type=consensus_fraction,
        default=1.0,
        help=(
            "Minimum fraction of valid cluster accessions that must support an "
            "assignment (default: 1.0, equivalent to strict LCA)"
        ),
    )
    return parser.parse_args(argv)


def load_taxonomy(taxdump):
    return taxidTools.read_taxdump(
        os.path.join(taxdump, "nodes.dmp"),
        os.path.join(taxdump, "rankedlineage.dmp"),
        os.path.join(taxdump, "merged.dmp")
    )


def parse_cluster_taxonomy(input_file):
    """Parse cluster_taxonomy.tsv (cluster_id, accession, taxid, phylum, class, order, family, genus, species)"""
    clusters_taxid = {}    # cluster_id -> list of taxids
    rows = []              # keep original rows for output

    with open(input_file, 'r') as fi:
        for line in fi:
            elements = line.rstrip('\n').split('\t')
            cluster_id = elements[0]
            taxid = elements[2]

            clusters_taxid.setdefault(cluster_id, [])
            clusters_taxid[cluster_id].append(taxid)

            rows.append(elements)

    return clusters_taxid, rows


def cluster_consensus(tax, taxids, min_consensus):
    """Return name, taxid, rank, and input-taxon names for one cluster."""
    unique_taxids = list(dict.fromkeys(taxids))
    disambiguation = ";".join(
        name
        for name in (tax.getName(taxid) for taxid in unique_taxids)
        if name is not None
    )
    valid_taxids = [taxid for taxid in taxids if tax.getName(taxid) is not None]

    if not valid_taxids:
        return ["Unclassified", "Unknown", "no rank", disambiguation]

    consensus = tax.consensus(
        valid_taxids,
        min_consensus=min_consensus,
        ignore_missing=True,
    )
    if consensus is None:
        return ["Unclassified", "Unknown", "no rank", disambiguation]

    return [
        consensus.name,
        str(consensus.taxid),
        consensus.rank,
        disambiguation,
    ]


def main(input_file, taxdump, output, min_consensus=1.0):
    tax = load_taxonomy(taxdump)

    # Parse input
    clusters_taxid, rows = parse_cluster_taxonomy(input_file)

    # Get consensus, rank, and disambiguation per cluster
    clusters_cons = {}
    for cluster_id, taxids in clusters_taxid.items():
        clusters_cons[cluster_id] = cluster_consensus(
            tax,
            taxids,
            min_consensus,
        )

    # Dump: original row + accession_name + LCA columns appended
    with open(output, 'w') as fo:
        for elements in rows:
            cluster_id = elements[0]
            accession_taxid = elements[2]
            node = tax.get(accession_taxid)
            accession_name = node.name if node else "Unknown"
            
            line = "\t".join(elements + [accession_name] + clusters_cons[cluster_id])
            fo.write(line + "\n")


if __name__ == "__main__":
    args = parse_args()
    main(
        args.input,
        args.taxdump,
        args.output,
        args.min_consensus,
    )

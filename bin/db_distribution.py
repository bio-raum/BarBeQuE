#!/usr/bin/env python3
"""
db_distribution.py — BarBeQuE DB composition table
Resolves taxids from the reference DB to their full lineage.

Usage:
    db_distribution.py --input <taxids_counts.tsv> --output <out.tsv> --taxdump <taxdump_dir>
"""

import argparse
import os

import taxidTools

RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

def read_taxid_counts(path):
    """Read taxid\tcount TSV, return {taxid: count}."""
    db_taxids = {}
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 2 or parts[0].lower() == "taxid":
                continue
            try:
                db_taxids[int(parts[0])] = int(parts[1])
            except ValueError:
                continue
    return db_taxids


def write_tsv(db_taxids, all_lineages, all_ranks, all_names, output_path):
    """Write one row per taxid with full lineage and sequence count."""
    with open(output_path, "w") as out:
        out.write("taxid\tcount\tresolved_rank\t" + "\t".join(RANKS) + "\n")

        for taxid, count in sorted(db_taxids.items()):
            lin = all_lineages.get(taxid, [])

            ranked_lin = [
                (tid, all_ranks.get(tid, ""), all_names.get(tid, "?"))
                for tid in lin
                if all_ranks.get(tid) in RANKS
            ]

            rank_map = {rank: name for _, rank, name in ranked_lin}
            my_rank  = all_ranks.get(taxid, "no rank")

            cols = (
                [str(taxid), str(count), my_rank]
                + [rank_map.get(r, "") for r in RANKS]
            )
            out.write("\t".join(cols) + "\n")


def resolve_lineages(tax, taxid_list):
    all_lineages = {}
    all_ranks = {}
    all_names = {}
    for tid in taxid_list:
        lin = []
        curr = tax.get(str(tid))
        if curr is None:
            all_lineages[tid] = []
            continue
            
        while curr:
            curr_tid = int(curr.taxid)
            lin.append(curr_tid)
            all_ranks[curr_tid] = curr.rank if curr.rank else "no rank"
            all_names[curr_tid] = curr.name if curr.name else str(curr_tid)
            
            if curr.parent and curr.parent.taxid != curr.taxid:
                curr = curr.parent
            else:
                curr = None
                
        all_lineages[tid] = lin[::-1] 
    
    return all_lineages, all_ranks, all_names


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",  required=True,  help="taxid counts TSV from BUILD_DB_TAXIDS")
    parser.add_argument("--output", required=True,  help="output TSV path")
    parser.add_argument("--taxdump", required=True, help="Path to the taxdump folder")
    args = parser.parse_args()

    tax = taxidTools.read_taxdump(
        os.path.join(args.taxdump, "nodes.dmp"),
        os.path.join(args.taxdump, "rankedlineage.dmp"),
        os.path.join(args.taxdump, "merged.dmp")
    )

    db_taxids = read_taxid_counts(args.input)
    all_lineages, all_ranks, all_names = resolve_lineages(tax, list(db_taxids.keys()))
    write_tsv(db_taxids, all_lineages, all_ranks, all_names, args.output)


if __name__ == "__main__":
    main()
#!/usr/bin/env python3

import argparse
import csv
from collections import Counter

from ete3 import NCBITaxa


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate taxonomic coverage and create MultiQC tables."
    )
    parser.add_argument("--taxon", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sample", required=True)
    return parser.parse_args()


def read_taxids(filename, column=None):
    taxids = set()

    with open(filename) as handle:
        for line in handle:
            if not line.strip():
                continue

            if column is None:
                taxid = line.strip()
            else:
                fields = line.rstrip("\n").split("\t")
                if len(fields) <= column:
                    continue
                taxid = fields[column].strip()

            taxids.add(taxid)

    return taxids


def main():
    args = parse_args()
    ncbi = NCBITaxa()

    matches = ncbi.get_name_translator([args.taxon])

    if args.taxon not in matches:
        raise ValueError(
            f"Taxon '{args.taxon}' was not found in the NCBI taxonomy."
        )

    root_taxid = matches[args.taxon][0]

    amplified_taxids = read_taxids(args.report, column=2)
    database_taxids = read_taxids(args.reference)

    tree = ncbi.get_descendant_taxa(
        root_taxid,
        collapse_subspecies=False,
        return_tree=True,
    )

    colors = {
        "PASS": "#7ee076",
        "FAIL": "#dfc2b1",
        "NOT IN DATABASE": "#eeeeee",
    }

    coverage = []
    species_nodes = []

    for node in tree.traverse():
        if node.rank != "species":
            continue

        taxid = str(node.name)
        taxon = node.sci_name

        if taxid in amplified_taxids:
            status = "PASS"
        elif taxid in database_taxids:
            status = "FAIL"
        else:
            status = "NOT IN DATABASE"

        coverage.append((taxon, status, taxid, colors[status]))
        node.name = taxon
        species_nodes.append(node)

    tree.prune(species_nodes)
    tree.write(outfile=args.output + ".nwk", format=1)

    with open(args.output + ".tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Taxon", "Status", "TaxID", "Color"])
        writer.writerows(coverage)

    counts = Counter(
        status for _taxon, status, _taxid, _color in coverage
    )

    with open(args.output + "_summary_mqc.tsv", "w", newline="") as handle:
        handle.write('# id: "taxonomic_coverage_summary"\n')
        handle.write('# parent_id: "taxonomic_coverage"\n')
        handle.write('# parent_name: "Taxonomic Coverage"\n')
        handle.write('# section_name: "Summary"\n')
        handle.write('# plot_type: "table"\n')

        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "Primer / Database",
            "Pass",
            "Fail",
            "Not in database",
        ])
        writer.writerow([
            args.sample,
            counts["PASS"],
            counts["FAIL"],
            counts["NOT IN DATABASE"],
        ])

    with open(args.output + "_details_mqc.tsv", "w", newline="") as handle:
        handle.write('# id: "taxonomic_coverage_details"\n')
        handle.write('# parent_id: "taxonomic_coverage"\n')
        handle.write('# parent_name: "Taxonomic Coverage"\n')
        handle.write('# section_name: "Taxon Details"\n')
        handle.write('# plot_type: "table"\n')
        handle.write('# pconfig:\n')
        handle.write('#   rows_are_samples: false\n')
        handle.write('#   flat_if_very_large: true\n')
        handle.write('# headers:\n')
        handle.write('#   TaxID:\n')
        handle.write('#     scale: false\n')

        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "Record",
            "Primer / Database",
            "Taxon",
            "TaxID",
            "Status",
        ])

        for taxon, status, taxid, _color in coverage:
            writer.writerow([
                f"{args.sample} | {taxid}",
                args.sample,
                taxon,
                taxid,
                status,
            ])

    print(f"Sample: {args.sample}")
    print(f"Taxon: {args.taxon}")
    print(f"Total species: {len(coverage)}")
    print(f"Pass: {counts['PASS']}")
    print(f"Fail: {counts['FAIL']}")
    print(f"Not in database: {counts['NOT IN DATABASE']}")


if __name__ == "__main__":
    main()

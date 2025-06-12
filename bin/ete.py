#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from ete3 import NCBITaxa
import re
import sys

parser = argparse.ArgumentParser(description="Script options")
parser.add_argument("--taxon", help="Taxonomy to use as root")
parser.add_argument("--reference", help="Taxonomy lists from the DB")
parser.add_argument("--report", help="The consensus report")
parser.add_argument("--output")
args = parser.parse_args()

ncbi = NCBITaxa()


def main(taxname, refs, report, output):

    # We pass a taxon name, and must translate it to an NCBI id
    name2taxid = ncbi.get_name_translator([taxname])

    # the first entry, of which the first value is taken
    taxid = list(name2taxid.items())[0][1][0]

    # the consensus report
    r = open(report, "r")
    lines = r.readlines()
    bucket = {}
    lines.pop()  # remove the header

    for line in lines:
        elements = line.split("\t")
        sci = elements[1]
        tax = elements[2]
        bucket[tax] = sci

    r.close()

    # the list taxids included in the BlastDB
    with open(refs, "r") as taxids:
        blast_tax = taxids.readlines()
        blast_tax = [tax.rstrip('\n') for tax in blast_tax]

    ok = "#7ee076"
    fail = "#dfc2b1"
    missing = "#eeeeee"
    tree = ncbi.get_descendant_taxa(taxid, collapse_subspecies=True, return_tree=True)

    data = {}

    nodes = []
    for n in tree.traverse():
        if n.is_leaf:
            # Make sure we only look at terminal nodes
            if n.rank == "species":
                tid = n.name
                n.name = n.sci_name

                # NCBI taxonomy is full of non-species level terminal leafs
                # we skip all leafs not matching the 'Genus species' pattern
                if re.match(r'^[A-Z][a-z]*\s[a-z]*$', n.sci_name):
                    # This taxon is in the blast db and was found
                    if tid in bucket:
                        data[n.sci_name] = f"OK\t{tid}\t{ok}"
                    # This taxon is in the blast db and was not found
                    elif tid in blast_tax:
                        data[n.sci_name] = f"FAIL\t{tid}\t{fail}"
                    # this taxon was not in the blast db
                    else:
                        data[n.sci_name] = f"NO_DATA\t{tid}\t{missing}"
                    nodes.append(n)
                else:
                    sys.stderr.write(f"Deleting {n.name}\n")
                    n.detach()

    tree.prune(nodes)

    print(tree.write(outfile=output + ".nwk", format=1))

    f = open(output + ".tsv", "w")
    f.write("Taxon\tStatus\tTaxid\tColor\n")
    for taxon in data:
        status = data[taxon]
        f.write(f"{taxon}\t{status}\n")

    f.close()


if __name__ == '__main__':
    main(args.taxon, args.reference, args.report, args.output)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import json

parser = argparse.ArgumentParser(description="Script options")
parser.add_argument("--input", help="A BarBeQue consensus file in tab format")
parser.add_argument("--primer", help="The primer name")
parser.add_argument("--database", help="The database name")
parser.add_argument("--output")
args = parser.parse_args()


def main(input, primer, db, output):

    data = {"id": "amplicon_size_histo_groups",
            "section_name": "Amplicon length by taxonomic group",
            "description": "Shows the lengths distribution of amplicons per taxonomic group",
            "plot_type": "linegraph",
            "pconfig": {"id": "amplicon_histo", "title": "Amplicon sizes", "ylab": "Frequency", "xlab": "Amplicon length"},
            "data": {}
            }

    histo = {}  # The basic histogram
    histo_normalized  = {}  # The histogram, each taxonomic class normalized to 100%
    all_lengths = []

    with open(input) as tsv:
        tsvreader = csv.DictReader(tsv, delimiter="\t")

        for entry in tsvreader:
            amlen = len(entry["amplicon"])
            tax_class = entry["class"]
            all_lengths.append(amlen)

            if tax_class in histo:

                if amlen in histo[tax_class]:
                    histo[tax_class][amlen] += 1
                else:
                    histo[tax_class][amlen] = 1

            else:

                histo[tax_class] = {amlen: 1}

    
    sorted(list(set(all_lengths)))
    print(all_lengths)

    for tclass, counts in histo.items():

        total_counts = 0

        for ampl_len, ampl_count in counts.items():
            all_lengths.append(ampl_len)

            total_counts += int(ampl_count)

        histo_normalized[tclass] = {}

        for this_len in all_lengths:

            if this_len in counts:
                ampl_count = counts[this_len]
            else:
                ampl_count = 0.0

            this_perc = 100 * round(int(ampl_count) / total_counts, 2)
            histo_normalized[tclass][this_len] = float(this_perc)

    data["data"] = histo

    with open(output, "w") as json_out:
        json.dump(data, json_out, indent=4, sort_keys=True)


if __name__ == '__main__':

    main(args.input, args.primer, args.database, args.output)

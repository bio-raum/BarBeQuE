#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import json
import statistics

parser = argparse.ArgumentParser(description="Script options")
parser.add_argument("--input", help="A BarBeQue consensus file in tab format")
parser.add_argument("--output")
args = parser.parse_args()


def main(input, output):

    data = {"id": "amplicon_mean_size",
            "section_name": "Amplicon mean length by taxonomic group",
            "description": "Mean length and std-dev of amplicons per taxonomic group",
            "plot_type": "table",
            "pconfig": {"id": "amp_mean_size", "col1_header": "Taxon"},
            "data": {}
            }

    histo = {}

    with open(input) as tsv:
        tsvreader = csv.DictReader(tsv, delimiter="\t")

        for entry in tsvreader:
            amlen = len(entry["amplicon"])
            tax_class = entry["class"]

            if tax_class in histo:

                if amlen in histo[tax_class]:
                    histo[tax_class][amlen] += 1
                else:
                    histo[tax_class][amlen] = 1

            else:
                histo[tax_class] = {amlen: 1}

    for tgroup, counts in histo.items():

        values = []
        for this_len, this_count in counts.items():
            for i in range(this_count):
                values.append(this_len)

        if len(values) > 10:
            mean = statistics.mean(values)
            stddev = statistics.stdev(values)

            data["data"][tgroup] = {"Mean": mean, "Stddev": stddev}

    with open(output, "w") as json_out:
        json.dump(data, json_out, indent=4, sort_keys=True)


if __name__ == '__main__':

    main(args.input, args.output)

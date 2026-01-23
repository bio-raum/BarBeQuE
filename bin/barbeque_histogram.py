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

    data = {"id": "amplicon_size_histo",
            "section_name": "Amplicon length frequency",
            "description": "Shows the frequency of amplicon lengths for the chosen primer and database",
            "plot_type": "linegraph",
            "pconfig": {"id": "amplicon_histo", "title": "Amplicon sizes", "ylab": "Frequency", "xlab": "Amplicon length"},
            "data": {}
            }

    histo = {}

    with open(input) as tsv:
        tsvreader = csv.DictReader(tsv, delimiter="\t")

        for entry in tsvreader:
            amlen = len(entry["amplicon"])

            if amlen in histo:
                histo[amlen] += 1
            else:
                histo[amlen] = 1

    data["data"][f"{primer}_{db}"] = histo

    with open(output, "w") as json_out:
        json.dump(data, json_out, indent=4, sort_keys=True)


if __name__ == '__main__':

    main(args.input, args.primer, args.database, args.output)

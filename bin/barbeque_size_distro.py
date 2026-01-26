#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import json
import statistics
import subprocess

parser = argparse.ArgumentParser(description="Script options")
parser.add_argument("--input", help="A BarBeQue consensus file in tab format")
parser.add_argument("--database", help="The Crabs database")
parser.add_argument("--output")
args = parser.parse_args()


def main(input, database, output):

    data = {"id": "amplicon_mean_size",
            "section_name": "Amplicon mean length by taxonomic group",
            "description": "Mean length and std-dev of amplicons per taxonomic group",
            "plot_type": "table",
            "pconfig": {"id": "amp_mean_size", "col1_header": "Taxon"},
            "data": {}
            }

    histo = {}
    tax_count = {}

    with open(input) as tsv:
        tsvreader = csv.DictReader(tsv, delimiter="\t")

        for entry in tsvreader:
            amlen = len(entry["amplicon"])
            tax_class = entry["class"]
            species = entry["species"]

            if tax_class in histo:

                if amlen in histo[tax_class]:
                    histo[tax_class][amlen] += 1
                else:
                    histo[tax_class][amlen] = 1

            else:
                histo[tax_class] = {amlen: 1}

            if tax_class in tax_count:
                tax_count[tax_class].append(species)
            else:
                tax_count[tax_class] = [species]

    for tgroup, counts in histo.items():

        # Read the database to get all possible species entries for this class
        known = []
        with open(database, "r") as db:
            for l_no, this_line in enumerate(db):
                if tgroup in this_line:
                    species = this_line.split("\t")[2]
                    known.append(species)

        total_species = len(list(set(known)))
        seen_species = len(list(set(tax_count[tgroup])))

        values = []
        for this_len, this_count in counts.items():

            for i in range(this_count):
                values.append(this_len)

        if len(values) > 10:
            mean = round(statistics.mean(values), 0)
            stddev = statistics.stdev(values)

            # We assume that our counts correspond to number of species 
            # so we can compare this number to the total number of species for this tax group in the database
            tax_cov = 100 * round(seen_species / total_species, 2)

            data["data"][tgroup] = {"Mean amplicon length": mean, "Stddev": stddev, "Taxonomic coverage": tax_cov}

    with open(output, "w") as json_out:
        json.dump(data, json_out, indent=4, sort_keys=True)


if __name__ == '__main__':

    main(args.input, args.database, args.output)

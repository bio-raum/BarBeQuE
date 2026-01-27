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
            "section_name": "Amplicon metrics across taxonomic groups",
            "description": "Mean length and std-dev of amplicons per taxonomic group",
            "plot_type": "table",
            "pconfig": {"id": "amp_mean_size", "col1_header": "Taxon", "col2_header": "Mean amplicon length" },
            "data": {}
            }

    histo = {}  # the amplicon size distribution by taxonomic group
    tax_count = {}  # The list of species per taxonomic group

    # Parse the consensus file and collect data
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

    # Iterate the amplicon distribution
    for tgroup, counts in histo.items():

        # Read the database to get all possible species entries for this taxonomic class
        known = []
        with open(database, "r") as db:
            for l_no, this_line in enumerate(db):
                if tgroup in this_line:
                    species = this_line.split("\t")[1]
                    known.append(species)

        total_species = len(list(set(known)))  # unique species in database
        seen_species = len(list(set(tax_count[tgroup])))  # unique species in consensus file

        amplicon_lengths = []
        # amplicon length vs observed counts
        for this_len, this_count in counts.items():

            for i in range(this_count):
                amplicon_lengths.append(this_len)

        # We ommit all taxa for which less than 10 measurements are present
        if len(amplicon_lengths) > 10:

            mean = round(statistics.mean(amplicon_lengths), 0)
            stddev = round(statistics.stdev(amplicon_lengths), 2)

            # We assume that our counts correspond to number of species 
            # so we can compare this number to the total number of species for this tax group in the database
            tax_cov = 100 * round(seen_species / total_species, 2)

            data["data"][tgroup] = {"Mean amplicon length": mean, "Stddev": stddev, "Database coverage": tax_cov}

    with open(output, "w") as json_out:
        json.dump(data, json_out, indent=4, sort_keys=True)


if __name__ == '__main__':

    main(args.input, args.database, args.output)

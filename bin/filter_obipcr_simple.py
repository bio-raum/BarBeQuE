#!/usr/bin/env python3

import argparse
import json
import sys
from pathlib import Path

iupac = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT"
}

# Hardcoded settings for testing by hand (Nextflow passes these as arguments)
# AMPLICONS_FILE = "/home/saj/barbeque-alpha/BarBeQuE/its2-test-hamming/raw/obipcr/ITS2_collapsed_euphyllophyta_raw.fasta"
# PRIMER_INPUT_FILE = "/home/saj/barbeque-alpha/samplesheet/ITS2_collapsed.fasta"
# PRIMER_NAME = "ITS2_collapsed"
# DATABASE_NAME = "euphyllophyta"
# MISMATCHES = 2
# OUTPUT_FILE = "/home/saj/barbeque-alpha/BarBeQuE/ITS2_collapsed_filtered.fasta"
# STATS_FILE = "/home/saj/barbeque-alpha/BarBeQuE/ITS2_collapsed_filter_stats.tsv"

# Settings from the command line (given by the FILTER_OBIPCR module)
parser = argparse.ArgumentParser()
parser.add_argument("--amplicons", required=True)
parser.add_argument("--primer-input", required=True)
parser.add_argument("--primer", required=True)
parser.add_argument("--db", default="")
parser.add_argument("--mismatches", required=True, type=int)
parser.add_argument("--out", required=True)
parser.add_argument("--stats", required=True)
args = parser.parse_args()

AMPLICONS_FILE = args.amplicons
PRIMER_INPUT_FILE = args.primer_input
PRIMER_NAME = args.primer
DATABASE_NAME = args.db
MISMATCHES = args.mismatches
OUTPUT_FILE = args.out
STATS_FILE = args.stats

# Primer input can be the file itself, or a folder with <primer>.fasta inside
if Path(PRIMER_INPUT_FILE).is_dir():
    PRIMER_INPUT_FILE = Path(PRIMER_INPUT_FILE) / (PRIMER_NAME + ".fasta")

# Read original primers
primer_sequence = []
with open(PRIMER_INPUT_FILE, "r") as primer_file:
    sequence = ""
    for line in primer_file:
        line = line.strip()
        if line.startswith(">"):
            if sequence != "":
                primer_sequence.append(sequence.upper())
            sequence = ""
        else:
            sequence += line
    if sequence != "":
        primer_sequence.append(sequence.upper())

total_count = 0
kept_count = 0
forward_failed_count = 0
reverse_failed_count = 0
keep_record = False

with open(AMPLICONS_FILE, "r") as obi_pcr_fasta, open(OUTPUT_FILE, "w") as output_fasta:
    for line in obi_pcr_fasta:
        line = line.strip()

        if line.startswith(">"):
            total_count += 1
            keep_record = False

            sequence_id, metadata_text = line[1:].split(" ", 1)
            header = json.loads(metadata_text)

            forward_match = header["forward_match"].upper()
            reverse_match = header["reverse_match"].upper()
            forward_primer = header["forward_primer"]
            reverse_primer = header["reverse_primer"]

            # Find # positions in forward primer
            forward_fixed_positions = []
            position = -1
            for character in forward_primer:
                if character == "#":
                    forward_fixed_positions.append(position)
                else:
                    position += 1

            # Find # positions in reverse primer
            reverse_fixed_positions = []
            position = -1
            for character in reverse_primer:
                if character == "#":
                    reverse_fixed_positions.append(position)
                else:
                    position += 1

            # Find best original forward primer
            best_forward_primer = None
            best_forward_distance = None
            forward_same_length = False

            for primer in primer_sequence:
                if len(primer) != len(forward_match):
                    continue

                forward_same_length = True
                distance = 0
                fixed_mismatch = False

                for i in range(len(primer)):
                    primer_base = primer[i]
                    binding_base = forward_match[i]

                    if binding_base not in iupac[primer_base]:
                        distance += 1
                        if i in forward_fixed_positions:
                            fixed_mismatch = True

                if fixed_mismatch:
                    continue

                if best_forward_distance is None or distance < best_forward_distance:
                    best_forward_primer = primer
                    best_forward_distance = distance

            # Find best original reverse primer
            best_reverse_primer = None
            best_reverse_distance = None
            reverse_same_length = False

            for primer in primer_sequence:
                if len(primer) != len(reverse_match):
                    continue

                reverse_same_length = True
                distance = 0
                fixed_mismatch = False

                for i in range(len(primer)):
                    primer_base = primer[i]
                    binding_base = reverse_match[i]

                    if binding_base not in iupac[primer_base]:
                        distance += 1
                        if i in reverse_fixed_positions:
                            fixed_mismatch = True

                if fixed_mismatch:
                    continue

                if best_reverse_distance is None or distance < best_reverse_distance:
                    best_reverse_primer = primer
                    best_reverse_distance = distance

            # Wrong primer file / wrong primer lengths
            if not forward_same_length:
                sys.exit(f"ERROR: No primer has the same length as forward match ({len(forward_match)} bases).")

            if not reverse_same_length:
                sys.exit(f"ERROR: No primer has the same length as reverse match ({len(reverse_match)} bases).")

            # Check forward failure
            forward_failed = best_forward_distance is None or best_forward_distance > MISMATCHES
            if forward_failed:
                forward_failed_count += 1

            # Check reverse failure
            reverse_failed = best_reverse_distance is None or best_reverse_distance > MISMATCHES
            if reverse_failed:
                reverse_failed_count += 1

            # Drop whole amplicon if either side fails
            if forward_failed or reverse_failed:
                continue

            # Put # positions back into best forward primer
            corrected_forward_primer = ""
            for i in range(len(best_forward_primer)):
                corrected_forward_primer += best_forward_primer[i]
                if i in forward_fixed_positions:
                    corrected_forward_primer += "#"

            # Put # positions back into best reverse primer
            corrected_reverse_primer = ""
            for i in range(len(best_reverse_primer)):
                corrected_reverse_primer += best_reverse_primer[i]
                if i in reverse_fixed_positions:
                    corrected_reverse_primer += "#"

            # Update metadata for PARSE_OBIPCR / TSV
            header["forward_primer"] = corrected_forward_primer
            header["forward_error"] = best_forward_distance
            header["reverse_primer"] = corrected_reverse_primer
            header["reverse_error"] = best_reverse_distance

            output_header = ">" + sequence_id + " " + json.dumps(header, separators=(",", ":"))
            output_fasta.write(output_header + "\n")

            keep_record = True
            kept_count += 1

        else:
            if keep_record:
                output_fasta.write(line + "\n")

discarded_count = total_count - kept_count

with open(STATS_FILE, "w") as stats:
    stats.write("primer\tdb\ttotal\tkept\tdiscarded\tforward_unexplained\treverse_unexplained\n")
    stats.write(
        f"{PRIMER_NAME}\t{DATABASE_NAME}\t{total_count}\t{kept_count}\t"
        f"{discarded_count}\t{forward_failed_count}\t{reverse_failed_count}\n"
    )

print("Total:", total_count)
print("Kept:", kept_count)
print("Discarded:", discarded_count)
print("Forward failed:", forward_failed_count)
print("Reverse failed:", reverse_failed_count)
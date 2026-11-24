#!/usr/bin/env python3

import argparse
import csv
import json

FIELDS = [
    "Sequence_ID",
    "Amplicon_Length",
    "Forward_Primer",
    "Forward_Match",
    "Forward_Errors",
    "Forward_Mismatch_Positions_Primer",
    "Forward_Primer_GC",
    "Forward_Match_GC",
    "Reverse_Primer",
    "Reverse_Match",
    "Reverse_Errors",
    "Reverse_Mismatch_Positions_Primer",
    "Reverse_Primer_GC",
    "Reverse_Match_GC",
    "Amplicon_GC",
    "Amplicon_Sequence"
]


IUPAC = {
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


# Calculate GC percentage
def calculate_gc(sequence):

    if sequence == "":
        raise ValueError("Sequence is empty")

    gc = 0

    for base in sequence.upper():

        if base not in IUPAC:
            raise ValueError(f"Unknown DNA character: {base}")

        bases = IUPAC[base]

        for possible_base in bases:
            if possible_base == "G" or possible_base == "C":
                gc += 1 / len(bases)

    return round(gc / len(sequence) * 100, 2)


# Find mismatch positions
def find_mismatches(primer, match):

    if len(primer) != len(match):
        raise ValueError(
            f"Primer and match have different lengths: "
            f"{len(primer)} vs {len(match)}"
        )

    positions = []
    positions_3prime = []

    for i in range(len(primer)):

        primer_base = primer[i].upper()
        match_base = match[i].upper()

        if primer_base not in IUPAC:
            raise ValueError(f"Unknown primer character: {primer_base}")

        if match_base not in IUPAC:
            raise ValueError(f"Unknown match character: {match_base}")

        if match_base not in IUPAC[primer_base]:
            positions.append(str(i + 1))
            positions_3prime.append(str(len(primer) - i))

    if len(positions) == 0:
        return "None", "None"

    return ",".join(positions), ",".join(positions_3prime)


def parse_fasta(input_fasta, output_tsv):

    with open(input_fasta, "r") as fasta, open(output_tsv, "w", newline="") as out:

        writer = csv.writer(out, delimiter="\t")
        writer.writerow(FIELDS)

        sequence = ""
        data = None

        for line in fasta:
            line = line.strip()

            # Header line
            if line.startswith(">"):

                # Write previous record
                if data is not None:
                    data.append(calculate_gc(sequence))
                    data.append(sequence)
                    writer.writerow(data)

                # First word of the header
                first_word = line.split(" ")[0]

                if "_sub[" not in first_word or ".." not in first_word:
                    raise ValueError(f"Malformed FASTA header: {line}")

                # Remove >
                first_word = first_word.replace(">", "")
                # Sequence ID
                sequence_id = first_word.split("_sub[")[0]
                # Get positions
                positions = first_word.split("_sub[")[1]
                positions = positions.replace("]", "")

                start = int(positions.split("..")[0])
                end = int(positions.split("..")[1])

                # Calculate amplicon length
                amplicon_length = end - start + 1

                # Everything after the first space is metadata
                metadata_text = line.split(" ", 1)[1]
                metadata = json.loads(metadata_text)

                #remove the # if the primers are clamped
                forward_primer = metadata["forward_primer"].replace("#", "")
                forward_match = metadata["forward_match"]
                forward_errors = metadata["forward_error"]

                reverse_primer = metadata["reverse_primer"].replace("#", "")
                reverse_match = metadata["reverse_match"]
                reverse_errors = metadata["reverse_error"]

                # Calculate mismatch positions (counted from the 5' end)
                forward_mismatch, _ = find_mismatches(
                    forward_primer,
                    forward_match
                )

                reverse_mismatch, _ = find_mismatches(
                    reverse_primer,
                    reverse_match
                )

                # Create TSV row
                data = [
                    sequence_id,
                    amplicon_length,

                    forward_primer,
                    forward_match,
                    forward_errors,
                    forward_mismatch,
                    calculate_gc(forward_primer),
                    calculate_gc(forward_match),

                    reverse_primer,
                    reverse_match,
                    reverse_errors,
                    reverse_mismatch,
                    calculate_gc(reverse_primer),
                    calculate_gc(reverse_match)
                ]

                sequence = ""

            else:
                sequence += line

        # Write final record
        if data is not None:
            data.append(calculate_gc(sequence))
            data.append(sequence)
            writer.writerow(data)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input",
        help="OBI-PCR FASTA file"
    )
    parser.add_argument(
        "output",
        help="Output TSV file"
    )
    args = parser.parse_args()

    parse_fasta(
        args.input,
        args.output
    )


if __name__ == "__main__":
    main()

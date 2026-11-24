#!/usr/bin/env python3

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from mask import read_fasta
from parse_obipcr import IUPAC
from parse_primers import find_input_files, make_name_safe, split_directions
from parse_primers import read_fasta as read_primer_fasta


# Read the obipcr JSON annotations after the first space of a header
def read_annotations(header):
    return json.loads(header.split(" ", 1)[1])


# Load real primer variants
def load_primers(folder, name):
    for path in find_input_files(folder):
        if make_name_safe(path.stem) == name:
            records = read_primer_fasta(path)
            groups, error = split_directions(records, path.name)
            if error:
                sys.exit(error)
            return [p.replace("U", "T") for p in groups["fwd"] + groups["rev"]]
    sys.exit(f"Primer set '{name}' not found")


# Get fixed 3' positions
def clamps(pattern):
    fixed, pos = set(), -1
    for base in pattern:
        if base == "#":
            fixed.add(pos)
        else:
            pos += 1
    return fixed


# Check whether any real primer matches
def passes(site, primers, max_mm, fixed):
    for primer in primers:
        if len(primer) != len(site):
            continue

        mm = {
            i for i, (p, s) in enumerate(zip(primer, site))
            if s not in IUPAC[p]
        }

        if len(mm) <= max_mm and not mm & fixed:
            return True

    return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--amplicons", required=True)
    parser.add_argument("--primer-input", required=True)
    parser.add_argument("--primer", required=True)
    parser.add_argument("--db", default="")
    parser.add_argument("--mismatches", required=True, type=int)
    parser.add_argument("--out", required=True)
    parser.add_argument("--stats", required=True)
    args = parser.parse_args()

    primers = load_primers(args.primer_input, args.primer)
    total = kept = f_fail = r_fail = 0

    # Filter OBIPCR amplicons
    with open(args.out, "w") as out:
        for header, seq in read_fasta(args.amplicons):
            total += 1
            a = read_annotations(header)

            f_ok = passes(
                a["forward_match"].upper(),
                primers,
                args.mismatches,
                clamps(a["forward_primer"]),
            )

            r_ok = passes(
                a["reverse_match"].upper(),
                primers,
                args.mismatches,
                clamps(a["reverse_primer"]),
            )

            f_fail += not f_ok
            r_fail += not r_ok

            if f_ok and r_ok:
                out.write(f">{header}\n{seq}\n")
                kept += 1

    # Save summary
    discarded = total - kept

    with open(args.stats, "w") as out:
        out.write(
            "primer\tdb\ttotal\tkept\tdiscarded\t"
            "forward_unexplained\treverse_unexplained\n"
        )
        out.write(
            f"{args.primer}\t{args.db}\t{total}\t{kept}\t{discarded}\t"
            f"{f_fail}\t{r_fail}\n"
        )


if __name__ == "__main__":
    main()

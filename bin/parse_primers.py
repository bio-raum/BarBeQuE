#!/usr/bin/env python3

"""
Convert primer FASTA files into an OBIPCR-ready TSV samplesheet.

One FASTA file is one primer set:
- Forward and reverse primers are identified from the FASTA headers.
- All forward primers in the file are combined into a single IUPAC-degenerate
  forward primer, and all reverse primers into a single reverse primer.
- The primer set is named after the file, so a directory holding three FASTA
  files produces three primer sets - one per file.
- Every forward primer must have the same length, and so must every reverse
  primer. Primers of differing lengths cannot be combined and are an error.
- Existing IUPAC bases such as Y, R, and N are supported.
- Invalid sequences, missing forward/reverse primers, and unrecognized
  directions are reported as errors.

Best practices for primer FASTA files:
- Keep one primer set per file, and name the file after the primer set.
- Include the primer direction in every header using:
    fwd, forward, rev, or reverse

  Good examples:
      >ITS2_fwd1
      >ITS2_fwd2
      >ITS2_rev1
      >ITS2_rev2

- Variants of the same primer are simply extra records with the same direction.
  Where the variant tag sits does not matter, so all of these are reverse
  primers of the same set:

      >MA_REV
      >MA_REV2
      >MA_ALT_REV

- A file holding exactly two records and no direction tokens at all is read as
  forward first, reverse second:

      >mlCOIintF
      >jgHCO2198

- Use only valid DNA/IUPAC nucleotide symbols:
    A, C, G, T, R, Y, S, W, K, M, B, D, H, V, N
- Avoid spaces, gaps, or other non-nucleotide characters inside sequences.
- Primer sequences may span multiple FASTA lines, but one line per primer is
  easier to read and maintain.
- Keep filenames descriptive, because the filename becomes the primer name:
      ITS2.fasta
      rbcL.fasta
      16S_ASU184.fasta

--input takes any number of primer FASTAs, directories of them, or a mix:

    parse_primers.py --input ITS2.fasta rbcL.fasta --min 100 --max 500 --out primers.tsv
    parse_primers.py --input primers/ --min 100 --max 500 --out primers.tsv

The output TSV contains:
    primer    fwd    rev    min    max
"""

import argparse
import re
import sys
from pathlib import Path

# 1. IUPAC CODES

# What each IUPAC symbol means.
IUPAC_TO_BASES = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "U": {"T"},

    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"C", "G"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},

    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},

    "N": {"A", "C", "G", "T"},
}


# Reverse lookup:
#
# {"C", "T"} -> Y
# {"A", "G"} -> R
# etc.
BASES_TO_IUPAC = {
    frozenset({"A"}): "A",
    frozenset({"C"}): "C",
    frozenset({"G"}): "G",
    frozenset({"T"}): "T",

    frozenset({"A", "G"}): "R",
    frozenset({"C", "T"}): "Y",
    frozenset({"C", "G"}): "S",
    frozenset({"A", "T"}): "W",
    frozenset({"G", "T"}): "K",
    frozenset({"A", "C"}): "M",

    frozenset({"C", "G", "T"}): "B",
    frozenset({"A", "G", "T"}): "D",
    frozenset({"A", "C", "T"}): "H",
    frozenset({"A", "C", "G"}): "V",

    frozenset({"A", "C", "G", "T"}): "N",
}


VALID_BASES = set(IUPAC_TO_BASES.keys())


# 2. SETTINGS

FASTA_SUFFIXES = (".fa", ".fasta", ".fna")

TSV_HEADER = [
    "primer",
    "fwd",
    "rev",
    "min",
    "max"
]


# Recognizes:
#
# MA_FWD
# MA_FORWARD
# MA_ALT_REV
# ITS2-short_fwd1.1
# ITS2-short_rev3.2
#
# but does not simply match "fwd" or "rev" anywhere inside a word.
DIRECTION_RE = re.compile(
    r"(?i)"
    r"(?:^|[\s_.-])"
    r"(fwd|forward|rev|reverse)"
    r"(?=\d|[\s_.-]|$)"
)


# 3. READ FASTA

def read_fasta(path):
    """
    Read a FASTA file.

    Returns:
        [
            (header, sequence),
            (header, sequence),
            ...
        ]
    """

    records = []

    current_header = None
    current_sequence = []

    with open(path) as file:
        for line in file:
            line = line.strip()

            # Ignore empty lines
            if not line:
                continue

            # New FASTA record
            if line.startswith(">"):
                # Save previous record
                if current_header is not None:
                    sequence = "".join(current_sequence)

                    records.append((current_header, sequence))

                current_header = line[1:].strip()
                current_sequence = []

            else:
                current_sequence.append(line)

    # Save final record
    if current_header is not None:
        sequence = "".join(current_sequence)

        records.append((current_header, sequence))

    return records


# 4. FIND THE PRIMER DIRECTION

def find_direction(header):
    """
    Return "fwd", "rev", or None when the header carries no direction.

    Only the direction matters - everything else in the header is a name for
    humans, so MA_REV, MA_REV2 and MA_ALT_REV are all reverse primers.
    """

    match = DIRECTION_RE.search(header)

    if match is None:
        return None

    if match.group(1).lower() in ("fwd", "forward"):
        return "fwd"

    return "rev"


# 5. VALIDATE A PRIMER SEQUENCE

def validate_sequence(header, sequence):
    """
    Return an error message if the primer sequence is invalid.
    Otherwise return None.
    """

    if sequence == "":
        return f"record '{header}' has an empty sequence"

    sequence = sequence.upper()

    invalid_characters = set()

    for character in sequence:
        if character not in VALID_BASES:
            invalid_characters.add(character)

    if invalid_characters:
        invalid_text = ", ".join(sorted(invalid_characters))

        return (
            f"record '{header}' contains invalid character(s): "
            f"{invalid_text}"
        )

    return None


# 6. SPLIT ONE FILE INTO FORWARD AND REVERSE PRIMERS

def split_directions(records, filename):
    """
    Sort every record of one FASTA into forward and reverse primers.

    Example:

        MA_FWD
        MA_REV
        MA_ALT_REV

    becomes:

        fwd = [MA_FWD]
        rev = [MA_REV, MA_ALT_REV]

    Returns (primers, error).
    """

    primers = {
        "fwd": [],
        "rev": []
    }

    unlabelled_headers = []

    for header, sequence in records:
        direction = find_direction(header)

        if direction is None:
            unlabelled_headers.append(header)

            continue

        primers[direction].append(sequence.upper())

    # --------------------------------------------------------
    # Special case:
    #
    # Exactly two primers and neither has FWD/REV labels.
    #
    # Assume:
    # first = FWD
    # second = REV
    # --------------------------------------------------------

    if len(unlabelled_headers) == 2 and len(records) == 2:
        return {
            "fwd": [records[0][1].upper()],
            "rev": [records[1][1].upper()]
        }, None

    # --------------------------------------------------------
    # If some primer directions cannot be identified, fail.
    # --------------------------------------------------------

    if unlabelled_headers:
        return None, (
            f"{filename}: could not determine FWD/REV "
            f"for these primer(s): {unlabelled_headers}"
        )

    # --------------------------------------------------------
    # A primer set needs both directions.
    # --------------------------------------------------------

    if not primers["fwd"] or not primers["rev"]:
        return None, (
            f"{filename}: found "
            f"{len(primers['fwd'])} fwd and "
            f"{len(primers['rev'])} rev - "
            f"a primer set needs both."
        )

    return primers, None


# 7. COMBINE THE PRIMERS OF ONE DIRECTION

def collapse_iupac(sequences):
    """
    Collapse primer variants into one IUPAC primer.

    Example:

        ACGT
        ACCT

    position 3 contains:
        G + C

    G/C = S

    result:
        ACST
    """

    sequence_length = len(sequences[0])

    collapsed_sequence = ""

    # Go through each nucleotide position
    for position in range(sequence_length):
        bases_at_position = set()

        # Look at this position in every primer
        for sequence in sequences:
            symbol = sequence[position].upper()

            # Expand IUPAC symbols.
            #
            # Example:
            # Y -> {C, T}
            # R -> {A, G}
            # A -> {A}
            possible_bases = IUPAC_TO_BASES[symbol]

            bases_at_position.update(possible_bases)

        # Convert the combined bases back to IUPAC.
        iupac_code = BASES_TO_IUPAC[
            frozenset(bases_at_position)
        ]

        collapsed_sequence += iupac_code

    return collapsed_sequence


def collapse_direction(sequences, direction, filename):
    """
    Combine every primer of one direction into a single sequence.

    Primers can only be combined position by position, so they all have to be
    the same length. Differing lengths mean the file holds more than one primer
    set, which this script does not support.

    Returns (sequence, error).
    """

    lengths = sorted({len(sequence) for sequence in sequences})

    if len(lengths) > 1:
        length_text = ", ".join(str(length) for length in lengths)

        return None, (
            f"{filename}: {direction} primers have different lengths "
            f"({length_text}) - a primer FASTA must contain one primer set; "
            f"split them into separate files."
        )

    return collapse_iupac(sequences), None


# 8. CONVERT ONE FASTA INTO A TSV ROW

def process_fasta(path, min_length, max_length):
    """
    Process one primer FASTA into the single primer set it describes.

    Returns:

        rows
        warnings
        errors
    """

    rows = []
    warnings = []
    errors = []

    # --------------------------------------------------------
    # Read FASTA
    # --------------------------------------------------------

    records = read_fasta(path)

    if len(records) == 0:
        errors.append(f"{path.name}: no FASTA records found")

        return rows, warnings, errors

    # --------------------------------------------------------
    # Validate every primer sequence
    # --------------------------------------------------------

    for header, sequence in records:
        error = validate_sequence(header, sequence)

        if error is not None:
            errors.append(f"{path.name}: {error}")

    # Stop processing this file if sequences are invalid.
    if errors:
        return rows, warnings, errors

    # --------------------------------------------------------
    # Separate FWD and REV primers
    # --------------------------------------------------------

    primers, error = split_directions(records, path.name)

    if error is not None:
        errors.append(error)

        return rows, warnings, errors

    # --------------------------------------------------------
    # Combine each direction into one primer
    # --------------------------------------------------------

    forward, forward_error = collapse_direction(primers["fwd"], "fwd", path.name)

    reverse, reverse_error = collapse_direction(primers["rev"], "rev", path.name)

    for error in (forward_error, reverse_error):
        if error is not None:
            errors.append(error)

    if errors:
        return rows, warnings, errors

    # --------------------------------------------------------
    # Combining variants puts a sequence in the samplesheet that appears in no
    # input file, so record that it happened.
    # --------------------------------------------------------

    if len(primers["fwd"]) > 1 or len(primers["rev"]) > 1:
        warnings.append(
            f"{path.name}: combined "
            f"{len(primers['fwd'])} fwd and "
            f"{len(primers['rev'])} rev primer(s) "
            f"into one degenerate primer pair."
        )

    rows.append({
        "primer": make_name_safe(path.stem),
        "fwd": forward,
        "rev": reverse,
        "min": min_length,
        "max": max_length,
    })

    return rows, warnings, errors


# 9. FIND INPUT FASTA FILES

# Characters that are unsafe in a file or directory name. Primer names
# become output paths later, so they are replaced with "_".
UNSAFE_NAME_RE = re.compile(r'[\s/\\:*?"<>|]')


def make_name_safe(name):
    """Replace path-unsafe characters in a primer or file name with '_'."""

    return UNSAFE_NAME_RE.sub("_", name)


def is_fasta(path):
    """
    True when the first non-blank line starts with ">".

    Used to tell a primer FASTA apart from a TSV samplesheet.
    """

    with open(path) as file:
        for line in file:
            if line.strip():
                return line.startswith(">")

    return False


def find_input_files(inputs):
    """
    --input may be any number of:

        primer FASTA files

    and:

        directories containing primer FASTA files
    """

    # A single path is also accepted, so callers do not have to wrap it.
    if isinstance(inputs, (str, Path)):
        inputs = [inputs]

    fasta_files = []

    for raw_input in inputs:
        input_path = Path(raw_input)

        # Single file
        if input_path.is_file():
            fasta_files.append(input_path)

            continue

        # Directory
        if input_path.is_dir():
            found_files = []

            for path in sorted(input_path.iterdir()):
                if path.suffix.lower() in FASTA_SUFFIXES:
                    found_files.append(path)

            if len(found_files) == 0:
                sys.exit(
                    f"ERROR: No FASTA files found in "
                    f"{input_path}"
                )

            fasta_files.extend(found_files)

            continue

        sys.exit(
            f"ERROR: Input does not exist: "
            f"{input_path}"
        )

    return fasta_files


# 10. PROCESS EVERY FASTA

def collect_rows(paths, min_length, max_length):
    """
    Process every primer FASTA and gather the results.

    Every file is checked before anything is written, so all problems are
    reported in one go. If any file failed, no rows are returned at all -
    a typo can never silently drop a primer set from the comparison.
    """

    all_rows = []
    all_warnings = []
    all_errors = []

    # Primer names become output paths, so two files cannot share one name.
    files_by_name = {}

    for path in paths:
        path = Path(path)

        rows, warnings, errors = process_fasta(
            path,
            min_length,
            max_length
        )

        for row in rows:
            previous_file = files_by_name.get(row["primer"])

            if previous_file is None:
                files_by_name[row["primer"]] = path

            else:
                all_errors.append(
                    f"{path}: primer name '{row['primer']}' is already used "
                    f"by {previous_file} - every primer FASTA needs its own name."
                )

        all_rows.extend(rows)
        all_warnings.extend(warnings)
        all_errors.extend(errors)

    if all_errors:
        return [], all_warnings, all_errors

    return all_rows, all_warnings, all_errors


# 11. WRITE TSV

def write_tsv(rows, output_file):
    with open(output_file, "w") as file:
        file.write("\t".join(TSV_HEADER) + "\n")

        for row in rows:
            values = []

            for column in TSV_HEADER:
                values.append(str(row[column]))

            file.write("\t".join(values) + "\n")


# 12. MAIN

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "--input",
        required=True,
        nargs="+",
        help="Primer FASTA files and/or directories of primer FASTAs"
    )

    parser.add_argument(
        "--min",
        required=True,
        type=int,
        help="Minimum amplicon length"
    )

    parser.add_argument(
        "--max",
        required=True,
        type=int,
        help="Maximum amplicon length"
    )

    parser.add_argument("--out", required=True, help="Output TSV samplesheet")

    parser.add_argument("--warnings", help="Optional warnings file")

    args = parser.parse_args()

    # --------------------------------------------------------
    # Find FASTA files
    # --------------------------------------------------------

    fasta_files = find_input_files(args.input)

    # --------------------------------------------------------
    # Process every FASTA
    # --------------------------------------------------------

    all_rows, all_warnings, all_errors = collect_rows(fasta_files, args.min, args.max)

    # --------------------------------------------------------
    # Print warnings
    # --------------------------------------------------------

    for warning in all_warnings:
        print("WARNING:", warning, file=sys.stderr)

    # Optional warning file
    if args.warnings:
        with open(args.warnings, "w") as file:
            for warning in all_warnings:
                file.write(warning + "\n")

    # --------------------------------------------------------
    # If any FASTA failed, write NO samplesheet.
    # --------------------------------------------------------

    if all_errors:
        print("\nERROR: Primer parsing failed.", file=sys.stderr)

        for error in all_errors:
            print(
                f"  - {error}",
                file=sys.stderr
            )

        print("\nNo primer samplesheet was created.", file=sys.stderr)

        sys.exit(1)

    # --------------------------------------------------------
    # Everything passed -> write output
    # --------------------------------------------------------

    write_tsv(all_rows, args.out)

    print(
        f"Created {args.out} "
        f"with {len(all_rows)} primer set(s)."
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse


def read_fasta(path):
    """Yield (header, sequence) pairs, joining multi-line sequences."""
    header = None
    seq_parts = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            yield header, "".join(seq_parts)

#Replace the middle of one sequence with 20 Ns if it exceeds 2*read_length, keeping both ends.
def mask_paired_end(sequence, read_length):
    length = len(sequence)
    if length <= 2 * read_length:
        return sequence
    return sequence[:read_length] + "N" * 20 + sequence[-read_length:]

#No N's added for single end
def mask_single_end(sequence, read_length):
    length = len(sequence)
    if length <= read_length:
        return sequence
    return sequence[:read_length]


def main():
    parser = argparse.ArgumentParser(
        description="Mask the center of amplicons to simulate a paired-end gap."
    )
    parser.add_argument("--input", required=True, help="input amplicon FASTA")
    parser.add_argument("--output", required=True, help="masked FASTA to write")
    parser.add_argument(
        "--read-length",
        type=int,
        default=None,
        help="bases kept at each end (the R1/R2 read length) or max sequence length for single-end. Defaults to 400 for single-end and 150 for paired-end.",
    )
    parser.add_argument(
        "--single-end",
        action="store_true",
        help="Use single-end masking logic (truncate to read-length, no Ns)",
    )
    args = parser.parse_args()

    if args.read_length is None:
        args.read_length = 400 if args.single_end else 150

    masked = 0
    untouched = 0
    with open(args.output, "w") as out:
        for header, sequence in read_fasta(args.input):
            if args.single_end:
                new_sequence = mask_single_end(sequence, args.read_length)
            else:
                new_sequence = mask_paired_end(sequence, args.read_length)
            if new_sequence != sequence:
                masked += 1
            else:
                untouched += 1
            out.write(">" + header + "\n" + new_sequence + "\n")

    print("masked   :", masked)
    print("untouched:", untouched)


if __name__ == "__main__":
    main()
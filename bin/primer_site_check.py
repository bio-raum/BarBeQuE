#!/usr/bin/env python3
"""
primer_site_check.py - why did a primer not bind in a record that has the barcode but gave no amplicon?

The best BLAST hit shows where the barcode is in the record. The centroids have no primers,
so each primer site lies right outside the barcode. For each record and primer the result is:
    PRIMER_CUTOFF  the site lies (partly) outside the record: never sequenced, a database problem
    UNSEQUENCED    the site is there, but has N bases: not sequenced, a database problem
    MISMATCH       the site is recognizable, but has more mismatches than obipcr allows
    PRIMER_SITE_NOT_FOUND
                   the best candidate has more mismatches than the configured site-recognition limit
    BOUND          the site is there and the primer binds: the cause is something else
"""

import argparse
import csv

import pandas as pd
from Bio import SeqIO
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Seq import reverse_complement


def accession(name):
    # BLAST drops a trailing ';' from the record name, so compare only the part before the first ';'
    # Example: AF218533;tax=k:Viridiplantae,...; -> AF218533
    return name.split(";")[0]


def check_site(seq, primer, start, turn, allowed, site_not_found_mismatches=4, window=100):
    """start: predicted first base of the site (1 = first base of the record).
    turn: on the record the site reads as the reverse complement of the primer.
    window: how many bases on each side of the predicted place to search. The predicted place is a guess
    when BLAST did not align the end of the barcode: the record's unaligned part can be longer or shorter
    than the centroid's (seen up to 25 bp in ITS2)."""
    end = start + len(primer) - 1

    # Sliding window: move the primer one base at a time and count the mismatches at each start
    # No gaps: an extra or missing base in the site shows up as several mismatches
    sites = []
    for s in range(max(1, start - window), min(len(seq) - len(primer) + 1, start + window) + 1):
        # Step 1: cut out a piece of the record as long as the primer
        site = seq[s - 1:s - 1 + len(primer)]
        # Reverse primer: turn the piece around so it reads like the primer
        if turn:
            site = reverse_complement(site)
        # Step 2: count the mismatches
        # Mismatch: the record base does not fit the primer base (the primer may have IUPAC codes, R = A or G)
        # Positions from the 5' end of the primer; N counts as a mismatch, as in obipcr
        mismatches = [i for i, (code, base) in enumerate(zip(primer, site), start=1)
                      if base not in ambiguous_dna_values[code]]
        # Step 3: does the primer bind here?
        binds = len(mismatches) <= allowed and "N" not in site
        sites.append({"start": s, "site": site, "mismatches": mismatches, "binds": binds})

    # Site (partly) outside the record and no primer-like site inside it: never sequenced
    if (start < 1 or end > len(seq)) and not any(len(s["mismatches"]) <= site_not_found_mismatches for s in sites):
        return ["PRIMER_CUTOFF", "NA", "NA", start, end, "NA"]

    # Step 4: best start: binds first, then the fewest mismatches, then the fewest N (ties: the leftmost start)
    best = min(sites, key=lambda site: (not site["binds"], len(site["mismatches"]), site["site"].count("N")))

    # Step 5: give the best site a status
    if len(best["mismatches"]) > site_not_found_mismatches:
        status = "PRIMER_SITE_NOT_FOUND"
    elif len(best["mismatches"]) > allowed:
        status = "MISMATCH"
    elif "N" in best["site"]:
        status = "UNSEQUENCED"
    else:
        status = "BOUND"

    positions = ",".join(map(str, best["mismatches"])) or "None"
    return [status, len(best["mismatches"]), positions, best["start"], best["start"] + len(primer) - 1, best["site"]]


def read_parsed_primers(path):
    primers = {"forward": set(), "reverse": set()}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            primers["forward"].add(row["Forward_Primer"].upper())
            primers["reverse"].add(row["Reverse_Primer"].upper())

    for direction, sequences in primers.items():
        if not sequences:
            raise ValueError(
                f"parsed OBIPCR table contains no {direction} primers: {path}"
            )
        lengths = {len(sequence) for sequence in sequences}
        if len(lengths) != 1:
            raise ValueError(
                f"parsed OBIPCR table has different-length {direction} primers: "
                f"{sorted(lengths)}"
            )

    return {
        direction: sorted(sequences)
        for direction, sequences in primers.items()
    }


def check_primer_variants(
    seq,
    primers,
    start,
    turn,
    allowed,
    site_not_found_mismatches,
    window,
):
    candidates = []
    for primer in primers:
        result = check_site(
            seq,
            primer,
            start,
            turn,
            allowed,
            site_not_found_mismatches,
            window,
        )
        if result[0] == "PRIMER_CUTOFF":
            key = (True, float("inf"), float("inf"))
        else:
            key = (
                result[0] != "BOUND",
                int(result[1]),
                result[5].count("N"),
            )
        candidates.append((key, primer, result))

    _, best_primer, best_result = min(candidates)
    return best_primer, best_result


def write_mismatch_fasta(table, primers, prefix):
    """One FASTA per primer with the MISMATCH sites, to see why the primer did not bind.
    The original primers come first. Each site reads in primer orientation, so all sites line up under
    the primer without aligning; bases that do not fit the primer are lowercase.
    Sites with more mismatches than --site-not-found-mismatches are PRIMER_SITE_NOT_FOUND and left out.
    Header: >accession|species|mismatches|positions=...|primer=best original primer"""
    for primer_name in ("forward", "reverse"):
        rows = table[(table.primer == primer_name) & (table.status == "MISMATCH")]
        rows = rows.sort_values(["mismatches", "site"])
        with open(f"{prefix}.{primer_name}_mismatch_sites.fasta", "w") as out:
            out.writelines(f">PRIMER_{primer_name}_{i}\n{primer}\n" for i, primer in enumerate(primers[primer_name], start=1))
            for row in rows.itertuples():
                # "Curcuma longa internal transcribed spacer 2, ..." -> Curcuma_longa
                # A hybrid keeps its third word: "Mentha x piperita ..." -> Mentha_x_piperita
                words = row.description.split()
                species = "_".join(words[:3] if len(words) > 2 and words[1] == "x" else words[:2])
                site = "".join(base if base in ambiguous_dna_values[code] else base.lower()
                               for code, base in zip(row.best_original_primer, row.site))
                out.write(f">{row.record}|{species}|{row.mismatches}|positions={row.mismatch_positions}"
                          f"|primer={row.best_original_primer}\n{site}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--blast", required=True)
    parser.add_argument("--obipcr-tsv", required=True)
    parser.add_argument("--mismatches", type=int, required=True)
    parser.add_argument(
        "--site-not-found-mismatches",
        type=int,
        default=4,
        help=(
            "maximum mismatches for a candidate to remain classified as "
            "MISMATCH; candidates above this are PRIMER_SITE_NOT_FOUND "
            "(default: 4)"
        ),
    )
    parser.add_argument(
        "--window",
        type=int,
        default=100,
        help=(
            "bases searched on each side of the primer site predicted from "
            "the BLAST hit (default: 100)"
        ),
    )
    parser.add_argument("--name", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--fasta-prefix",
        help=(
            "also write <prefix>.forward_mismatch_sites.fasta and "
            "<prefix>.reverse_mismatch_sites.fasta with the MISMATCH sites"
        ),
    )
    args = parser.parse_args()

    if args.site_not_found_mismatches < args.mismatches:
        parser.error(
            "--site-not-found-mismatches must be greater than or equal to "
            "--mismatches"
        )

    primers = read_parsed_primers(args.obipcr_tsv)

    # Best hit per record: the highest bitscore (the most similar and longest alignment)
    # dtype=str: an accession like 1234 must stay text, not become a number
    blast = pd.read_csv(args.blast, sep="\t", dtype={"record": str})
    blast["record"] = blast["record"].map(accession)
    best_hits = blast.sort_values("bitscore", ascending=False).drop_duplicates("record").set_index("record")

    rows = []
    for record in SeqIO.parse(args.fasta, "fasta"):
        hit = best_hits.loc[accession(record.id)]
        seq = str(record.seq).upper()

        # Centroid bases that did not align, before and after the alignment
        before = hit.centroid_from - 1
        after = hit.centroid_len - hit.centroid_to

        if hit.record_from < hit.record_to:
            # + strand: forward site before the barcode, reverse site after it (as reverse complement)
            fwd_len = len(primers["forward"][0])
            rev_len = len(primers["reverse"][0])
            sites = {
                "forward": (hit.record_from - before - fwd_len, False),
                "reverse": (hit.record_to + after + 1, True),
            }
        else:
            # - strand: the record is turned around, so the forward site is after the barcode
            # (as reverse complement) and the reverse site before it
            fwd_len = len(primers["forward"][0])
            rev_len = len(primers["reverse"][0])
            sites = {
                "forward": (hit.record_from + before + 1, True),
                "reverse": (hit.record_to - after - rev_len, False),
            }

        # The record title without the accession, e.g. "Ilex paraguariensis psbA-trnH ..., partial sequence"
        description = record.description[len(record.id):].strip()

        for primer_name, (start, turn) in sites.items():
            best_primer, result = check_primer_variants(
                seq,
                primers[primer_name],
                start,
                turn,
                args.mismatches,
                args.site_not_found_mismatches,
                args.window,
            )
            rows.append(
                [args.name, record.id, description, len(seq), primer_name]
                + result
                + [best_primer]
            )

    columns = ["sample", "record", "description", "record_len", "primer",
               "status", "mismatches", "mismatch_positions", "site_start", "site_end", "site",
               "best_original_primer"]
    table = pd.DataFrame(rows, columns=columns)
    table.to_csv(args.output, sep="\t", index=False)

    if args.fasta_prefix:
        write_mismatch_fasta(table, primers, args.fasta_prefix)


if __name__ == "__main__":
    main()

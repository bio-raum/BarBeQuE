#!/usr/bin/env python3
"""
taxid_db_filter.py — keep only sequences belonging to a taxid (or its descendants) in a reference db.

Usage:
    taxid_db_filter.py --fasta db.fasta --accession-taxid accession_taxid.tsv \\
        --taxids taxids.txt --prefix <db_id>

--accession-taxid is a plain accession<TAB>taxid mapping (e.g. NCBI accession2taxid/genbank2taxid).
--taxids is one taxid per line (the requested taxid plus every descendant, from `taxonkit list`).
Matching is version-insensitive (MK123456 matches MK123456.1).
"""

import argparse
import re

from Bio import SeqIO


def strip_version(accession):
    """'FM163243;tax=k:...;' -> 'FM163243', 'MK123456.1' -> 'MK123456'.

    Drops the SINTAX annotation first so sintax-formatted reference dbs normalise to the same
    bare accession as plain GenBank headers, then the version/range suffix.
    """
    return re.sub(r"(\.[0-9]+)+$", "", accession.split(";")[0])


def load_taxids(path):
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def load_matching_accessions(path, keep_taxids):
    accessions = set()
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            # BarBeQuE's generated lookup is accession<TAB>taxid, whereas
            # NCBI nucl_gb.accession2taxid is accession<TAB>accession.version
            # <TAB>taxid<TAB>gi. Support both forms.
            taxid_index = 2 if len(fields) >= 3 else 1
            if len(fields) > taxid_index and fields[taxid_index] in keep_taxids:
                accessions.add(strip_version(fields[0]))
    return accessions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--accession-taxid", required=True)
    parser.add_argument("--taxids", required=True)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    keep_taxids = load_taxids(args.taxids)
    keep_accessions = load_matching_accessions(args.accession_taxid, keep_taxids)

    kept = [
        record
        for record in SeqIO.parse(args.fasta, "fasta")
        if strip_version(record.id) in keep_accessions
    ]
    SeqIO.write(kept, f"{args.prefix}.taxid_filtered.fasta", "fasta")


if __name__ == "__main__":
    main()

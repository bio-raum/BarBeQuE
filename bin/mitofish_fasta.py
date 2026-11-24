#!/usr/bin/env python3

"""Decompress the downloaded MitoFish FASTA file."""

import gzip
import shutil
import sys
from pathlib import Path

input_gz, output_fasta = sys.argv[1:3]

# Decompress the FASTA and check its first header.
with gzip.open(input_gz, "rb") as source:
    with open(output_fasta + ".partial", "wb") as output:
        first_line = source.readline()
        if not first_line.startswith(b">"):
            raise ValueError("MitoFish file does not start with a FASTA header")
        output.write(first_line)
        shutil.copyfileobj(source, output)

# Rename only after successful decompression.
Path(output_fasta + ".partial").replace(output_fasta)

print(f"MitoFish FASTA created: {output_fasta}")

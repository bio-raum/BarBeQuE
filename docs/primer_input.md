# Primer Input

Choose exactly one primer input:

| Input | Parameter | Amplicon bounds |
| --- | --- | --- |
| Samplesheet | `--input primers.tsv` | `min`/`max` columns, per row |
| Primer FASTA | `--input primers.fasta` or `--input primer_dir/` | `--primer_min`/`--primer_max`, for all sets |
| Named sets | `--primer_set COI,Fish16S` | FooDMe2 catalog, per set |

`--input` is read as FASTA when it is a directory, or when the file's first non-blank line starts
with `>`. Anything else is treated as a samplesheet. The file extension doesn't matter.

## Primer FASTA rules

- **One file is one primer set.** The set is named after the file (`ITS2.fasta` becomes `ITS2`).
  A directory of three FASTAs gives three sets.
- **Directories:** only `.fa`, `.fasta` and `.fna` files directly inside are read. Subdirectories are
  ignored.
- **Direction comes from the header.** Each header needs a `fwd`, `forward`, `rev` or `reverse` token
  (case-insensitive), separated from the rest of the name by `_`, `-`, `.` or a space. The rest of
  the header is ignored:

  ```text
  >ITS2_fwd1
  >ITS2_rev1
  >ITS2_ALT_REV
  ```

  A file with exactly two records and no direction tokens is read as forward, then reverse.

## Combining variants

obipcr takes one forward and one reverse primer, so all forward primers in a file are merged into
one IUPAC-degenerate primer, and likewise for reverse primers:

```text
>MA_fwd_1  ACGT
>MA_fwd_2  ATGT   ->  fwd AYGT
>MA_rev    TGCA       rev TGCA
```

- Primers of the same direction must all be the same length. Differing lengths are an error.
- Two different markers in one file are merged into one set. To keep them separate, use separate
  files.
- A file with one forward and one reverse primer passes through unchanged.
- Every file whose variants were merged is listed in `primers/*.primer_warnings.txt`.

## Filtering merged-primer artefacts

A merged primer also matches base combinations that none of the original primers have
(`ACGTA` + `ACCTT` → `ACSTW`, which also matches `ACGTT`). For FASTA input, `FILTER_OBIPCR`
re-checks each amplicon's binding sites against the original primers, using the same mismatch limit
and 3' clamp. It keeps an amplicon only if both sites match an original primer.

- Counts are written to `filtered_obipcr/<primer>_<db>_filter_stats.tsv`.
- Disable the filter with `--filter_collapsed_primers false`.
- The filter is not applied to samplesheet or `--primer_set` input.

## A file is rejected when

- it has no FASTA records, or a record has an empty sequence
- a sequence contains anything other than `ACGTU` and IUPAC codes
- it lacks either forward or reverse primers
- primers of one direction differ in length
- a header has no direction token (unless the file is a plain two-record pair)
- two files resolve to the same primer name

All files are checked before anything runs, and every error is reported together. If any file fails,
nothing is benchmarked.

## FASTA folder or samplesheet?

A FASTA folder uses one `--primer_min`/`--primer_max` for every set. If your markers have different
amplicon lengths, use a samplesheet with per-row bounds instead:

```text
primer  fwd                         rev                   min  max
ITS2    CGAGTYTTTGAAYGCAAGTTG       YCCCGYCTGAYCTGRGGT    200  500
rbcL    ATGTCACCACAAACAGAGACTAAAGC  GTAAAATCAAGTCCACCRCG  500  800
```

The samplesheet uses the primer strings as written and doesn't merge variants. To use merged primers
in a samplesheet, run the FASTA through the pipeline once and copy the sequences from
`primers/primers.tsv`.

The primers that were actually benchmarked are always listed in the TSV under `primers/` in your
results.

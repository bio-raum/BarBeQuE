#!/usr/bin/env python3
"""Pairwise sequence-distance analysis for one consensus cluster.

The module joins a cluster's accessions to parsed OBI-PCR amplicon sequences,
collapses exact duplicates, calculates normalized global edit distances, and
projects them with classical PCoA. It is importable by Streamlit and also
provides a small command-line interface that writes reviewable TSV files.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

CONSENSUS_COLUMNS = [
    "cluster_id",
    "accession",
    "accession_taxid",
    "accession_name",
    "assigned_name",
    "assigned_taxid",
    "assigned_rank",
    "disambiguation",
]


def accession_key(value: object) -> str:
    return pd.Series([str(value)]).str.replace(r"\.\d+$", "", regex=True).iloc[0]


def normalized_edit_distance(left: str, right: str) -> float:
    """Global Levenshtein distance divided by the longer sequence length."""
    left = str(left).upper()
    right = str(right).upper()
    if left == right:
        return 0.0
    if not left or not right:
        return 1.0
    if len(left) < len(right):
        left, right = right, left

    previous = list(range(len(right) + 1))
    for row_index, left_base in enumerate(left, start=1):
        current = [row_index]
        for column_index, right_base in enumerate(right, start=1):
            current.append(
                min(
                    current[-1] + 1,
                    previous[column_index] + 1,
                    previous[column_index - 1] + (left_base != right_base),
                )
            )
        previous = current
    return previous[-1] / max(len(left), len(right))


def load_cluster_sequences(
    consensus_path: Path,
    parsed_path: Path,
    cluster_id: int,
    focal_taxid: str,
    max_unique_sequences: int = 300,
) -> tuple[pd.DataFrame, np.ndarray]:
    consensus = pd.read_csv(
        consensus_path,
        sep="\t",
        header=None,
        names=CONSENSUS_COLUMNS,
        dtype={"accession_taxid": "string", "assigned_taxid": "string"},
    )
    cluster = consensus[consensus["cluster_id"].eq(int(cluster_id))].copy()
    if cluster.empty:
        raise ValueError(f"Cluster {cluster_id} was not found in {consensus_path.name}")

    parsed = pd.read_csv(
        parsed_path,
        sep="\t",
        usecols=["Sequence_ID", "Amplicon_Length", "Amplicon_Sequence"],
    ).rename(
        columns={
            "Sequence_ID": "sequence_id",
            "Amplicon_Length": "amplicon_length",
            "Amplicon_Sequence": "sequence",
        }
    )
    cluster["accession_key"] = cluster["accession"].map(accession_key)
    parsed["accession_key"] = parsed["sequence_id"].map(accession_key)
    joined = cluster.merge(parsed, on="accession_key", how="left")
    joined = joined.dropna(subset=["sequence"]).copy()
    if joined.empty:
        raise ValueError("No parsed amplicon sequences matched the cluster accessions")

    joined["accession_taxid"] = joined["accession_taxid"].astype(str).str.strip()
    rows = []
    for sequence_number, (sequence, members) in enumerate(
        joined.groupby("sequence", sort=False), start=1
    ):
        accessions = sorted(members["accession"].astype(str).unique())
        taxids = sorted(members["accession_taxid"].unique())
        taxon_names = sorted(
            name
            for name in members["accession_name"].dropna().astype(str).unique()
            if name and name.lower() != "unknown"
        )
        contains_focal = str(focal_taxid).strip() in taxids
        label = (
            f"S{sequence_number} | n={len(accessions)} | "
            f"taxa={';'.join(taxon_names) if taxon_names else ','.join(taxids)}"
        )
        rows.append(
            {
                "sequence_id": f"S{sequence_number}",
                "label": label,
                "sequence": sequence,
                "length": len(sequence),
                "record_count": len(accessions),
                "accessions": ";".join(accessions),
                "taxids": ";".join(taxids),
                "taxon_names": ";".join(taxon_names),
                "contains_focal_taxid": contains_focal,
            }
        )
    metadata = pd.DataFrame(rows)
    if len(metadata) > max_unique_sequences:
        raise ValueError(
            f"Cluster {cluster_id} has {len(metadata):,} unique amplicons; "
            f"the interactive safety limit is {max_unique_sequences:,}. "
            "Dereplicate more strictly or raise --max-unique for an offline run."
        )

    count = len(metadata)
    distances = np.zeros((count, count), dtype=float)
    for left_index in range(count):
        for right_index in range(left_index + 1, count):
            distance = normalized_edit_distance(
                metadata.iloc[left_index]["sequence"],
                metadata.iloc[right_index]["sequence"],
            )
            distances[left_index, right_index] = distances[right_index, left_index] = distance
    return metadata, distances


def order_by_similarity(coordinates: np.ndarray) -> list[int]:
    """Return a display order that keeps similar sequences adjacent.

    Sorting by the first PCoA axis is a cheap seriation: the axis captures the
    dominant gradient in the distance matrix, so neighbouring sequences end up
    next to each other and a heatmap ordered this way still shows its blocks.
    """
    if len(coordinates) < 2:
        return list(range(len(coordinates)))
    return list(np.argsort(coordinates[:, 0]))


def classical_pcoa(distances: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return 2-D classical PCoA coordinates and explained positive eigenvalues."""
    count = len(distances)
    if count < 2:
        return np.zeros((count, 2)), np.array([1.0, 0.0])
    centering = np.eye(count) - np.ones((count, count)) / count
    gram = -0.5 * centering @ (distances**2) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    positive = np.clip(eigenvalues[:2], 0, None)
    coordinates = eigenvectors[:, :2] * np.sqrt(positive)
    positive_total = np.clip(eigenvalues, 0, None).sum()
    explained = positive / positive_total if positive_total else np.array([0.0, 0.0])
    return coordinates, explained


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--consensus", required=True, type=Path)
    parser.add_argument("--parsed", required=True, type=Path)
    parser.add_argument("--cluster-id", required=True, type=int)
    parser.add_argument("--focal-taxid", required=True)
    parser.add_argument("--max-unique", type=int, default=300)
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args()

    metadata, distances = load_cluster_sequences(
        args.consensus,
        args.parsed,
        args.cluster_id,
        args.focal_taxid,
        args.max_unique,
    )
    coordinates, _ = classical_pcoa(distances)
    metadata.drop(columns=["sequence"]).to_csv(
        f"{args.output_prefix}.metadata.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        distances, index=metadata["sequence_id"], columns=metadata["sequence_id"]
    ).to_csv(f"{args.output_prefix}.distances.tsv", sep="\t")
    pd.DataFrame(
        coordinates, index=metadata["sequence_id"], columns=["pcoa1", "pcoa2"]
    ).to_csv(f"{args.output_prefix}.pcoa.tsv", sep="\t")


if __name__ == "__main__":
    main()

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cluster_pairwise_analysis import (
    classical_pcoa,
    load_cluster_sequences,
    order_by_similarity,
)

st.set_page_config(page_title="Tea List Analysis", page_icon="📊", layout="wide")


def get_results_dir() -> Path | None:
    """Return the results directory supplied to the Streamlit application."""
    if len(sys.argv) < 2:
        return None
    return Path(sys.argv[1]).expanduser().resolve()


def get_data_path(results_dir: Path) -> Path | None:
    """Return the supplied target list, or one published with this run."""
    if len(sys.argv) >= 3:
        supplied_path = Path(sys.argv[2]).expanduser().resolve()
        if supplied_path.is_file():
            return supplied_path

    candidates = [results_dir / "teeliste.tsv"]
    return next((path for path in candidates if path.is_file()), None)


@st.cache_data
def load_tea_list(path: Path, modified_time_ns: int) -> pd.DataFrame:
    # modified_time_ns is intentionally part of the cache key so edits reload.
    source_columns = ["German name", "lat.", "taxid", "taxonomic_rank", "Grundlage"]
    data = pd.read_csv(
        path,
        sep="\t",
        header=0,
        names=source_columns,
        usecols=range(5),
        dtype={"taxid": "string"},
        engine="python",
    )
    data = data.rename(
        columns={
            "German name": "german_name",
            "lat.": "latin_name",
            "taxonomic_rank": "rank",
            "Grundlage": "source",
        }
    )
    required = {"german_name", "latin_name", "taxid", "rank", "source"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing columns: {', '.join(sorted(missing))}")

    for column in required:
        data[column] = data[column].fillna("Unknown").astype(str).str.strip()
    data["rank"] = data["rank"].replace("", "unknown").str.lower()
    return data


def get_available_runs(results_dir: Path) -> list[dict[str, str]]:
    """Return primer/database combinations represented by consensus files."""
    db_dir = results_dir / "build_db_taxids"
    known_dbs = {
        path.name.removesuffix(".db_taxids_counts.tsv")
        for path in db_dir.glob("*.db_taxids_counts.tsv")
    }
    runs = []
    for path in sorted((results_dir / "consensus").glob("*.cluster_consensus.tsv")):
        stem = path.name.removesuffix(".cluster_consensus.tsv")
        database = "unknown"
        primer = stem
        for candidate in sorted(known_dbs, key=len, reverse=True):
            suffix = f"_{candidate}"
            if stem.endswith(suffix):
                primer = stem[: -len(suffix)]
                database = candidate
                break
        runs.append({"stem": stem, "primer": primer, "database": database})
    return runs


def load_primer_counts(
    consensus_path: Path, modified_time_ns: int, run_stem: str
) -> pd.Series:
    """Count amplified accessions for the selected run without cross-run caching."""
    consensus_columns = [
        "cluster_id",
        "accession",
        "accession_taxid",
        "accession_name",
        "assigned_name",
        "assigned_taxid",
        "assigned_rank",
        "disambiguation",
    ]
    consensus = pd.read_csv(
        consensus_path,
        sep="\t",
        header=None,
        names=consensus_columns,
        usecols=["accession", "accession_taxid"],
        dtype={"accession_taxid": "string"},
    )
    return consensus["accession_taxid"].str.strip().value_counts()


@st.cache_data
def load_database_counts(database_path: Path, modified_time_ns: int) -> pd.Series:
    database = pd.read_csv(
        database_path,
        sep="\t",
        header=None,
        names=["taxid", "count"],
        dtype={"taxid": "string"},
    )
    database["taxid"] = database["taxid"].str.strip()
    database["count"] = pd.to_numeric(database["count"], errors="coerce").fillna(0)
    return database.groupby("taxid")["count"].sum()


@st.cache_data
def load_amplified_sequence_details(
    consensus_path: Path,
    consensus_modified_time_ns: int,
    parsed_path: Path,
    parsed_modified_time_ns: int,
    taxid: str,
) -> pd.DataFrame:
    """Join matching consensus accessions to their OBI-PCR amplicon sequences."""
    consensus_columns = [
        "cluster_id",
        "accession",
        "accession_taxid",
        "accession_name",
        "assigned_name",
        "assigned_taxid",
        "assigned_rank",
        "disambiguation",
    ]
    consensus = pd.read_csv(
        consensus_path,
        sep="\t",
        header=None,
        names=consensus_columns,
        dtype={"accession_taxid": "string", "assigned_taxid": "string"},
    )
    matching = consensus[
        consensus["accession_taxid"].str.strip().eq(str(taxid).strip())
    ].copy()
    if matching.empty:
        return pd.DataFrame()

    matching["accession_key"] = matching["accession"].astype(str).str.replace(
        r"\.\d+$", "", regex=True
    )
    if not parsed_path.is_file():
        matching["sequence_id"] = pd.NA
        matching["amplicon_length"] = pd.NA
        matching["amplicon_sequence"] = pd.NA
        return matching

    parsed_header = pd.read_csv(parsed_path, sep="\t", nrows=0).columns
    wanted_columns = [
        column
        for column in ["Sequence_ID", "Amplicon_Length", "Amplicon_Sequence"]
        if column in parsed_header
    ]
    if "Sequence_ID" not in wanted_columns:
        return matching

    parsed = pd.read_csv(parsed_path, sep="\t", usecols=wanted_columns)
    parsed["accession_key"] = parsed["Sequence_ID"].astype(str).str.replace(
        r"\.\d+$", "", regex=True
    )
    parsed = parsed.rename(
        columns={
            "Sequence_ID": "sequence_id",
            "Amplicon_Length": "amplicon_length",
            "Amplicon_Sequence": "amplicon_sequence",
        }
    )
    return matching.merge(parsed, on="accession_key", how="left")


@st.cache_data
def analyze_consensus_cluster(
    consensus_path: Path,
    consensus_modified_time_ns: int,
    parsed_path: Path,
    parsed_modified_time_ns: int,
    cluster_id: int,
    focal_taxid: str,
    cache_version: int = 3,
):
    metadata, distances = load_cluster_sequences(
        consensus_path, parsed_path, cluster_id, focal_taxid
    )
    coordinates, explained = classical_pcoa(distances)
    order = order_by_similarity(coordinates)
    return metadata, distances, order, coordinates, explained


st.title("Tea List Analysis")
st.caption("Explore the species, genera, varieties, and other taxonomic ranks in teeliste.tsv.")

results_dir = get_results_dir()
if results_dir is None or not results_dir.is_dir():
    st.error(
        "No valid results directory was supplied. Start the dashboard with "
        "`streamlit run app.py -- <results_dir> [teeliste.tsv]`."
    )
    st.stop()

data_path = get_data_path(results_dir)
if data_path is None:
    st.error(
        "No target-list TSV was supplied and teeliste.tsv was not found in the "
        f"selected results directory: {results_dir}"
    )
    st.stop()

try:
    tea_df = load_tea_list(data_path, data_path.stat().st_mtime_ns)
except (OSError, ValueError, pd.errors.ParserError) as exc:
    st.error(f"Could not load {data_path}: {exc}")
    st.stop()

rank_counts = (
    tea_df.groupby("rank", as_index=False)
    .size()
    .rename(columns={"size": "entries"})
    .sort_values("entries", ascending=False)
)

metric_columns = st.columns(4)
metric_columns[0].metric("Total entries", f"{len(tea_df):,}")
metric_columns[1].metric("Species", f"{int((tea_df['rank'] == 'species').sum()):,}")
metric_columns[2].metric("Genera", f"{int((tea_df['rank'] == 'genus').sum()):,}")
metric_columns[3].metric("Varieties", f"{int(tea_df['rank'].isin(['varietas', 'variety']).sum()):,}")

chart_col, summary_col = st.columns([2, 1])
with chart_col:
    rank_pie = px.pie(
        rank_counts,
        names="rank",
        values="entries",
        hole=0.35,
        title="Taxonomic rank distribution",
        color_discrete_sequence=px.colors.qualitative.Set2,
    )
    rank_pie.update_traces(textposition="inside", textinfo="percent+label")
    rank_pie.update_layout(legend_title_text="Rank", margin=dict(l=20, r=20, t=60, b=20))
    st.plotly_chart(rank_pie, width="stretch")

with summary_col:
    st.subheader("Rank counts")
    st.dataframe(rank_counts, width="stretch", hide_index=True)
    st.caption(f"Data source: {data_path}")

st.subheader("Database Representation vs Primer Amplification")
available_runs = get_available_runs(results_dir)
if not available_runs:
    st.info(f"No consensus results are available in {results_dir} for primer comparison.")
else:
    runs_df = pd.DataFrame(available_runs)
    database_col, selector_col = st.columns(2)
    with database_col:
        selected_database = st.selectbox(
            "Database",
            sorted(runs_df["database"].unique()),
            key="analysis_database",
        )
    database_runs = runs_df[runs_df["database"] == selected_database]
    with selector_col:
        selected_primer = st.selectbox(
            "Primer",
            sorted(database_runs["primer"].unique()),
            key="analysis_primer",
        )

    selected_stem = database_runs.loc[
        database_runs["primer"] == selected_primer, "stem"
    ].iloc[0]
    consensus_path = results_dir / "consensus" / f"{selected_stem}.cluster_consensus.tsv"
    database_path = (
        results_dir
        / "build_db_taxids"
        / f"{selected_database}.db_taxids_counts.tsv"
    )

    if not database_path.is_file():
        st.warning(f"Database count file not found: {database_path.name}")
    else:
        primer_counts = load_primer_counts(
            consensus_path,
            consensus_path.stat().st_mtime_ns,
            selected_stem,
        )
        database_counts = load_database_counts(
            database_path, database_path.stat().st_mtime_ns
        )
        comparison_df = tea_df.copy()
        comparison_df["database_sequences"] = (
            comparison_df["taxid"].map(database_counts).fillna(0).astype(int)
        )
        comparison_df["amplified_sequences"] = (
            comparison_df["taxid"].map(primer_counts).fillna(0).astype(int)
        )
        comparison_df["database_log"] = np.log10(
            comparison_df["database_sequences"] + 1
        )
        comparison_df["amplified_log"] = np.log10(
            comparison_df["amplified_sequences"] + 1
        )
        comparison_df["diagnosis"] = np.select(
            [
                comparison_df["amplified_sequences"] > 0,
                comparison_df["database_sequences"] == 0,
            ],
            ["Amplified", "Missing from database"],
            default="In database, not amplified",
        )
        amplified_target_count = int(
            (comparison_df["amplified_sequences"] > 0).sum()
        )

        diagnosis_colors = {
            "Amplified": "#2ca02c",
            "In database, not amplified": "#ff7f0e",
            "Missing from database": "#d62728",
        }
        comparison_chart = px.scatter(
            comparison_df,
            x="database_log",
            y="amplified_log",
            color="diagnosis",
            color_discrete_map=diagnosis_colors,
            hover_name="latin_name",
            hover_data={
                "german_name": True,
                "taxid": True,
                "database_sequences": ":,",
                "amplified_sequences": ":,",
                "database_log": False,
                "amplified_log": False,
            },
            labels={
                "database_log": "Database sequences (log10(count + 1))",
                "amplified_log": "Primer-amplified sequences (log10(count + 1))",
                "diagnosis": "Result",
            },
            title=f"{selected_primer} / {selected_database}",
        )
        comparison_chart.update_traces(marker={"size": 11, "opacity": 0.8})
        comparison_chart.update_layout(height=520)
        st.info(
            f"Current run: {selected_stem} · Amplified target taxa: "
            f"{amplified_target_count:,} of {len(comparison_df):,}"
        )
        st.plotly_chart(
            comparison_chart,
            width="stretch",
            key=f"database_primer_scatter_{selected_stem}",
        )
        st.caption(
            "Red: absent from the database. Orange: present in the database but not "
            "amplified. Green: amplified by the selected primer."
        )

        st.subheader("Amplified Sequence Explorer")
        taxon_options = comparison_df.sort_values(
            ["german_name", "latin_name"]
        ).to_dict("records")
        selected_taxon = st.selectbox(
            "Target taxon",
            taxon_options,
            index=0,
            format_func=lambda row: (
                f"{row['german_name']} — {row['latin_name']} (taxid {row['taxid']})"
            ),
            key="sequence_taxon",
        )
        parsed_path = results_dir / "parsed_obipcr" / f"{selected_stem}.tsv"
        sequence_details = load_amplified_sequence_details(
            consensus_path,
            consensus_path.stat().st_mtime_ns,
            parsed_path,
            parsed_path.stat().st_mtime_ns if parsed_path.is_file() else 0,
            str(selected_taxon["taxid"]),
        )

        if sequence_details.empty:
            st.info(
                f"No sequences from taxid {selected_taxon['taxid']} were amplified "
                f"by {selected_primer}."
            )
        else:
            display_columns = [
                "cluster_id",
                "accession",
                "sequence_id",
                "assigned_name",
                "assigned_taxid",
                "assigned_rank",
                "amplicon_length",
                "amplicon_sequence",
            ]
            display_columns = [
                column for column in display_columns if column in sequence_details.columns
            ]
            unique_accessions = sequence_details["accession"].nunique()
            sequences_available = (
                sequence_details.get("amplicon_sequence", pd.Series(dtype="string"))
                .notna()
                .sum()
            )
            count_col, sequence_col = st.columns(2)
            count_col.metric("Unique amplified accessions", f"{unique_accessions:,}")
            sequence_col.metric(
                "Amplicon sequences available", f"{int(sequences_available):,}"
            )
            st.dataframe(
                sequence_details[display_columns].drop_duplicates(),
                width="stretch",
                hide_index=True,
                height=420,
            )
            st.download_button(
                "Download amplified sequences",
                data=sequence_details[display_columns]
                .drop_duplicates()
                .to_csv(sep="\t", index=False)
                .encode("utf-8"),
                file_name=(
                    f"taxid_{selected_taxon['taxid']}_{selected_primer}_amplicons.tsv"
                ),
                mime="text/tab-separated-values",
            )

            st.subheader("Pairwise Cluster Analysis")
            st.caption(
                "Compares every unique amplicon in the selected consensus cluster. "
                "Exact duplicate sequences are collapsed, while all contributing "
                "accessions and taxids remain available in the metadata."
            )
            cluster_options = sorted(
                sequence_details["cluster_id"].dropna().astype(int).unique().tolist()
            )
            selected_cluster = st.selectbox(
                "Consensus cluster",
                cluster_options,
                key="pairwise_cluster_id",
            )
            try:
                (
                    cluster_metadata,
                    cluster_distances,
                    cluster_order,
                    cluster_coordinates,
                    cluster_explained,
                ) = analyze_consensus_cluster(
                    consensus_path,
                    consensus_path.stat().st_mtime_ns,
                    parsed_path,
                    parsed_path.stat().st_mtime_ns,
                    selected_cluster,
                    str(selected_taxon["taxid"]),
                    3,
                )
            except (OSError, ValueError, pd.errors.ParserError) as exc:
                st.error(f"Could not analyze cluster {selected_cluster}: {exc}")
            else:
                unique_count = len(cluster_metadata)
                focal_unique_count = int(
                    cluster_metadata["contains_focal_taxid"].sum()
                )
                other_taxids = set()
                for values in cluster_metadata["taxids"]:
                    other_taxids.update(str(values).split(";"))
                other_taxids.discard(str(selected_taxon["taxid"]))

                metric_one, metric_two, metric_three = st.columns(3)
                metric_one.metric("Unique amplicons", f"{unique_count:,}")
                metric_two.metric(
                    "Containing focal taxid", f"{focal_unique_count:,}"
                )
                metric_three.metric("Other taxids", f"{len(other_taxids):,}")

                if unique_count < 2:
                    st.info("This cluster contains only one unique amplicon sequence.")
                else:
                    plot_labels = [
                        ("★ " if row.contains_focal_taxid else "")
                        + f"{row.sequence_id} | n={row.record_count} | taxids={row.taxids}"
                        for row in cluster_metadata.itertuples()
                    ]
                    st.caption(
                        "Each row is one unique amplicon: the sequence ID, the number "
                        "of source records it collapses, and the taxids and accessions "
                        "it came from."
                    )
                    st.dataframe(
                        cluster_metadata[
                            [
                                "sequence_id",
                                "contains_focal_taxid",
                                "record_count",
                                "length",
                                "taxon_names",
                                "taxids",
                                "accessions",
                            ]
                        ],
                        width="stretch",
                        hide_index=True,
                        height=280,
                    )
                    heatmap_tab, pcoa_tab = st.tabs(["Distance heatmap", "PCoA"])
                    with heatmap_tab:
                        ordered_indices = list(cluster_order)
                        ordered_labels = [plot_labels[index] for index in ordered_indices]
                        ordered_distances = cluster_distances[np.ix_(
                            ordered_indices, ordered_indices
                        )]
                        heatmap = go.Figure(
                            go.Heatmap(
                                z=ordered_distances,
                                x=ordered_labels,
                                y=ordered_labels,
                                colorscale="Blues",
                                colorbar={"title": "Distance"},
                                zmin=0,
                                customdata=1 - ordered_distances,
                                hovertemplate=(
                                    "%{y} vs %{x}<br>Distance: %{z:.4f}"
                                    "<br>Approx. identity: %{customdata:.2%}<extra></extra>"
                                ),
                            )
                        )
                        heatmap.update_layout(
                            title="Pairwise normalized global edit distances",
                            height=max(520, 24 * unique_count),
                            xaxis={"tickangle": -70},
                            margin={"l": 90, "r": 20, "t": 70, "b": 130},
                        )
                        st.plotly_chart(heatmap, width="stretch")
                        st.caption(
                            "Rows and columns are ordered along the first PCoA axis, so "
                            "similar sequences sit next to each other and blocks of low "
                            "distance stand out. ★ marks a sequence shared with the "
                            "selected focal taxid."
                        )
                    with pcoa_tab:
                        pcoa_df = cluster_metadata.drop(columns=["sequence"]).copy()
                        pcoa_df["PCoA1"] = cluster_coordinates[:, 0]
                        pcoa_df["PCoA2"] = cluster_coordinates[:, 1]
                        pcoa_df["group"] = np.where(
                            pcoa_df["contains_focal_taxid"],
                            f"Contains taxid {selected_taxon['taxid']}",
                            "Other taxids",
                        )
                        pcoa = px.scatter(
                            pcoa_df,
                            x="PCoA1",
                            y="PCoA2",
                            color="group",
                            symbol="group",
                            size="record_count",
                            hover_name="sequence_id",
                            hover_data=[
                                "taxids",
                                "accessions",
                                "length",
                                "record_count",
                            ],
                            color_discrete_map={
                                f"Contains taxid {selected_taxon['taxid']}": "#d95f02",
                                "Other taxids": "#355f8a",
                            },
                            labels={
                                "PCoA1": f"PCoA 1 ({cluster_explained[0]:.1%})",
                                "PCoA2": f"PCoA 2 ({cluster_explained[1]:.1%})",
                                "group": "Sequence membership",
                            },
                            title="Two-dimensional view of pairwise sequence distances",
                        )
                        pcoa.update_traces(marker={"line": {"width": 1, "color": "#222"}})
                        pcoa.update_layout(height=560)
                        st.plotly_chart(pcoa, width="stretch")

                    metadata_download = cluster_metadata.drop(columns=["sequence"])
                    distance_download = pd.DataFrame(
                        cluster_distances,
                        index=cluster_metadata["sequence_id"],
                        columns=cluster_metadata["sequence_id"],
                    )
                    download_col_one, download_col_two = st.columns(2)
                    download_col_one.download_button(
                        "Download cluster metadata",
                        metadata_download.to_csv(sep="\t", index=False).encode("utf-8"),
                        file_name=f"cluster_{selected_cluster}_metadata.tsv",
                        mime="text/tab-separated-values",
                    )
                    download_col_two.download_button(
                        "Download distance matrix",
                        distance_download.to_csv(sep="\t").encode("utf-8"),
                        file_name=f"cluster_{selected_cluster}_distances.tsv",
                        mime="text/tab-separated-values",
                    )

                    st.markdown("**Copy or download amplicons for BLAST**")
                    all_fasta_records = []
                    for row in cluster_metadata.itertuples():
                        for accession in str(row.accessions).split(";"):
                            header = (
                                f">{accession}|cluster={selected_cluster}"
                                f"|sequence_group={row.sequence_id}|taxids={row.taxids}"
                            )
                            all_fasta_records.append(f"{header}\n{row.sequence}")
                    all_cluster_fasta = "\n".join(all_fasta_records) + "\n"

                    all_fasta_option = "__all_cluster_accessions__"
                    blast_options = [all_fasta_option] + cluster_metadata[
                        "sequence_id"
                    ].tolist()

                    def format_blast_option(sequence_id):
                        if sequence_id == all_fasta_option:
                            return f"All accessions in cluster {selected_cluster}"
                        return cluster_metadata.loc[
                            cluster_metadata["sequence_id"].eq(sequence_id), "label"
                        ].iloc[0]

                    blast_sequence_id = st.selectbox(
                        "Amplicon selection",
                        blast_options,
                        format_func=format_blast_option,
                        key="blast_sequence_id",
                    )
                    if blast_sequence_id == all_fasta_option:
                        blast_fasta = all_cluster_fasta
                        selected_fasta_filename = (
                            f"cluster_{selected_cluster}_all_accessions.fasta"
                        )
                    else:
                        blast_row = cluster_metadata.loc[
                            cluster_metadata["sequence_id"].eq(blast_sequence_id)
                        ].iloc[0]
                        selected_fasta_records = []
                        for accession in str(blast_row["accessions"]).split(";"):
                            blast_header = (
                                f">{accession}|cluster={selected_cluster}"
                                f"|sequence_group={blast_row['sequence_id']}"
                                f"|taxids={blast_row['taxids']}"
                            )
                            selected_fasta_records.append(
                                f"{blast_header}\n{blast_row['sequence']}"
                            )
                        blast_fasta = "\n".join(selected_fasta_records) + "\n"
                        selected_fasta_filename = (
                            f"cluster_{selected_cluster}_{blast_sequence_id}.fasta"
                        )
                    st.text_area(
                        "BLAST-ready FASTA (select all and copy)",
                        value=blast_fasta,
                        height=180,
                        key="blast_fasta_text",
                    )

                    fasta_col_one, fasta_col_two = st.columns(2)
                    fasta_col_one.download_button(
                        "Download current FASTA selection",
                        data=blast_fasta.encode("utf-8"),
                        file_name=selected_fasta_filename,
                        mime="text/plain",
                    )
                    fasta_col_two.download_button(
                        "Download all cluster accessions FASTA",
                        data=all_cluster_fasta.encode("utf-8"),
                        file_name=f"cluster_{selected_cluster}_all_accessions.fasta",
                        mime="text/plain",
                    )

st.subheader("Explore entries")
filter_col, search_col = st.columns([1, 2])
with filter_col:
    selected_ranks = st.multiselect(
        "Taxonomic rank",
        options=rank_counts["rank"].tolist(),
        default=rank_counts["rank"].tolist(),
    )
with search_col:
    search = st.text_input("Search German name, Latin name, taxid, or source")

filtered_df = tea_df[tea_df["rank"].isin(selected_ranks)].copy()
if search.strip():
    search_mask = filtered_df.astype(str).apply(
        lambda column: column.str.contains(search.strip(), case=False, na=False, regex=False)
    ).any(axis=1)
    filtered_df = filtered_df[search_mask]

st.dataframe(filtered_df, width="stretch", hide_index=True, height=440)
st.download_button(
    "Download filtered data",
    data=filtered_df.to_csv(sep="\t", index=False).encode("utf-8"),
    file_name="teeliste_filtered.tsv",
    mime="text/tab-separated-values",
)

st.page_link("app.py", label="Back to the taxonomy dashboard", icon="↩️")

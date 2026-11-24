import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from venn import venn

st.set_page_config(page_title="Primer Evaluation (Multi-Primer Comparison)", layout="wide")

if len(sys.argv) < 2:
    st.error("Usage: streamlit run app.py -- <path_to_results_dir>")
    st.stop()

DATA_DIR = Path(sys.argv[1])
if not DATA_DIR.exists():
    st.error(f"Directory not found: {DATA_DIR}")
    st.stop()

TEELISTE_PATH = (
    Path(sys.argv[2]) if len(sys.argv) >= 3 else DATA_DIR / "teeliste.tsv"
)

# Helper functions
def get_known_dbs(outdir):
    known = set()
    for subdir, suffix in [("db_distribution", ".db_distribution.tsv"), ("build_db_taxids", ".db_taxids_counts.tsv")]:
        d = Path(outdir) / subdir
        if d.exists():
            known.update(f.name.replace(suffix, "") for f in d.glob(f"*{suffix}"))
    return known

def split_primer_db(stem, known_dbs):
    for db in sorted(known_dbs, key=len, reverse=True):
        suffix = f"_{db}"
        if stem.endswith(suffix) and len(stem) > len(suffix):
            return stem[: -len(suffix)], db
    return stem, "unknown"

def get_available_runs(outdir):
    consensus_dir = Path(outdir) / "consensus"
    if not consensus_dir.exists():
        return []
    known_dbs = get_known_dbs(outdir)
    stems = sorted(f.name.replace(".cluster_consensus.tsv", "") for f in consensus_dir.glob("*.cluster_consensus.tsv"))
    runs = []
    for stem in stems:
        primer, db = split_primer_db(stem, known_dbs)
        runs.append({"stem": stem, "primer": primer, "db": db})
    return runs

def compact_plot(fig, height=360):
    fig.update_layout(
        height=height,
        margin=dict(l=20, r=20, t=60, b=90),
        xaxis_tickangle=-45
    )
    return style_plot_text(fig)

def style_plot_text(fig):
    fig.update_layout(
        font=dict(color="#111111", size=14),
        title_font=dict(color="#111111", size=20),
        legend=dict(font=dict(color="#111111", size=13)),
    )
    fig.update_xaxes(
        title_font=dict(color="#111111", size=15),
        tickfont=dict(color="#111111", size=13),
    )
    fig.update_yaxes(
        title_font=dict(color="#111111", size=15),
        tickfont=dict(color="#111111", size=13),
    )
    fig.update_traces(textfont=dict(color="#111111", size=13), selector=dict(type="bar"))
    fig.update_traces(textfont=dict(color="#111111", size=13), selector=dict(type="heatmap"))
    return fig

def plot_download_config(filename, width=4000, height=2600, scale=4):
    return {
        "toImageButtonOptions": {
            "format": "png",
            "filename": filename,
            "width": width,
            "height": height,
            "scale": scale,
        }
    }

def normalize_rank(rank):
    rank = str(rank).strip()
    return rank if rank else "unclassified"

RANK_ORDER = [
    "forma",
    "varietas",
    "subspecies",
    "species",
    "subsection",
    "section",
    "subgenus",
    "genus",
    "subtribe",
    "tribe",
    "subfamily",
    "family",
    "suborder",
    "order",
    "subclass",
    "class",
    "phylum",
    "kingdom",
    "superkingdom",
    "domain",
    "clade",
    "no rank",
    "unclassified",
]

def rank_sort_key(rank):
    rank = normalize_rank(rank)
    if rank in RANK_ORDER:
        return (RANK_ORDER.index(rank), "")
    return (len(RANK_ORDER), rank)

def get_rank_threshold_options(consensus_ranks):
    ordered = sorted({normalize_rank(rank) for rank in consensus_ranks}, key=rank_sort_key)
    if ordered:
        return ordered
    return ["species"]

def get_ranks_at_or_finer_than(rank_cutoff, consensus_ranks):
    cutoff_key = rank_sort_key(rank_cutoff)
    return {
        rank
        for rank in {normalize_rank(r) for r in consensus_ranks}
        if rank_sort_key(rank) <= cutoff_key
    }

def format_rank_list(ranks):
    ordered = sorted(ranks, key=rank_sort_key)
    return ", ".join(ordered[:8]) + ("..." if len(ordered) > 8 else "")

@st.cache_data
def get_consensus_ranks_for_db(data_dir, runs_list):
    ranks = set()
    consensus_dir = Path(data_dir) / "consensus"
    for run in runs_list:
        filepath = consensus_dir / f"{run['stem']}.cluster_consensus.tsv"
        if not filepath.exists():
            continue

        with open(filepath, "r", encoding="utf-8") as f:
            for line in f:
                columns = line.rstrip("\n").split("\t")
                if len(columns) > 6:
                    ranks.add(normalize_rank(columns[6]))

    return sorted(ranks)

@st.cache_data
def load_teeliste_data(t_path, db_path):
    # dtype="string" for taxid: rows with a blank taxid would otherwise upcast
    # the whole column to float64, turning e.g. "4442" into "4442.0" and
    # breaking every taxid match against accession/assigned_taxid downstream.
    t_df = pd.read_csv(t_path, sep="\t", header=0, dtype={"taxid": "string"})
    t_df = t_df.rename(columns={"German name": "german_name", "lat.": "latin_name"})
    t_df["rank"] = t_df.get("taxonomic_rank", "Unknown")
    t_df = t_df[["german_name", "latin_name", "taxid", "rank"]]

    if db_path.exists():
        db_df = pd.read_csv(db_path, sep="\t", header=None)
        db_df.columns = ["taxid", "count"]
    else:
        db_df = pd.DataFrame(columns=["taxid", "count"])
        
    return t_df, db_df

def accession_without_version(accession):
    """Match an OBI-PCR sequence id to the unversioned accession taxid table."""
    accession = str(accession).strip()
    base, separator, version = accession.rpartition(".")
    return base if separator and version.isdigit() else accession

@st.cache_data
def get_trimmed_sequence_counts(data_dir, runs_list, target_taxids, accession_taxid_path):
    """
    Count distinct database accessions that produced an in-silico amplicon for
    each primer and target taxid.  This is deliberately based on parsed OBI-PCR
    output, rather than consensus clusters: it answers how many source database
    sequences were trimmed/amplified, even when a later cluster has a broader
    taxonomic assignment.
    """
    accession_taxid_path = Path(accession_taxid_path)
    if not accession_taxid_path.exists():
        return {}

    accession_taxids = pd.read_csv(
        accession_taxid_path,
        sep="\t",
        header=None,
        names=["accession", "taxid"],
        dtype="string",
    )
    accession_to_taxid = dict(zip(
        accession_taxids["accession"].str.strip(),
        accession_taxids["taxid"].str.strip(),
    ))
    target_taxids = set(target_taxids)
    trimmed_counts = {}

    for run in runs_list:
        parsed_path = Path(data_dir) / "parsed_obipcr" / f"{run['stem']}.tsv"
        if not parsed_path.exists():
            continue

        sequence_ids = pd.read_csv(
            parsed_path,
            sep="\t",
            usecols=["Sequence_ID"],
            dtype="string",
        )["Sequence_ID"].dropna().unique()
        taxids = pd.Series(sequence_ids, dtype="string").map(
            lambda accession: accession_to_taxid.get(accession_without_version(accession))
        )
        counts = taxids[taxids.isin(target_taxids)].value_counts().astype(int).to_dict()
        trimmed_counts[run["primer"]] = counts

    return trimmed_counts

@st.cache_data
def parse_all_consensus_files_for_db(data_dir, runs_list, target_taxids):
    primer_data = {}
    consensus_dir = Path(data_dir) / "consensus"
    for run in runs_list:
        stem = run["stem"]
        primer = run["primer"]
        
        filepath = consensus_dir / f"{stem}.cluster_consensus.tsv"
        if not filepath.exists():
            continue
            
        records = []
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                columns = line.split('\t')
                if len(columns) > 6:
                    accession = columns[1].strip()
                    accession_taxid = columns[2].strip()
                    assigned_name = columns[4].strip() if len(columns) > 4 else "Unknown"
                    assigned_taxid = columns[5].strip()
                    rank = normalize_rank(columns[6] if len(columns) > 6 else "")
                    
                    if accession_taxid in target_taxids or assigned_taxid in target_taxids:
                        records.append({
                            "seqid": accession,
                            "acc_taxid": accession_taxid,
                            "ass_taxid": assigned_taxid,
                            "ass_name": assigned_name,
                            "rank": rank
                        })
        if records:
            primer_data[primer] = records
            
    return primer_data

def get_resolved_sequence_counts(primer_data, target_taxids):
    """
    Count accessions whose own taxid differs from a target taxid but whose
    consensus cluster is assigned to that target.  For example, sequences
    labelled as C. sinensis varieties that resolve to C. sinensis are counted
    here, separately from direct accession-taxid matches.
    """
    target_taxids = set(target_taxids)
    resolved_counts = {}

    for primer, records in primer_data.items():
        resolved_accessions = {}
        for record in records:
            assigned_taxid = record["ass_taxid"]
            if (
                assigned_taxid in target_taxids
                and record["acc_taxid"] != assigned_taxid
            ):
                resolved_accessions.setdefault(assigned_taxid, set()).add(record["seqid"])
        resolved_counts[primer] = {
            taxid: len(accessions)
            for taxid, accessions in resolved_accessions.items()
        }

    return resolved_counts

def filter_venn_data(primer_data, teeliste_taxids, selected_primers, selected_ranks, metric):
    venn_data = {}
    selected_rank_set = {normalize_rank(rank) for rank in selected_ranks}
    for primer in selected_primers:
        if primer not in primer_data:
            continue
            
        covered_taxids = set()
        for record in primer_data[primer]:
            if "All Ranks" not in selected_ranks and record["rank"] not in selected_rank_set:
                continue
                
            if ("Maximum Coverage" in metric or "Accession TaxIDs" in metric) and record["acc_taxid"] in teeliste_taxids:
                covered_taxids.add(record["acc_taxid"])
                
            if ("Maximum Coverage" in metric or "Assigned TaxIDs" in metric) and record["ass_taxid"] in teeliste_taxids:
                covered_taxids.add(record["ass_taxid"])
                
        venn_data[primer] = covered_taxids
    return venn_data

def rank_primers_by_coverage(primer_data, primers, target_taxids, selected_ranks, metric):
    rows = []
    total_targets = len(target_taxids)
    selected_rank_set = {normalize_rank(rank) for rank in selected_ranks}

    for primer in primers:
        covered_taxids = set()
        matching_records = 0
        rank_counts = {}

        for record in primer_data.get(primer, []):
            rank = normalize_rank(record["rank"])
            if rank not in selected_rank_set:
                continue

            matched = False
            if ("Maximum Coverage" in metric or "Accession TaxIDs" in metric) and record["acc_taxid"] in target_taxids:
                covered_taxids.add(record["acc_taxid"])
                matched = True
            if ("Maximum Coverage" in metric or "Assigned TaxIDs" in metric) and record["ass_taxid"] in target_taxids:
                covered_taxids.add(record["ass_taxid"])
                matched = True

            if matched:
                matching_records += 1
                rank_counts[rank] = rank_counts.get(rank, 0) + 1

        best_rank = min(rank_counts, key=rank_sort_key) if rank_counts else "none"
        covered_count = len(covered_taxids)
        rows.append({
            "Primer": primer,
            "Covered Target TaxIDs": covered_count,
            "Coverage %": round((covered_count / total_targets * 100) if total_targets else 0, 1),
            "Matching Records": matching_records,
            "Best Rank Seen": best_rank,
        })

    return (
        pd.DataFrame(rows)
        .sort_values(["Covered Target TaxIDs", "Matching Records", "Primer"], ascending=[False, False, True])
        .reset_index(drop=True)
    )

def get_resolution_ranks(consensus_ranks):
    return sorted({normalize_rank(rank) for rank in consensus_ranks}, key=rank_sort_key)

def choose_best_rank(current_rank, new_rank):
    if current_rank is None:
        return new_rank
    return new_rank if rank_sort_key(new_rank) < rank_sort_key(current_rank) else current_rank

def primer_resolution_table(primer_data, primers, target_taxids, selected_ranks, metric, resolution_ranks):
    selected_rank_set = {normalize_rank(rank) for rank in selected_ranks}
    resolution_ranks = list(resolution_ranks)
    rows = []

    for primer in primers:
        best_rank_by_taxid = {}

        for record in primer_data.get(primer, []):
            rank = normalize_rank(record["rank"])
            if rank not in selected_rank_set:
                continue

            matched_taxids = []
            if ("Maximum Coverage" in metric or "Accession TaxIDs" in metric) and record["acc_taxid"] in target_taxids:
                matched_taxids.append(record["acc_taxid"])
            if ("Maximum Coverage" in metric or "Assigned TaxIDs" in metric) and record["ass_taxid"] in target_taxids:
                matched_taxids.append(record["ass_taxid"])

            for taxid in matched_taxids:
                best_rank_by_taxid[taxid] = choose_best_rank(best_rank_by_taxid.get(taxid), rank)

        rank_counts = {rank: 0 for rank in resolution_ranks}
        broader_or_other = 0
        for rank in best_rank_by_taxid.values():
            if rank in rank_counts:
                rank_counts[rank] += 1
            else:
                broader_or_other += 1

        covered = len(best_rank_by_taxid)
        species_or_finer = sum(
            count for rank, count in rank_counts.items()
            if rank_sort_key(rank) <= rank_sort_key("species")
        )
        genus_or_finer = sum(
            count for rank, count in rank_counts.items()
            if rank_sort_key(rank) <= rank_sort_key("genus")
        )
        missing = len(target_taxids) - covered

        rows.append({
            "Primer": primer,
            **rank_counts,
            "broader_or_other": broader_or_other,
            "missing": missing,
            "covered": covered,
            "species_or_finer": species_or_finer,
            "genus_or_finer": genus_or_finer,
            "coverage_pct": round((covered / len(target_taxids) * 100) if target_taxids else 0, 1),
            "species_or_finer_pct": round((species_or_finer / len(target_taxids) * 100) if target_taxids else 0, 1),
        })

    sort_cols = ["genus_or_finer", "species_or_finer", "covered", "missing", "Primer"]
    return pd.DataFrame(rows).sort_values(sort_cols, ascending=[False, False, False, True, True]).reset_index(drop=True)

def get_nonzero_resolution_columns(resolution_df, resolution_ranks):
    return [
        col for col in list(resolution_ranks) + ["broader_or_other", "missing"]
        if col in resolution_df.columns and resolution_df[col].sum() > 0
    ]

def plot_resolution_heatmap(resolution_df, resolution_ranks):
    value_cols = get_nonzero_resolution_columns(resolution_df, resolution_ranks)
    if not value_cols:
        return None
    heatmap_df = resolution_df.set_index("Primer")[value_cols]
    zmax = heatmap_df.to_numpy().max()
    text_threshold = zmax * 0.55 if zmax else 0

    fig = go.Figure(
        data=go.Heatmap(
            z=heatmap_df.values,
            x=heatmap_df.columns,
            y=heatmap_df.index,
            colorscale="YlGnBu",
            colorbar=dict(title="Target taxa"),
        )
    )

    annotations = []
    for row_index, primer in enumerate(heatmap_df.index):
        for col_index, rank in enumerate(heatmap_df.columns):
            value = heatmap_df.iloc[row_index, col_index]
            annotations.append(
                dict(
                    x=rank,
                    y=primer,
                    text=str(value),
                    showarrow=False,
                    font=dict(
                        color="#ffffff" if value >= text_threshold and value > 0 else "#111111",
                        size=13,
                    ),
                )
            )

    fig.update_layout(
        title="Primer Resolution Matrix",
        height=max(420, 42 * len(resolution_df) + 160),
        margin=dict(l=20, r=20, t=60, b=60),
        annotations=annotations,
    )
    fig.update_xaxes(title_text="Resolution rank")
    fig.update_yaxes(title_text="Primer", autorange="reversed")
    return style_plot_text(fig)

def plot_resolution_stacked_bar(resolution_df, resolution_ranks):
    value_cols = get_nonzero_resolution_columns(resolution_df, resolution_ranks)
    if not value_cols:
        return None
    long_df = resolution_df.melt(
        id_vars=["Primer"],
        value_vars=value_cols,
        var_name="Rank",
        value_name="Target TaxIDs",
    )
    primer_totals = resolution_df.set_index("Primer")[value_cols].sum(axis=1)
    long_df["Primer Total"] = long_df["Primer"].map(primer_totals)
    long_df["Percentage"] = (
        long_df["Target TaxIDs"]
        .div(long_df["Primer Total"].replace(0, pd.NA))
        .mul(100)
        .fillna(0)
    )
    long_df["Label"] = long_df.apply(
        lambda row: f'{int(row["Target TaxIDs"]):,} ({row["Percentage"]:.1f}%)',
        axis=1,
    )
    fig = px.bar(
        long_df,
        y="Primer",
        x="Percentage",
        color="Rank",
        text="Label",
        orientation="h",
        category_orders={"Primer": resolution_df["Primer"].tolist()[::-1], "Rank": value_cols},
        title="Primer Resolution Efficiency by Best Rank",
        labels={"Percentage": "% of target taxa"},
        hover_data={
            "Target TaxIDs": ":,",
            "Primer Total": ":,",
            "Percentage": ":.1f",
            "Label": False,
        },
    )
    fig.update_traces(textposition="inside")
    fig.update_xaxes(range=[0, 100], ticksuffix="%")
    fig.update_layout(height=max(420, 38 * len(resolution_df) + 140), margin=dict(l=20, r=20, t=60, b=40))
    return style_plot_text(fig)

st.title("Primer Evaluation (Multi-Primer Comparison)")
st.markdown("Visually compare the performance of multiple primers using a Venn diagram and Coverage Matrix.")

all_runs = get_available_runs(DATA_DIR)
if not all_runs:
    st.warning("No consensus files found in the selected folder.")
    st.stop()

runs_df = pd.DataFrame(all_runs)

available_dbs = sorted(runs_df["db"].unique())
if not available_dbs:
    st.warning("No databases found.")
    st.stop()
selected_db = st.sidebar.selectbox("Select a database", available_dbs)

teeliste_path = TEELISTE_PATH
taxids_db_path = DATA_DIR / "build_db_taxids" / f"{selected_db}.db_taxids_counts.tsv"
accession_taxid_path = DATA_DIR / "build_db_taxids" / f"{selected_db}.accession_taxid.tsv"

try:
    db_runs = [r for r in all_runs if r["db"] == selected_db]

    if db_runs and teeliste_path.exists():
        teeliste_df, taxids_db_df = load_teeliste_data(teeliste_path, taxids_db_path)
        target_taxids = teeliste_df["taxid"].astype(str).str.strip().tolist()
        target_taxid_set = set(target_taxids)

        all_db_primers = sorted(list(set(r["primer"] for r in db_runs)))
        trimmed_sequence_counts = get_trimmed_sequence_counts(
            DATA_DIR,
            db_runs,
            target_taxid_set,
            accession_taxid_path,
        )
        
        st.markdown(f"**Filter Settings for {selected_db}**")
        col1, col2, col3 = st.columns(3)
        
        with col1:
            metric_options = ["Maximum Coverage (Accession + Assigned TaxIDs)", "Assigned TaxIDs Only", "Accession TaxIDs Only"]
            metric = st.radio("Coverage Metric:", metric_options, key="venn_metric")
            
        with col2:
            consensus_ranks = get_consensus_ranks_for_db(DATA_DIR, db_runs)
            rank_threshold_options = get_rank_threshold_options(consensus_ranks)
            default_rank = "species" if "species" in rank_threshold_options else rank_threshold_options[0]
            rank_cutoff = st.select_slider(
                "Rank cutoff:",
                options=rank_threshold_options,
                value=default_rank,
                key="rank_cutoff"
            )
            selected_ranks = get_ranks_at_or_finer_than(rank_cutoff, consensus_ranks)
            resolution_ranks = get_resolution_ranks(consensus_ranks)
            st.caption(f"Including {format_rank_list(selected_ranks)}")
            
        with col3:
            preferred_default_primers = [
                "ITS2_collapsed",
                "intron-region-of-a-transfer-RNA-gene-200-1000",
                "PCR_psbA-trnH",
            ]
            default_primers = [
                primer for primer in preferred_default_primers
                if primer in all_db_primers
            ]
            if not default_primers:
                default_primers = all_db_primers[:3]
            selected_eval_primers = st.multiselect("Select Primers (Max 6):", all_db_primers, default=default_primers, max_selections=6, key="venn_primers")

        primer_data = parse_all_consensus_files_for_db(DATA_DIR, db_runs, target_taxid_set)
        resolved_sequence_counts = get_resolved_sequence_counts(
            primer_data,
            target_taxid_set,
        )

        st.markdown("### Primer Resolution Detail")
        st.caption("Resolution detail has its own rank slider, independent of the Venn/Best Primers rank cutoff above.")
        default_detail_rank = "family" if "family" in resolution_ranks else resolution_ranks[-1]
        detail_rank_cutoff = st.select_slider(
            "Resolution detail broadest rank:",
            options=resolution_ranks,
            value=default_detail_rank,
            key="resolution_rank_cutoff",
        )
        detail_resolution_ranks = sorted(
            get_ranks_at_or_finer_than(detail_rank_cutoff, resolution_ranks),
            key=rank_sort_key,
        )
        st.caption(f"Resolution detail including {format_rank_list(detail_resolution_ranks)}")
        selected_resolution_primers = st.multiselect(
            "Select primers for resolution detail:",
            all_db_primers,
            default=all_db_primers,
            key="resolution_primers",
        )
        show_zero_resolution_ranks = st.checkbox(
            "Show ranks with zero coverage",
            value=False,
            key="show_zero_resolution_ranks",
        )

        if selected_resolution_primers:
            resolution_detail = primer_resolution_table(
                primer_data,
                selected_resolution_primers,
                target_taxid_set,
                set(detail_resolution_ranks),
                metric,
                detail_resolution_ranks,
            )
            plot_ranks = detail_resolution_ranks if show_zero_resolution_ranks else [
                col for col in detail_resolution_ranks
                if col in resolution_detail.columns and resolution_detail[col].sum() > 0
            ]
            fig_resolution_heatmap = plot_resolution_heatmap(resolution_detail, plot_ranks)
            fig_resolution_bar = plot_resolution_stacked_bar(resolution_detail, plot_ranks)
            if fig_resolution_heatmap is not None:
                st.plotly_chart(
                    fig_resolution_heatmap,
                    width="stretch",
                    config=plot_download_config("primer_resolution_matrix"),
                )
            if fig_resolution_bar is not None:
                st.plotly_chart(
                    fig_resolution_bar,
                    width="stretch",
                    config=plot_download_config("primer_resolution_stacked_bar"),
                )

            if show_zero_resolution_ranks:
                rank_cols = [
                    col for col in detail_resolution_ranks + ["broader_or_other", "missing"]
                    if col in resolution_detail.columns
                ]
            else:
                rank_cols = get_nonzero_resolution_columns(resolution_detail, detail_resolution_ranks)
            summary_cols = [
                "covered",
                "species_or_finer",
                "genus_or_finer",
                "coverage_pct",
                "species_or_finer_pct",
            ]
            display_cols = ["Primer"] + rank_cols + [
                col for col in summary_cols if col in resolution_detail.columns
            ]
            st.dataframe(resolution_detail[display_cols], width="stretch", hide_index=True)
        else:
            st.info("Select at least one primer to show resolution detail.")

        st.markdown("### Best Primers")
        primer_rankings = rank_primers_by_coverage(
            primer_data,
            all_db_primers,
            target_taxid_set,
            selected_ranks,
            metric
        )
        top_primers = primer_rankings.head(20)
        fig_best = px.bar(
            top_primers,
            x="Primer",
            y="Coverage %",
            hover_data=["Covered Target TaxIDs", "Matching Records", "Best Rank Seen"],
            title=f"Top primers at {rank_cutoff}-level or finer resolution",
            labels={"Coverage %": "Target coverage (%)"}
        )
        fig_best = compact_plot(fig_best)
        st.plotly_chart(
            fig_best,
            width="stretch",
            config=plot_download_config("best_primers", height=1200),
        )
        st.dataframe(primer_rankings, width="stretch", hide_index=True)

        if len(selected_eval_primers) > 0:
            venn_data = filter_venn_data(primer_data, target_taxid_set, selected_eval_primers, selected_ranks, metric)
            
            st.markdown("### Venn Diagram")
            if len(venn_data) < 2:
                st.info("Select at least 2 primers with data to show the Venn diagram.")
            elif not any(len(s) > 0 for s in venn_data.values()):
                st.info("None of the selected primers amplified any target taxa matching the filter criteria.")
            else:
                fig_venn, ax_venn = plt.subplots(figsize=(10, 6))
                try:
                    venn(venn_data, ax=ax_venn)
                    legend = ax_venn.get_legend()
                    if legend:
                        legend.set_bbox_to_anchor((1.05, 0.5))
                        legend._loc = 6
                    fig_venn.tight_layout()
                    st.pyplot(fig_venn, width="content")
                except Exception as e:
                    st.error(f"Could not draw Venn Diagram: {e}")
                    
                all_covered_taxids = set()
                for t_ids in venn_data.values():
                    all_covered_taxids = all_covered_taxids.union(t_ids)
                    
                covered_count = len(all_covered_taxids)
                omitted_count = len(target_taxids) - covered_count
                st.info(f"**Target Species:** {len(target_taxids)} | **Covered in Venn:** {covered_count} | **Omitted:** {omitted_count}")

            st.markdown("### Coverage Matrix")
            st.caption(
                "Each primer cell shows whether the taxon meets the selected coverage "
                "filter, then direct trimmed sequences / original database sequences. "
                "The 🔄 count is a different set of accessions whose consensus taxonomy "
                "resolved to the target taxid (for example, a variety resolving to its species)."
            )
            taxid_names = dict(zip(teeliste_df["taxid"].astype(str).str.strip(), teeliste_df["german_name"] + " (" + teeliste_df["latin_name"] + ")"))
            taxid_ranks = dict(zip(teeliste_df["taxid"].astype(str).str.strip(), teeliste_df["rank"].fillna("Unknown")))
            db_counts = dict(zip(taxids_db_df["taxid"].astype(str).str.strip(), taxids_db_df["count"]))
            primer_column_labels = {
                primer: f"{primer} trimmed / DB"
                for primer in selected_eval_primers
            }
            
            matrix_data = []
            for taxid in target_taxids:
                row = {
                    "TaxID": taxid,
                    "Name": taxid_names.get(taxid, "Unknown"),
                    "Rank": taxid_ranks.get(taxid, "Unknown"),
                    "DB Count": db_counts.get(taxid, 0)
                }
                total_covered = 0
                for primer in selected_eval_primers:
                    trimmed_count = trimmed_sequence_counts.get(primer, {}).get(taxid, 0)
                    resolved_count = resolved_sequence_counts.get(primer, {}).get(taxid, 0)
                    db_count = pd.to_numeric(db_counts.get(taxid, 0), errors="coerce")
                    db_count = 0 if pd.isna(db_count) else int(db_count)
                    if taxid in venn_data.get(primer, set()):
                        coverage_marker = "✅"
                        total_covered += 1
                    else:
                        coverage_marker = "❌"
                    resolved_suffix = f" + 🔄 {resolved_count:,}" if resolved_count else ""
                    row[primer_column_labels[primer]] = (
                        f"{coverage_marker} {trimmed_count:,} / {db_count:,}{resolved_suffix}"
                        if db_count else f"{coverage_marker} {trimmed_count:,}{resolved_suffix}"
                    )
                row["Total"] = total_covered
                matrix_data.append(row)
                
            matrix_df = pd.DataFrame(matrix_data).sort_values(by="Total", ascending=False).reset_index(drop=True)
            
            total_db_count = sum(pd.to_numeric(pd.Series(list(db_counts.values())), errors="coerce").fillna(0))
            total_row = {"TaxID": "TOTAL", "Name": "Total Species Covered", "Rank": "", "DB Count": total_db_count}
            
            for primer in selected_eval_primers:
                covered = len(venn_data.get(primer, set()))
                trimmed_total = sum(trimmed_sequence_counts.get(primer, {}).values())
                resolved_total = sum(resolved_sequence_counts.get(primer, {}).values())
                total_row[primer_column_labels[primer]] = (
                    f"{covered} / {len(target_taxids)} taxa | "
                    f"{trimmed_total:,} direct + 🔄 {resolved_total:,} resolved"
                )
            total_row["Total"] = None
            
            matrix_df = pd.concat([matrix_df, pd.DataFrame([total_row])], ignore_index=True)
            matrix_df["Total"] = matrix_df["Total"].astype(str)
            st.dataframe(matrix_df, width="stretch")

            st.markdown("### Amplicon Length Distribution of Target Taxa")
            # Real amplicon lengths come from the parsed OBI-PCR tables, as in app.py
            length_dir = Path(DATA_DIR) / "parsed_obipcr"
            if length_dir.exists():
                all_merged_data = []
                for primer in selected_eval_primers:
                    if primer not in primer_data:
                        continue
                    valid_records = []
                    selected_rank_set = {normalize_rank(rank) for rank in selected_ranks}
                    for record in primer_data[primer]:
                        if record["rank"] not in selected_rank_set:
                            continue
                        is_valid = False
                        if ("Maximum Coverage" in metric or "Accession TaxIDs" in metric) and record["acc_taxid"] in target_taxid_set:
                            is_valid = True
                        if ("Maximum Coverage" in metric or "Assigned TaxIDs" in metric) and record["ass_taxid"] in target_taxid_set:
                            is_valid = True
                        if is_valid:
                            valid_records.append(record)
                            
                    if not valid_records:
                        continue
                        
                    df_records = pd.DataFrame(valid_records)
                    
                    run_stems = [r["stem"] for r in db_runs if r["primer"] == primer]
                    if not run_stems:
                        continue
                    primer_stem = run_stems[0]
                    
                    length_file = length_dir / f"{primer_stem}.tsv"
                    if not length_file.exists():
                        continue
                        
                    len_df = pd.read_csv(length_file, sep="\t", usecols=["Sequence_ID", "Amplicon_Length"])
                    len_df = len_df.rename(columns={"Sequence_ID": "header", "Amplicon_Length": "length"})
                    len_df["seqid"] = len_df["header"].str.split().str[0]
                    len_df["match_id"] = len_df["seqid"].str.replace(r'\.\d+$', '', regex=True)
                    df_records["match_id"] = df_records["seqid"].str.replace(r'\.\d+$', '', regex=True)
                    
                    merged = pd.merge(len_df, df_records, on="match_id", how="inner")
                    merged["primer"] = primer
                    if not merged.empty:
                        all_merged_data.append(merged)
                        
                if all_merged_data:
                    final_df = pd.concat(all_merged_data, ignore_index=True)
                    fig_lengths = px.box(
                        final_df,
                        x="ass_name",
                        y="length",
                        color="primer",
                        title="Amplicon Length Distribution of Target Taxa",
                        labels={"length": "Amplicon length (bp)", "ass_name": "Target Taxon"}
                    )
                    fig_lengths.update_layout(xaxis_tickangle=-45)
                    fig_lengths = compact_plot(fig_lengths)
                    st.plotly_chart(
                        fig_lengths,
                        width="stretch",
                        config=plot_download_config("target_amplicon_lengths", height=1200),
                    )
                else:
                    st.info("No amplicon lengths found for the filtered target taxa.")
            else:
                st.info("Amplicon length directory not found.")
    else:
        if not teeliste_path.exists():
            st.info("No `teeliste.tsv` file found in the results directory.")
except Exception as e:
    st.error(f"Error rendering Multi-Primer Comparison: {e}")

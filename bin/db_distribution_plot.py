#!/usr/bin/env python3
"""
db_distribution_plot.py — BarBeQuE DB composition treemap
Reads the TSV from db_distribution.py and generates a standalone HTML treemap.

Usage:
    db_distribution_plot.py --input <db_distribution.tsv> --output <out.html>
"""

import argparse

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio

RANKS = ["class", "order", "family", "genus", "species"]



# I/O

def read_distribution(path):
    """Read db_distribution TSV, return dataframe."""

    df = pd.read_csv(path, sep="\t", dtype=str)
    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(0).astype(int)

    return df


def write_html(fig, path):
    """Write self-contained HTML (no CDN needed)."""

    pio.write_html(fig, file=path, full_html=True, include_plotlyjs=True, auto_open=False)

    print(f"[db_distribution_plot] HTML -> {path}")


# Figure

def build_treemap(df):
    """Build plotly treemap figure from distribution dataframe."""

    # drop rows with no class resolution
    plot_df = df[df["class"].notna() & (df["class"].str.strip() != "")].copy()

    # Plotly treemap requires that all entries in the hierarchy path are not None/NaN
    path_cols = ["class", "order", "family", "genus"]
    plot_df[path_cols] = plot_df[path_cols].fillna("Unknown")
    for col in path_cols:
        plot_df.loc[plot_df[col].str.strip() == "", col] = "Unknown"

    total_seqs = df["count"].sum()
    n_species  = df[df["resolved_rank"] == "species"].shape[0]
    n_genera   = plot_df["genus"].dropna().nunique()
    n_families = plot_df["family"].dropna().nunique()

    fig = px.treemap(
        plot_df,
        path=["class", "order", "family", "genus"],
        values="count",
        color="count",
        color_continuous_scale="Viridis",
        template="plotly_dark",
    )

    fig.update_traces(
        textinfo="label+value",
        hovertemplate=(
            "<b>%{label}</b><br>"
            "sequences: %{value:,}<br>"
            "%{percentRoot:.2%} of total"
            "<extra></extra>"
        ),
    )

    fig.update_layout(
        margin=dict(t=60, l=10, r=10, b=10),
        paper_bgcolor="#0d0d1a",
        coloraxis_showscale=False,
        title=dict(
            text=(
                f"BarBeQuE — DB Distribution&nbsp;&nbsp;"
                f"<span style='font-size:13px; color:#556'>"
                f"{total_seqs:,} sequences &nbsp;·&nbsp; "
                f"{n_species:,} species &nbsp;·&nbsp; "
                f"{n_genera:,} genus &nbsp;·&nbsp; "
                f"{n_families:,} families"
                f"</span>"
            ),
            font=dict(size=20, color="#ffffff", family="Courier New"),
            x=0.01,
        ),
    )

    return fig


def build_species_histogram(df):
    """Sequences-per-species histogram, log-spaced bins."""

    species = df[df["resolved_rank"] == "species"]["count"].values

    if len(species) == 0:
        return None

    bins        = np.logspace(np.log10(max(1, species.min())),
                              np.log10(species.max()), 60)
    hist, edges = np.histogram(species, bins=bins)
    singletons  = (species == 1).sum()
    pct         = 100 * singletons / len(species)

    fig = go.Figure(go.Bar(
        x=edges[:-1],
        y=hist,
        width=np.diff(edges),
        marker_color="steelblue",
        marker_line_width=0,
    ))

    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor="#0d0d1a",
        title=dict(
            text=(
                f"Sequences per species &nbsp;·&nbsp; "
                f"{len(species):,} species &nbsp;·&nbsp; "
                f"{singletons:,} singletons ({pct:.1f}%)"
            ),
            font=dict(size=16, color="#ffffff", family="Courier New"),
            x=0.01,
        ),
        xaxis=dict(title="Sequences per species", type="log"),
        yaxis=dict(title="Number of species"),
        bargap=0,
        margin=dict(t=60, l=60, r=20, b=60),
    )

    return fig


# Main

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--input",  required=True, help="db_distribution TSV")
    parser.add_argument("--output", required=True, help="output HTML path")
    args = parser.parse_args()

    df  = read_distribution(args.input)
    fig = build_treemap(df)
    write_html(fig, args.output)

    hist = build_species_histogram(df)
    if hist:
        write_html(hist, args.output.replace(".html", "_histogram.html"))


if __name__ == "__main__":
    main()
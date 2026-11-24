import csv
import json
import sys
from collections import Counter, defaultdict

#headers of cluster consensus tsv
headers = [
    "cluster_id",
    "accession",
    "accession_taxid",
    "accession_name",
    "consensus_name",
    "consensus_taxid",
    "consensus_rank",
    "disambiguation"
]
# Taxonomic ranks grouped by reporting resolution
species_ranks = {
    "strain",
    "isolate",
    "forma specialis",
    "forma",
    "varietas",
    "subspecies",
    "species",
}

genus_ranks = {
    "species subgroup",
    "species group",
    "series",
    "section",
    "subgenus",
    "genus",
}

family_ranks = {
    "subtribe",
    "tribe",
    "subfamily",
    "family",
}
#condensing ranks for better visualization  
order_ranks = {
    "superfamily",
    "parvorder",
    "infraorder",
    "suborder",
    "order",
}

class_ranks = {
    "superorder",
    "subcohort",
    "cohort",
    "infraclass",
    "subclass",
    "class",
}

phylum_ranks = {
    "superclass",
    "subphylum",
    "phylum",
}

kingdom_ranks = {
    "superphylum",
    "subkingdom",
    "kingdom",
    "superkingdom",
}

# Containers for results
amplified_taxa = set()
cluster_ranks = {}

# Track amplified taxa and consensus taxa separately
amplified_clusters = defaultdict(set)
consensus_clusters = defaultdict(set)

with open(sys.argv[1]) as file:
    reader = csv.reader(file, delimiter="\t")

    for row in reader:
        if len(row) != len(headers):
            raise ValueError(
                f"Expected {len(headers)} columns, but found {len(row)}."
            )
        data = dict(zip(headers, row))  
        # Store unique amplified taxa
        amplified_taxa.add((data["accession_taxid"], data["accession_name"]))

        # Store one consensus rank per cluster
        cluster_ranks[data["cluster_id"]] = data["consensus_rank"].lower()

        # Store where each amplified taxon occurs eg : where each Camellia sinensis sequnce went to 
        amplified_clusters[data["accession_name"]].add(
    (
        data["accession_taxid"],
        data["accession"],
        data["cluster_id"],
        data["consensus_taxid"],
        data["consensus_name"],
        data["consensus_rank"],
    )
)

        # Store where each consensus taxon occurs eg: Camellia sinensis was resolved to this cluster , 4442 might have varietas and it resolved to Camellia siensis
        consensus_clusters[data["consensus_name"]].add(
            (
                data["consensus_taxid"],
                data["cluster_id"],
                data["accession"],
                data["accession_taxid"],
                data["accession_name"],
                data["consensus_rank"],
            )
        )



       

#resolution bar chart 
def calculate_taxonomic_resolution(cluster_ranks):
    resolution_counts = Counter()

    for rank in cluster_ranks.values():

        if rank in species_ranks:
            resolution_counts["Species"] += 1

        elif rank in genus_ranks:
            resolution_counts["Genus"] += 1

        elif rank in family_ranks:
            resolution_counts["Family"] += 1

        elif rank in order_ranks:
            resolution_counts["Order"] += 1

        elif rank in class_ranks:
            resolution_counts["Class"] += 1

        elif rank in phylum_ranks:
            resolution_counts["Phylum"] += 1

        elif rank in kingdom_ranks:
            resolution_counts["Kingdom"] += 1

        elif rank in {"no rank", "unranked"}:
            resolution_counts["Unclassified"] += 1

        else:
            resolution_counts["Other"] += 1

    return resolution_counts


resolution_counts = calculate_taxonomic_resolution(cluster_ranks)

#exporting

# MultiQC table: unique amplified taxa
with open("amplified_taxa_mqc.tsv", "w", encoding="utf-8", newline="") as file:
    file.write('# id: "amplified_taxa"\n')
    file.write('# parent_id: "cluster_consensus"\n')
    file.write('# parent_name: "Cluster Consensus"\n')
    file.write('# section_name: "Amplified Taxa"\n')
    file.write('# plot_type: "table"\n')
    file.write('# pconfig:\n')
    file.write('#   no_violin: true\n')

    writer = csv.writer(file, delimiter="\t")

    writer.writerow([
        "Taxon",
        "TaxID",
    ])

    for taxid, taxon in sorted(amplified_taxa):
        writer.writerow([
            taxon,
            taxid,
        ])

# Duplicate taxonomic-resolution export intentionally disabled.
# The main reporting script already creates the wired, prefixed MultiQC file.
# with open("taxonomic_resolution_mqc.tsv", "w", encoding="utf-8", newline="") as file:
#     file.write('# id: "taxonomic_resolution"\n')
#     file.write('# parent_id: "cluster_consensus"\n')
#     file.write('# parent_name: "Cluster Consensus"\n')
#     file.write('# section_name: "Taxonomic Resolution"\n')
#     file.write('# plot_type: "bargraph"\n')
#     writer = csv.writer(file, delimiter="\t")
#     for category in [
#         "Species", "Genus", "Family", "Order", "Class", "Phylum",
#         "Kingdom", "Unclassified", "Other",
#     ]:
#         writer.writerow([category, resolution_counts[category]])

# Prepare information for searching taxa by scientific name or TaxID
taxon_search = defaultdict(
    lambda: {
        "amplified_clusters": {},
        "consensus_clusters": {},
    }
)
taxid_aliases = {}


# Store where sequences belonging to each amplified taxon went
for taxon_name, records in amplified_clusters.items():

    for (
        taxon_taxid,
        accession,
        cluster_id,
        consensus_taxid,
        consensus_name,
        consensus_rank,
    ) in records:

        # Search names are stored in lowercase
        search_entry = taxon_search[taxon_name.lower()]

        cluster = search_entry["amplified_clusters"].setdefault(
            cluster_id,
            {
                "consensus_name": consensus_name,
                "consensus_taxid": consensus_taxid,
                "consensus_rank": consensus_rank,
                "sequences": 0,
            }
        )

        cluster["sequences"] += 1

        # The TaxID points to the same search result
        taxid_aliases[taxon_taxid] = taxon_name.lower()


# Store which clusters resolved to each consensus taxon
# and which amplified sequences are inside those clusters
for consensus_name, records in consensus_clusters.items():

    for (
        consensus_taxid,
        cluster_id,
        accession,
        accession_taxid,
        accession_name,
        consensus_rank,
    ) in records:

        # Search names are stored in lowercase
        search_entry = taxon_search[consensus_name.lower()]

        cluster = search_entry["consensus_clusters"].setdefault(
            cluster_id,
            []
        )

        cluster.append(
            {
                "accession": accession,
                "accession_taxid": accession_taxid,
                "accession_name": accession_name,
                "consensus_rank": consensus_rank,
            }
        )

        # The TaxID points to the same search result
        taxid_aliases[consensus_taxid] = consensus_name.lower()

# Convert the prepared lookup to JSON for the MultiQC search
# Convert search data to JSON
search_json = json.dumps(taxon_search).replace("</", "<\\/")
taxid_aliases_json = json.dumps(taxid_aliases).replace("</", "<\\/")

# Create interactive MultiQC search section
html = f"""
<!--
id: "taxon_search"
parent_id: "cluster_consensus"
parent_name: "Cluster Consensus"
section_name: "Taxon Search"
description: "Search by scientific name or NCBI TaxID."
-->

<p><a href="https://www.ncbi.nlm.nih.gov/datasets/taxonomy/browser/" target="_blank" rel="noopener noreferrer">Open NCBI Taxonomy Browser</a></p>
<input id="taxon_search_query" placeholder="Enter taxon name or TaxID">
<button onclick="searchTaxon()">Search</button>

<h4>Where did the amplified sequences go?</h4>
<div id="amplified_results"></div>

<h4>Which clusters resolved to this taxon?</h4>
<div id="consensus_results"></div>

<script>

const data = {search_json};
const taxidAliases = {taxid_aliases_json};

function searchTaxon() {{

    // Read search term
    const query = document
        .getElementById("taxon_search_query")
        .value
        .trim()
        .toLowerCase();

    // Find matching taxon
    const result = data[query] || data[taxidAliases[query]];

    const amplified = document.getElementById("amplified_results");
    const consensus = document.getElementById("consensus_results");

    // Clear previous results
    amplified.innerHTML = "";
    consensus.innerHTML = "";

    // Stop if taxon is not found
    if (!result) {{
        amplified.innerHTML = "Taxon not found.";
        return;
    }}

    // Count resolution ranks for amplified sequences
    const rankCounts = {{}};

    // Show where amplified sequences went
    for (const [cluster, info] of Object.entries(result.amplified_clusters)) {{

        amplified.innerHTML +=
            "Cluster " + cluster +
            " | " + info.consensus_name +
            " | " + info.consensus_rank +
            " | " + info.sequences + " sequence(s)<br>";

        rankCounts[info.consensus_rank] =
            (rankCounts[info.consensus_rank] || 0) + 1;
    }}

    // Show resolution summary
    amplified.innerHTML += "<br><b>Resolution summary:</b><br>";

    for (const [rank, count] of Object.entries(rankCounts)) {{
        amplified.innerHTML +=
            rank + ": " + count + " cluster(s)<br>";
    }}

    // Count clusters resolved to searched taxon
    const resolvedClusters =
        Object.keys(result.consensus_clusters);

    consensus.innerHTML +=
        resolvedClusters.length +
        " cluster(s) resolved to this taxon.<br><br>";

    // Show contents of resolved clusters
    for (const [cluster, members] of Object.entries(
        result.consensus_clusters
    )) {{

        consensus.innerHTML +=
            "<b>Cluster " + cluster + "</b><br>";

        members.forEach(member => {{

            consensus.innerHTML +=
                member.accession +
                " | " +
                member.accession_name +
                " | TaxID " +
                member.accession_taxid +
                "<br>";
        }});

        consensus.innerHTML += "<br>";
    }}
}}

</script>
"""


# Write MultiQC HTML file
with open(
    "taxon_search_mqc.html",
    "w",
    encoding="utf-8"
) as file:
    file.write(html)

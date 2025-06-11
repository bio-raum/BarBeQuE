#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import argparse
import taxidTools


parser = argparse.ArgumentParser(description="Script options")
parser.add_argument("--clusters", help="Path to VSearch clusters table")
parser.add_argument("--table", help="Path to CRABS amplicon table")
parser.add_argument("--taxdump", help="Path to tthe taxdump folder")
parser.add_argument("--output", help="Path to output table")
args = parser.parse_args()


def load_taxonomy(taxdump):
    return taxidTools.read_taxdump(
        os.path.join(taxdump, "nodes.dmp"),
        os.path.join(taxdump, "rankedlineage.dmp"),
        os.path.join(taxdump, "merged.dmp")
    )


def parse_crabs(table):
    """parsing crabs in something more efficient than pandas"""
    crabs = {}
    with open(table, 'r') as fi:
        for line in fi.readlines():
            l = line.split('\t')
            crabs[l[0]] = l[1:]
    return crabs


def main(clusters, table, taxdump, output):
    tax = load_taxonomy(taxdump)

    # Parse clusters and lookup taxid for seqid in crabs table, link cluster id to taxid in hashmap
    crabs = parse_crabs(table)

    # Parse clusters
    clusters_seqid, clusters_taxid = {}, {}
    with open(clusters, 'r') as fi:
        for line in fi.readlines():
            l = line.split('\t')
            if l[0] == 'C':
                continue
            clusters_seqid[l[8]] = l[1]  # SeqID -> ClusterID
            clusters_taxid.setdefault(l[1], [])
            clusters_taxid[l[1]].append(crabs[l[8]][1])  # ClusterID -> taxid list

    # Get consensus, rank, and disambiguation
    clusters_cons = {}
    for k, v in clusters_taxid.items():
        disambiguation = ";".join([tax.getName(t) for t in v])
        cons = tax.lca(v, ignore_missing=True)
        clusters_cons[k] = [
            cons.name,
            cons.taxid,
            cons.rank,
            disambiguation,
        ]
    
    # Dump
    tab = '\t'
    with open(output, 'w') as fo:
        for k, v in crabs.items():
            line = "\t".join([k, '\t'.join(v[:-1]), '\t'.join(clusters_cons[clusters_seqid[k]]), v[-1]])
            fo.write(
                line
            )


if __name__ == "__main__":
    main(
        args.clusters,
        args.table,
        args.taxdump,
        args.output
    )

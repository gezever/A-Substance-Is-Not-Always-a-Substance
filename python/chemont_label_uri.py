#!/usr/bin/env python3
"""
Extract skos:prefLabel → concept URI mapping from ChemOnt SKOS TTL.

Usage:
  python chemont_label_uri.py <skos_ttl> <output_csv>

Output CSV columns: label, uri
"""

import csv
import sys
import rdflib

SKOS = rdflib.Namespace("http://www.w3.org/2004/02/skos/core#")


def main():
    if len(sys.argv) != 3:
        print("Usage: python chemont_label_uri.py <skos_ttl> <output_csv>")
        sys.exit(1)

    skos_ttl, output_csv = sys.argv[1], sys.argv[2]

    g = rdflib.Graph()
    g.parse(skos_ttl, format="turtle")

    rows = [
        (str(label), str(concept))
        for concept, label in g.subject_objects(SKOS.prefLabel)
    ]
    rows.sort(key=lambda x: x[0])

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["label", "uri"])
        writer.writerows(rows)

    print(f"Extracted {len(rows)} labels → {output_csv}")


if __name__ == "__main__":
    main()

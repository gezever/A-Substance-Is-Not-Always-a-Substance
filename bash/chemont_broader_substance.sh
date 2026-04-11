#!/bin/bash

set -e

OUTPUT="data/processed/chemont_ground_truth.csv"

/opt/apache-jena-5.6.0/bin/sparql \
    --results=CSV \
    --data=data/processed/rdf/substances_taxonomy.ttl \
    --query=bash/chemont_broader_substance.rq \
    > "$OUTPUT"

echo "Written to $OUTPUT ($(wc -l < "$OUTPUT") rows)"

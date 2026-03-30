#!/bin/bash


/opt/apache-jena-5.6.0/bin/riot --formatted=TURTLE data/processed/rdf/substances.ttl data/source/chemont/ChemOnt_2_1_skos.ttl > /tmp/merge.ttl
/opt/apache-jena-5.6.0/bin/sparql --results=TURTLE --data=/tmp/merge.ttl --query bash/merge.rq > /tmp/merge2.ttl
/opt/apache-jena-5.6.0/bin/sparql --results=TURTLE --data=/tmp/merge2.ttl --query bash/chemont_to_xkos.rq > /tmp/merge3.ttl
/opt/apache-jena-5.6.0/bin/riot --formatted=TURTLE data/processed/rdf/substances.ttl /tmp/merge2.ttl /tmp/merge3.ttl > data/processed/rdf/substances_taxonomy.ttl
/opt/apache-jena-5.6.0/bin/riot --formatted=TURTLE data/processed/rdf/substances_taxonomy.ttl  data/processed/rdf/clusters.ttl  data/processed/rdf/linkability_taxonomy.ttl data/processed/rdf/annotations.ttl > /tmp/merge5.ttl

cp /tmp/merge5.ttl data/processed/rdf/substances_taxonomy.ttl


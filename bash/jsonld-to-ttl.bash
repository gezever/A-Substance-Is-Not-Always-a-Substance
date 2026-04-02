#!/bin/bash


/opt/apache-jena-5.6.0/bin/riot --formatted=TURTLE data/processed/rdf/classyfire.jsonld > data/processed/rdf/classyfire.ttl
/opt/apache-jena-5.6.0/bin/riot --formatted=JSONLD data/processed/rdf/classyfire.jsonld > /tmp/classyfire.jsonld
mv /tmp/classyfire.jsonld data/processed/rdf/classyfire.jsonld


/opt/apache-jena-5.6.0/bin/riot --formatted=TURTLE data/processed/rdf/substances.jsonld > data/processed/rdf/substances.ttl
/opt/apache-jena-5.6.0/bin/riot --formatted=JSONLD data/processed/rdf/substances.jsonld > /tmp/substances.jsonld
mv /tmp/substances.jsonld data/processed/rdf/substances.jsonld
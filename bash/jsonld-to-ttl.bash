#!/bin/bash

cd data/processed/rdf
/opt/apache-jena-5.6.0/bin/riot --formatted=TURTLE substances.jsonld > substances.ttl
/opt/apache-jena-5.6.0/bin/riot --formatted=TURTLE classyfire.jsonld > classyfire.ttl


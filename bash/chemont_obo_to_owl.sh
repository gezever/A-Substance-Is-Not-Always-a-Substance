#!/bin/bash

# this script uses robot, a java tool https://github.com/ontodev/robot
# git clone and mvn clean install before use
# sparql is a jena cli tool
# transforms a obo file in http://classyfire.wishartlab.com/system/downloads/1_0/chemont/ChemOnt_2_1.obo.zip to csv
curl --request GET -sL \
     --url 'http://classyfire.wishartlab.com/system/downloads/1_0/chemont/ChemOnt_2_1.obo.zip'\
     --output '../data/source/chemont/ChemOnt_2_1.obo.zip'
unzip  ../data/source/chemont/ChemOnt_2_1.obo.zip -d ../data/source/chemont/

rm  ../data/source/chemont/ChemOnt_2_1.obo.zip

echo "convert obo to owl-model"

robot convert -i '../data/source/chemont/ChemOnt_2_1.obo'  --format ttl -o '../data/source/chemont/ChemOnt_2_1.ttl'

echo "convert ttl to csv"

sparql --results=CSV --data=../data/source/chemont/ChemOnt_2_1.ttl  --query chemont_to_table.rq  > '../data/source/chemont/ChemOnt_2_1.csv'

sparql --results=TURTLE --data=../data/source/chemont/ChemOnt_2_1.ttl  --query chemont_to_skos.rq > '../data/source/chemont/ChemOnt_2_1_skos.ttl'

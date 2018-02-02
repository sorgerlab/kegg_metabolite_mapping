Introduction
============

This tools maps metabolite names from a CSV file in a format produced by a
particular assay and tries to map them to KEGG Compound IDs via the Human
Metabolite Database and some further manual curation.


Setup
=====

Download the "All Metabolites" file from HMDB at http://www.hmdb.ca/downloads .
Unzip it into the `input` directory to obtain `hmdb_metabolites.xml` (1.6GB).


Usage
=====

python map_metabolites.py input/data.csv > output/mapping.csv


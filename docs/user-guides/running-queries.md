# Running Queries: Overview

TODO

## Run the VCF dimensions command

For performance reasons and to ensure query feasibility, key metadata from the VCF files must first be cached in DivBase.
For each VCF in the DivBase project, the system extracts file name, number and name of all samples, number and name of all scaffolds, number of variants,

DivBase will use this whenever a user submits a query to the project. For instance, the user might make a sample metadata filtering that results in only certain samples. The system knows which file names each requested sample are located in, and will ensure that only those files will be transferred to the worker.

example
dimensions show

## Side car sample metadata queries

for sample metadata linked VCF queries, it can be good to do a dry run first

how to write the sample metadata TSV (+template)

how to use dimensions show to get all samples in the project

how to write query

## VCF queries

can be run with or without sample metadata filtering

how to write a query

how to optimize a query (see separate markdown)

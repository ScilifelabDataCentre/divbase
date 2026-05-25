# How to create efficient DivBase queries

TODO add more tips here

## General advice

Split large VCF files into smaller files. For instance by chromosome. If the assembly does not have chromosome level contiguity, we suggest to store a range of scaffolds per file.

If you want to check out all samples across all files, it will be more efficient to download the files than to run a query

## bcftools command pipe order

For example, it is typically faster to first subset on genomic range (=drop rows in the VCF that do not fulfil the range filter) and then subset on samples than to do the reverse. This is because sample subset will require writes to each row in the VCF; by first reducing the number of rows, there will be fewer write operations and thus a faster operation.

## Read next

- [Tutorial: Running a query on a public dataset](tutorial-query-on-public-data.md)

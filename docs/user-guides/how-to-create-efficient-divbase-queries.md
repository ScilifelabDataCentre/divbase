# How to create efficient DivBase queries

Split large VCF files into smaller files. For instance by chromosome. If the assembly does not have chromosome level contiguity, we suggest to store a range of scaffolds per file.

divbase bcftools "pipes"

If subsets based on variant range are included, run then first before any sample subsets. subsetting on samples requres inspecting and potentially updating every row in the VCF (N+1 problem); subsetting on samples will be faster on smaller datasets and therefore better to run after all row-based subsets have been done.

if you want all samples across all files, it will be more efficient to download the files than to run a query

## bcftools command pipe order

For example, it is typically faster to first subset on genomic range (=drop rows in the VCF that do not fullfil the range filter) and then subset on samples than to do the reverse. This is because sample subset will require writes to each row in the VCF; by first reducing the number of rows, there will be fewer write operations and thus a faster operation.

## Read next

- [Tutorial: Running a query on a public dataset](tutorial-query-on-public-data.md)

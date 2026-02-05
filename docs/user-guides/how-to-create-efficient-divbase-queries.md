# How to create efficient DivBase queries

Split large VCF files into smaller files. For instance by chromosome. If the assembly does not have chromosome level contiguity, we suggest to store a range of scaffolds per file.

divbase bcftools "pipes"

If subsets based on variant range are included, run then first before any sample subsets. subsetting on samples requres inspecting and potentially updating every row in the VCF (N+1 problem); subsetting on samples will be faster on smaller datasets and therefore better to run after all row-based subsets have been done.

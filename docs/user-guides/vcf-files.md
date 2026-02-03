# VCF Files

TODO

- what format? VCF.GZ
  - can it be compressed with gzip, or do we recommend bcftools

- does it need to have an index file
  - no, divbase creates these during the queries. DivBase uses CSI indexed instead of TBI indexes since the former can accomate larger genome assemblies (TODO: find the reference for the memory size the two indexes can handle)

- file size limitation?
  - each DivBsae project has file quota. ask the staff if you have questions about this
  - comment that splitting large VCFs by chromosome will likely be more efficient for DivBase.
    - we can provide suggestions on how to split
    - example script from mouse VCF benchmarking to split a larger VCF into 20 smaller by chromosome: <https://github.com/ScilifelabDataCentre/divbase/blob/per-task-cpu-ram-metrics/scripts/benchmarking/split_mouse_vcf_per_scaffold.sh>

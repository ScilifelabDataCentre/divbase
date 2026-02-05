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

## Managing VCF files in DivBase

DivBase is built to help users manage their VCF data. For this to work well, users need to consider the following when uploading data

DivBase is built on with good data management priciples in mind (TODO: add an external ref for this? RDMkit, Turingway Etc.)

- Data should not be duplicated in a DivBase project.

  - The same sample or the same variant can be present in multiple files, but the same sample and variant cannot!

    Toy examples. For the sake of demonstration, let's assume that we have two VCF files that only contain four samples and one variant.

    The following example is compatible with DivBase: the same samples are present in both files but the files describe different variants. Furhtermore, DivBase is built to handle that big VCF files are be split into smaller files by chromosome or by scafold/contig; in fact, it is often more performant to store the data in that way.

    ```
    # >zcat file_1.vcf.gz | grep -v "^#" | head
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE001   SAMPLE002   SAMPLE003   SAMPLE004
    1      12345        1_12345     T       C       .       .       PR;AN=10;AC=0   GT      0/0     1/0     0/0     0/

    # >zcat file_2.vcf.gz | grep -v "^#" | head
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE001   SAMPLE002   SAMPLE003   SAMPLE004
    1      56789        1_56789     T       C       .       .       PR;AN=10;AC=0   GT      0/0     0/0     0/0     0/0
    ```

    The following is also compatible with DivBase: the same variant is present in both files but the samples it describes are different between the two files.

    ```
    # >zcat file_1.vcf.gz | grep -v "^#" | head
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE001   SAMPLE002   SAMPLE003   SAMPLE004
    1      12345        1_12345     T       C       .       .       PR;AN=10;AC=0   GT      0/0     1/0     0/0     0/

    # >zcat file_2.vcf.gz | grep -v "^#" | head
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE005   SAMPLE006   SAMPLE007   SAMPLE008
    1      12345        1_12345     T       C       .       .       PR;AN=10;AC=0   GT      0/0     1/0     0/0     0/
    ```

    THE FOLLOWING IS NOT COMPATIBLE with DivBase: here one sample (SAMPLE001) occurs for the exact same variant in both files. DivBase will not be able to process this.

    ```
    # >zcat file_1.vcf.gz | grep -v "^#" | head
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE001   SAMPLE002   SAMPLE003   SAMPLE004
    1      12345        1_12345     T       C       .       .       PR;AN=10;AC=0   GT      0/0     1/0     0/0     0/

    # >zcat file_2.vcf.gz | grep -v "^#" | head
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE001   SAMPLE005   SAMPLE006   SAMPLE007
    1      12345        1_12345     T       C       .       .       PR;AN=10;AC=0   GT      0/0     1/0     0/0     0/
    ```

- VCF and sidecar sample metadata files are versioned in the DivBase project

  - The system will use the latest version for all queries.
    - TODO can file version be set with a flag?

  -

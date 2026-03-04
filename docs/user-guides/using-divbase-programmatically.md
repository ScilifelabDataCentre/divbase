# Tips for using DivBase programmatically

TODO - Should be a set of tips for users who want to use DivBase in scripts etc..

E.G. how to wait for a query to be complete, and download the query results programmatically.

## Parse divbase-cli files ls/info commands programmatically

1. You can make the output of the `divbase-cli files info` and `divbase-cli files ls` commands in TSV format for easier parsing. Use the `--tsv` flag:

    ```bash
    divbase-cli files info --tsv
    divbase-cli files ls --tsv
    ```

2. Rather than first downloading a file, you can stream a file from the command line and pipe it into other tools for processing directly without saving it to disk.

    ```bash
    divbase-cli files stream my_file.vcf.gz | zcat | less
    ```

    !!! Info
        BCFTools accepts stdin as input, so you can also pipe a VCF file directly into BCFTools without saving it first:

        ```bash
        divbase-cli files stream my_file.vcf.gz | bcftools view -h -
        ```

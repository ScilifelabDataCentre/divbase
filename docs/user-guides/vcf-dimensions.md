# VCF Dimensions caching

TODO

This a prerequisite step before running any queries

The dimensions update job will:

- For each VCF files in the project data store:
  - store filename and data store file version ID
  - store timestamp of last update [TODO IS THIS LAST TIME DIMENSIONS WAS RUN OR LAST TIME FILE WAS UPDATED?]
  - Fetche and store all sample names from the VCF file
  - Fetch and store all  chromosome/scaffold/contig names
  - count the number of samples and variants in the VCF file
- not consider results files produced by DivBase
  - by filename (`result_of_job_` prefix)
  - as an additional check, looks for bcftools annotate stamp in header (all DivBase results files are branded like this)

Run `divbase-cli dimensions update` whenever:

- A new VCF file is uploaded to the project
- A new version of an existing VCF file is uploaded to the project

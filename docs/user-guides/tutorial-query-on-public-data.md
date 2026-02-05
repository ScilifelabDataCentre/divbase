# Tutorial: Running a DivBase query on a public dataset

This tutorial assumes that you have an account on DivBase and have a membership in a DivBase project with at least an EDIT role (i.e. can upload files and run queries).

We will use a mouse (_Mus musculus_) data set availalbe on the European Nucleotide Archive: <https://www.ebi.ac.uk/ena/browser/view/ERZ022025>. It a 5.5 Gb VCF.gz file that contains 18 samples and 66,007,044 variants.

## 1. Obtain the data and upload it to your DivBase project

First we need download the files from EVA to our local computer. We can do that from the ENA homepage with a web browser, or from the terminal with a tool like `curl` or `wget`:

```bash
curl -o mgp.v3.snps.rsIDdbSNPv137.vcf.gz ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ022/ERZ022025/mgp.v3.snps.rsIDdbSNPv137.vcf.gz
```

On the terminal, log in to DivBase:

```bash
divbase-cli auth login <USER_EMAIL>
```

We then need to upload the file to the DivBase project. If the project you want to upload to is your default project, you can skip `--project <YOUR_DIVBASE_PROJECT_NAME>`, but otherwise you will need to specify that.

```bash
divbase-cli files upload mgp.v3.snps.rsIDdbSNPv137.vcf.gz --project <YOUR_DIVBASE_PROJECT_NAME>
```

For the sake of demonstration, we can create a sidecar metadata file for the 18 samples in this dataset so that we can use that for the query. In fact, a file prepared for this demo can be downloaded from the DivBase repo with:

```bash
curl -o tutorial_mock_metadata_mgpv3snps.tsv https://raw.githubusercontent.com/ScilifelabDataCentre/divbase/refs/heads/docs-on-queries/tests/fixtures/tutorial_mock_metadata_mgpv3snps.tsv
```

TODO update the TSV url when it has been merged to main.

And then upload it to the DivBase project with:

```bash
divbase-cli files upload tutorial_mock_metadata_mgpv3snps.tsv --project <YOUR_DIVBASE_PROJECT_NAME>
```

With this, we are set in terms of the files we need in the DivBase project to be able to run the query.

### A comment on the mock sidecar metadata TSV

If we look at the first few lines of this TSV, we see that the first column has the SAMPLE IDs from the VCF. There is a column for population grouping, and a column for the geographical area the sample comes from. All of this is made up for the tutorial, just to illustrate what a sidecar metadata file can look like. We will use the Area column later on to subset on samples.

```bash
head tests/fixtures/tutorial_mock_metadata_mgpv3snps.tsv

#Sample_ID      Population      Area
129P2   1       North
129S1   2       East
129S5   3       South
AJ      4       West
```

For comparison, we can inspect all the samples in the VCF with:

```bash
# On MacOS, use gzcat instead of zcat
zcat mgp.v3.snps.rsIDdbSNPv137.vcf.gz | grep "#CHROM"

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  129P2   129S1   129S5   AJ      AKRJ    BALBcJ  C3HHeJ  C57BL6NJ        CASTEiJ CBAJ    DBA2J   FVBNJ   LPJ     NODShiLtJ       NZOHlLtJ   PWKPhJ  SPRETEiJ        WSBEiJ
```

## 2. Submit a job to cache the VCF dimensions in the DivBase backend

DivBase needs this to make calculations. It can be seen as a pre-processing step that needs to be run every time a new version of a VCF is uploaded to DivBase.

```bash
divbase-cli dimensions update --project <YOUR_DIVBASE_PROJECT_NAME>
```

The task should now have been submitted. The terminal prints a DivBase Task ID with a message like this:

```
# Example with Task ID 102
Job submitted successfully with task id: 102
```

Make note of the Task ID integer for for now; we will use it later in the tutorial for a variable named `<THE_JOB_ID_OF_THE_QUERY>`.

We need to wait until this has finished before we can send the actual DivBase query.

```bash
divbase-cli task-history user
```

## 3. Submit a query job

```bash
divbase-cli query bcftools-pipe --tsv-filter 'Area:North,East' --command 'view -s SAMPLES; view -r 1:15000000-25000000' --metadata-tsv-name tutorial_mock_metadata_mgpv3snps.tsv --project <YOUR_DIVBASE_PROJECT_NAME>
```

note that we are specifically using the `tutorial_mock_metadata_mgpv3snps.tsv` sample metadata. If this is not specified, the query will default to `sample_metadata.tsv`.

Depending on queue length, this job will take some time to run. Check that it has started, and then leave it running for half an hour.

```bash
divbase-cli task-history user
```

## 4. Download the results file

```bash
divbase-cli file download result_of_job_<THE_JOB_ID_OF_THE_QUERY>.vcf.gz--project <YOUR_DIVBASE_PROJECT_NAME>
```

We can now run some quick sanity-checks on the result file.

```bash
# On MacOS, use gzcat instead of zcat
zcat result_of_job_<THE_JOB_ID_OF_THE_QUERY>.vcf.gz | grep -v "^#" |wc -l

# Expected terminal output:
297415
```

If you want, you can compare this to the original file:

```bash
# Note! This will take a little time since this file has many rows
zzcat mgp.v3.snps.rsIDdbSNPv137.vcf.gz | grep -v "^#" |wc -l

# Expected terminal output:
66007044
```

As for the samples, the expected result is the following:

```bash
zcat mgp.v3.snps.rsIDdbSNPv137.vcf.gz | grep "#CHROM"
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  129P2   129S1   129S5   AJ      AKRJ    BALBcJ  C3HHeJ  C57BL6NJ        CASTEiJ CBAJ    DBA2J   FVBNJ   LPJ     NODShiLtJ       NZOHlLtJ   PWKPhJ  SPRETEiJ        WSBEiJ

zcat result_of_job_<THE_JOB_ID_OF_THE_QUERY>.vcf.gz | grep "#CHROM"
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  129P2   129S1   AKRJ    BALBcJ  CASTEiJ CBAJ    LPJ     NODShiLtJ       SPRETEiJ        WSBEiJ
```

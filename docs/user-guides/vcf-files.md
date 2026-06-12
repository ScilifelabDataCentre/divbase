# Working with VCF Files in DivBase

This page describes how to prepare and organize VCF files in a DivBase project to ensure that DivBase queries will run reliably.

If you are looking for the VCF query syntax, see the [DivBase VCF query syntax](vcf-query-syntax.md) guide.

## 1. VCF file requirements for DivBase

The VCF files need to be sorted by coordinate and compressed with [BGZF](https://en.wikipedia.org/wiki/BGZF) and have the extension (`.vcf.gz`). This can be achieved with the `bgzip` tool which comes bundled with `bcftools`.

This can be achieved with:

```bash
bcftools sort input.vcf | bgzip -c > output.vcf.gz
```

The data in the VCF files also need follow the considerations for [how to organize multiple VCF files in one DivBase project](#2-how-to-organize-multiple-vcf-files-in-a-divbase-project) discussed below.

The DivBase server will handle creation of CSI indexes during VCF processing, so users do not need to upload index files manually.

!!! note
    All VCF files in a DivBase project should descibe the same version of the same reference genome. If you need to upload data for another reference genome, please contact the DivBase staff about creating a new DivBase project.

## 2. How to organize multiple VCF files in a DivBase project

A core idea behind DivBase is that the files in a project are the single source-of-truth for the project’s VCF data. It probably goes without saying that data should not be duplicated between VCF files in the project, and this section will describe what that means.

### 2.1. Performance considerations: monolithic VCF files versus split files VCF

For large VCF data, it will often be faster for the DivBase server to run queries when data is split into multiple VCF files by scaffold/chromosome rather than one very large monolithic file. The reason for this is that the server will need to fetch a copy of each VCF file from the project's data storage, and if the data is split across more files chances are that not all files are needed by the queries and thus are smaller and faster to download for the server.

There is no file size limit per VCF file, but smaller files will be faster to process due to how the DivBase server handles files. However, each DivBase project has storage quota limits for the total file content in the project. Contact the DivBase staff if needed.

!!! tip
    Split your VCF files by chromsome or by a small set of scaffolds to improve the chances of faster DivBase queries. Ensure that the Sample column order is the same in all split files and that no exact row (variant and samples) is duplicated across the files.

### 2.2. Sample set overlap - important for DivBase compatibility

The main constraints come from `bcftools concat` and `bcftools merge`, which DivBase uses internally when a query needs data from multiple files. This has to do with what combination of samples are found in each VCF files. We will refer to this as _sample sets_ (see also the [bcftools merge manual](https://samtools.github.io/bcftools/bcftools.html#merge)).

To illustrate this, let's consider an example of few VCF files that each contain certain sample sets in a fixed order based on the column order in the `#CHROM heading` of the VCF files:

- file_A.vcf.gz samples: `{S1, S2, S3}`

- file_B.vcf.gz samples: `{S1,S2,S3}`

- file_C.vcf.gz samples: `{S4,S5}`

- file_D.vcf.gz samples: `{S3,S4,S5}`

From this example, we can delineate four difference sample set overlap cases, three of which are supported by DivBase/`bcftools`:

1. **Identical sets** (same sample IDs, same sample column order)
   Example: A=`{S1,S2,S3}`, B=`{S1,S2,S3}`
   This is `bcftools concat`-compatible. The main use case for this is when the VCF data is split by chromosome/scaffold into several VCF files.

2. **Non-overlapping sets** (no shared sample IDs)
   Example: A=`{S1,S2,S3}`, C=`{S4,S5}`
   This is `bcftools merge`-compatible.

3. **Partial overlap between sets** (some shared sample IDs, but not all)
   Example: A=`{S1,S2,S3}`, D=`{S3,S4,S5}`
   This is not supported by DivBase queries since there is a non-identical overlap of the sample sets. Please do not organize your VCF data like this.

4. **Mixed group of overlapping and non-overlaping sets**
   Example: A=`{S1,S2,S3}`, B=`{S1,S2,S3}`, C=`{S4,S5}`
   This is supported by DivBase. First `bcftools concat` is applied to the identical-set groups (`A+B`), and then `bcftools merge` for the two resulting groups (`(A+B)+C`) since they have non-overlapping sample sets (`{S1,S2,S3}` and `{S4,S5}`).

Links to relevant `bcftools` documentation:

- `bcftools concat` (same sample columns in the same order; input order matters): <https://samtools.github.io/bcftools/bcftools#bcftools-concat-options-file1-file2>
- `bcftools merge` (designed for non-overlapping sample sets; duplicate sample names require force mode): <https://samtools.github.io/bcftools/bcftools#bcftools-merge-options-avcfgz-bvcfgz>

### 2.3. Additional examples of compatible and incompatible VCF data overlaps

In addition to the sample set overlaps described above, the combination of variant and sample columns also need to be unique across the files. Below are a few examples of compatible and incompatible variants.

#### 2.3.1 Compatible: same samples, different variants

This matches chromosome/region split workflows and is typically concat-compatible.

```text
# file_1.vcf.gz
#CHROM  POS     ID      REF ALT ... FORMAT SAMPLE001 SAMPLE002 SAMPLE003 SAMPLE004
1       12345   1_12345 T   C   ... GT     0/0       1/0       0/0       0/0

# file_2.vcf.gz
#CHROM  POS     ID      REF ALT ... FORMAT SAMPLE001 SAMPLE002 SAMPLE003 SAMPLE004
1       56789   1_56789 T   C   ... GT     0/0       0/0       0/0       0/0
```

#### 2.3.2. Compatible: same variant, different samples

This is typically merge-compatible.

```text
# file_1.vcf.gz
#CHROM  POS     ID      REF ALT ... FORMAT SAMPLE001 SAMPLE002 SAMPLE003 SAMPLE004
1       12345   1_12345 T   C   ... GT     0/0       1/0       0/0       0/0

# file_2.vcf.gz
#CHROM  POS     ID      REF ALT ... FORMAT SAMPLE005 SAMPLE006 SAMPLE007 SAMPLE008
1       12345   1_12345 T   C   ... GT     0/0       1/0       0/0       0/0
```

#### 2.3.3. Incompatible: duplicated sample+variant across files

This is a data duplication issue, and is not supported in DivBase:

```text
# file_1.vcf.gz
#CHROM  POS     ID      REF ALT ... FORMAT SAMPLE001 SAMPLE002 SAMPLE003 SAMPLE004
1       12345   1_12345 T   C   ... GT     0/0       1/0       0/0       0/0

# file_2.vcf.gz
#CHROM  POS     ID      REF ALT ... FORMAT SAMPLE001 SAMPLE005 SAMPLE006 SAMPLE007
1       12345   1_12345 T   C   ... GT     0/0       1/0       0/0       0/0
```

#### 2.3.4 Sample order matters for split-VCF files

If the VCF data has been split per chromosome/scaffolds in multiple VCF files (which is generally recommended for DivBase), please ensure that the the sample column order is identical order across the files.

Supported (identical sample sets):

- `chr1.vcf.gz`: `S1,S2,S3,S4`
- `chr2.vcf.gz`: `S1,S2,S3,S4`

Not supported (non-identical sample sets):

- `chr1.vcf.gz`: `S1,S2,S3,S4`
- `chr2.vcf.gz`: `S2,S1,S3,S4`

This follows `bcftools concat` requirements:
<https://samtools.github.io/bcftools/bcftools#concat>

## 3. After uploading the files: update dimensions before querying

After adding a new VCF file or uploading a new version of an existing file, please ensure that the DivBase project's [VCF dimensions cache](vcf-dimensions.md) is up-to-date by running:

```bash
divbase-cli dimensions update --project <PROJECT_NAME>
```

## 4. DivBase versioning of files

VCF files (and sidecar metadata files) are versioned in the DivBase project's data storage. The DivBase query workflows use the latest file versions of each VCF file, as indexed in the VCF dimensions cache. This is why it is important to ensure that the dimensions cache is up-to-date.

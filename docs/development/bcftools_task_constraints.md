# Bcftools Celery Task Constraints

DivBase uses `bcftools` to perform queries that subset VCF files based on user input. Thus, DivBase needs to adhere to the syntax and requirements for the different `bcftools` commands it uses as described in the [`bcftools` manual](https://samtools.github.io/bcftools/bcftools.html).

What sets DivBase apart from regular `bcftools` scripting is that the DivBase server dynamically orchestrates `bcftools` workflows based on user queries, and will take into account multiple VCF files from the DivBase project's data storage if needed ([as described in e.g. the VCF query user guide](../user-guides/vcf-query-syntax.md/#53-how-does-divbase-process-the-vcf-files)). There are two ways to combine VCF files with `bcftools`: `bcftools merge` and `bcftools concat`. The choice of command depends on whether or not the sample columns overlap or not in the VCF files that are to be combined and are described in [Section 2](#2-bcftools-merge) and [Section 3](#3-bcftools-concat), respectively.

Several of the findings from this document has been implemented in the `_check_if_samples_can_be_combined_with_bcftools` helper function in `tasks.py`. It uses the VCF dimensions index for the DivBase project to analyse whether or not the input VCFs of a given query fullfil the rules for processing by the bcftools logic in `tasks.bcftools_query`. If not, the task exits to save on compute resources. The exit happens before the transfer of VCF files to the worker containers from S3 is started.

These rules are considered a system contract for the design of the DivBase `bcftools` orchestration logic. Be very careful when refactoring code that touches upon these rules. Several regression tests (prefixed with `test_regression_*`) have been implemented to guard the parts of the codebase that rely on these rules.

## 1. General rules

### 1.1. VCF files need to be sorted by position

Sorted the variants is a major assumption for working VCF files, so it is not very likely that users will come to DivBase with unsorted files. Nevertheless, the system need to ensure that files are sorted. Sorting can be done with `bcftools sort`, but it should be up to the users to ensure that their files are sorted before upload to DivBase.

The variants need to be sorted at least per scaffold: it is possible to have the scaffolds in different order if the VCF files are indexed, as discussed below.

If attempting to subset a single VCF file that has unsorted positions, the `bcftools` will ask for an index:

```bash
bcftools view -r 21 tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_scramble_coordinates.vcf.gz

# [E::idx_find_and_load] Could not retrieve index file for 'tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_scramble_coordinates.vcf.gz'
```

but indexing is actually not possible in this case, since it requires sorted positions:

```bash
bcftools index tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_scramble_coordinates.vcf.gz

# [E::hts_idx_push] Unsorted positions on sequence #1: 22053057 followed by 17504018
# index: failed to create index for "tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_scramble_coordinates.vcf.gz"
```

**Implementation:**

- Check that all VCF files are sorted before running `bcftools` on the files. Currently implemented in the module-level `ensure_csi_index` function in `services/queries.py`, shared by both `BcftoolsQueryManager` and `VCFDimensionCalculator`. This is a two birds in one stone-case: files should preferrably be indexed, and unsorted files will error upon indexing (as illustrated above). If the files are not sorted, a message is returned to the user that tells them to sort the files with `bcftools sort` and upload them to their DivBase project and try the query again.

**Ideas for future improvements:**

- Implement a VCF specification checking job to ensure that user-uploaded VCFs adhere to the [VCF specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf).

- Sorting of VCF files could possibly be done by jobs triggered when a VCF file is uploaded to a bucket? Could be tied to VCF specification checks: if file pass all other checks but is unsorted, run it through `bcftools sort` and upload that as a new version?

### 1.2. VCF files need to be indexed

Most - but not all -  `bcftools` commands requires that the VCF files are indexed either with the TBI or CSI index format. `bcftools view` requires an index and since that operation is the core of the DivBase bcftools orchestration, index files are thus required in DivBase.

`bcftools merge` will for instance raise this error if one of the files that are to be combined do not have an index in the same folder as the VCF file; there is a flag `--force-no-index`, but as mentioned earlier, indexes are needed for other bcftools commands.

`[E::idx_find_and_load] Could not retrieve index file for '<VCF_FILENAME>'`

For the commands in this document to work, the following files needs to be indexed:

```bash
bcftools index tests/fixtures/HOM_20ind_17SNPs_last_10_samples.vcf.gz
bcftools index tests/fixtures/HOM_20ind_17SNPs_first_10_samples.vcf.gz
bcftools index tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_break_merge.vcf.gz
bcftools index tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_scramble_coordinates.vcf.gz
bcftools index tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_have_different_scaffolds.vcf.gz
bcftools index tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_have_different_variants.vcf.gz
```

**Implementation:**

- The module-level `ensure_csi_index` function in `services/queries.py` is called for all VCF files operated on in a given query (including temp files) and also during VCF dimension calculation. `BcftoolsQueryManager.ensure_csi_index` delegates to it, passing `self.run_bcftools` so Docker-exec routing is preserved. The index files are not uploaded to the DivBase project (at the time of writing), and thus they are recalculated each time a VCF file is subject to a query.

**Ideas for future improvements:**

- An alternative solution would be to store index files in the bucket along with the VCF files. The challenge there is to ensure that the index files are kept up to date with the VCF files. There is risk for drift if the VCF files are updated but not their associated CSI index files.

### 1.3. There cannot be duplicate sample names cannot recur in a single VCF file

This should not come as a surprise, but if a VCF file is modified so that the exact same sample name ocurrs in two sample columns in the same VCF file, it cannot even be bgzipped:

```bash
bcftools view -Oz -o tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_duplicate_a_sample_name.vcf.gz tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_duplicate_a_sample_name.vcf

# [E::bcf_hdr_add_sample_len] Duplicated sample name '8_HOM-E57'
# Failed to read from tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_duplicate_a_sample_name.vcf: could not parse header
```

**Solution:**

- It unlikely that users will upload with VCF files to DivBase with this error. But it would be easy to implement a function that checks that this does not occur.

### 1.4. What is required in the header?

- All scaffolds that are contained in the VCF have to be defined in the header

`[W::vcf_parse] Contig '51' is not defined in the header. (Quick workaround: index the file with tabix.)`

adding it to the header as ##contig=<ID=51,length=22053058> will resolve this.

- TODO: investigate the header requirements further

- TODO: what happens if there are coordinate that are outside the is within the stated scaffold length in the header?

## 2. bcftools merge

As stated in the [bcftools manual](https://www.htslib.org/doc/1.1/bcftools.html#merge), `bcftools merge` is used to:

> Merge multiple VCF/BCF files from non-overlapping sample sets to create one multi-sample file. For example, when merging file A.vcf.gz containing samples S1, S2 and S3 and file B.vcf.gz containing samples S3 and S4, the output file will contain four samples named S1, S2, S3, 2:S3 and S4.

Note that it is responsibility of the user to ensure that the sample names are unique across all files. If they are not, the program will exit with an error unless the option --force-samples is given. Note that sample names can be also given explicitly using the --print-header and --use-header options.

### 2.1. Only records from different VCF files can merged

there is a command called `bcftools norm` that can do melt operations within a file, but that is not needed for DivBase.

### 2.2. Sample names must be unique across all input files

`--force-samples` can be used to override this, but it gives undesireable results for DivBase

This is allowed since the VCF files have different samples names:

```bash
bcftools merge -Ov -o test1.vcf tests/fixtures/HOM_20ind_17SNPs_first_10_samples.vcf.gz tests/fixtures/HOM_20ind_17SNPs_last_10_samples.vcf.gz
```

This is not allowed since the second VCF file contains one sample name that is also in the first VCF:

```bash
bcftools merge -Ov -o test2.vcf tests/fixtures/HOM_20ind_17SNPs_first_10_samples.vcf.gz tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_break_merge.vcf.gz
```

It gives the error: `Error: Duplicate sample names (5a_HOM-I7), use --force-samples to proceed anyway.`

**Solution:**

- The VCF dimensions file could store all sample names for all VCF files in the bucket, and then when a query is submitted, a function should check if there are overlapping samples in a manner that cannot be handled by `bcftools concat` either.

- A VCF with sample set "S1,S2,S3" cannot be combined with `bcftools merge` or `bcftools concat` with a VCF that contains "S3,S4,S5". It cannot use `merge` since there are overlapping samples. It cannot use `concat` since there is only a partial overlap and not a complete overlap.

### 2.3. What about the order of the scaffolds?

from `bcftools merge -help`:

```text
--no-index  Merge unindexed files, the same chromosomal order is required and -r/-R are not allowed
```

This is not allowed:

```bash
bcftools merge --no-index  -Ov -o test3.vcf tests/fixtures/HOM_20ind_17SNPs_first_10_samples.vcf.gz tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_scramble_scaffold_order.vcf.gz
```

and gives this warning from using `--no-index`:

```text
[W::bcf_sr_add_reader] Using multiple unindexed files may produce errors, make sure chromosomes are in the same order!
```

and this error from the incorrect scaffold order:

```text
[E::_reader_fill_buffer] Sequences out of order, cannot stream multiple unindexed files: tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_scramble_scaffold_order.vcf.gz
```

If we ensure that the files are indexed and run the command without the `--no-index` flag, it works, however:

````bash
bcftools merge -Ov -o test3.vcf tests/fixtures/HOM_20ind_17SNPs_first_10_samples.vcf.gz tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_scramble_scaffold_order.vcf.gz
````

**Solution:**

- As long as the VCF files are indexed and that the `--no-index`  is not used, the order of the scaffolds is not a problem.

### 2.4. Can there be different variants in the VCF files?

Yes! This was probably to be expected, but this example shows it:

```bash
bcftools merge -Ov -o test5.vcf tests/fixtures/HOM_20ind_17SNPs_first_10_samples.vcf.gz tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_have_different_variants.vcf.gz
```

**Solution:**

- Not much to solve here, so this is more of a conclusion: as long as the sample names do not overlap, the VCF files can have different variant positions and still combine fine.

### 2.5. Can there be different scaffolds in the VCF files?

Yes! `bcftools merge` works even if there are no overlapping scaffolds:

```bash
bcftools merge -Ov -o test4.vcf tests/fixtures/HOM_20ind_17SNPs_first_10_samples.vcf.gz tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_have_different_scaffolds.vcf.gz
```

**Solution:**

- Again, more of a conclusion of a non-issue: as long as the sample names do not overlap, the VCF files do not need to share any scaffolds between them.

### 2.6. What about the headers? Do they need to be the same?

TODO: investigate this further

### 2.7. What if one file contains different INFO column values than the other?

from `bcftools merge -help`:

```text
-i, --info-rules TAG:METHOD,..    Rules for merging INFO fields (method is one of sum,avg,min,max,join) or "-" to turn off the default [DP:sum,DP4:sum]
```

TODO: investigate this further

## 3. bcftools concat

As stated in the [bcftools manual](https://www.htslib.org/doc/1.1/bcftools.html#concat), `bcftools concat` is used to:

> Concatenate or combine VCF/BCF files. All source files must have the same sample columns appearing in the same order. Can be used, for example, to concatenate chromosome VCFs into one VCF, or combine a SNP VCF and an indel VCF into one. The input files must be sorted by chr and position. The files must be given in the correct order to produce sorted VCF on output unless the -a, --allow-overlaps option is specified.

### 3.1. Sample names must be identical across all files

This is allowed:

```bash
bcftools concat -Ov -o test6.vcf tests/fixtures/HOM_20ind_17SNPs.1.vcf.gz tests/fixtures/HOM_20ind_17SNPs.4.vcf.gz
```

But this is not allowed:

```bash
bcftools concat -Ov -o test6.vcf tests/fixtures/HOM_20ind_17SNPs.1.vcf.gz tests/fixtures/HOM_20ind_17SNPs_last_10_samples.vcf.gz

# Checking the headers and starting positions of 2 files
# Different number of samples in tests/fixtures/HOM_20ind_17SNPs_last_10_samples.vcf.gz. Perhaps "bcftools merge" is what you are looking for?
```

**Solution:**

- There is logic in the `BcftoolsQueryManager` class checks the sample sets across the temp files generated during the pipeline and determined if some files should be processed with `concat` or `merge`. It can also handle the more complex  mixed concat-merge case where there are two groups for VCF files: files that contain completely overlapping samples (concat); and files with samples that do not overlap whatsoever with any other file in the query (merge). Thus, concat can first be performed and the resulting temp file is merged together with the other files.

- However, anything that deviates from the requirements of just `merge`, just `concat`, or the specific mixed concat-merge case described in the previous bullet, should raise an error to the user to ask them to change their query or even how they have arranged their samples in their bucket.

### 3.2. Sample names must be in the same order

This is not allowed:

```bash
bcftools concat -Ov -o test6.vcf tests/fixtures/HOM_20ind_17SNPs.1.vcf.gz tests/fixtures/HOM_20ind_17SNPs.4_with_edits_to_change_sample_order.vcf

# Checking the headers and starting positions of 2 files
# Different sample names in tests/fixtures/HOM_20ind_17SNPs.4_with_edits_to_change_sample_order.vcf. Perhaps "bcftools merge" is what you are looking for?
```

**Solution:**

- As proposed earlier in the document, the VCF dimensions file could store all sample names for all VCF files in the bucket, but it should also store them in immutable order so that it is possible to check not only that they sample names overlap, but that that the sample name order is preserved.

### 3.3. Does the bcftools concat input files be sorted by chr and position?

For `bcftools concat`, yes.

> Concatenate or combine VCF/BCF files. All source files must have the same sample columns appearing in the same order. Can be used, for example, to concatenate chromosome VCFs into one VCF, or combine a SNP VCF and an indel VCF into one. The input files must be sorted by chr and position. [bcftools concat manual](https://samtools.github.io/bcftools/bcftools.html#concat)

### 3.4. Does the files stritcly be given in the correct order to produce sorted VCFs with bcftools concat?

Yes, if the intent is to get sorted output from `bcftools concat`, the files need to be given in genomic/chromosomal order.

TODO: investigate if incorrect order, e.g. if the next file “goes backward” relative to the previous one can actually break the concat operation!

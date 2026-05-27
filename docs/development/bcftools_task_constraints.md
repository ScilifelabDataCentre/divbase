# Bcftools Celery Task Constraints

DivBase uses `bcftools` for several VCF parsing steps, including VCF dimensions caching and VCF queries. This document outlines the rules and constraints of different `bcftools` commands used in the DivBase backend. The main reference for this is the [`bcftools` manual](https://samtools.github.io/bcftools/bcftools.html). The commands were run using `bcftools` v1.21.

!!! note
    There is intentional overlap between this page and [Key design decisions related to VCF files](key_design_decisions_related_to_vcf_files.md). This page goes deep into the details of `bcftools`, whereas the other document describe high level design choices.

What sets DivBase apart from regular `bcftools` scripting is that the DivBase server dynamically orchestrates `bcftools` workflows based on user queries, and will take into account multiple VCF files from the DivBase project's data storage if needed ([as described in e.g. the VCF query user guide](../user-guides/vcf-query-syntax.md/#53-how-does-divbase-process-the-vcf-files)). There are two ways to combine VCF files with `bcftools`: `bcftools merge` and `bcftools concat`. The choice of command depends on whether or not the sample columns overlap or not in the VCF files that are to be combined. An overview of the possible sample column overlap cases and are found in [Section 2](#2-how-divbase-chooses-between-bcftools-merge-and-concat) and described in detail in [Section 3 (`bcftools merge`)](#3-bcftools-merge) and [Section 4 (`bcftools concat`)](#4-bcftools-concat).

!!! Important
    The rules described in this document are considered a system contract for the design of the DivBase `bcftools` orchestration logic. Be very careful when refactoring code that touches upon these rules. Several regression tests (prefixed with `test_regression_*`) have been implemented to guard the parts of the codebase that rely on these rules.

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

- Check that all VCF files are sorted before running `bcftools` on the files. Currently implemented in the module-level `ensure_csi_index` function in `services/bcftools_helpers.py`, shared by both `BcftoolsQueryManager` and `VCFDimensionCalculator`. This is a two birds in one stone-case: files should preferrably be indexed, and unsorted files will error upon indexing (as illustrated above). If the files are not sorted, a message is returned to the user that tells them to sort the files with `bcftools sort` and upload them to their DivBase project and try the query again.

**Regression test coverage:**

- `test_regression_ensure_csi_index_raises_task_user_error_on_unsorted_positions` in `tests/unit/divbase_api/test_vcf_query_task.py`.

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

- The module-level `ensure_csi_index` function in `services/bcftools_helpers.py` is called for all VCF files operated on in a given query (including temp files) and also during VCF dimension calculation. Docker-exec routing is preserved because `ensure_csi_index` calls `run_bcftools`, which handles the Docker/k8s dispatch. The index files are not uploaded to the DivBase project (at the time of writing), and thus they are recalculated each time a VCF file is subject to a query.

**Regression test coverage:**

- `test_regression_ensure_csi_index_raises_task_user_error_on_unsorted_positions` in `tests/unit/divbase_api/test_vcf_query_task.py`.

**Ideas for future improvements:**

- An alternative solution would be to store index files in the bucket along with the VCF files. The challenge there is to ensure that the index files are kept up to date with the VCF files. There is risk for drift if the VCF files are updated but not their associated CSI index files.

### 1.3. Duplicate sample names cannot recur in a single VCF file

If a VCF file is modified so that the same sample name occurs in two sample columns in the same file, it cannot be processed by `bcftools`:

```bash
bcftools view -Oz -o tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_duplicate_a_sample_name.vcf.gz tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_duplicate_a_sample_name.vcf

# [E::bcf_hdr_add_sample_len] Duplicated sample name '8_HOM-E57'
# Failed to read from tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_duplicate_a_sample_name.vcf: could not parse header
```

**Implementation:**

- This is guarded against during the VCF dimensions caching. `services/bcftools_helpers.py` contains shared bcftools stderr classification logic in `_raise_task_user_error_from_bcftools_stderr`. `services/vcf_dimension_indexing.py` calls that shared helper when header parsing and bgzip/index-related steps fail. If stderr contains `Duplicated sample name`, DivBase raises a user-facing `TaskUserError` with guidance to fix duplicate sample IDs and re-run dimensions update.

**Regression test coverage:**

- `test_regression_update_dimensions_fails_for_vcf_with_duplicate_sample_ids_in_header` in `tests/e2e_integration/cli_commands/test_dimensions_cli.py`.

### 1.4. What is required in the header?

There is no single "one-line rule" for VCF headers across all `bcftools` commands, but for DivBase dimensions caching and query orchestration the following header-related errors need to be handled:

- Invalid/unparseable header content (for example `unknown file type` or `could not parse header`)
- Duplicate sample IDs in header (`Duplicated sample name ...`)

**Implementation:**

- `services/vcf_dimension_indexing.py` performs header parsing (`bcftools view --header-only`).
- Any stderr from this step is routed through shared classifier logic in `services/bcftools_helpers.py::_raise_task_user_error_from_bcftools_stderr`.
- This gives explicit `TaskUserError` messages for header/content failures instead of opaque bcftools errors.

**Regression test coverage:**

- `test_regression_update_dimensions_fails_for_vcf_with_invalid_header_content` in `tests/e2e_integration/cli_commands/test_dimensions_cli.py`.
- `test_regression_vcf_query_fails_with_descriptive_error_for_malformed_sample_columns` in `tests/e2e_integration/cli_commands/test_query_cli.py`.
- `test_regression_update_dimensions_fails_for_non_vcf_file_disguised_as_vcf_gz` in `tests/e2e_integration/cli_commands/test_dimensions_cli.py`.

**Ideas for future improvements:**

- Missing contig definitions are often warning-level diagnostics in bcftools, for example:
  `[W::vcf_parse] Contig '51' is not defined in the header. (Quick workaround: index the file with tabix.)`
  Add explicit DivBase policy checks if we want this to be a strict upload/dimensions error.
- Add explicit checks for coordinates outside declared contig lengths and decide whether to treat this as a strict DivBase validation failure.
- Add broader VCF specification validation as a preflight/upload-time job.
- Add a validation step (e.g. at file upload or at dimensions update) to check for header-data mismatches, e.g. sample-column count mismatches (`Number of columns ... does not match the number of samples`). A full test would require passing over every variant line in the VCF file, which is an O(n) process. Another idea would be to validate only the first e.g. 1000 lines and assume that potential systematic formatting errors caused by the upstream VCF generation workflow would be caught in that subsubset of variants.

## 2. How DivBase chooses between bcftools merge and concat

`bcfools merge` and `concat` are designed to consider so-called "sample sets", i.e. ["sample columns appearing in the same order"](https://samtools.github.io/bcftools/bcftools.html#concat) in the `#CHROM` header of the VCF file. For example, if `S1`- `S6` are names of different samples, a sample set is an ordered combination of samples in a given VCF file, such as `S1,S3,S4`.

The relationship between sample sets across the input files to a DivBase query determines which `bcftools` command (`merge` or `concat`) — or whether the query is rejected before bcftools is called at all.

DivBase handles five cases of sample set overlap:

- **Completely non-overlapping sample sets**

    E.g. VCF file 1: `S1,S2,S3` vs VCF file 2: `S4,S5,S6`. Handled with `bcftools merge` (see [Section 3](#3-bcftools-merge)).

- **Identical sample sets (same sample column order)**

    E.g. VCF file 1: `S1,S2,S3` vs VCF file 2: `S1,S2,S3`. Handled with `bcftools concat` (see [Section 4](#4-bcftools-concat)).

- **Mixed identical-sample set and non-overlapping sample set**

    E.g. VCF file 1: `S1,S2,S3`, VCF file 2: `S1,S2,S3`, and VCF file 3: `S4,S5,S6`. The identical-sample set files (`S1,S2,S3`) are concatenated first, then the output file from the concatenation (`S1,S2,S3`) is merged with the non-overlapping file `S4,S5,S6`. Handled with `bcftools concat` followed by `bcftools merge`.

- **Same samples but different sample column order**

    E.g. VCF file 1: `S1,S2,S3` vs VCF file 2: `S2,S1,S3`. Violates `bctools` contstrains, rejected with a `TaskUserError` before the VCF processing of the VCF query starts.

- **Partially overlapping sample sets**

    E.g. VCF file 1: `S1,S2,S3` vs VCF file 2: `S3,S4,S5` (`S3` is in both sample sets). Violates `bctools` contstrains, rejected with a `TaskUserError` before the VCF processing of the VCF query starts.

This classification is performed by `_check_if_samples_can_be_combined_with_bcftools` in `worker/tasks.py`, using sample names stored in the dimensions cache to determine whether the input VCFs fulfil the rules for the bcftools orchestration in `tasks.bcftools_query`. If they do not, the task exits before VCF files are downloaded from S3, saving compute resources.

## 3. bcftools merge

This section covers a few specific requirements for `bcftools merge` and add details to the sample set overlap cases described in [Section 2](#2-how-divbase-chooses-between-bcftools-merge-and-concat).

As stated in the [bcftools manual](https://www.htslib.org/doc/1.1/bcftools.html#merge), `bcftools merge` is used to:

> Merge multiple VCF/BCF files from non-overlapping sample sets to create one multi-sample file. For example, when merging file A.vcf.gz containing samples S1, S2 and S3 and file B.vcf.gz containing samples S3 and S4, the output file will contain four samples named S1, S2, S3, 2:S3 and S4.

!!! Note
    The bash examples throughout Sections 3 and 4 use `-Ov -o *.vcf` (uncompressed VCF text) so that output is human-readable when running the commands manually. `BcftoolsQueryManager` in `services/vcf_queries.py` uses `-Ou -o *.bcf` (uncompressed BCF binary) for all intermediate temp files, which avoids compression and decompression overhead between pipeline steps. Only the final `bcftools sort` step writes `-Oz -o *.vcf.gz` (bgzipped VCF) as the query result file.

### 3.1. Only records from different VCF files can be merged

`bcftools merge` combines records across different input files — it does not transform or reshape records within a single file. If the query only involves a single source VCF file, `merge_or_concat_bcftools_temp_files` will receive a list of one temp file and skip the merge/concat step entirely, renaming that file to the output path instead (`len(output_temp_files) == 1` branch in `services/vcf_queries.py`).

For intra-file normalization such as splitting multiallelic variants into biallelic records, `bcftools norm` is the relevant command, but that is not part of the current scope of DivBase. Users should perform such file operations prior to uploading files to a DivBase project.

**Implementation:**

No real need to guard or enforce this, much thanks to the single source VCF file handling described above.

**Ideas for future improvements:**

TODO: exact duplicate variants in to VCF files with different filename are currently not guarded against. That is a potentially big refactoring of the VCF dimensions caching. (marked as TODO since it is important)

### 3.2. Sample names must be unique across all input files

This covers the partially overlapping case in [Section 2](#2-how-divbase-chooses-between-bcftools-merge-and-concat). `bcftools merge` has a `--force-samples` flag that overrides the duplicate-sample-name error, but DivBase deliberately does not use it: when triggered, it silently renames conflicting samples by prepending the file index (e.g. `S3` becomes `2:S3`), which would corrupt the sample columns in the query results without any visible error. DivBase instead rejects incompatible sample sets before the merge step runs.

The following fixtures can be merged the VCF files have different samples names:

```bash
bcftools merge -Ov -o test1.vcf tests/fixtures/HOM_20ind_17SNPs_first_10_samples.vcf.gz tests/fixtures/HOM_20ind_17SNPs_last_10_samples.vcf.gz
```

The following fixtures cannot be merged since the second VCF file contains one sample name that is also in the first VCF:

```bash
bcftools merge -Ov -o test2.vcf tests/fixtures/HOM_20ind_17SNPs_first_10_samples.vcf.gz tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_break_merge.vcf.gz
```

It gives the error: `Error: Duplicate sample names (5a_HOM-I7), use --force-samples to proceed anyway.`

**Implementation:**

- The VCF dimensions cache stores all sample names for all VCF files in the bucket. When a query is submitted, `_check_if_samples_can_be_combined_with_bcftools` in `worker/tasks.py` reads these from the dimensions cache and checks whether the sample sets across the requested files can be combined. Partial overlap is one of the blocked cases (see [Section 2](#2-how-divbase-chooses-between-bcftools-merge-and-concat)).

**Regression test coverage:**

- `test_regression_query_fails_for_vcf_with_partly_overlapping_sample_sets` in `tests/e2e_integration/cli_commands/test_query_cli.py`.

### 3.3. What about the order of the scaffolds?

As long as the VCF files have a CSI index, `bcftools merge` will use the it to parse each file by region and does not therefore not require variants to to appear in any particular order.

It is possible to attempt to merge VCF files that do not have CSI indexes with `bcftools merge --no-index`. This is not used in DivBase (since indexing and sorting of files is a prerequisite, as described in [Section 1](#1-general-rules)).

!!! Example
    This is a small example that illustrate what happens when trying to use `--no-index` and unindexed VCFs:

    Docstring from `bcftools merge -help`:

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

    ```bash
    bcftools merge -Ov -o test3.vcf tests/fixtures/HOM_20ind_17SNPs_first_10_samples.vcf.gz tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_scramble_scaffold_order.vcf.gz
    ```

**Implementation:**

- DivBase never passes `--no-index` to any `bcftools merge` call. All files are indexed with a CSI index before the merge step runs, via the same `ensure_csi_index` function in `services/bcftools_helpers.py` described in [Section 1.2](#12-vcf-files-need-to-be-indexed). The scaffold record order of the input files is therefore irrelevant to the merge.

**Regression test coverage:**

- No dedicated test for this case. The indexing guarantee that makes it a non-issue is covered by `test_regression_ensure_csi_index_raises_task_user_error_on_unsorted_positions` in `tests/unit/divbase_api/test_vcf_query_task.py` (see [Section 1.2](#12-vcf-files-need-to-be-indexed)).

### 3.4. How does merge handle overlapping and non-overlapping variant positions and scaffolds?

`bcftools merge` operates at the variant level: for each unique CHROM/POS across all input files it produces one output record, combining genotype columns from all samples. This is the key behavioral difference from `bcftools concat`, which simply stacks records from each file without any position-level merging.

The following three examples use VCF files with non-overlapping sample sets (see [Section 2](#2-how-divbase-chooses-between-bcftools-merge-and-concat)) and show how `bcftools merge` handles varying degrees of variant and scaffold overlap.

**Example 1 — Same variant coordinates in both files (the common case):**

`HOM_20ind_17SNPs_first_10_samples.vcf.gz` and `HOM_20ind_17SNPs_last_10_samples.vcf.gz` share all 17 variant positions but have non-overlapping sample sets. Merging them produces 17 output records (not 34), with genotype data from all 20 samples present in each row:

```bash
bcftools merge -Ov -o test5a.vcf tests/fixtures/HOM_20ind_17SNPs_first_10_samples.vcf.gz tests/fixtures/HOM_20ind_17SNPs_last_10_samples.vcf.gz
```

**Example 2 — Non-overlapping variant coordinates, same scaffolds:**

`HOM_20ind_17SNPs_last_10_samples_with_edit_to_have_different_variants.vcf.gz` has all coordinates shifted by +1 relative to `HOM_20ind_17SNPs_first_10_samples.vcf.gz`, so no positions are shared. Merging them produces 34 output records (17 from each file). For each record, samples from the file that does not contain that position receive missing genotypes (`./.`).

`bcftools concat` cannot be used to avoid the missing genotypes here: concat requires identical sample sets in the same column order, which is the opposite precondition from merge. Because these two files have non-overlapping sample sets they must be merged, and the missing genotypes are the correct biological output — those samples have no data at those positions.

```bash
bcftools merge -Ov -o test5b.vcf tests/fixtures/HOM_20ind_17SNPs_first_10_samples.vcf.gz tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_have_different_variants.vcf.gz
```

**Example 3 — Completely different scaffold names:**

`HOM_20ind_17SNPs_last_10_samples_with_edit_to_have_different_scaffolds.vcf.gz` has all scaffold IDs prefixed with `5` (e.g. chromosome `1` → `51`, `4` → `54`, `21` → `521`), so there is no scaffold overlap. This is the same mechanism as Example 2: because no CHROM/POS in one file matches any CHROM/POS in the other, the merge produces 34 output records (17 per file) and non-represented samples receive missing genotypes (`./.`).

```bash
bcftools merge -Ov -o test5c.vcf tests/fixtures/HOM_20ind_17SNPs_first_10_samples.vcf.gz tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_have_different_scaffolds.vcf.gz
```

**Implementation:**

- More of a conclusion than a problem: no special handling is needed in DivBase. `merge_or_concat_bcftools_temp_files` in `services/vcf_queries.py` passes the temp files directly to `bcftools merge` without any variant-level flags. The position-merging and missing-genotype-filling is `bcftools` default behaviour: shared positions are merged into single records; positions unique to one file get missing genotypes for samples from the other file. This applies regardless of whether the position differences arise from shifted coordinates or entirely different scaffold names.

**Regression test coverage:**

- No dedicated test for this case. It is implicitly covered by the checksum-based query result tests in `tests/e2e_integration/cli_commands/test_query_cli.py`, which assert on stable output across multi-file merge queries.

### 3.5. What about the headers? Do they need to be the same?

No, the headers do not need to be identical. `bcftools merge` builds the output header by computing a union of all input file headers. Deduplication is performed bason on header lines (`FILTER`, `INFO`, `FORMAT`, and `contig`) whose `ID` does not yet appear in the output header; lines whose `ID` is already present are silently skipped (first-seen wins tie-breaks).

No error or warning is emitted for structural differences such as different sets of `contig` lines or different `FILTER` definitions. The only header-related warning is for `INFO` or `FORMAT` fields that appear in two files with a different `Number` or `Type` attribute: `bcftools` emits a warning to stderr but continues the merge, keeping the first file's definition. There seem to be no `bcftools`option to control this behaviour.

There is a `--use-header FILE` option to bypass the header union logic and load the output header from a user-supplied file. DivBase does not use this option.

!!! Example
    The `--print-header` flag prints the merged header without processing any records, which is useful for inspection. For example, merging files with completely different scaffold names (see Example 3 in [Section 3.4](#34-how-does-merge-handle-overlapping-and-non-overlapping-variant-positions-and-scaffolds)) produces a header that is the union of both files' `contig` lines:

    ```bash
    bcftools merge --print-header tests/fixtures/HOM_20ind_17SNPs_first_10_samples.vcf.gz tests/fixtures/HOM_20ind_17SNPs_last_10_samples_with_edit_to_have_different_scaffolds.vcf.gz
    ```

    The merged header contains `contig` entries from the first file (`1`, `4`, `5`, `6`, `7`, `8`, `13`, `18`, `20`, `21`, `22`, `24`) followed by entries from the second file (`51`, `54`, `55`, `56`, `57`, `58`, `513`, `518`, `520`, `521`, `522`, `524`). The `INFO`, `FORMAT`, and `FILTER` definitions are identical in these two fixtures (both derive from the same source VCF), so the deduplication is a no-op for those lines.

**Implementation:**

- No special handling is needed in DivBase. `merge_or_concat_bcftools_temp_files` in `services/vcf_queries.py` does not pass any header-related flags to `bcftools merge`.

**Regression test coverage:**

- No dedicated test for this case. The default union header behaviour is exercised implicitly by every multi-file merge query in `tests/e2e_integration/cli_commands/test_query_cli.py`.

### 3.6. What if one file contains different INFO column values than the other?

When merging files with non-overlapping sample sets, the `INFO` column values at a shared position typically differ between files: each file's `AC`, `AN`, `DP`, etc. reflect only the samples in that file. `bcftools merge` handles this differently depending on whether not the merged values can be recalculated from the values in the input VCFs:

**Standard fields derivable from genotype data (e.g. `AC`, `AN`):**

`bcftools` recalculates these from the merged genotype columns rather than combining the input values. For example, the `AC` and `AN` values for all positions in the merged output correctly reflect all 20 samples, i.e. it took the 10 samples in each input file into account:

```bash
bcftools query -f '%POS\t%INFO/AC\t%INFO/AN\n' tests/fixtures/HOM_20ind_17SNPs_first_10_samples.vcf.gz | head -3
# 17504018    7    20
# 22053057    10   20
# 13086614    1    20

bcftools query -f '%POS\t%INFO/AC\t%INFO/AN\n' tests/fixtures/HOM_20ind_17SNPs_last_10_samples.vcf.gz | head -3
# 17504018    0    20
# 22053057    10   18
# 13086614    0    20

bcftools merge -Ov -o test_info.vcf tests/fixtures/HOM_20ind_17SNPs_first_10_samples.vcf.gz tests/fixtures/HOM_20ind_17SNPs_last_10_samples.vcf.gz
bcftools query -f '%POS\t%INFO/AC\t%INFO/AN\n' test_info.vcf | head -3
# 17504018    7    40
# 22053057    20   38
# 13086614    1    40
```

**Non-derivable fields (e.g. `DP` for read depth, custom annotations):**

For fields that cannot be computed from genotype data, `bcftools` applies the `-i, --info-rules` mechanism:

```text
-i, --info-rules TAG:METHOD,..    Rules for merging INFO fields (method is one of sum,avg,min,max,join) or "-" to turn off the default [DP:sum,DP4:sum]
```

DP and DP4 are the only two fields with default sum rules. Every other field not covered by an explicit `-i` rule takes its value from the first input file at the shared positions (same CHROM/POS in both files). This has the following implications:

- **Both files have the field but with different values:** the merged record uses the first file's value; the second file's value is silently discarded.
- **Only the second (or later) file has the field:** the merged record takes the first file's "value" — which is absent — so the field is missing in the output even though the second file had data there.

Both cases apply only at positions shared by both files. At positions unique to one file (as in Examples 2 and 3 in [Section 3.4](#34-how-does-merge-handle-overlapping-and-non-overlapping-variant-positions-and-scaffolds), where variant coordinates do not overlap), there is no conflict and the single file's `INFO` values pass through unchanged.

**Implementation:**

- DivBase does not pass `-i` to `bcftools merge`. This means that the default `bcftools` behavior is used. Currently, users who rely on custom `INFO` fields being correctly combined across files need to post-process the query output themselves, which might not be optimal.

**Regression test coverage:**

- No dedicated test for this case. The standard `AC`/`AN` recalculation behaviour is implicitly included in every multi-file merge query in `tests/e2e_integration/cli_commands/test_query_cli.py`.

**Ideas for future improvements:**

Gather feedback from users regarding the `-i` rules to learn if this is a pain point. One solution to this would be to allow users to control the merge rules in `divbase-cli query vcf`

## 4. bcftools concat

This section covers a few specific requirements for `bcftools concat` and add details to the sample set overlap cases described in [Section 2](#2-how-divbase-chooses-between-bcftools-merge-and-concat).

As stated in the [bcftools manual](https://samtools.github.io/bcftools/bcftools.html#concat), `bcftools concat` is used to:

> Concatenate or combine VCF/BCF files. All source files must have the same sample columns appearing in the same order. Can be used, for example, to concatenate chromosome VCFs into one VCF, or combine a SNP VCF and an indel VCF into one. The input files must be sorted by chr and position. The files must be given in the correct order to produce sorted VCF on output unless the `-a, --allow-overlaps` option is specified.

!!! Note
    As with Section 3, the bash examples below use `-Ov -o *.vcf` for readability. DivBase uses `-Ou -o *.bcf` for the actual concat temp files (see the note on output formats at the start of [Section 3](#3-bcftools-merge)).

### 4.1. Sample sets must be identical across all VCF files (sample names and order)

As the name implies, the `concat` command concatenates rows from the VCF files that are to be combined, and thus requires that the sample sets are identical between the files.

A typical use case for `concat` is variant data for a sample set that has been separated by chromosome or scaffold. For example, these fixtures describe the same exact sample set, but the variants in each VCF describe chromosome 1 and 4, respectively:

```bash
bcftools concat -Ov -o test6.vcf tests/fixtures/HOM_20ind_17SNPs.1.vcf.gz tests/fixtures/HOM_20ind_17SNPs.4.vcf.gz
```

However, trying to combine two files with different sample sets is not allowed.
In this example: the two VCF files have partially overlapping sample sets. Note that the error message suggests to try `bcftools merge` instead:

```bash
bcftools concat -Ov -o test6.vcf tests/fixtures/HOM_20ind_17SNPs.1.vcf.gz tests/fixtures/HOM_20ind_17SNPs_last_10_samples.vcf.gz

# Checking the headers and starting positions of 2 files
# Different number of samples in tests/fixtures/HOM_20ind_17SNPs_last_10_samples.vcf.gz. Perhaps "bcftools merge" is what you are looking for?
```

Another example of non-identical sample sets: sample names must be in the same order. This is not allowed, since the sample order in the chromosome 4 VCF file has been scrambled relative to the order in the file for chromosome 1:

```bash
bcftools concat -Ov -o test6.vcf tests/fixtures/HOM_20ind_17SNPs.1.vcf.gz tests/fixtures/HOM_20ind_17SNPs.4_with_edits_to_change_sample_order.vcf

# Checking the headers and starting positions of 2 files
# Different sample names in tests/fixtures/HOM_20ind_17SNPs.4_with_edits_to_change_sample_order.vcf. Perhaps "bcftools merge" is what you are looking for?
```

**Implementation:**

- `_check_if_samples_can_be_combined_with_bcftools` in `worker/tasks.py` groups VCF files by their sample set and only passes files with completely identical sample sets (same names, same order) to `bcftools concat`. This is also used to handle the mixed concat-merge case (see case 5 in [Section 2](#2-how-divbase-chooses-between-bcftools-merge-and-concat)).

**Regression test coverage:**

- `test_regression_sample_order_preservation_allows_concat_compatible_sample_sets` in `tests/unit/divbase_api/test_worker_crud_dimensions.py`.
- `test_merge_or_concat_raises_on_step_failure` (concat case) in `tests/unit/divbase_api/test_vcf_query_task.py`.
- `test_vcf_query_result_file_by_headerless_checksum` (Case 2) in `tests/e2e_integration/cli_commands/test_query_cli.py`.
- `test_regression_sample_order_preservation_enforces_bcftools_concat_guard` in `tests/unit/divbase_api/test_worker_crud_dimensions.py`.

**Ideas for future improvements:**

For case 4 in [Section 2](#2-how-divbase-chooses-between-bcftools-merge-and-concat), i.e. sample sets with same sample names but different order, it would technically be possible to sort the files to have identical sample sets. The [manual states](https://samtools.github.io/bcftools/bcftools.html#common_options) that when using `--samples` or `--samples-file`, "The sample order is updated to reflect that given in the input file".

Something like the following could be used to ensure that the sample order in a VCF file is sorted based on the order of another VCF file (assuming that they contain the same sample names):

```bash
bcftools query -l input_file_1.vcf > order_of_samples_in_file1.txt
bcftools view --samples-file order_of_samples_in_file1.txt input_file_2.vcf -o input_file_2_sorted.vcf
```

### 4.2. What about the order that the input files are given to bcftools concat?

`bcftools concat` requires that a) each VCF file is sorted by chromosome and position (same requirement as for `bcftools merge`), and b) that the input files are supplied to the `concat` command such that the same chromosome does not appear non-contiguously across files.

The requirement in a) is handled by the `ensure_csi_index` function in `services/bcftools_helpers.py` (see [Section 1.1](#11-vcf-files-need-to-be-sorted-by-position) and [Section 1.2](#12-vcf-files-need-to-be-indexed)).

For requirement b): if a chromosome that has already been processed in an earlier file reappears in a later file after a different chromosome has appeared in between, `bcftools` exits with a hard error:

```text
The chromosome block <chr> is not contiguous, consider running with -a.
```

See [Section 4.3](#43-what-happens-if-input-files-have-overlapping-or-duplicate-variant-positions) for a full discussion of all overlapping and duplicate position cases, including what `--allow-overlaps` does and its caveats.

**Implementation:**

- DivBase does not pass `--allow-overlaps` to `bcftools concat`. Within-file sort order is guaranteed by `ensure_csi_index` (see [Section 1.2](#12-vcf-files-need-to-be-indexed)). The across-file order is determined by the sequence in which VCF files appear in the query processing; DivBase does currently not explicitly sort the concat input list by chromosomal start coordinate. See Ideas for future improvements below.

**Regression test coverage:**

- No dedicated test for the across-file ordering case. Within-file sort order is covered by `test_regression_ensure_csi_index_raises_task_user_error_on_unsorted_positions` in `tests/unit/divbase_api/test_vcf_query_task.py` (see [Section 1.2](#12-vcf-files-need-to-be-indexed)).

**Ideas for future improvements:**

- TODO: Consider sorting the `concat` input file list by the first chromosome/position declared in each file's header before passing them to `bcftools concat`. Alternatively, use `concat --allow-overlaps --remove-duplicates` followed by `bcftools sort`; see [Section 4.3](#43-what-happens-if-input-files-have-overlapping-or-duplicate-variant-positions) for caveats about `--allow-overlaps`. (marked as TODO since it is important)

### 4.3. What happens if input files have overlapping or duplicate variant positions?

Overlapping or duplicate variant positions across `concat` input files are handled differently depending on the specific overlap pattern. Manual testing with the DivBase VCF test fixtures identified three different overlap cases and how they are handled differently by the default `bcftools concat` command with any special options.

**Case 1 — Same chromosome reappears non-contiguously across files (hard error):**

If a chromosome ID reappears in a later input file after a *different* chromosome has already been processed, `bcftools` exits with a non-zero exit code:

```bash
# HOM_20ind_17SNPs.1.vcf.gz: chr1 (positions 17504018, 22053057)
# HOM_20ind_17SNPs.4.vcf.gz: chr4 (position 13086614)
bcftools concat -Ov -o test_noncontig.vcf \
  tests/fixtures/HOM_20ind_17SNPs.1.vcf.gz \
  tests/fixtures/HOM_20ind_17SNPs.4.vcf.gz \
  tests/fixtures/HOM_20ind_17SNPs.1.vcf.gz

# The chromosome block 1 is not contiguous, consider running with -a.
# exit code: 255
```

**Case 2 — Exact duplicate records (no error, silent duplicate rows in output):**

If two files share the same chromosome and contain records at the same positions, `bcftools concat` does NOT error. All records are written to the output in the given file order. Can for instance be demonstrated by concatenating `HOM_20ind_17SNPs.1.vcf.gz` with itself:

```bash
# HOM_20ind_17SNPs.1.vcf.gz: 2 records (chr1:17504018, chr1:22053057)
bcftools concat -Ov -o test_dups.vcf \
  tests/fixtures/HOM_20ind_17SNPs.1.vcf.gz \
  tests/fixtures/HOM_20ind_17SNPs.1.vcf.gz

# No error. Exit code: 0.
# Output: 4 records — chr1:17504018 and chr1:22053057 appear twice each.
```

**Case 3 — Positions go backward across a file boundary (no error, unsorted output):**

If the first position of file N+1 is lower than the last position of file N on the same chromosome, `bcftools concat` also does NOT error. The records are written in the given file order, producing unsorted output with no warning. Verified empirically with bcftools 1.21 by constructing a file with chr1:17504018 followed by a file with chr1:17000000 — exit code was 0 and the output was chr1:17504018 then chr1:17000000 (unsorted).

**Conclusion:** Only case 1 is caught by `bcftools`. Cases 2 and 3 produce silent bad output. The manual's statement that "the files must be given in the correct order to produce sorted VCF on output" is a user responsibility, not a bcftools enforcement.

**Handling overlapping positions with `--allow-overlaps` and `-d`:**

The `--allow-overlaps` (`-a`) flag for `concat` mentioned in the error message of case 1 can be used to override requirement b) of [Section 4.2](#42-what-about-the-order-that-the-input-files-are-given-to-bcftools-concat) by allowing that the ["First coordinate of the next file can precede last record of the current file."](https://samtools.github.io/bcftools/bcftools.html#concat). When using `--allow-overlaps`, deduplication can be performed with `-d/--rm-dups`.

The difference between the options `snps|indels|both|all|exact` of `-d` are not that well explained in the manual, but [other](https://github.com/samtools/bcftools/issues/917) [sources](https://open.bioqueue.org/home/knowledge/showKnowledge/sig/bcftools-concat) [online](https://www.biocomputix.com/post/combining-multiple-vcf-bcf-files-using-bcftools-concat) seem to converge on that `exact` is the argument to remove only identical duplicates (same CROM, POS, REF and ALT allele) and that this is the reason that `-D, --remove-duplicates` is a short form for `-d exact`. The arguments `snps`, `indels`, `both` seem to removes duplicates of the specified variant type by dropping the ALT alleles. `all` seem to drop all duplicates completely.

!!! Note
    `--allow-overlaps` requires all input files to be bgzip-compressed and indexed. This is not documented in the manual but is confirmed in [GitHub issue #340](https://github.com/samtools/bcftools/issues/340).

**Implementation:**

- DivBase does not check for overlapping or duplicate variant positions before calling `bcftools concat`, and does not pass `--allow-overlaps` or `-d`. Case 1 propagates back as a task failure since `bcftools` exits with a non-zero return code. Cases 2 and 3 are silent failure modes that produce malformed output (duplicate rows or unsorted positions) without raising a task failure. Currently, users are responsible for ensuring their chromosome-split or region-split VCF files have non-overlapping, non-duplicate coordinates, but this will not be sustainable as the complexity of the VCF data in a DivBase project grows...

**Regression test coverage:**

- No dedicated test for any of these cases.

**Ideas for future improvements:**

- TODO: Handle overlapping variants for the concat case. Perhaps `--allow-overlaps --remove-duplicates` could be used? (marked as TODO since it is important)

### 4.4. What about headers and INFO field values?

The topics covered in [Section 3.5](#35-what-about-the-headers-do-they-need-to-be-the-same) and [Section 3.6](#36-what-if-one-file-contains-different-info-column-values-than-the-other) for `bcftools merge` are more straightforward for `bcftools concat`:

**Headers:**

The `bcftools concat` manual is silent on header handling in standard mode. Based on source code analysis, `bcftools concat` builds the output header using the same union logic as `bcftools merge` (see [Section 3.5](#35-what-about-the-headers-do-they-need-to-be-the-same)): new `INFO`, `FORMAT`, `FILTER`, and `contig` lines from later files are added to the output header, and lines whose `ID` already exists are silently skipped (first-seen wins). In practice, files used with `concat` have identical sample sets and often originate from the same upstream pipeline, so the likelihood is high that the headers are identical and the union is a no-op.

!!! Note
    There is an option `bcftools concat --naive` for fast concatenation without recompression (raw compressed blocks are streamed directly). The manual states it requires all files to be of the same type (all VCF or all BCF) and have the same headers, and that "A header compatibility check is performed and the program throws an error if it is not safe to use the option." `--naive-force` bypasses that check; the manual warns "Dangerous, use with caution." Neither flag is used in DivBase.

**INFO field values:**

Unlike `bcftools merge`, `concat` does not join sample column genotypes — it stacks variant rows from different files without combining any values. Each row in the output comes from exactly one input file and its `INFO` values are copied unchanged. There is no equivalent of the `--info-rules` mechanism from `merge` to consider.

**Implementation:**

- No special handling is needed in DivBase for headers or INFO columns. `merge_or_concat_bcftools_temp_files` in `services/vcf_queries.py` passes the temp files directly to `concat -Ou -o <output> <inputs>` without any header or INFO flags. The union header behavior and INFO column pass-through are both `bcftools` defaults.

**Regression test coverage:**

- No dedicated test for this case. Header union and INFO pass-through are exercised implicitly by the concat path in `test_vcf_query_result_file_by_headerless_checksum` (Case 2) in `tests/e2e_integration/cli_commands/test_query_cli.py`.

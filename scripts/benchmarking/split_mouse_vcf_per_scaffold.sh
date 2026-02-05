#!/bin/bash

# Companion script to the other benchmarking scripts that use the public mouse VCF file (https://www.ebi.ac.uk/ena/browser/view/ERZ022025).
# This script splits the mouse VCF by scaffold to benchmark how DivBase performs when the data is in multiple VCF files instead of a single VCF.
# This assumes that all samples are present in each scaffold VCF file, but that no variants are duplicated across the scaffold-split VCF files.

VCF="mgp.v3.snps.rsIDdbSNPv137.vcf.gz"
VCF_URL="ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ022/ERZ022025/mgp.v3.snps.rsIDdbSNPv137.vcf.gz"

# 1. Check if VCF file exists in current directory, otherwise download
if [ ! -f "$VCF" ]; then
  echo "VCF file $VCF not found. Downloading..."
  wget "$VCF_URL" || { echo "Download failed!"; exit 1; }
else
  echo "VCF file $VCF found."
fi

# 2. Check if .csi index exists in current directory, otherwise create it
if [ ! -f "$VCF.csi" ]; then
  echo "Index file $VCF.csi not found. Creating CSI index..."
  bcftools index -c "$VCF"
else
  echo "Index file $VCF.csi found."
fi

# 3. Extract scaffold names from the index
echo "Extracting scaffold names..."
scaffolds=$(bcftools index -s "$VCF" | grep -v '^#' | awk '{print $1}')
scaffold_list=$(echo $scaffolds)
echo "Scaffolds found: $scaffold_list"

# 4. Loop over scaffolds and split VCF if output does not exist for that scaffold
echo "Splitting VCF by scaffold..."
for scaffold in $scaffolds; do
  outvcf="mgp.v3.snps.rsIDdbSNPv137.$scaffold.vcf.gz"
  if [ -f "$outvcf" ]; then
    echo "Output $outvcf already exists. Skipping."
  else
    echo "Processing scaffold: $scaffold"
    bcftools view -r "$scaffold" -Oz -o "$outvcf" "$VCF"
  fi
done


# 5. Create a file list for divbase-cli upload
filelist="split_scaffold_files.txt"
rm -f "$filelist"
for scaffold in $scaffolds; do
  outvcf="mgp.v3.snps.rsIDdbSNPv137.$scaffold.vcf.gz"
  echo "$outvcf" >> "$filelist"
done

echo "File list for upload written to $filelist"
echo "Done."

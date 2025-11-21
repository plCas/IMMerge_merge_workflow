# Dependencies
[cyvcf2 v0.31.1](https://github.com/brentp/cyvcf2) <br>
[IMMerge](https://github.com/belowlab/IMMerge) <br>


# Process
  1. Remain only genotypes to reduce file size and generate .info.gz file by IMMerge. [bcftools_remain_gt_Info.sh](scripts/bcftools_remain_gt_Info.sh) <br>
  2. Calculate chunks for each chromosome (default: 3000000 variants/chunk; can adjust the number of variants according to the memory usage). [chunk_vcf_by_counts.py](scripts/chunk_vcf_by_counts.py) <br>
  3. Split vcf file and .info.gz according to pre-calculated chunks (generate from step 2). [bcftools_chunk_infos.sh](scripts/bcftools_chunk_infos.sh) <br>
  4. Merge all batches for each chunk. [IMMerge_chunks.sh](scripts/IMMerge_chunks.sh) <br>
  5. Concat all chunks together for each chromosome. [bcftools_concat.sh](scripts/bcftools_concat.sh) <br>

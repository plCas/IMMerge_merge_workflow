#!/usr/bin/env python3
import argparse
import sys
from typing import Iterator, Optional, Tuple

try:
    from cyvcf2 import VCF
except ImportError:
    sys.stderr.write("cyvcf2 is required. Install with: pip install cyvcf2\n")
    raise


def generate_regions_by_variant_count(vcf_path: str, chunk_size: int) -> Iterator[Tuple[str, int]]:
    """
    Iterate through a VCF and yield (region, count) pairs where each region is
    a string "chrom:start-end" that contains up to `chunk_size` variants.
    The count is the number of variants included in that region.
    """
    vcf = VCF(vcf_path)

    current_chrom: Optional[str] = None
    region_start: Optional[int] = None
    last_pos: Optional[int] = None
    current_count: int = 0
    next_region_start: Optional[int] = None

    for variant in vcf:
        chrom = variant.CHROM
        pos = int(variant.POS)
        print(f"{chrom}:{pos}\t{current_count}")
        if current_chrom is None:
            current_chrom = chrom
            next_region_start = 1

        # Chromosome change: flush any partial region
        if chrom != current_chrom:
            if current_count > 0 and region_start is not None and last_pos is not None:
                print(f"{current_chrom}:{region_start}-{last_pos}\t{current_count}")
                yield f"{current_chrom}:{region_start}-{last_pos}", current_count
            # reset for new chromosome
            current_chrom = chrom
            region_start = None
            last_pos = None
            current_count = 0
            next_region_start = 1

        # Start a new region if needed
        if region_start is None:
            region_start = next_region_start if next_region_start is not None else pos

        # Count this variant
        current_count += 1
        last_pos = pos

        # If we reached the target count, emit the region and reset
        if current_count == chunk_size:
            yield f"{current_chrom}:{region_start}-{last_pos}", current_count
            print(f"{current_chrom}:{region_start}-{last_pos}\t{current_count}")
            next_region_start = last_pos + 1
            region_start = None
            last_pos = None
            current_count = 0

    # Flush tail region at end of file
    if current_count > 0 and current_chrom is not None and region_start is not None and last_pos is not None:
        yield f"{current_chrom}:{region_start}-", current_count
        print(f"{current_chrom}:{region_start}-\t{current_count}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Output genomic regions each containing N variants from a VCF (count-based chunking).",
    )
    parser.add_argument("vcf", help="Input VCF/BCF (.vcf, .vcf.gz, or .bcf)")
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=3000000,
        help="Number of variants per region (default: 3,000,000)",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output file path (default: stdout)",
    )

    args = parser.parse_args()

    out_stream = sys.stdout
    close_stream = False
    total_variants = 0
    try:
        if args.output:
            out_stream = open(args.output, "w", encoding="utf-8")
            close_stream = True

        for region, region_count in generate_regions_by_variant_count(args.vcf, args.chunk_size):
            out_stream.write(region + "\n")
            total_variants += region_count
    finally:
        if close_stream:
            out_stream.close()

    # Write summary to stderr so it does not pollute stdout results
    sys.stderr.write(f"total_variants: {total_variants}\n")


if __name__ == "__main__":
    # usage : python3 chunk_vcf_by_count.py your.vcf.gz > regions.txt
    main()





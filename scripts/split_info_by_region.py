#!/usr/bin/env python3
"""
Filter a tab-delimited .info.gz file by one or more genomic regions.

Input file is expected to be gzipped with a header line followed by rows where
the first field contains a variant identifier including chromosome and position,
such as "chr9:245:TCAC:T" or "9:245:A:G". Only chromosome and position are used
for filtering.

Usage:
  python split_info_by_region.py -i input.info.gz -r chr9:200-3000 chr12:12000-
  python split_info_by_region.py -i input.info.gz -o prefix_ -r chr9:200-3000 chr12:12000-
  python split_info_by_region.py -i input.info.gz --name-style chunk -r chr12:1-100 chr12:200-300
  python split_info_by_region.py -i input.info.gz --name-style chunk -r chr12:1-100 chr12:200-300 -c chunk0 chunk5

For each provided region, an output file will be written containing the header and
matching rows. Naming styles:
 - region (default): chr<chrom>_<start>_<end|end>.gz
 - chunk:            chr<chrom>.chunk<index>.inf.gz (index follows provided region order)
"""

import argparse
import sys
import os
import gzip
import re
from typing import Optional, Tuple

def normalize_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return chrom

def parse_region(region_str: str) -> Tuple[str, int, Optional[int]]:
    """
    Accepts forms like:
      chr9:200-3000
      9:200-3000
      chrX:1-100
      chr12:12000-        (open-ended; through end of chromosome)
    Returns (chrom, start, end) with chrom normalized (no 'chr'), and ints for coords.
    If 'end' is omitted, returns end=None to indicate open-ended.
    """
    m = re.match(r'^(chr)?(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)?$', region_str, flags=re.IGNORECASE)
    if not m:
        raise ValueError(f"Invalid region format: {region_str}. Expected like 'chr9:200-3000' or 'chr12:12000-'.")
    chrom = normalize_chrom(m.group('chrom'))
    start = int(m.group('start'))
    end_str = m.group('end')
    end = int(end_str) if end_str is not None else None
    if end is not None and start > end:
        start, end = end, start
    return chrom, start, end

def region_output_name(region_str: str) -> str:
    chrom, start, end = parse_region(region_str)
    # Keep 'chr' prefix in filename for clarity
    end_part = "end" if end is None else str(end)
    return f"chr{chrom}_{start}_{end_part}.gz"

def record_in_region(snp_field: str, target_chrom: str, start: int, end: Optional[int]) -> bool:
    """
    snp_field examples:
      'chr9:245:TCAC:T'
      '9:245:A:G'
      'chr9:60930:G:C'
    We only need chromosome and position (the 1st and 2nd colon-separated tokens).
    """
    parts = re.split(r'[:_]', snp_field.strip())
    if len(parts) < 2:
        return False
    chrom_raw = parts[0]
    pos_str = parts[1]
    chrom = normalize_chrom(chrom_raw)
    if chrom != target_chrom:
        return False
    try:
        pos = int(pos_str)
    except ValueError:
        return False
    if end is None:
        return pos >= start
    return start <= pos <= end

def filter_by_region(input_gz: str, region_str: str, output_gz: Optional[str] = None) -> int:
    tgt_chrom, start, end = parse_region(region_str)
    if output_gz is None:
        output_gz = region_output_name(region_str)

    with gzip.open(input_gz, "rt") as fin:
        header = fin.readline()
        if not header:
            raise ValueError("Input file appears empty or missing header.")
        # Detect Rsq and Genotyped column indices
        # - Rsq: case-insensitive match
        # - Genotyped: case-insensitive match on column name; header itself is not modified
        header_cols = header.rstrip("\n").split("\t")
        rsq_idx = -1
        for i, col_name in enumerate(header_cols):
            if col_name.strip().lower() == "rsq":
                rsq_idx = i
                break
        genotyped_idx = -1
        for i, col_name in enumerate(header_cols):
            if col_name.strip().lower() == "genotyped":
                genotyped_idx = i
                break

        kept = 0
        with gzip.open(output_gz, "wt") as fout:
            fout.write(header)
            for line in fin:
                if not line.strip():
                    continue
                first_field = line.split("\t", 1)[0]
                if record_in_region(first_field, tgt_chrom, start, end):
                    parts = line.rstrip("\n").split("\t")
                    # Normalize SNP/variant id (first column) to 'chr<chrom>:<pos>:<ref>:<alt>'
                    if parts:
                        snp = parts[0].strip()
                        toks = re.split(r'[:_]', snp)
                        if len(toks) >= 4:
                            chrom_with_chr = f"chr{normalize_chrom(toks[0])}"
                            pos_token = toks[1]
                            try:
                                pos_token = str(int(pos_token))
                            except ValueError:
                                pass
                            ref = toks[2]
                            alt = toks[3]
                            parts[0] = f"{chrom_with_chr}:{pos_token}:{ref}:{alt}"
                        else:
                            # Fallback: replace underscores with colons verbatim
                            parts[0] = snp.replace("_", ":")
                    # Normalize Rsq '-' to '1'
                    if rsq_idx >= 0 and rsq_idx < len(parts) and parts[rsq_idx] == "-":
                        parts[rsq_idx] = "1"
                    # Normalize Genotyped values
                    if genotyped_idx >= 0 and genotyped_idx < len(parts):
                        gv = parts[genotyped_idx]
                        if gv == "Imputed":
                            parts[genotyped_idx] = "IMPUTED"
                        elif gv == "Genotyped":
                            parts[genotyped_idx] = "TYPED"
                        elif gv == "Typed_Only":
                            parts[genotyped_idx] = "TYPED_ONLY"
                    line = "\t".join(parts) + "\n"
                    fout.write(line)
                    kept += 1

    end_str = "" if end is None else str(end)
    print(f"Wrote {output_gz} with {kept} rows (region chr{tgt_chrom}:{start}-{end_str}).")
    return kept

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write one filtered .info.gz per region from an input .info.gz file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        help="Optional output filename prefix to prepend to each generated file",
    )
    parser.add_argument(
        "--name-style",
        choices=["region", "chunk"],
        default="region",
        help="Output naming style: 'region' (chr<chrom>_<start>_<end|end>.gz) or 'chunk' (chr<chrom>.chunk<index>.inf.gz)",
    )
    parser.add_argument(
        "-c",
        "--chunk-index",
        nargs="+",
        dest="chunk_indices",
        help="Optional chunk indices matching -r regions, e.g., 0 5 or chunk0 chunk5",
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        dest="input_flag",
        help="Path to input .info.gz file",
    )
    parser.add_argument(
        "-r",
        "--region",
        action="append",
        nargs="+",
        dest="regions",
        required=True,
        help="One or more regions like chr9:200-3000 chr12:12000-. Can repeat -r",
    )
    return parser.parse_args()

def main() -> None:
    args = parse_args()

    input_gz = args.input_flag
    # Flatten list of lists from action=append with nargs="+"
    regions = [r for group in args.regions for r in group]
    output_prefix = args.output or ""
    name_style = args.name_style
    raw_indices = args.chunk_indices or []

    def parse_chunk_token(tok: str) -> Optional[int]:
        t = tok.strip().lower()
        if t.startswith("chunk"):
            t = t[5:]
        try:
            return int(t)
        except ValueError:
            return None

    explicit_indices = []
    for tok in raw_indices:
        val = parse_chunk_token(tok)
        if val is None:
            sys.stderr.write(f"Error: invalid chunk index token: {tok}\n")
            sys.exit(1)
        explicit_indices.append(val)
    if name_style == "chunk" and explicit_indices and len(explicit_indices) != len(regions):
        sys.stderr.write("Error: number of --chunk-index values must match number of regions when using chunk naming.\n")
        sys.exit(1)

    if not os.path.exists(input_gz):
        sys.stderr.write(f"Error: input file not found: {input_gz}\n")
        sys.exit(1)

    for idx, region in enumerate(regions):
        try:
            if name_style == "chunk":
                chrom, _, _ = parse_region(region)
                use_index = explicit_indices[idx] if explicit_indices else idx
                base_name = f"chr{chrom}.chunk{use_index}.info.gz"
            else:
                base_name = region_output_name(region)
            if output_prefix and (output_prefix.endswith(os.sep) or output_prefix.endswith("/")):
                outname = os.path.join(output_prefix, base_name)
            else:
                outname = f"{output_prefix}{base_name}"
            filter_by_region(input_gz, region, outname)
        except Exception as e:
            sys.stderr.write(f"Failed for region {region}: {e}\n")
            sys.exit(1)

if __name__ == "__main__":
    main()



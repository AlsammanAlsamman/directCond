#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract GWAS and reference-panel subset using explicit chr:start-end bounds"
    )
    parser.add_argument("--gwas-file", required=True)
    parser.add_argument("--ref-prefix", required=True)
    parser.add_argument("--chr", required=True)
    parser.add_argument("--start", type=int)
    parser.add_argument("--end", type=int)
    # Legacy compatibility: allow pos/window_bp if start/end are not provided.
    parser.add_argument("--pos", type=int)
    parser.add_argument("--window-bp", type=int)
    parser.add_argument("--out-gwas-window", required=True)
    parser.add_argument("--out-ref-prefix", required=True)
    parser.add_argument("--out-done", required=True)
    parser.add_argument("--plink-bin", default="plink")
    return parser.parse_args()


def find_column(fieldnames: list[str], candidates: list[str]) -> str:
    lowered = {name.lower(): name for name in fieldnames}
    for candidate in candidates:
        if candidate in lowered:
            return lowered[candidate]
    raise ValueError(f"None of columns found: {candidates}")


def extract_gwas_window(gwas_file: str, chrom: str, start_bp: int, end_bp: int, out_file: str) -> int:
    Path(out_file).parent.mkdir(parents=True, exist_ok=True)

    # Fast path for large GWAS files: delegate row filtering to awk.
    try:
        with open(gwas_file, "r", encoding="utf-8", newline="") as src:
            header_line = src.readline().rstrip("\n\r")
        if header_line:
            header_cols = header_line.split("\t")
            chr_col_name = find_column(header_cols, ["chrom", "chr"])
            pos_col_name = find_column(header_cols, ["pos", "position", "bp"])
            chr_idx = header_cols.index(chr_col_name) + 1
            pos_idx = header_cols.index(pos_col_name) + 1
            target_chr = chrom.strip().lower().replace("chr", "")

            awk_program = (
                "BEGIN { FS=\"\\t\"; OFS=\"\\t\"; n=0 } "
                "NR==1 { print; next } "
                "{ c=tolower($C); gsub(/^chr/, \"\", c); p=$P+0; "
                "if (c==T && p>=S && p<=E) { print; n++ } } "
                "END { print n > \"/dev/stderr\" }"
            )
            awk_program = awk_program.replace("$C", f"${chr_idx}").replace("$P", f"${pos_idx}")

            with open(out_file, "w", encoding="utf-8", newline="") as dst:
                proc = subprocess.run(
                    [
                        "awk",
                        "-v",
                        f"T={target_chr}",
                        "-v",
                        f"S={start_bp}",
                        "-v",
                        f"E={end_bp}",
                        awk_program,
                        gwas_file,
                    ],
                    stdout=dst,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )
            count_text = (proc.stderr or "").strip().splitlines()
            if count_text:
                try:
                    return int(count_text[-1].strip())
                except ValueError:
                    pass
            # Fallback if awk stderr count parsing failed.
            return sum(1 for _ in open(out_file, "r", encoding="utf-8")) - 1
    except FileNotFoundError:
        # awk not available; fallback to pure-Python parsing below.
        pass
    except Exception:
        # Any awk/header issue: fallback to pure-Python parsing below.
        pass

    rows = 0

    with open(gwas_file, "r", encoding="utf-8", newline="") as src, open(
        out_file, "w", encoding="utf-8", newline=""
    ) as dst:
        reader = csv.DictReader(src, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("GWAS file has no header")

        chr_col = find_column(reader.fieldnames, ["chrom", "chr"])
        pos_col = find_column(reader.fieldnames, ["pos", "position", "bp"])

        writer = csv.DictWriter(dst, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()

        for row in reader:
            row_chr = str(row[chr_col]).strip().lower().replace("chr", "")
            target_chr = chrom.strip().lower().replace("chr", "")
            try:
                row_pos = int(float(row[pos_col]))
            except (TypeError, ValueError):
                continue

            if row_chr == target_chr and start_bp <= row_pos <= end_bp:
                writer.writerow(row)
                rows += 1

    return rows


def run_plink_subset(plink_bin: str, ref_prefix: str, chrom: str, start_bp: int, end_bp: int, out_prefix: str) -> None:
    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        plink_bin,
        "--bfile",
        ref_prefix,
        "--chr",
        chrom,
        "--from-bp",
        str(start_bp),
        "--to-bp",
        str(end_bp),
        "--make-bed",
        "--out",
        out_prefix,
    ]
    subprocess.run(cmd, check=True)


def main() -> int:
    args = parse_args()
    if args.start is not None and args.end is not None:
        start_bp = args.start
        end_bp = args.end
    elif args.pos is not None and args.window_bp is not None:
        half_window = args.window_bp // 2
        start_bp = max(1, args.pos - half_window)
        end_bp = args.pos + half_window
    else:
        raise ValueError("Provide either --start/--end or legacy --pos/--window-bp")

    if start_bp > end_bp:
        raise ValueError(f"Invalid interval: start ({start_bp}) > end ({end_bp})")

    extracted = extract_gwas_window(args.gwas_file, args.chr, start_bp, end_bp, args.out_gwas_window)
    if extracted == 0:
        print("No GWAS rows extracted for requested region", file=sys.stderr)

    run_plink_subset(args.plink_bin, args.ref_prefix, args.chr, start_bp, end_bp, args.out_ref_prefix)

    Path(args.out_done).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_done).write_text(
        f"ok\tchr={args.chr}\tstart={start_bp}\tend={end_bp}\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

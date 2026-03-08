#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Harmonize GWAS window with reference subset BIM and add rsid."
    )
    parser.add_argument("--gwas-window", required=True, help="Input GWAS window TSV")
    parser.add_argument("--ref-bim", required=True, help="Reference subset BIM file")
    parser.add_argument("--out-gwas", required=True, help="Output harmonized GWAS TSV")
    parser.add_argument("--out-done", required=True, help="Done marker file")
    return parser.parse_args()


def norm_chr(value: str) -> str:
    text = str(value).strip().lower()
    if text.startswith("chr"):
        return text[3:]
    return text


def norm_allele(value: str) -> str:
    return str(value).strip().upper()


def find_column(fieldnames: List[str], candidates: List[str]) -> str:
    lowered = {name.lower(): name for name in fieldnames}
    for candidate in candidates:
        if candidate in lowered:
            return lowered[candidate]
    raise ValueError(f"Required column not found; expected one of: {candidates}")


def load_bim_index(bim_path: str) -> Dict[Tuple[str, int], List[Tuple[str, str, str]]]:
    index: Dict[Tuple[str, int], List[Tuple[str, str, str]]] = {}
    with open(bim_path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            chrom = norm_chr(parts[0])
            rsid = parts[1]
            pos = int(parts[3])
            allele1 = norm_allele(parts[4])
            allele2 = norm_allele(parts[5])
            key = (chrom, pos)
            index.setdefault(key, []).append((rsid, allele1, allele2))
    return index


def harmonize_gwas(gwas_path: str, bim_index: Dict[Tuple[str, int], List[Tuple[str, str, str]]], out_path: str) -> Tuple[int, int, int]:
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)

    with open(gwas_path, "r", encoding="utf-8", newline="") as src, open(
        out_path, "w", encoding="utf-8", newline=""
    ) as dst:
        reader = csv.DictReader(src, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("GWAS file has no header")

        chr_col = find_column(reader.fieldnames, ["chrom", "chr"])
        pos_col = find_column(reader.fieldnames, ["pos", "position", "bp"])
        ea_col = find_column(reader.fieldnames, ["ea", "effect_allele", "effectallele", "a1"])
        nea_col = find_column(reader.fieldnames, ["nea", "other_allele", "otherallele", "a2", "non_effect_allele"])
        beta_col = find_column(reader.fieldnames, ["beta", "b", "effect"])

        output_fields = list(reader.fieldnames)
        if "rsid" not in output_fields:
            output_fields.append("rsid")
        if "harmonization_status" not in output_fields:
            output_fields.append("harmonization_status")

        writer = csv.DictWriter(dst, fieldnames=output_fields, delimiter="\t")
        writer.writeheader()

        total = 0
        aligned = 0
        flipped = 0

        for row in reader:
            total += 1
            row_chr = norm_chr(row[chr_col])
            try:
                row_pos = int(float(row[pos_col]))
            except (TypeError, ValueError):
                row["harmonization_status"] = "invalid_pos"
                if "rsid" not in row or not row.get("rsid"):
                    row["rsid"] = ""
                writer.writerow(row)
                continue

            ea = norm_allele(row[ea_col])
            nea = norm_allele(row[nea_col])
            key = (row_chr, row_pos)
            candidates = bim_index.get(key, [])

            chosen_rsid = ""
            status = "no_ref_match"
            was_flipped = False

            for rsid, allele1, allele2 in candidates:
                if ea == allele2 and nea == allele1:
                    chosen_rsid = rsid
                    status = "aligned"
                    break
                if ea == allele1 and nea == allele2:
                    chosen_rsid = rsid
                    status = "flipped"
                    was_flipped = True
                    break

            if was_flipped:
                old_beta = row[beta_col]
                try:
                    row[beta_col] = str(-float(old_beta))
                except (TypeError, ValueError):
                    row[beta_col] = old_beta
                row[ea_col], row[nea_col] = row[nea_col], row[ea_col]
                flipped += 1
            elif status == "aligned":
                aligned += 1

            row["rsid"] = chosen_rsid
            row["harmonization_status"] = status
            writer.writerow(row)

    return total, aligned, flipped


def main() -> int:
    args = parse_args()
    bim_index = load_bim_index(args.ref_bim)
    total, aligned, flipped = harmonize_gwas(args.gwas_window, bim_index, args.out_gwas)

    Path(args.out_done).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_done).write_text(
        f"ok\ttotal={total}\taligned={aligned}\tflipped={flipped}\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

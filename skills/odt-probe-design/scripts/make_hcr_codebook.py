#!/usr/bin/env python3
"""Write an HCR codebook (one initiator per gene) for the ODT hcr_probe_designer.

Usage:
  make_hcr_codebook.py --genes genes.txt --initiators hcr_initiators.tsv --out hcr_codebook.tsv \
      [--assign GENE=bit_2 ...]

- genes.txt: the gene list used as `file_regions` (one gene ID per line).
- hcr_initiators.tsv: the initiator table used as `file_initiator_table` (CSV/TSV with a `bit` column).
- --assign: put a gene on a specific initiator; all other genes get the remaining initiators in
  table order.

The output has a `region_id` column and one `bit_*` column per initiator, with exactly one 1 per
row, which is what ODT v2.3.0 expects (rows with zero or several 1s are silently dropped by ODT).
Only uses the standard library.
"""

import argparse
import csv
import sys


def read_genes(file_genes: str) -> list[str]:
    genes: list[str] = []
    with open(file_genes) as handle:
        for line in handle:
            gene = line.strip()
            if gene and gene not in genes:
                genes.append(gene)
    if not genes:
        sys.exit(f"No gene IDs found in {file_genes}")
    return genes


def read_bits(file_initiators: str) -> list[str]:
    with open(file_initiators, newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,;")
        rows = list(csv.DictReader(handle, dialect=dialect))
    if not rows or "bit" not in rows[0]:
        sys.exit(f"{file_initiators} needs a 'bit' column (and initiator_L_sequence, initiator_R_sequence)")
    bits = [row["bit"].strip() for row in rows]
    invalid = [bit for bit in bits if not bit.startswith("bit_")]
    if invalid:
        sys.exit(f"Bit names must start with 'bit_', found: {invalid}")
    if len(set(bits)) != len(bits):
        sys.exit(f"Duplicate bit names in {file_initiators}: {bits}")
    return bits


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--genes", required=True, help="gene list file (one gene ID per line)")
    parser.add_argument("--initiators", required=True, help="initiator table with a 'bit' column")
    parser.add_argument("--out", required=True, help="output codebook TSV")
    parser.add_argument("--assign", action="append", default=[], metavar="GENE=BIT", help="fixed assignment")
    args = parser.parse_args()

    genes = read_genes(args.genes)
    bits = read_bits(args.initiators)

    assignment: dict[str, str] = {}
    for item in args.assign:
        if "=" not in item:
            sys.exit(f"--assign expects GENE=BIT, got '{item}'")
        gene, bit = (x.strip() for x in item.split("=", 1))
        if gene not in genes:
            sys.exit(f"--assign: gene '{gene}' is not in {args.genes}")
        if bit not in bits:
            sys.exit(f"--assign: '{bit}' is not in the initiator table (available: {', '.join(bits)})")
        if bit in assignment.values():
            sys.exit(f"--assign: '{bit}' is used for more than one gene")
        assignment[gene] = bit

    free_bits = [bit for bit in bits if bit not in assignment.values()]
    unassigned = [gene for gene in genes if gene not in assignment]
    if len(unassigned) > len(free_bits):
        sys.exit(
            f"{len(genes)} genes but only {len(bits)} initiators ({len(free_bits)} still free for "
            f"{len(unassigned)} genes). Add initiators to the table or reduce the gene list."
        )
    assignment.update(zip(unassigned, free_bits))

    with open(args.out, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["region_id"] + bits)
        for gene in genes:
            writer.writerow([gene] + [1 if bit == assignment[gene] else 0 for bit in bits])

    print(f"Wrote {args.out}:")
    for gene in genes:
        print(f"  {gene}\t{assignment[gene]}")


if __name__ == "__main__":
    main()

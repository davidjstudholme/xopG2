import sys
from Bio import SeqIO
from pathlib import Path

def find_gene_feature(record, gene_name):
    """
    Find the first gene/CDS feature matching gene_name.
    Returns (start, end) in 0-based coordinates, or None.
    """
    for feature in record.features:
        if feature.type not in {"gene", "CDS"}:
            continue

        genes = feature.qualifiers.get("gene", [])
        if gene_name in genes:
            return int(feature.location.start), int(feature.location.end)

    return None

def main(gene_a, gene_b):
    for gb_file in Path(".").glob("*.gb*"):
        print(f"Processing {gb_file}")

        records = list(SeqIO.parse(gb_file, "genbank"))
        if not records:
            print("  ⚠ No records found")
            continue

        for record in records:
            loc_a = find_gene_feature(record, gene_a)
            loc_b = find_gene_feature(record, gene_b)

            if not loc_a or not loc_b:
                print(f"  ✖ Missing {gene_a} or {gene_b}")
                continue

            start = min(loc_a[0], loc_b[0])
            end   = max(loc_a[1], loc_b[1])

            subrecord = record[start:end]

            subrecord.id = f"{record.id}_{gene_a}_to_{gene_b}"
            subrecord.name = subrecord.id
            subrecord.description = f"Region from {gene_a} to {gene_b}"

            out_file = gb_file.with_suffix(f".{gene_a}_to_{gene_b}.gbk")
            SeqIO.write(subrecord, out_file, "genbank")

            print(f"  ✔ Wrote {out_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_region.py <geneA> <geneB>")
        sys.exit(1)

    geneA = sys.argv[1]
    geneB = sys.argv[2]

    main(geneA, geneB)


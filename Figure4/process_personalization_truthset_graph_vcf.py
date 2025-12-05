#!/usr/bin/env python3
import argparse
import gzip
import csv

def compare_genotypes(gt1, gt2):
    """
    Compare two haploid genotypes using the following encoding:
      0: both missing ('.')
      1: same genotype
      2: one missing
      3: different genotype
    """
    if gt1 == '.' and gt2 == '.':
        return '0'
    if gt1 == '.' or gt2 == '.':
        return '2'
    if gt1 == gt2:
        return '1'
    return '3'

def main():
    parser = argparse.ArgumentParser(
        description="Process a gzipped VCF to compute pairwise haplotype genotype comparisons."
    )
    parser.add_argument("input", help="Input gzipped VCF file")
    parser.add_argument("output", help="Output tab-delimited file")
    args = parser.parse_args()

    with gzip.open(args.input, 'rt') as infile, open(args.output, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        # Write header for the output file.
        out_header = [
            "CHROM", "POS", "REF.Length", "ALT.Length", "DELTA.Length",
            "h1_h2", "h1_h3", "h1_h4", "h2_h3", "h2_h4", "h3_h4"
        ]
        writer.writerow(out_header)

        sample_indices = None  # will hold indices for the 4 sample columns

        for line in infile:
            line = line.strip()
            if not line:
                continue

            # Skip meta-information lines
            if line.startswith("##"):
                continue

            # Header line with column names
            if line.startswith("#CHROM"):
                parts = line.split("\t")
                # Assuming the first 9 columns are fixed VCF columns,
                # the remaining columns are the samples.
                sample_indices = list(range(9, len(parts)))
                if len(sample_indices) != 4:
                    raise ValueError("Expected exactly 4 sample columns, but found {}.".format(len(sample_indices)))
                continue

            # Process variant lines.
            parts = line.split("\t")
            chrom = parts[0]
            pos = parts[1]
            ref = parts[3]
            alt = parts[4]

            ref_len = len(ref)

            # Handle multiple ALT alleles
            alt_alleles = alt.split(',')
            alt_lengths = [len(a) for a in alt_alleles]
            alt_lengths_str = ",".join(str(l) for l in alt_lengths)
            max_alt_length = max(alt_lengths)
            delta_len = max_alt_length - ref_len

            # Extract the genotype for each haplotype.
            # Here we assume that the genotype is the first colon-separated field.
            haplotypes = []
            for idx in sample_indices:
                sample_field = parts[idx]
                gt = sample_field.split(":")[0]
                haplotypes.append(gt)

            # Map the haplotypes to h1, h2, h3, h4 in the order they appear.
            h1, h2, h3, h4 = haplotypes

            # Compute pairwise comparisons.
            comp_h1_h2 = compare_genotypes(h1, h2)
            comp_h1_h3 = compare_genotypes(h1, h3)
            comp_h1_h4 = compare_genotypes(h1, h4)
            comp_h2_h3 = compare_genotypes(h2, h3)
            comp_h2_h4 = compare_genotypes(h2, h4)
            comp_h3_h4 = compare_genotypes(h3, h4)

            # Write the computed values along with variant information.
            writer.writerow([
                chrom, pos, ref_len, alt_lengths_str, delta_len,
                comp_h1_h2, comp_h1_h3, comp_h1_h4, comp_h2_h3, comp_h2_h4, comp_h3_h4
            ])

if __name__ == '__main__':
    main()

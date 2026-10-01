#!/usr/bin/env python3
"""
在 SNV 位点统计 tumor 和 normal BAM/CRAM 中的 coverage、reads 平均 MAPQ，
以及 ref / alt supporting reads 数量。

输入位点文件支持:
  - VCF / VCF.gz：取 CHROM、POS、REF、ALT（多等位 ALT 拆成多行输出）
  - BED（或类 BED）：取第 3 列作为 1-based 位置
      标准 BED (start=pos-1, end=pos) 和 start==end==pos 的写法都适用
      REF/ALT 来源（按优先级）：
        1) --ref-col / --alt-col 指定的列（1-based 列号）
        2) REF 从 -r 参考 fasta 取；ALT 取 tumor 中出现最多的非 REF 碱基

输出 TSV 列:
  chrom pos ref alt
  tumor_depth  tumor_mean_mapq  tumor_mapq0_frac  tumor_ref_reads  tumor_alt_reads  tumor_vaf
  normal_depth normal_mean_mapq normal_mapq0_frac normal_ref_reads normal_alt_reads normal_vaf

用法:
  python snv_mapq_depth.py -s snv.vcf.gz -t tumor.bam -n normal.bam -o out.tsv
  python snv_mapq_depth.py -s snv.bed -t tumor.bam -n normal.bam -r ref.fa -o out.tsv
  python snv_mapq_depth.py -s snv.bed --ref-col 4 --alt-col 5 -t T.bam -n N.bam -o out.tsv
"""
import argparse
import gzip
import sys
from collections import Counter

import pysam

BASES = "ACGT"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-s", "--sites", required=True, help="SNV 位点文件 (BED 或 VCF[.gz])")
    p.add_argument("-t", "--tumor", required=True, help="tumor BAM/CRAM（需要索引）")
    p.add_argument("-n", "--normal", required=True, help="normal BAM/CRAM（需要索引）")
    p.add_argument("-r", "--reference", default=None,
                   help="参考 fasta（CRAM 必需；BED 无 REF 列时用于取 REF 碱基，需 .fai）")
    p.add_argument("-o", "--output", default="-", help="输出 TSV，默认 stdout")
    p.add_argument("--ref-col", type=int, default=None, help="BED 中 REF 碱基所在列（1-based）")
    p.add_argument("--alt-col", type=int, default=None, help="BED 中 ALT 碱基所在列（1-based）")
    p.add_argument("--min-mapq", type=int, default=0,
                   help="只统计 MAPQ >= 此值的 reads（默认 0）")
    p.add_argument("--min-bq", type=int, default=0,
                   help="该位点碱基质量 >= 此值才计数（默认 0；samtools mpileup 默认 13）")
    p.add_argument("--exclude-flags", type=lambda x: int(x, 0), default=0x704,
                   help="排除的 SAM flag，默认 0x704 = UNMAP|SECONDARY|QCFAIL|DUP")
    p.add_argument("--chrom-map", default=None,
                   help="可选：两列文件 <位点文件中的染色体名> <BAM 中的染色体名>")
    return p.parse_args()


def load_chrom_map(path):
    m = {}
    if path:
        with open(path) as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    a, b = line.split()[:2]
                    m[a] = b
    return m


def read_sites(path, ref_col, alt_col):
    """产出 (chrom, pos_1based, ref 或 None, alt 或 None)，去重。"""
    is_vcf = ".vcf" in path
    opener = gzip.open if path.endswith(".gz") else open
    seen = set()
    with opener(path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            c = line.rstrip("\n").split("\t")
            chrom = c[0]
            if is_vcf:
                pos, ref, alts = int(c[1]), c[3].upper(), c[4].upper().split(",")
            else:
                pos = int(c[2])
                ref = c[ref_col - 1].upper() if ref_col else None
                alts = [c[alt_col - 1].upper()] if alt_col else [None]
            for alt in alts:
                key = (chrom, pos, alt)
                if key in seen:
                    continue
                seen.add(key)
                yield chrom, pos, ref, alt


def pileup_site(bam, chrom, pos, args):
    """返回 (mapq 列表, 碱基计数 Counter)；染色体不存在时返回 None。"""
    if chrom not in bam.references:
        return None
    mapqs, bases = [], Counter()
    for col in bam.pileup(chrom, pos - 1, pos, truncate=True,
                          stepper="nofilter", min_base_quality=0,
                          ignore_overlaps=False, ignore_orphans=False,
                          max_depth=10_000_000):
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            aln = pr.alignment
            if aln.flag & args.exclude_flags:
                continue
            if aln.mapping_quality < args.min_mapq:
                continue
            qp = pr.query_position
            if args.min_bq > 0:
                quals = aln.query_qualities
                if quals is not None and quals[qp] < args.min_bq:
                    continue
            mapqs.append(aln.mapping_quality)
            bases[aln.query_sequence[qp].upper()] += 1
    return mapqs, bases


def summarize(res, ref, alt):
    if res is None:
        return ["NA"] * 6
    mapqs, bases = res
    depth = len(mapqs)
    if depth == 0:
        return ["0", "NA", "NA", "0", "0", "NA"]
    ref_n = bases.get(ref, 0) if ref else 0
    alt_n = bases.get(alt, 0) if alt else 0
    return [str(depth),
            f"{sum(mapqs) / depth:.2f}",
            f"{sum(1 for q in mapqs if q == 0) / depth:.4f}",
            str(ref_n) if ref else "NA",
            str(alt_n) if alt else "NA",
            f"{alt_n / depth:.4f}" if alt else "NA"]


def main():
    args = parse_args()
    cmap = load_chrom_map(args.chrom_map)
    is_vcf = ".vcf" in args.sites
    if not is_vcf and not args.ref_col and not args.reference:
        sys.exit("[error] BED 输入没有 REF 列时需要 -r 参考 fasta 来获取 REF 碱基")

    tumor = pysam.AlignmentFile(args.tumor, reference_filename=args.reference)
    normal = pysam.AlignmentFile(args.normal, reference_filename=args.reference)
    fasta = pysam.FastaFile(args.reference) if args.reference else None

    out = sys.stdout if args.output == "-" else open(args.output, "w")
    header = ["chrom", "pos", "ref", "alt"]
    for s in ("tumor", "normal"):
        header += [f"{s}_depth", f"{s}_mean_mapq", f"{s}_mapq0_frac",
                   f"{s}_ref_reads", f"{s}_alt_reads", f"{s}_vaf"]
    out.write("\t".join(header) + "\n")

    warned = set()
    for chrom, pos, ref, alt in read_sites(args.sites, args.ref_col, args.alt_col):
        bam_chrom = cmap.get(chrom, chrom)
        t = pileup_site(tumor, bam_chrom, pos, args)
        n = pileup_site(normal, bam_chrom, pos, args)
        if (t is None or n is None) and bam_chrom not in warned:
            print(f"[warn] 染色体 {bam_chrom} 不在 BAM header 中，输出 NA（可用 --chrom-map 转换）",
                  file=sys.stderr)
            warned.add(bam_chrom)

        # BED 没有 REF：从 fasta 取
        if ref is None and fasta is not None and bam_chrom in fasta.references:
            ref = fasta.fetch(bam_chrom, pos - 1, pos).upper()
        # 没有 ALT：取 tumor 中最多的非 REF 碱基
        if alt is None and t is not None:
            cand = [(cnt, b) for b, cnt in t[1].items() if b in BASES and b != ref]
            alt = max(cand)[1] if cand else None

        out.write("\t".join([chrom, str(pos), ref or "NA", alt or "NA"]
                            + summarize(t, ref, alt) + summarize(n, ref, alt)) + "\n")

    if out is not sys.stdout:
        out.close()
    tumor.close()
    normal.close()
    if fasta:
        fasta.close()


if __name__ == "__main__":
    main()

#!/bin/bash
lifted_vcf="fp.vcf.gz"
original_vcf="/scratch/qfu/COLO829T_fiberseq/COLO829BL_per_R1_20x_2_new_matched/COLO829T_COLO829BL_per_20x_new_2_PacBio_deepsomatic_pass_norm_snvs_pileup_filtered_alt2.vcf.gz"

# 提取 OriginalContig 和 OriginalStart
bcftools query -f '%INFO/OriginalContig\t%INFO/OriginalStart\n' $lifted_vcf \
  | awk 'NF==2 {print $1"\t"$2-1"\t"$2}' > tmp_original_coords.bed

# 从原始 VCF 提取对应位置的变异
bcftools view -R tmp_original_coords.bed -Oz -o fp_20x_2_unique_on_20x_2.vcf.gz $original_vcf

# 索引输出 VCF
bcftools index -t fp_20x_2_unique_on_20x_2.vcf.gz
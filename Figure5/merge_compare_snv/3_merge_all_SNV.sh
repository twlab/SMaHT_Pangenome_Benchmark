module load bcftools

bcftools concat -a -O z -o all_inferred_DSA_SNV.vcf.gz tp-baseline_shared.vcf.gz fn.vcf.gz fp_20x_2_unique_on_20x_2.vcf.gz /scratch/qfu/COLO829T_fiberseq/COLO829BL_per_R1_20x_2_new_matched/liftover_COLO829BL_R1_chm13_clipped_20x_new_1/source_to_target-rejected_pass_snvs_alt2_w_origin.vcf.gz
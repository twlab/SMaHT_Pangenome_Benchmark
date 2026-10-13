#!/bin/bash
#SBATCH --job-name=rtg_liftover_hg38
#SBATCH --partition=general
#SBATCH --cpus-per-task=10
#SBATCH --time=6-3:00:00
#SBATCH --output=./rtg_liftover_hg38.log
#SBATCH --error=./rtg_liftover_hg38.err
#SBATCH --mem=60G


ref="hg38"
threads=10
template="/scratch/qfu/reference/hg38_rtg"
baseline="/scratch/qfu/COLO829T_fiberseq/hg38_matched/COLO829T_hg38_PacBio_deepsomatic_pass_norm_snvs_pileup_filtered_alt2.vcf.gz"
calls="/scratch/qfu/COLO829T_fiberseq/COLO829BL_per_R1_20x_2_new_matched/liftover_COLO829BL_R1_chm13_clipped_20x_new_1/rtg_COLO829BL_R1_chm13_clipped_20x_new_1/rtg_1_2_alt2/liftover_merged_hg38/all_inferred_DSA_SNV_lift_hg38.vcf.gz"
output_dir="./rtg_liftover_hg38"

/scratch/qfu/rtg-tools-3.12.1/rtg vcfeval -T "${threads}" \
    --template "${template}" \
    --baseline "${baseline}" \
    --calls "${calls}" \
    --sample ALT,ALT \
    --squash-ploidy \
    --output "${output_dir}"

#!/bin/bash
#SBATCH --job-name=rtg_1_2_alt2
#SBATCH --partition=general
#SBATCH --cpus-per-task=10
#SBATCH --time=6-3:00:00
#SBATCH --output=./rtg_1_2_alt2.log
#SBATCH --error=./rtg_1_2_alt2.err
#SBATCH --mem=60G


ref="20x_1"
threads=10
template="/scratch/qfu/COLO829BL_WGS/personalized/final_personalization/rtg/hs_hprc-v1.1-mc-chm13_20x_1/"
baseline="/scratch/qfu/COLO829T_fiberseq/COLO829BL_per_R1_20x_1_new_matched/COLO829T_COLO829BL_per_20x_new_1_PacBio_deepsomatic_pass_norm_snvs_pileup_filtered_alt2.vcf.gz"
calls="/scratch/qfu/COLO829T_fiberseq/COLO829BL_per_R1_20x_2_new_matched/liftover_COLO829BL_R1_chm13_clipped_20x_new_1/COLO829T_COLO829BL_per_20x_new_2_PacBio_deepsomatic_alt2_lift_COLO829BL_per_R1_20x_2_pass_snvs_w_origin.vcf.gz"
output_dir="./rtg_1_2_alt2"

/scratch/qfu/rtg-tools-3.12.1/rtg vcfeval -T "${threads}" \
    --template "${template}" \
    --baseline "${baseline}" \
    --calls "${calls}" \
    --sample ALT,ALT \
    --squash-ploidy \
    --output "${output_dir}"

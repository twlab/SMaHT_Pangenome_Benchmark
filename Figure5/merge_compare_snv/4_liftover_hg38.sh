#!/bin/bash
#SBATCH --job-name=liftover
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=6-3:00:00
#SBATCH --output=./liftover_pass_snvs.log
#SBATCH --error=./liftover_pass_snvs.err
#SBATCH --mem=50G


java -jar /scratch/qfu/picard.jar LiftoverVcf \
    I=../all_inferred_DSA_SNV.vcf.gz \
    O=all_inferred_DSA_SNV_lift_hg38.vcf.gz \
    CHAIN=/scratch/qfu/COLO829BL_WGS/personalized/final_personalization/chain_file/hprc-v1.1-mc-chm13_20x_1_2_hg38.chain \
    REFERENCE_SEQUENCE=/scratch/qfu/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    REJECT=source_to_target-rejected_pass_snvs.vcf \
    TMP_DIR=/scratch/qfu/COLO829T_fiberseq/COLO829BL_per_R1_20x_2_new_matched/liftover_COLO829BL_R1_chm13_clipped_20x_new_1/rtg_COLO829BL_R1_chm13_clipped_20x_new_1/rtg_1_2_alt2/liftover_merged_hg38/tmp \
    RECOVER_SWAPPED_REF_ALT=true \
    WRITE_ORIGINAL_POSITION=true \
    WRITE_ORIGINAL_ALLELES=true

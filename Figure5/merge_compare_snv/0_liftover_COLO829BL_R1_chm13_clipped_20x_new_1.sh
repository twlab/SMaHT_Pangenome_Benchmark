#!/bin/bash
#SBATCH --job-name=liftover
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=6-3:00:00
#SBATCH --output=./liftover_pass_snvs_alt2_w_origin.log
#SBATCH --error=./liftover_pass_snvs_alt2_w_origin.err
#SBATCH --mem=100G

java -jar /scratch/qfu/picard.jar LiftoverVcf -I ../COLO829T_COLO829BL_per_20x_new_2_PacBio_deepsomatic_pass_norm_snvs_pileup_filtered_alt2.vcf.gz \
-O COLO829T_COLO829BL_per_20x_new_2_PacBio_deepsomatic_alt2_lift_COLO829BL_per_R1_20x_2_pass_snvs_w_origin.vcf.gz -C /scratch/qfu/COLO829BL_WGS/personalized/final_personalization/chain_file/COLO829BL_R1_chm13_clipped_20x_new_2_COLO829BL_R1_chm13_clipped_20x_new_1/chain_hprc-v1.1-mc-chm13_20x_2_hprc-v1.1-mc-chm13_20x_1/chainnet/liftover.chain \
-R /scratch/qfu/COLO829BL_WGS/personalized/final_personalization/hs_hprc-v1.1-mc-chm13_20x_1.fasta --RECOVER_SWAPPED_REF_ALT true \
--REJECT source_to_target-rejected_pass_snvs_alt2_w_origin.vcf --TMP_DIR /scratch/qfu/COLO829T_fiberseq/COLO829BL_per_R1_20x_2_new_matched/liftover_COLO829BL_R1_chm13_clipped_20x_new_1/tmp/  --WRITE_ORIGINAL_POSITION true --WRITE_ORIGINAL_ALLELES true

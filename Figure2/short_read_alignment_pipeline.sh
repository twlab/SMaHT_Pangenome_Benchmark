#!/bin/bash

# SLURM Job Directives
#SBATCH --job-name=variant_calling_pipeline
#SBATCH --partition=general
#SBATCH --array=0-7 # Adjust based on the number of _trimed_1.fastq files (8 in your example)
#SBATCH --cpus-per-task=10
#SBATCH --time=14-00:00:00 # Increased time for comprehensive pipeline
#SBATCH --output=./slurm_logs/pipeline_%A_%a.log
#SBATCH --error=./slurm_logs/pipeline_%A_%a.err
#SBATCH --mem=300G # Increased memory for all steps

# Define common directories and paths
export PATH="/ref/twlab/software/zxin/bwa/bwa-0.7.17:$PATH"
export PATH="/ref/twlab/software/zxin/samtools/bin:$PATH"
export PATH="/ref/twlab/software/zxin/java/jdk1.8.0_92/bin:$PATH"
export PATH="/ref/twlab/software/zxin/gatk/gatk-4.1.0.0:$PATH"

# Reference Genome and Known Sites Paths
GENOME_HG38="/ref/twlab/data/qfu/references/mgi_non_HG005_per_1/mgi_non_HG005_ragtag_hg005_1_sorted_normalized.fasta"
DBSNP_VCF="/ref/twlab/data/qfu/references/mgi_non_HG005_per_1/dbsnp_indels/Homo_sapiens_assembly38.dbsnp138_lift_recomb1.vcf.gz"
MILLS_VCF="/ref/twlab/data/qfu/references/mgi_non_HG005_per_1/dbsnp_indels/Mills_and_1000G_gold_standard.indels.hg38_lift_recomb1.vcf.gz"

# Pipeline Specifics
REF_GENOME_TAG="mgi_non_HG005_per_1" # Tag for the reference genome

# Base output directory where all processed data will be stored
BASE_OUTPUT_DIR="/scratch/twlab/qfu/real_mixture/mgi_non_HG005_per_1_short_bwa"
mkdir -p "${BASE_OUTPUT_DIR}/slurm_logs" # Create log directory upfront

# Temporary directory setup
CUSTOM_TMPDIR="/scratch/twlab/qfu/tmp"

# Create a unique temporary directory for this job (will be deleted on exit)
TMPDIR_GLOBAL=$(mktemp -d "${CUSTOM_TMPDIR}/pipeline_tmp.XXXXXXX")
export TMPDIR="${TMPDIR_GLOBAL}" # Set TMPDIR for tools that respect it

function cleanup() {
    echo "Cleaning up temporary directory: ${TMPDIR_GLOBAL}"
    rm -rf "${TMPDIR_GLOBAL}"
}
trap cleanup EXIT # Ensure cleanup function runs on script exit

# Dynamically get R1 files from your input data directory
# Ensure this path is correct and accessible
INPUT_DATA_DIR="/ref/twlab/data/zxin/real_mixture_SMaHT/data/short"
declare -a FILES=($(ls "${INPUT_DATA_DIR}"/*_trimed_1.fastq))

# Get SLURM array task ID
ARRAY_ID=${SLURM_ARRAY_TASK_ID}
CURRENT_FASTQ_R1_PATH="${FILES[${ARRAY_ID}]}"

# Extract just the R1 filename (e.g., LIB027514_223KKHLT4_S13_L006_trimed_1.fastq)
CURRENT_FASTQ_R1_FILENAME=$(basename "${CURRENT_FASTQ_R1_PATH}")

# Extract the base name for alignment (e.g., LIB027514_223KKHLT4_S13_L006_trimed)
BWA_FILE_PREFIX=$(echo "${CURRENT_FASTQ_R1_FILENAME}" | sed 's/_1.fastq$//')

# Extract the primary sample name (e.g., LIB027514) for output BAM prefixes
# This will be the main identifier for sample-specific folders
CURRENT_SAMPLE=$(echo "${BWA_FILE_PREFIX}" | cut -d'_' -f1)

# Full paths for R1 and R2
R1_FULL_PATH="${INPUT_DATA_DIR}/${BWA_FILE_PREFIX}_1.fastq"
R2_FULL_PATH="${INPUT_DATA_DIR}/${BWA_FILE_PREFIX}_2.fastq"

echo "Processing sample: ${CURRENT_SAMPLE} (from file: ${CURRENT_FASTQ_R1_FILENAME})"
echo "Temporary directory for this job: ${TMPDIR_GLOBAL}"
echo "----------------------------------------------------"

# --- Define Output Directories for Each Step ---
BWA_OUTPUT_DIR="${BASE_OUTPUT_DIR}/1_bwa_align/${CURRENT_SAMPLE}"
ARG_OUTPUT_DIR="${BASE_OUTPUT_DIR}/2_add_read_groups/${CURRENT_SAMPLE}"
MD_OUTPUT_DIR="${BASE_OUTPUT_DIR}/3_mark_duplicates/${CURRENT_SAMPLE}"
IR_OUTPUT_DIR="${BASE_OUTPUT_DIR}/4_indel_realign/${CURRENT_SAMPLE}"
BQSR_OUTPUT_DIR="${BASE_OUTPUT_DIR}/5_bqsr/${CURRENT_SAMPLE}"
MERGE_OUTPUT_DIR="${BASE_OUTPUT_DIR}/6_merged_bam" # This will be created once at the end

# Create directories for the current sample's output
mkdir -p "${BWA_OUTPUT_DIR}"
mkdir -p "${ARG_OUTPUT_DIR}"
mkdir -p "${MD_OUTPUT_DIR}"
mkdir -p "${IR_OUTPUT_DIR}"
mkdir -p "${BQSR_OUTPUT_DIR}"

# --- Step 1: BWA Alignment ---
echo "Starting BWA Alignment for sample: ${CURRENT_SAMPLE}"

BWA_THREADS=${SLURM_CPUS_PER_TASK} # Use allocated CPUs for BWA

# Create symbolic links to the input FASTQ files in the current working directory
# BWA will be run in the BWA_OUTPUT_DIR, so links are made there.
(cd "${BWA_OUTPUT_DIR}" && ln -s "${R1_FULL_PATH}" "${BWA_FILE_PREFIX}_1.fastq")
(cd "${BWA_OUTPUT_DIR}" && ln -s "${R2_FULL_PATH}" "${BWA_FILE_PREFIX}_2.fastq")

# Perform alignment and sorting
# Output is directed to the BWA_OUTPUT_DIR
bwa mem -K 10000000 -t ${BWA_THREADS} "${GENOME_HG38}" "${BWA_OUTPUT_DIR}/${BWA_FILE_PREFIX}_1.fastq" "${BWA_OUTPUT_DIR}/${BWA_FILE_PREFIX}_2.fastq" | \
samtools sort --no-PG -@ ${BWA_THREADS} -o "${BWA_OUTPUT_DIR}/${CURRENT_SAMPLE}_sorted.bam" -

# Index the sorted BAM file
samtools index -@ ${BWA_THREADS} "${BWA_OUTPUT_DIR}/${CURRENT_SAMPLE}_sorted.bam"

# Clean up symbolic links after alignment
rm "${BWA_OUTPUT_DIR}/${BWA_FILE_PREFIX}_1.fastq"
rm "${BWA_OUTPUT_DIR}/${BWA_FILE_PREFIX}_2.fastq"

INPUT_BAM_FOR_ARG="${BWA_OUTPUT_DIR}/${CURRENT_SAMPLE}_sorted.bam"

# echo "BWA Alignment finished for sample: ${CURRENT_SAMPLE}"
# echo "----------------------------------------------------"

# --- Step 2: Add Read Groups ---
echo "Starting Add Read Groups for sample: ${CURRENT_SAMPLE}"

ARG_THREADS=5 # As defined in 2_run_arg_hg38_rep1.sh

OUTPUT_BAM_AFTER_ARG="${ARG_OUTPUT_DIR}/${CURRENT_SAMPLE}_sorted_rg.bam"

# Call the Python script directly. Specify output directory using -d
# Assumes AddReadGroups_ris.py is in the current directory or PATH.
python AddReadGroups_ris.py \
    -i "${INPUT_BAM_FOR_ARG}" \
    -s "HapMap_Mixture" \
    -l "${CURRENT_SAMPLE}" \
    -p ILLUMINA \
    -t ${ARG_THREADS} \
    -x \
    -d "${ARG_OUTPUT_DIR}"

echo "Add Read Groups finished for sample: ${CURRENT_SAMPLE}"
echo "----------------------------------------------------"

# --- Step 3: Mark Duplicates ---
echo "Starting Mark Duplicates for sample: ${CURRENT_SAMPLE}"

INPUT_BAM_FOR_MD="${ARG_OUTPUT_DIR}/${CURRENT_SAMPLE}_sorted_rg.bam"
OUTPUT_BAM_AFTER_MD="${MD_OUTPUT_DIR}/${CURRENT_SAMPLE}_marked_duplicates.bam"
METRICS_FILE="${MD_OUTPUT_DIR}/${CURRENT_SAMPLE}_markduplicates_metrics.txt"

# Using GATK MarkDuplicates. Ensure GATK is in your PATH.
gatk MarkDuplicates \
    -I "${INPUT_BAM_FOR_MD}" \
    -O "${OUTPUT_BAM_AFTER_MD}" \
    -M "${METRICS_FILE}" \
    --CREATE_INDEX true \
    --TMP_DIR "${TMPDIR_GLOBAL}"

echo "Mark Duplicates finished for sample: ${CURRENT_SAMPLE}"
echo "----------------------------------------------------"

# --- Step 4: Indel Realignment ---
echo "Starting Indel Realignment for sample: ${CURRENT_SAMPLE}"

INPUT_BAM_FOR_IR="${MD_OUTPUT_DIR}/${CURRENT_SAMPLE}_marked_duplicates.bam"
# The indel realignment tools output files based on the basename of the input BAM.
# We'll run the commands from the IR_OUTPUT_DIR to ensure all output files land there.
# The GATK 3.x JAR path is hardcoded as per your provided indelrealignment.sh
GATK3_JAR="/ref/twlab/software/zxin/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"

# The output names from GATK 3.x based on input base name
IR_OUTPUT_BAM_BASENAME=$(basename "${INPUT_BAM_FOR_IR}" _marked_duplicates.bam) # e.g., LIB027514
IR_TARGET_LIST="${IR_OUTPUT_DIR}/${IR_OUTPUT_BAM_BASENAME}_target.list"
OUTPUT_BAM_AFTER_IR="${IR_OUTPUT_DIR}/${IR_OUTPUT_BAM_BASENAME}_realigned.bam"

# Step 4a: RealignerTargetCreator
java -jar "${GATK3_JAR}" -T RealignerTargetCreator \
	-nt ${SLURM_CPUS_PER_TASK} \
	-known "${MILLS_VCF}" \
	-R  "${GENOME_HG38}" \
	-I "${INPUT_BAM_FOR_IR}" \
	-o "${IR_TARGET_LIST}" \
    -l INFO # Add log level for clearer output

# Step 4b: IndelRealigner
java -jar "${GATK3_JAR}" -T IndelRealigner \
	-R "${GENOME_HG38}" \
	-I "${INPUT_BAM_FOR_IR}" \
	-targetIntervals "${IR_TARGET_LIST}" \
	-known "${MILLS_VCF}" \
	-o "${OUTPUT_BAM_AFTER_IR}" \
    -l INFO # Add log level for clearer output

# Ensure the realigned BAM is indexed by samtools (GATK 3.x IndelRealigner doesn't always index)
samtools index -@ ${SLURM_CPUS_PER_TASK} "${OUTPUT_BAM_AFTER_IR}"

echo "Indel Realignment finished for sample: ${CURRENT_SAMPLE}"
echo "----------------------------------------------------"


# --- Step 5: Base Quality Score Recalibration (BQSR) ---
echo "Starting BQSR for sample: ${CURRENT_SAMPLE}"

INPUT_BAM_FOR_BQSR="${IR_OUTPUT_DIR}/${CURRENT_SAMPLE}_realigned.bam" # Uses realigned BAM
OUTPUT_TABLE_BQSR="${BQSR_OUTPUT_DIR}/${CURRENT_SAMPLE}_BQSR.table"
OUTPUT_BAM_AFTER_BQSR="${BQSR_OUTPUT_DIR}/${CURRENT_SAMPLE}_recalibrated.bam"

# GATK BaseRecalibrator
gatk BaseRecalibrator \
        -R "${GENOME_HG38}" \
        -I "${INPUT_BAM_FOR_BQSR}" \
        --enable-baq \
        --known-sites "${DBSNP_VCF}" \
        --known-sites "${MILLS_VCF}" \
        -O "${OUTPUT_TABLE_BQSR}" \
        --tmp-dir "${TMPDIR_GLOBAL}" && \
# GATK ApplyBQSR
gatk ApplyBQSR \
        -R "${GENOME_HG38}" \
        -I "${INPUT_BAM_FOR_BQSR}" \
        -bqsr "${OUTPUT_TABLE_BQSR}" \
        -O "${OUTPUT_BAM_AFTER_BQSR}" \
        --tmp-dir "${TMPDIR_GLOBAL}"

echo "BQSR finished for sample: ${CURRENT_SAMPLE}"
echo "----------------------------------------------------"


echo "Pipeline finished for sample: ${CURRENT_SAMPLE}"
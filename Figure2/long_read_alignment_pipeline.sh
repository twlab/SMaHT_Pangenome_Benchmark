#!/bin/bash
#SBATCH --job-name=pbmm2_alignment
#SBATCH --partition=general
#SBATCH --array=0-2 # Adjust based on the number of .fastq.gz files (3 in your example: 0, 1, 2)
#SBATCH --cpus-per-task=20
#SBATCH --time=6-3:00:00
#SBATCH --output=./slurm_logs/pbmm2_%A_%a.log
#SBATCH --error=./slurm_logs/pbmm2_%A_%a.err
#SBATCH --mem=50G

# Set PATH for pbmm2
export PATH="/ref/twlab/software/zxin/pbmm2:$PATH"

# Reference Genome Path
GENOME_HG38="/ref/twlab/data/zxin/genome/new_graph_dual_ref/mgi_non_HG005/ris/personalization/HG005_short/ragtag/mgi_non_HG005_ragtag_hg005_1_sorted_normalized.fasta"

# Base output directory for alignments
BASE_OUTPUT_DIR="/scratch/twlab/qfu/real_mixture/mgi_non_HG005_per_1_long/pbmm2"
mkdir -p "${BASE_OUTPUT_DIR}/slurm_logs" # Create log directory upfront

# Temporary directory setup
custom_tmpdir=/scratch/twlab/qfu/tmp
tmpdir=$(mktemp -d "${custom_tmpdir}/tmp.XXXXXXX")
export TMPDIR="${tmpdir}"

function cleanup() {
    rm -rf "${tmpdir}"
}
trap cleanup EXIT # Ensure cleanup function runs on script exit

# Dynamically get fastq.gz files from your input data directory
INPUT_DATA_DIR="/ref/twlab/data/zxin/real_mixture_SMaHT/data/long"
declare -a FILES=($(ls "${INPUT_DATA_DIR}"/*.fastq.gz)) # Get all .fastq.gz files

# Get SLURM array task ID
ARRAY_ID=${SLURM_ARRAY_TASK_ID}
CURRENT_FASTQ_PATH="${FILES[${ARRAY_ID}]}" # Full path to the current fastq file
CURRENT_FASTQ_FILENAME=$(basename "${CURRENT_FASTQ_PATH}") # Just the filename

# Extract sample name for output. Assuming filename like m84081_240301_s2.hifi_reads.bc2057.fastq.gz
# We'll use the full filename without .fastq.gz as the output prefix.
# Example: m84081_240301_s2.hifi_reads.bc2057
OUTPUT_PREFIX=$(echo "${CURRENT_FASTQ_FILENAME}" | sed 's/\.fastq\.gz$//')

# Define output directory for this specific sample
SAMPLE_OUTPUT_DIR="${BASE_OUTPUT_DIR}/${OUTPUT_PREFIX}"
mkdir -p "${SAMPLE_OUTPUT_DIR}"

echo "Processing file: ${CURRENT_FASTQ_FILENAME}"
echo "Output will be in: ${SAMPLE_OUTPUT_DIR}"
echo "Temporary directory for this job: ${TMPDIR}"
echo "----------------------------------------------------"

# --- PBMM2 Alignment ---
echo "Starting pbmm2 alignment for ${CURRENT_FASTQ_FILENAME}"

PBMM2_THREADS=${SLURM_CPUS_PER_TASK} # Use allocated CPUs for pbmm2

# Create a symbolic link to the fastq.gz file in the output directory for pbmm2 to access
# This ensures all outputs are contained within SAMPLE_OUTPUT_DIR
(cd "${SAMPLE_OUTPUT_DIR}" && ln -s "${CURRENT_FASTQ_PATH}" "${CURRENT_FASTQ_FILENAME}")

# Run pbmm2 alignment
# Output BAM will be named based on the fastq filename prefix
pbmm2 align \
    -j "${PBMM2_THREADS}" \
    --sort \
    --strip \
    --unmapped \
    "${GENOME_HG38}" \
    "${SAMPLE_OUTPUT_DIR}/${CURRENT_FASTQ_FILENAME}" \
    "${SAMPLE_OUTPUT_DIR}/${OUTPUT_PREFIX}_sorted.bam"

# Clean up symbolic link
rm "${SAMPLE_OUTPUT_DIR}/${CURRENT_FASTQ_FILENAME}"

echo "pbmm2 alignment finished for ${CURRENT_FASTQ_FILENAME}"
echo "----------------------------------------------------"

echo "Pipeline finished for ${CURRENT_FASTQ_FILENAME}"
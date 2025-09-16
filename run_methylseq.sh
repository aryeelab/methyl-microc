#!/bin/bash

# Script to run nf-core/methylseq pipeline
# Usage: bash run_methylseq.sh [nextflow_args...]
# Example: bash run_methylseq.sh --input tests/samplesheet.csv --fasta references/chr22/chr22.fa

set -e  # Exit on any error

echo "Setting up environment for nf-core/methylseq pipeline..."

# Set Java environment
export JAVA_HOME="/homes9/martin/java/jdk-24.0.2"
export PATH="$JAVA_HOME/bin:$PATH"

# Add Nextflow to PATH
export PATH="$PWD:$PATH"

# Activate conda environment
source /homes9/martin/miniforge3/etc/profile.d/conda.sh
conda activate methyl-microc

# Set Singularity cache directory (if using containers)
export NXF_SINGULARITY_CACHEDIR="/aryeelab/singularity"

# Verify installations
echo "Java version:"
java -version

echo -e "\nNextflow version:"
./nextflow -version

echo -e "\nnf-core version:"
nf-core --version

# Parse command line arguments to find input file
INPUT_FILE="samplesheet.csv"  # default
ARGS=()
i=1
while [[ $i -le $# ]]; do
    arg="${!i}"
    if [[ "$arg" == "--input" ]]; then
        ((i++))
        if [[ $i -le $# ]]; then
            INPUT_FILE="${!i}"
        fi
    elif [[ "$arg" == --input=* ]]; then
        INPUT_FILE="${arg#*=}"
    fi
    ((i++))
done

# Store all arguments for passing to nextflow
ARGS=("$@")

# Check if input files exist
echo -e "\nChecking input files..."
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "ERROR: Input file $INPUT_FILE not found!"
    exit 1
fi

if [[ ! -f "methylseq.config" ]]; then
    echo "ERROR: methylseq.config not found!"
    exit 1
fi

# Verify FASTQ files exist
echo "Verifying FASTQ files..."
while IFS=, read -r sample fastq1 fastq2; do
    if [[ "$sample" != "sample" ]]; then  # Skip header
        if [[ ! -f "$fastq1" ]]; then
            echo "ERROR: $fastq1 not found!"
            exit 1
        fi
        if [[ ! -f "$fastq2" ]]; then
            echo "ERROR: $fastq2 not found!"
            exit 1
        fi
        echo "✓ Found: $fastq1"
        echo "✓ Found: $fastq2"
    fi
done < "$INPUT_FILE"

echo -e "\nAll input files verified successfully!"

# Create results directory
mkdir -p results

echo -e "\nStarting nf-core/methylseq pipeline..."
echo "Pipeline: nf-core-methylseq_4.1.0/4_1_0"
echo "Input: $INPUT_FILE"
echo "Output: results/"

# Run the pipeline with all provided arguments and local genome files
./nextflow run nf-core-methylseq_4.1.0/4_1_0 \
    -c methylseq.config \
    -resume \
    --fasta_index "$PWD/references/chr22/chr22.fa.fai" \
    --bwameth_index "$PWD/references/bwameth_index/" \
    "${ARGS[@]}"

echo -e "\nPipeline execution completed!"
echo "Results are available in the 'results' directory"
echo "Check the following files for pipeline reports:"
echo "  - results/pipeline_info/execution_report.html"
echo "  - results/pipeline_info/execution_timeline.html"
echo "  - results/pipeline_info/execution_trace.txt"
echo "  - results/multiqc/multiqc_report.html"

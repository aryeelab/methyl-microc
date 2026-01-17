#!/bin/bash

# Script to run nf-core/methylseq pipeline
# Usage: bash run_methylseq.sh [nextflow_args...]
# Example: bash run_methylseq.sh --input tests/samplesheet.csv --fasta references/chr22/chr22.fa

set -e  # Exit on any error

echo "Setting up environment for nf-core/methylseq pipeline..."

# --- Java setup (portable) ---
# The README requires Java 17+ for Nextflow.
# Do not hard-code JAVA_HOME to a machine-specific path; instead, prefer the user's current setup.
if [[ -z "${JAVA_HOME:-}" ]]; then
    if command -v /usr/libexec/java_home >/dev/null 2>&1; then
        # macOS helper: prefer Java 17+ if available
        JAVA_HOME="$(/usr/libexec/java_home -v 17 2>/dev/null || /usr/libexec/java_home 2>/dev/null || true)"
        export JAVA_HOME
    fi
fi
if [[ -n "${JAVA_HOME:-}" ]]; then
    export PATH="$JAVA_HOME/bin:$PATH"
fi

# --- Conda setup (portable) ---
# If the user already activated the env, do not require conda.sh.
if [[ "${CONDA_DEFAULT_ENV:-}" != "methyl-microc" ]]; then
    if ! command -v conda >/dev/null 2>&1; then
        echo "ERROR: conda not found on PATH. Please install Miniconda/Miniforge and/or run 'conda init'." >&2
        exit 1
    fi

    CONDA_BASE="$(conda info --base 2>/dev/null || true)"
    if [[ -z "$CONDA_BASE" || ! -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
        echo "ERROR: Could not locate conda.sh (conda base='$CONDA_BASE')." >&2
        echo "       Please run: conda init <your-shell>  (or activate methyl-microc before running this script)." >&2
        exit 1
    fi

    # shellcheck disable=SC1090
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate methyl-microc
fi

# Set Singularity cache directory (if using containers). For conda profile this is typically unused.
if [[ -d "/aryeelab/singularity" ]]; then
    export NXF_SINGULARITY_CACHEDIR="/aryeelab/singularity"
fi

# Verify installations
echo "Java version:"
java -version

echo -e "\nNextflow version:"
nextflow -version

echo -e "\nnf-core version:"
nf-core --version

# Parse command line arguments to find input file
INPUT_FILE="samplesheet.csv"  # default
OUTDIR="results"              # default
FASTA_FILE=""                 # optional
USER_PROVIDED_FASTA_INDEX=0
USER_PROVIDED_BWAMETH_INDEX=0
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
    elif [[ "$arg" == "--outdir" ]]; then
        ((i++))
        if [[ $i -le $# ]]; then
            OUTDIR="${!i}"
        fi
    elif [[ "$arg" == --outdir=* ]]; then
        OUTDIR="${arg#*=}"
    elif [[ "$arg" == "--fasta" ]]; then
        ((i++))
        if [[ $i -le $# ]]; then
            FASTA_FILE="${!i}"
        fi
    elif [[ "$arg" == --fasta=* ]]; then
        FASTA_FILE="${arg#*=}"
    elif [[ "$arg" == "--fasta_index" || "$arg" == --fasta_index=* ]]; then
        USER_PROVIDED_FASTA_INDEX=1
    elif [[ "$arg" == "--bwameth_index" || "$arg" == --bwameth_index=* ]]; then
        USER_PROVIDED_BWAMETH_INDEX=1
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
mkdir -p "$OUTDIR"

echo -e "\nStarting nf-core/methylseq pipeline..."
echo "Pipeline: nf-core/methylseq (revision 4.1.0)"
echo "Input: $INPUT_FILE"
echo "Output: $OUTDIR"

# Add optional index arguments only if the user didn't specify them.
EXTRA_ARGS=()
if [[ $USER_PROVIDED_FASTA_INDEX -eq 0 && -n "$FASTA_FILE" && -f "${FASTA_FILE}.fai" ]]; then
    EXTRA_ARGS+=("--fasta_index" "${FASTA_FILE}.fai")
fi

# IMPORTANT: Do not force --bwameth_index unless you actually have a prebuilt index directory.
# If it's missing, the pipeline will create an index automatically from --fasta.
if [[ $USER_PROVIDED_BWAMETH_INDEX -eq 0 && -d "$PWD/references/bwameth_index" ]]; then
    EXTRA_ARGS+=("--bwameth_index" "$PWD/references/bwameth_index")
fi

# Run the pipeline with all provided arguments and local genome files
nextflow run nf-core/methylseq \
    -r 4.1.0 \
    -c methylseq.config \
    -resume \
    "${EXTRA_ARGS[@]}" \
    "${ARGS[@]}"

echo -e "\nPipeline execution completed!"
echo "Results are available in the '$OUTDIR' directory"
echo "Check the following files for pipeline reports:"
echo "  - $OUTDIR/pipeline_info/execution_report.html"
echo "  - $OUTDIR/pipeline_info/execution_timeline.html"
echo "  - $OUTDIR/pipeline_info/execution_trace.txt"
echo "  - $OUTDIR/multiqc/multiqc_report.html"

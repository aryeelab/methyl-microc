#!/bin/bash

set -e

echo "Setting up environment for nf-core/methylseq pipeline..."

if [[ -z "${JAVA_HOME:-}" ]]; then
    if command -v /usr/libexec/java_home >/dev/null 2>&1; then
        JAVA_HOME="$(/usr/libexec/java_home -v 17 2>/dev/null || /usr/libexec/java_home 2>/dev/null || true)"
        export JAVA_HOME
    fi
fi
if [[ -n "${JAVA_HOME:-}" ]]; then
    export PATH="$JAVA_HOME/bin:$PATH"
fi

if [[ "${CONDA_DEFAULT_ENV:-}" != "methyl-microc" ]]; then
    if ! command -v conda >/dev/null 2>&1; then
        echo "ERROR: conda not found on PATH." >&2
        exit 1
    fi

    CONDA_BASE="$(conda info --base 2>/dev/null || true)"
    if [[ -z "$CONDA_BASE" || ! -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
        echo "ERROR: Could not locate conda.sh (conda base='$CONDA_BASE')." >&2
        exit 1
    fi

    # shellcheck disable=SC1090
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate methyl-microc
fi

if [[ -d "/aryeelab/singularity" ]]; then
    export NXF_SINGULARITY_CACHEDIR="/aryeelab/singularity"
fi

echo "Java version:"
java -version

echo -e "\nNextflow version:"
nextflow -version

echo -e "\nnf-core version:"
nf-core --version

INPUT_FILE="samplesheet.csv"
OUTDIR="results"
FASTA_FILE=""
USER_PROVIDED_FASTA_INDEX=0
USER_PROVIDED_BWAMETH_INDEX=0

i=1
while [[ $i -le $# ]]; do
    arg="${!i}"
    if [[ "$arg" == "--input" ]]; then
        ((i++))
        if [[ $i -le $# ]]; then INPUT_FILE="${!i}"; fi
    elif [[ "$arg" == --input=* ]]; then
        INPUT_FILE="${arg#*=}"
    elif [[ "$arg" == "--outdir" ]]; then
        ((i++))
        if [[ $i -le $# ]]; then OUTDIR="${!i}"; fi
    elif [[ "$arg" == --outdir=* ]]; then
        OUTDIR="${arg#*=}"
    elif [[ "$arg" == "--fasta" ]]; then
        ((i++))
        if [[ $i -le $# ]]; then FASTA_FILE="${!i}"; fi
    elif [[ "$arg" == --fasta=* ]]; then
        FASTA_FILE="${arg#*=}"
    elif [[ "$arg" == "--fasta_index" || "$arg" == --fasta_index=* ]]; then
        USER_PROVIDED_FASTA_INDEX=1
    elif [[ "$arg" == "--bwameth_index" || "$arg" == --bwameth_index=* ]]; then
        USER_PROVIDED_BWAMETH_INDEX=1
    fi
    ((i++))
done

ARGS=("$@")

echo -e "\nChecking input files..."
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "ERROR: Input file $INPUT_FILE not found!" >&2
    exit 1
fi

if [[ ! -f "methylseq.config" ]]; then
    echo "ERROR: methylseq.config not found!" >&2
    exit 1
fi

echo "Verifying FASTQ files..."
while IFS=, read -r sample fastq1 fastq2; do
    if [[ "$sample" != "sample" ]]; then
        if [[ ! -f "$fastq1" ]]; then echo "ERROR: $fastq1 not found!" >&2; exit 1; fi
        if [[ ! -f "$fastq2" ]]; then echo "ERROR: $fastq2 not found!" >&2; exit 1; fi
        echo "Found: $fastq1"
        echo "Found: $fastq2"
    fi
done < "$INPUT_FILE"

mkdir -p "$OUTDIR"

echo -e "\nStarting nf-core/methylseq pipeline..."
echo "Pipeline: nf-core/methylseq (revision 4.1.0)"
echo "Input: $INPUT_FILE"
echo "Output: $OUTDIR"

EXTRA_ARGS=()
if [[ $USER_PROVIDED_FASTA_INDEX -eq 0 && -n "$FASTA_FILE" && -f "${FASTA_FILE}.fai" ]]; then
    EXTRA_ARGS+=("--fasta_index" "${FASTA_FILE}.fai")
fi

if [[ $USER_PROVIDED_BWAMETH_INDEX -eq 0 && -d "$PWD/references/bwameth_index" ]]; then
    EXTRA_ARGS+=("--bwameth_index" "$PWD/references/bwameth_index")
fi

nextflow run nf-core/methylseq \
    -r 4.1.0 \
    -c methylseq.config \
    --skip_fastqc \
    -resume \
    "${EXTRA_ARGS[@]}" \
    "${ARGS[@]}"

echo -e "\nPipeline execution completed!"

#!/bin/bash

set -e

echo "Setting up environment for nf-core/methylseq pipeline..."

# --- Java ---
if [[ -n "${JAVA_HOME:-}" ]]; then
    export PATH="$JAVA_HOME/bin:$PATH"
fi

# --- Conda check ---
if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "ERROR: No conda environment is active." >&2
    exit 1
fi

unset JAVA_CMD

# Use Java from the currently activated conda environment if available
if [[ -d "${CONDA_PREFIX}/lib/jvm" ]]; then
    export JAVA_HOME="${CONDA_PREFIX}/lib/jvm"
    export PATH="${JAVA_HOME}/bin:${PATH}"
fi

USE_PREBUILT="${USE_PREBUILT_METHYLSEQ_ENV:-false}"

if [[ "$USE_PREBUILT" == "true" ]]; then
    echo "[INFO] Using prebuilt environment (no nested conda)"
    export NXF_CONDA_ENABLED=false
    export CONDA_NO_PLUGINS=true
fi

echo "Java version:"
java -version

echo -e "\nNextflow version:"
nextflow -version

# --- parse args ---
INPUT_FILE="samplesheet.csv"
OUTDIR="results"
FASTA_FILE=""

ARGS=("$@")

for ((i=1; i<=$#; i++)); do
    arg="${!i}"
    if [[ "$arg" == "--input" ]]; then
        ((i++)); INPUT_FILE="${!i}"
    elif [[ "$arg" == --input=* ]]; then
        INPUT_FILE="${arg#*=}"
    elif [[ "$arg" == "--outdir" ]]; then
        ((i++)); OUTDIR="${!i}"
    elif [[ "$arg" == --outdir=* ]]; then
        OUTDIR="${arg#*=}"
    elif [[ "$arg" == "--fasta" ]]; then
        ((i++)); FASTA_FILE="${!i}"
    elif [[ "$arg" == --fasta=* ]]; then
        FASTA_FILE="${arg#*=}"
    fi
done

# --- checks ---
[[ -f "$INPUT_FILE" ]] || { echo "Missing $INPUT_FILE"; exit 1; }
[[ -f "methylseq.config" ]] || { echo "Missing methylseq.config"; exit 1; }

mkdir -p "$OUTDIR"

echo -e "\nRunning nf-core/methylseq..."

if [[ "$USE_PREBUILT" == "true" ]]; then

    echo "[INFO] Running nf-core/methylseq with prebuilt environment."

    cat > nfcore_prebuilt.config <<'EOF'
conda.enabled = false
docker.enabled = false
singularity.enabled = false

process {
    executor = 'local'
}
EOF

    nextflow run nf-core/methylseq \
        -r 4.1.0 \
        -c methylseq.config \
        -c nfcore_prebuilt.config \
        --skip_fastqc \
        --max_cpus 4 \
        -resume \
        "${ARGS[@]}"

else

    echo "[INFO] Running nf-core/methylseq with conda profile."

    export CONDA_NO_PLUGINS=true
    export NXF_CONDA_CACHEDIR="${NXF_CONDA_CACHEDIR:-$PWD/nfcore_conda_cache}"

    nextflow run nf-core/methylseq \
        -r 4.1.0 \
        -c methylseq.config \
        -profile conda \
        --skip_fastqc \
        --max_cpus 4 \
        -resume \
        "${ARGS[@]}"

fi

echo "Done."

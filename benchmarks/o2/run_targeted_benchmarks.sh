#!/usr/bin/env bash
set -euo pipefail

module load conda/miniforge3
conda activate methyl-microc

mkdir -p benchmarks/o2/logs benchmarks/o2/sacct results/bench work/bench benchmarks/o2/nextflow-conda-cache benchmarks/o2/outer-conda-cache

RUN_STAMP="${RUN_STAMP:-$(date +%Y%m%d_%H%M%S)}"
BWAMETH_INDEX="${BWAMETH_INDEX:-/home/mja8/projects/methyl-microc/references/chr22_bwameth_index}"

run_nf() {
    local main_script="$1"
    local label="$2"
    shift 2

    local run_id="${RUN_STAMP}_${label}"
    local log="benchmarks/o2/logs/${run_id}.log"
    local work_dir="work/bench/${run_id}"
    local out_dir="results/bench/${run_id}"

    echo "### START ${run_id} $(date -Is)" | tee "$log"

    NXF_ANSI_LOG=false nextflow run "$main_script" \
        --input test_input/samplesheet.csv \
        --fasta references/chr22.fa \
        --fai references/chr22.fa.fai \
        --outdir "$out_dir" \
        -profile cluster \
        -c o2_lab.config \
        -c benchmarks/o2/trace.config \
        -c benchmarks/o2/all_short.config \
        -c benchmarks/o2/fixed_inner_conda.config \
        "$@" \
        -work-dir "$work_dir" \
        2>&1 | tee -a "$log"

    echo "### END ${run_id} $(date -Is)" | tee -a "$log"

    if [[ -f "$out_dir/pipeline_info/trace.txt" ]]; then
        cp "$out_dir/pipeline_info/trace.txt" "benchmarks/o2/logs/${run_id}.trace.txt"
        local ids
        ids="$(awk 'NR>1{print $3}' "benchmarks/o2/logs/${run_id}.trace.txt" | paste -sd, -)"
        if [[ -n "$ids" ]]; then
            sacct -j "$ids" \
                --format=JobID,JobName%36,State,Elapsed,Submit,Start,End,AllocCPUS,ReqMem,MaxRSS,MaxVMSize,CPUTimeRAW,TotalCPU \
                -P > "benchmarks/o2/sacct/${run_id}.sacct.tsv" || true
        fi
    fi
}

run_nf main.nf fixed_cache_warmup
run_nf bench_index.nf index_skipmultiqc_cache \
    --bwameth_index "$BWAMETH_INDEX" \
    --skip_inner_multiqc true
run_nf bench_combined.nf combined_index_skipmultiqc_cache \
    -c benchmarks/o2/combined_short.config \
    --bwameth_index "$BWAMETH_INDEX" \
    --skip_inner_multiqc true

echo "Targeted benchmark run stamp: ${RUN_STAMP}"

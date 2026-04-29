#!/usr/bin/env bash
set -euo pipefail

module load conda/miniforge3
conda activate methyl-microc

mkdir -p benchmarks/o2/logs benchmarks/o2/sacct results/bench work/bench

RUN_STAMP="${RUN_STAMP:-$(date +%Y%m%d_%H%M%S)}"

run_outer() {
    local label="$1"
    shift

    local run_id="${RUN_STAMP}_${label}"
    local log="benchmarks/o2/logs/${run_id}.log"
    local work_dir="work/bench/${run_id}"
    local out_dir="results/bench/${run_id}"

    echo "### START ${run_id} $(date -Is)" | tee "$log"

    NXF_ANSI_LOG=false nextflow run main.nf \
        --input test_input/samplesheet.csv \
        --fasta references/chr22.fa \
        --fai references/chr22.fa.fai \
        --outdir "$out_dir" \
        -profile cluster \
        -c o2_lab.config \
        -c benchmarks/o2/trace.config \
        "$@" \
        -work-dir "$work_dir" \
        2>&1 | tee -a "$log"

    echo "### END ${run_id} $(date -Is)" | tee -a "$log"

    grep -E 'submitted process|Completed at|Duration|Succeeded|FAILED|ERROR|WARN' "$log" \
        > "benchmarks/o2/logs/${run_id}.summary.txt" || true

    if [[ -f "$out_dir/pipeline_info/trace.txt" ]]; then
        cp "$out_dir/pipeline_info/trace.txt" "benchmarks/o2/logs/${run_id}.trace.txt"
    fi

    local job_ids
    job_ids="$(grep -o 'jobId: [0-9]*' "$log" | awk '{print $2}' | paste -sd, - || true)"
    if [[ -n "$job_ids" ]]; then
        sacct -j "$job_ids" \
            --format=JobID,JobName%36,State,Elapsed,Submit,Start,End,AllocCPUS,ReqMem,MaxRSS,MaxVMSize,CPUTimeRAW,TotalCPU \
            -P > "benchmarks/o2/sacct/${run_id}.sacct.tsv" || true
    fi
}

submit_io_microbench() {
    local label="${RUN_STAMP}_io_microbench"
    echo "### START ${label} $(date -Is)" | tee "benchmarks/o2/logs/${label}.log"
    local job_id
    job_id="$(sbatch --wait --parsable benchmarks/o2/io_microbench.sbatch)"
    echo "job_id=${job_id}" | tee -a "benchmarks/o2/logs/${label}.log"
    sacct -j "$job_id" \
        --format=JobID,JobName%36,State,Elapsed,Submit,Start,End,AllocCPUS,ReqMem,MaxRSS,MaxVMSize,CPUTimeRAW,TotalCPU \
        -P > "benchmarks/o2/sacct/${label}.sacct.tsv" || true
    echo "### END ${label} $(date -Is)" | tee -a "benchmarks/o2/logs/${label}.log"
}

run_outer current
run_outer all_short -c benchmarks/o2/all_short.config
run_outer all_short_shared_conda -c benchmarks/o2/all_short.config -c benchmarks/o2/shared_inner_conda.config
submit_io_microbench

echo "Benchmark run stamp: ${RUN_STAMP}"
echo "Logs: benchmarks/o2/logs"
echo "Slurm accounting: benchmarks/o2/sacct"

# O2 Benchmark Handoff

Current date: 2026-04-28.

This handoff is for restarting the methyl-Micro-C O2 benchmark from `Martins-Mac-Studio.local`.

## Transfer The Working Tree

The current local repo is:

```bash
/Users/martin/projects/methyl-microc
```

From this laptop, transfer the working directory to the Mac Studio:

```bash
cd /Users/martin/projects
rsync -a --delete --exclude '.git/' methyl-microc/ Martins-Mac-Studio.local:~/projects/methyl-microc/
```

If `rsync` is not available, use `scp`:

```bash
scp -r /Users/martin/projects/methyl-microc Martins-Mac-Studio.local:~/projects/
```

Then start a new Codex session on `Martins-Mac-Studio.local` with:

```bash
cd ~/projects/methyl-microc
```

## O2 Login

Use a persistent SSH session to avoid repeated 2FA prompts:

```bash
ssh -t o2 'cd ~/projects/methyl-microc && exec ${SHELL:-/bin/bash} -l'
```

If `o2` alias is not configured:

```bash
ssh -t mja8@o2.hms.harvard.edu 'cd ~/projects/methyl-microc && exec ${SHELL:-/bin/bash} -l'
```

Avoid many one-off SSH commands because each can trigger 2FA.

## Current O2 State To Know

Two 30M-read whole-genome jobs were running as of 2026-04-28 09:07 EDT:

- Baseline/current-wrapper run: Slurm `38116211`, run id `20260428_054750_hct116_30m_combined_wg_writable_index`.
- Optimized comparison: Slurm `38121803`, run id `20260428_071300_hct116_30m_fasttrim_no_fastqc_16cpu`.

They are useful but too slow for iterative benchmarking. If they are still running and cluster usage matters, cancel them after collecting the partial evidence:

```bash
squeue -u mja8
scancel 38116211 38121803
```

Do not rely on the laptop foreground SSH session to keep running. The optimized run was launched with `nohup`; the baseline was launched in an interactive terminal and may die when the laptop disconnects.

## Existing 30M Data And References

Remote 30M chunk directory:

```bash
/n/data1/dfci/pathonc/johnstone/lab/martin/methyl_microc_30m_chunk
```

Files:

```bash
HCT116_30M_R1.fastq.gz
HCT116_30M_R2.fastq.gz
samplesheet.csv
```

Whole-genome reference:

```bash
FASTA=/n/data1/dfci/pathonc/johnstone/lab/egs/references/controls/GRCh38_with_controls.fa
FAI=/n/data1/dfci/pathonc/johnstone/lab/egs/references/controls/GRCh38_with_controls.fa.fai
```

BWA-meth index:

```bash
ORIG_INDEX=/n/data1/dfci/pathonc/johnstone/lab/egs/references/controls/bwameth_index
WRITABLE_INDEX=/n/data1/dfci/pathonc/johnstone/lab/martin/methyl_microc_30m_chunk/GRCh38_with_controls_bwameth_index_writable
```

The writable wrapper index is necessary because nf-core/methylseq creates a symlink inside the index directory. Passing the original shared index failed with permission denied.

## Findings So Far

The major slowdown is not just generic O2 IO.

Observed on the 30M dataset:

- Baseline `TRIMGALORE` used only `--cores 1` because nf-core computes paired-end Trim Galore cores as `task.cpus - 4`; with `task.cpus = 4`, it becomes 1.
- Baseline Trim Galore also ran embedded FastQC even though `--skip_fastqc` was passed. The nf-core Trim Galore module hard-codes `--fastqc` through `ext.args`.
- Baseline `TRIMGALORE` took `54m31s`, used only `139.2%` CPU, read `88.7 GB`, and wrote `85.3 GB`.
- Optimized Trim Galore with `cpus = 12`, `--cores 8`, and `ext.args = ''` took `6m27s`, used `1162.8%` CPU, read `126 GB`, and wrote `125.7 GB`.
- Baseline BWA-meth alignment used `-t 8`; optimized used `-t 12`.
- BWA-meth scaling was real but not linear: optimized processed about `848k` reads per `286-292s`; baseline processed about `564k` reads per `311-335s`.

Conclusion: a 5M-read test will save time. The 30M runs are still useful evidence, but for comparing configurations the 5M chunk should reduce alignment time roughly 6x and make iteration practical.

## Create A 5M Test Dataset On O2

Create it from the already-created 30M chunk to avoid waiting on source uploads:

```bash
base=/n/data1/dfci/pathonc/johnstone/lab/martin/methyl_microc_5m_chunk
src=/n/data1/dfci/pathonc/johnstone/lab/martin/methyl_microc_30m_chunk
mkdir -p "$base"

module load conda/miniforge3/24.11.3-0

zcat "$src/HCT116_30M_R1.fastq.gz" | head -n 20000000 | gzip -1 > "$base/HCT116_5M_R1.fastq.gz"
zcat "$src/HCT116_30M_R2.fastq.gz" | head -n 20000000 | gzip -1 > "$base/HCT116_5M_R2.fastq.gz"

cat > "$base/samplesheet.csv" <<EOF
sample,fastq_1,fastq_2
HCT116_5M,$base/HCT116_5M_R1.fastq.gz,$base/HCT116_5M_R2.fastq.gz
EOF
```

Validate read counts:

```bash
for f in "$base"/HCT116_5M_R*.fastq.gz; do
  lines=$(zcat "$f" | wc -l)
  echo "$f lines=$lines reads=$((lines / 4))"
done
```

Expected: `20000000` lines and `5000000` reads per mate.

## Restart Benchmark Plan On 5M

Use these local helper files now present in the repo:

- `bench_combined.nf`
- `benchmarks/o2/run_methylseq_index.nf`
- `benchmarks/o2/combined_postprocess.nf`
- `benchmarks/o2/trace.config`
- `benchmarks/o2/production_5m_16cpu.config`
- `benchmarks/o2/methylseq_fasttrim_no_fastqc_16cpu.config`
- `benchmarks/o2/run_methylseq_o2_uncapped.sh`

Baseline 5M run with current wrapper behavior:

```bash
run_id=20260428_5m_baseline_current
base=/n/data1/dfci/pathonc/johnstone/lab/martin/methyl_microc_5m_chunk
FASTA=/n/data1/dfci/pathonc/johnstone/lab/egs/references/controls/GRCh38_with_controls.fa
FAI=/n/data1/dfci/pathonc/johnstone/lab/egs/references/controls/GRCh38_with_controls.fa.fai
WRITABLE_INDEX=/n/data1/dfci/pathonc/johnstone/lab/martin/methyl_microc_30m_chunk/GRCh38_with_controls_bwameth_index_writable

NXF_ANSI_LOG=false nextflow run bench_combined.nf \
  --input "$base/samplesheet.csv" \
  --fasta "$FASTA" \
  --fai "$FAI" \
  --outdir "$base/results_$run_id" \
  --bwameth_index "$WRITABLE_INDEX" \
  --skip_inner_multiqc true \
  -profile cluster \
  -c o2_lab.config \
  -c benchmarks/o2/trace.config \
  -c benchmarks/o2/production_5m_16cpu.config \
  -work-dir "$base/work_$run_id" \
  2>&1 | tee "benchmarks/o2/logs/$run_id.log"
```

Optimized 5M run:

```bash
run_id=20260428_5m_fasttrim_no_fastqc_16cpu
base=/n/data1/dfci/pathonc/johnstone/lab/martin/methyl_microc_5m_chunk
FASTA=/n/data1/dfci/pathonc/johnstone/lab/egs/references/controls/GRCh38_with_controls.fa
FAI=/n/data1/dfci/pathonc/johnstone/lab/egs/references/controls/GRCh38_with_controls.fa.fai
WRITABLE_INDEX=/n/data1/dfci/pathonc/johnstone/lab/martin/methyl_microc_30m_chunk/GRCh38_with_controls_bwameth_index_writable

NXF_ANSI_LOG=false nextflow run bench_combined.nf \
  --input "$base/samplesheet.csv" \
  --fasta "$FASTA" \
  --fai "$FAI" \
  --outdir "$base/results_$run_id" \
  --bwameth_index "$WRITABLE_INDEX" \
  --skip_inner_multiqc true \
  --methylseq_runner "$PWD/benchmarks/o2/run_methylseq_o2_uncapped.sh" \
  --methylseq_config "$PWD/benchmarks/o2/methylseq_fasttrim_no_fastqc_16cpu.config" \
  -profile cluster \
  -c o2_lab.config \
  -c benchmarks/o2/trace.config \
  -c benchmarks/o2/production_5m_16cpu.config \
  -work-dir "$base/work_$run_id" \
  2>&1 | tee "benchmarks/o2/logs/$run_id.log"
```

## Metrics To Collect

For every run:

```bash
run_id=...
base=/n/data1/dfci/pathonc/johnstone/lab/martin/methyl_microc_5m_chunk

cat "benchmarks/o2/logs/$run_id.log"
cat "$base/results_$run_id/pipeline_info/execution_trace"*.txt 2>/dev/null || true
find "$base/work_$run_id" -name .command.trace -print -exec cat {} \;
find "$base/work_$run_id" -path '*/results/pipeline_info/execution_trace*.txt' -print -exec cat {} \;

sacct -u mja8 --starttime 2026-04-28 \
  --format=JobID,JobName,State,Elapsed,Submit,Start,End,AllocCPUS,ReqMem,MaxRSS,MaxVMSize,CPUTimeRAW,TotalCPU \
  > "benchmarks/o2/sacct/${run_id}.sacct.tsv"
```

Primary comparison points:

- Trim Galore duration, CPU percent, `rchar`, `wchar`.
- Whether `TRIMGALORE` command contains `--fastqc`.
- Whether `TRIMGALORE` command uses `--cores 1` or `--cores 8`.
- BWA-meth `-t` value and per-block `[M::mem_process_seqs]` throughput.
- Total wall time, queue wait time, and postprocess duration.

## Next Optimization Experiments

After the 5M baseline and optimized run:

1. Try local node `/tmp` scratch for Trim Galore and BWA-meth work directories.
2. Try BWA-meth `cpus = 16` or `20` if O2 `short` queue availability is acceptable.
3. Keep per-chunk MultiQC disabled and run final MultiQC once per sample.
4. Keep using the writable prebuilt BWA-meth index wrapper.
5. Keep combined `POSTPROCESS_PAIRS` to avoid second-scale Slurm jobs.

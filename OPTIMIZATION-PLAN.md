# O2 SLURM Optimization Plan

## Goal

Determine why the methyl-Micro-C pipeline is slower on O2 than on a laptop, then produce an evidence-based SLURM configuration for production runs.

Production assumptions:

- One sample per run.
- Whole-genome data.
- 500M-1000M reads per sample.
- FASTQs are split into approximately 30M-read chunks.
- O2 `/tmp` is local to each compute node and is not shared across nodes.

## Current Evidence

The completed O2 test run succeeded, but most of the wall time was not useful compute.

- Outer workflow completed in about 27 minutes.
- `RUN_METHYLSEQ` used about 11 minutes of wall time but only about 6 minutes of CPU.
- `PARSE_PAIRS`, `ANNOTATE_PAIRS`, and `VALIDATE_PAIRS` each used only seconds of compute.
- The outer work directory was on `/home`, reported by Nextflow as NFS.
- The nested nf-core/methylseq run created conda environments inside the task work directory.
- The nested nf-core/methylseq run rebuilt the bwa-meth index instead of using a prebuilt project index.
- The outer trace file was stale because `trace.overwrite` was not enabled.

Initial hypothesis:

1. Queue latency dominates second-scale tasks.
2. Nested Nextflow and per-run conda environment setup add avoidable overhead.
3. Rebuilding or repeatedly reading the whole-genome bwa-meth index is a major production risk.
4. NFS metadata and intermediate-file IO may be slower than node-local storage for high-churn stages.

## Measurement Rules

Every benchmark run must write to unique directories:

- `results/bench/<run_id>`
- `work/bench/<run_id>` or another explicit benchmark work directory

Every benchmark run must capture:

- Nextflow trace, timeline, and report.
- `.command.trace`, `.command.log`, `.command.out`, `.command.err`.
- Slurm `sacct` fields: `JobID`, `JobName`, `State`, `Elapsed`, `Submit`, `Start`, `End`, `AllocCPUS`, `ReqMem`, `MaxRSS`, `MaxVMSize`, `CPUTimeRAW`, `TotalCPU`.
- Queue wait time: `Start - Submit`.
- CPU efficiency: `TotalCPU / (Elapsed * AllocCPUS)`.
- Read and write bytes from `.command.trace`.

Benchmark configs must set:

```groovy
trace.enabled = true
trace.overwrite = true
timeline.enabled = true
report.enabled = true
```

## Benchmark Phases

### 1. Baseline Repeats

Run the current small chr22 test several times with a fixed, unique output directory.

Purpose:

- Confirm reproducibility.
- Quantify queue wait versus task runtime.
- Establish baseline behavior before optimization.

### 2. Partition Test

Compare current partition choices against O2's `short` partition.

Current test jobs are far below 12 hours, so `short` should be preferred for chunk-level work unless measured production chunks exceed that limit.

Variants:

- Current config.
- All outer processes on `short`.
- `RUN_METHYLSEQ` on `short`, postprocessing on `short`.

### 3. Representative 30M-Read Chunk Test

Use one real or representative whole-genome 30M-read chunk.

Variants:

- Current config.
- `short` partition.
- Prebuilt whole-genome bwa-meth index.
- Shared/prebuilt environments.
- MultiQC disabled during chunk processing.

Repeat enough times to separate queue variability from run-time behavior.

### 4. Local `/tmp` Scratch Test

Use node-local `/tmp` only inside an individual SLURM task. Do not use `/tmp` as the global Nextflow `workDir`, because different tasks may run on different nodes.

Each scratch-enabled task should use a pattern like:

```bash
SCRATCH=/tmp/${USER}/${SLURM_JOB_ID}
mkdir -p "$SCRATCH"
cp required inputs and indexes "$SCRATCH"/
cd "$SCRATCH"
run compute-heavy command
cp final outputs back to the Nextflow task work directory or publish directory
rm -rf "$SCRATCH"
```

Variants:

- No scratch: all work on `/home` NFS.
- Local `/tmp` for intermediates only.
- Local `/tmp` with FASTQs staged locally.
- Local `/tmp` with reference/index staged locally.
- Local `/tmp` with FASTQs and reference/index staged locally.

Measure whether copy-in/copy-out cost is outweighed by faster alignment, sorting, deduplication, and methylation extraction.

### 5. Index Reuse Test

Whole-genome production runs must not rebuild the bwa-meth index.

Variants:

- Current behavior.
- Absolute `--bwameth_index` pointing to the prebuilt shared index.
- Prebuilt index copied to node-local `/tmp` per SLURM job.
- Prebuilt index read directly from shared storage.

### 6. Environment Setup Test

The nested nf-core/methylseq run currently creates conda environments inside its task work directory.

Variants:

- Current nested conda behavior.
- Shared `NXF_CONDA_CACHEDIR`.
- Prebuilt conda environments.
- Apptainer/Singularity container if available and faster on O2.

### 7. Pipeline Shape Test

Compare orchestration designs:

- Current outer Nextflow process that launches inner nf-core/methylseq.
- Direct nf-core/methylseq submission per chunk.
- One SLURM job per chunk that runs methylseq plus parse, annotate, and validate.
- Combined postprocessing process instead of separate second-scale Slurm jobs.

### 8. Resource Scaling

For a representative 30M-read chunk, test:

- 4 CPUs.
- 8 CPUs.
- 12 CPUs.
- 16 CPUs.
- 20 CPUs.

Use appropriate memory levels and test with local `/tmp` scratch enabled and disabled.

O2 normal jobs are capped at 20 cores, so the useful scaling range ends at 20 CPUs.

## Expected Optimization Direction

The likely final SLURM strategy, pending benchmark results:

- Use `short` for 30M-read chunks when runtime is below 12 hours.
- Use `medium` only for measured jobs that exceed `short` limits.
- Prebuild and pass the whole-genome bwa-meth index with an absolute path.
- Disable per-chunk FastQC and MultiQC; run final MultiQC once after all chunks complete.
- Use shared or prebuilt environments to avoid per-run conda creation.
- Use local node `/tmp` inside each SLURM job for high-IO intermediates.
- Combine tiny postprocessing steps to avoid repeated queue latency.
- Tune CPU count per chunk based on measured CPU efficiency, likely in the 12-20 CPU range for whole-genome chunks.
- Optimize for throughput across many 30M-read chunks rather than a single monolithic job.

## Benchmark Results

Initial tests were run on O2 using the available chr22 test dataset. Those results are strongest for orchestration, queue, environment, index, and filesystem behavior.

On April 28, 2026, a representative HCT116 whole-genome 5M-read paired-end dataset was created on O2 from the existing 30M chunk:

- Source: `/n/data1/dfci/pathonc/johnstone/lab/martin/methyl_microc_30m_chunk`
- 5M dataset: `/n/data1/dfci/pathonc/johnstone/lab/martin/methyl_microc_5m_chunk`
- Validated counts: `20,000,000` FASTQ lines and `5,000,000` reads per mate.
- Both runs used the writable whole-genome bwa-meth index wrapper.

### 5M Whole-Genome Results

| Run | Run id | RUN_METHYLSEQ Slurm job | Queue wait | RUN_METHYLSEQ elapsed | Inner Trim Galore | Inner BWA-meth | Postprocess |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Baseline/current wrapper | `20260428_093826_5m_baseline_current` | `38131352` | 2m 50s | 2h 16m 39s | 10m 22s | 1h 59m 44s | 17m 34s |
| Optimized fast trim/no FastQC/16 CPU | `20260428_100318_5m_fasttrim_no_fastqc_16cpu` | `38133529` | 24m 50s | 1h 28m 59s | 1m 19s | 1h 20m 41s | 17m 30s |

The optimized run started later and overlapped the baseline BWA-meth job from 10:28:50 to 11:57:49, so BWA timings include concurrent cluster/filesystem load. The result is still useful because both jobs completed successfully on different nodes and the BWA block-level throughput was stable, but the BWA comparison should be repeated one-at-a-time before final CPU selection.

### 5M Wall-Clock Timelines

Legend:

- `Wait`: Slurm queue wait.
- `TG`: Trim Galore.
- `BWA`: BWA-meth alignment.
- `Post-align`: sort, index, stats, mark duplicates, MethylDackel.
- `PP wait`: Slurm queue wait for `POSTPROCESS_PAIRS`.
- `PP`: `POSTPROCESS_PAIRS`.

Baseline/current wrapper, `20260428_093826_5m_baseline_current`:

| Relative time | Clock time EDT | Event | Duration |
| ---: | --- | --- | ---: |
| 0:00:00 | 09:38:32 | Driver started | |
| 0:00:39 | 09:39:11 | `RUN_METHYLSEQ` submitted to Slurm | |
| 0:03:29 | 09:42:01 | `RUN_METHYLSEQ` started | 2m 50s wait |
| 0:04:04 | 09:42:36 | Inner nf-core tasks started | |
| 0:04:04 -> 0:14:28 | 09:42:36 -> 09:53:00 | `TRIMGALORE` | 10m 22s realtime |
| 0:14:28 -> 2:14:13 | 09:53:00 -> 11:52:46 | `BWAMETH_ALIGN` | 1h 59m 44s realtime |
| 2:14:14 -> 2:20:06 | 11:52:46 -> 11:58:38 | post-align methylseq tasks | about 5m 52s |
| 2:20:08 | 11:58:40 | `RUN_METHYLSEQ` Slurm job ended | 2h 16m 39s elapsed |
| 2:20:31 | 11:59:03 | `POSTPROCESS_PAIRS` submitted to Slurm | |
| 2:22:45 | 12:01:17 | `POSTPROCESS_PAIRS` started | 2m 14s wait |
| 2:40:25 | 12:18:57 | `POSTPROCESS_PAIRS` ended | 17m 40s elapsed |
| 2:41:16 | 12:19:48 | Driver ended | |

Optimized fast trim/no FastQC/16 CPU, `20260428_100318_5m_fasttrim_no_fastqc_16cpu`:

| Relative time | Clock time EDT | Event | Duration |
| ---: | --- | --- | ---: |
| 0:00:00 | 10:03:23 | Driver started | |
| 0:00:37 | 10:04:00 | `RUN_METHYLSEQ` submitted to Slurm | |
| 0:25:27 | 10:28:50 | `RUN_METHYLSEQ` started | 24m 50s wait |
| 0:26:25 | 10:29:48 | Inner nf-core tasks started | |
| 0:26:25 -> 0:27:44 | 10:29:48 -> 10:31:10 | `TRIMGALORE` | 1m 19s realtime |
| 0:27:47 -> 1:48:30 | 10:31:10 -> 11:51:53 | `BWAMETH_ALIGN` | 1h 20m 41s realtime |
| 1:48:30 -> 1:54:22 | 11:51:53 -> 11:57:45 | post-align methylseq tasks | about 5m 52s |
| 1:54:26 | 11:57:49 | `RUN_METHYLSEQ` Slurm job ended | 1h 28m 59s elapsed |
| 1:54:45 | 11:58:08 | `POSTPROCESS_PAIRS` submitted to Slurm | |
| 1:57:54 | 12:01:17 | `POSTPROCESS_PAIRS` started | 3m 10s wait |
| 2:15:30 | 12:18:53 | `POSTPROCESS_PAIRS` ended | 17m 36s elapsed |
| 2:16:15 | 12:19:38 | Driver ended | |

Approximate Gantt view, segment lengths rounded to about 5 minutes per character. Very short segments are widened to one character for visibility.

```text
Baseline  09:38 -> 12:19  |WTTBBBBBBBBBBBBBBBBBBBBBBBBDQPPPP|
Optimized 10:03 -> 12:19  |WWWWWTBBBBBBBBBBBBBBBBDQPPPP|

W = RUN_METHYLSEQ Slurm wait
T = Trim Galore
B = BWA-meth alignment
D = post-align methylseq work
Q = POSTPROCESS_PAIRS Slurm wait
P = POSTPROCESS_PAIRS
```

Queue wait was not consistent between the two configurations: the optimized `RUN_METHYLSEQ` job waited 24m 50s, while the baseline waited 2m 50s. For runtime comparisons, use Slurm elapsed/realtime for the task itself; for operational planning, include both because many chunked production runs will accumulate queue latency.

Trim Galore command evidence:

- Baseline command contained `--fastqc` and `--cores 1`.
- Optimized command removed `--fastqc` and used `--cores 8`.
- Baseline Trim Galore trace: 10m 22s realtime, 159.7% CPU, 15.3 GB `rchar`, 14.5 GB `wchar`.
- Optimized Trim Galore trace: 1m 19s realtime, 1066.8% CPU, 21.5 GB `rchar`, 21.2 GB `wchar`.

Conclusion: `--skip_fastqc` alone is not enough, because the nf-core Trim Galore module still hard-codes FastQC via `ext.args`. The optimized inner config must clear `ext.args` for `NFCORE_METHYLSEQ:METHYLSEQ:TRIMGALORE`.

BWA-meth command and batch evidence:

| Run | BWA threads | BWA batches | Validated reads | Sum of BWA batch real time | Throughput |
| --- | ---: | ---: | ---: | ---: | ---: |
| Baseline/current wrapper | `-t 8` | 18 | 9,936,056 | 7,152.199s | 1,389 reads/s |
| Optimized fast trim/no FastQC/16 CPU | `-t 12` | 12 | 9,936,056 | 4,797.180s | 2,071 reads/s |

BWA "batches" here are BWA-MEM internal input chunks reported by `[M::process]` and `[M::mem_process_seqs]`, not Nextflow or Slurm chunks. The command lines did not explicitly set BWA `-K`, so the observed batch sizes came from BWA/bwa-meth defaults: about 80M bp per batch for `-t 8`, and about 120M bp per batch for `-t 12`.

The baseline Slurm `RUN_METHYLSEQ` job used 16 CPUs for 2h 16m 39s and reported 16:35:19 total CPU, about 45.5% CPU efficiency. The optimized job used 16 CPUs for 1h 28m 59s and reported 16:35:31 total CPU, about 69.9% CPU efficiency. MaxRSS was about 15.3 GB baseline and 17.2 GB optimized, so 64 GB remains conservative for the tested 5M chunk.

Postprocessing on the 5M whole-genome BAM was no longer a second-scale task: both runs used 8 CPUs / 32 GB and took about 17.5 minutes elapsed, with about 2.1 GB MaxRSS and about 20 minutes total CPU. That request is memory-heavy and CPU-light for 5M; test 4 CPUs and 8-16 GB on 30M chunks.

### End-to-End Results

| Run | Wall time | Key changes |
| --- | ---: | --- |
| Current config | 29m 23s | `RUN_METHYLSEQ` on `medium`; postprocess split into 3 Slurm jobs |
| All `short` | 22m 50s | all outer tasks on `short` |
| Fixed shared conda warmup | 21m 50s | all `short`; shared inner `NXF_CONDA_CACHEDIR`; first run populated cache |
| Prebuilt index + warm conda cache | 14m 21s | all `short`; prebuilt chr22 bwa-meth index; warm shared inner conda cache |
| Combined postprocess + index/cache | 6m 52s | one postprocess Slurm job instead of parse/annotate/validate jobs |
| Combined postprocess + index/cache + literal skip MultiQC | 5m 02s | same as above, with inner nf-core MultiQC actually skipped |

### Queue and Runtime Observations

The current run spent a large fraction of wall time waiting in Slurm for tasks that executed in seconds.

Current config:

- `RUN_METHYLSEQ`: submitted 22:53:47, started 22:57:07, ended 23:11:12.
- `PARSE_PAIRS`: waited about 4m 07s, ran 9s.
- `ANNOTATE_PAIRS`: waited about 2m 31s, ran 4s.
- `VALIDATE_PAIRS`: waited about 2m 21s, ran 3s.

Best benchmark:

- `RUN_METHYLSEQ`: waited about 1m 22s, ran 42s.
- `POSTPROCESS_PAIRS`: waited about 1m 48s, ran 10s.

Conclusion: for small and moderate chunks, the number of Slurm submissions matters. Combining short postprocessing steps is a major win.

### Inner nf-core/methylseq Observations

In the original runs, nested nf-core/methylseq rebuilt the bwa-meth index and created conda environments under the task work directory.

Measured effects:

- Fixed shared `NXF_CONDA_CACHEDIR` worked. The warmup created inner envs under `benchmarks/o2/nextflow-conda-cache`; later runs logged `mamba found local env`.
- Passing an absolute prebuilt `--bwameth_index` removed the inner `BWAMETH_INDEX` process.
- Literal `--skip_multiqc` removed the inner `MULTIQC` process.
- With index reuse and warm conda cache, `RUN_METHYLSEQ` task IO dropped from roughly 8 GB read / 2 GB written to roughly 0.8 GB read / 0.03-0.06 GB written in the chr22 benchmark.

### Local `/tmp` Microbenchmark

One `short` Slurm job compared node-local `/tmp` on `compute-b-16-194` with project `/home` NFS.

| Test | node-local `/tmp` | `/home` NFS |
| --- | ---: | ---: |
| 2 GiB sequential write | 7.01s | 3.39s |
| 2 GiB sequential read | 0.38s | 3.24s |
| create 5000 small files | 0.16s | 18.70s |
| stat 5000 small files | 22.61s | 24.53s |

Conclusion: local `/tmp` is clearly better for repeated reads and small-file creation. It was not faster for this one sequential write test, so production scratch use should be targeted: use `/tmp` for high-churn intermediates and copied-in indexes/FASTQs, then copy only final outputs back.

## Recommended O2 SLURM Plan

1. Use `short` for 30M-read chunk jobs by default.

   The benchmarked tasks are far below 12 hours. Production 30M chunks should be measured, but `medium` should only be used for chunks proven to exceed `short` limits.

2. Prebuild the whole-genome bwa-meth index once and pass it with an absolute path.

   Do not rely on `$PWD/references/bwameth_index` inside a task work directory. That caused the inner nf-core pipeline to rebuild the index.

3. Use a shared inner Nextflow conda cache or prebuilt envs.

   Set this inside the `RUN_METHYLSEQ` Slurm job:

   ```bash
   export NXF_CONDA_CACHEDIR=/path/to/shared/nextflow-conda-cache
   ```

   The cache must be outside per-run `work/` directories.

4. Disable per-chunk FastQC and MultiQC.

   Pass a literal `--skip_multiqc` to the inner nf-core/methylseq run. Run MultiQC once after all chunks complete.

   For Trim Galore, `--skip_fastqc` is not sufficient in this workflow because nf-core/methylseq still injected `--fastqc` through module `ext.args`. Use an inner methylseq config equivalent to `benchmarks/o2/methylseq_fasttrim_no_fastqc_16cpu.config`: set `NFCORE_METHYLSEQ:METHYLSEQ:TRIMGALORE` to 12 CPUs, clear `ext.args`, and let Trim Galore use `--cores 8`.

5. Combine postprocessing.

   Replace separate `PARSE_PAIRS`, `ANNOTATE_PAIRS`, and `VALIDATE_PAIRS` Slurm tasks with one `POSTPROCESS_PAIRS` process where possible. The `pairtools` prebuilt env has the Python libraries needed by the annotation scripts in this test.

6. Use node-local `/tmp` inside each Slurm task, not as global Nextflow `workDir`.

   Pattern:

   ```bash
   SCRATCH=/tmp/${USER}/${SLURM_JOB_ID}
   mkdir -p "$SCRATCH"
   cp required inputs/indexes "$SCRATCH"/
   cd "$SCRATCH"
   run compute-heavy command
   cp final outputs back
   rm -rf "$SCRATCH"
   ```

   This is safe because each task uses only its own node-local `/tmp` and returns final outputs to shared storage.

7. Reduce over-requesting on postprocessing.

   The chr22 postprocess task ran in about 10 seconds with 4 CPUs and 8 GB. On the 5M whole-genome chunk, combined postprocess took about 17.5 minutes with 8 CPUs / 32 GB but used only about 2.1 GB MaxRSS and about 20 minutes total CPU. Production chunks will be larger, but 32 GB looks excessive for 5M; test 4 CPUs and 8-16 GB on 30M chunks.

8. Use the 16-CPU wrapper with 12-thread inner BWA-meth as the next whole-genome baseline.

   The 5M optimized run used outer `RUN_METHYLSEQ` `cpus = 16`, inner BWA-meth `-t 12`, and inner Trim Galore `--cores 8`. It reduced `RUN_METHYLSEQ` elapsed time from 2h 16m 39s to 1h 28m 59s on 5M reads. Treat this as the next baseline, not as the final production optimum.

9. Re-test with one real 30M-read whole-genome chunk.

   Before final production deployment, run the best configuration on one real 30M chunk and test 12, 16, and 20 outer CPUs, including one-at-a-time BWA-meth runs to avoid overlap. Track CPU efficiency, elapsed time, BWA block throughput, and IO bytes. Use that to choose the production CPU count.

## Concrete Target Configuration

Use this as the next production-oriented cluster profile baseline:

```groovy
process {
    executor = 'slurm'
    queue = 'short'
    shell = ['/bin/bash', '-euo', 'pipefail']

    withName: RUN_METHYLSEQ {
        queue = 'short'
        cpus = 16
        memory = 64.GB
        time = 12.h
        beforeScript = 'module load conda/miniforge3/24.11.3-0; export NXF_CONDA_CACHEDIR=/path/to/shared/nextflow-conda-cache'
    }

    withName: POSTPROCESS_PAIRS {
        queue = 'short'
        cpus = 4
        memory = 16.GB
        time = 4.h
        conda = '/path/to/prebuilt/pairtools-env'
    }
}

trace {
    enabled = true
    overwrite = true
}

timeline {
    enabled = true
    overwrite = true
}

report {
    enabled = true
    overwrite = true
}
```

For production runs, pass:

```bash
--bwameth_index /absolute/path/to/GRCh38/bwameth_index \
--methylseq_runner /path/to/run_methylseq_o2_uncapped.sh \
--methylseq_config /path/to/methylseq_fasttrim_no_fastqc_16cpu.config \
--skip_multiqc
```

Then run final MultiQC once after chunk aggregation.

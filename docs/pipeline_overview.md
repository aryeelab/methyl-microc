## Pipeline structure

```text
methyl-microc/
├── main.nf
├── nextflow.config
├── o2_lab.config
├── modules/
│   ├── read_qc.nf
│   ├── split_fastq.nf
│   ├── run_methylseq.nf
│   ├── parse_pairs.nf
│   ├── merge_dedup_pairs.nf
│   ├── pair_qc.nf
│   ├── annotate_pairs.nf
│   ├── methylation_qc.nf
│   └── validate_pairs.nf
├── bin/
│   ├── annotate_pairs_methylation.py
│   ├── validate_pairs_methylation.py
│   ├── pair_level_qc.py
│   └── methylation_qc.py
├── envs/
│   ├── methylseq.yml
│   ├── pairtools.yml
│   └── methyl.yml
└── prebuilt_envs/
    ├── methylseq/
    ├── pairtools/
    └── methyl/
```

## Pipeline components overview

### Workflow

```text
READ_QC                         (parallel with SPLIT_FASTQ)
→ SPLIT_FASTQ
→ RUN_METHYLSEQ per chunk      (parallel)
→ PARSE_PAIRS per chunk        (parse + sort only, no dedup)
→ MERGE_DEDUP_PAIRS            (global merge + dedup)
→ PAIR_QC                      (parallel with ANNOTATE_PAIRS)
→ ANNOTATE_PAIRS
→ METHYLATION_QC               (parallel with VALIDATE_PAIRS)
→ VALIDATE_PAIRS
```

### Top level files

* main.nf  
  Defines the overall workflow and connects all pipeline steps.  
  Handles input channels, parameter passing, and process chaining.

* nextflow.config  
  Specifies execution settings, including:
  - resource allocation (CPU, memory, time)
  - executor configuration (local / Slurm)
  - conda environments per process

* o2_lab.config  
  Cluster-specific configuration for O2: prebuilt environment usage


### modules/ 

* read_qc.nf  
  Generates read-level QC reports using FastQC and MultiQC.

* split_fastq.nf  
  Splits paired-end FASTQ files into chunks based on # reads_per_chunk for parallel processing.

* run_methylseq.nf  
  Aligns sequencing reads to the reference genome and generates BAM files.

* parse_pairs.nf  
  Converts BAM files into pairs format to extract chromatin contact information.

* merge_dedup_pairs.nf  
  Merges and deduplicates pairs files generated from multiple FASTQ chunks.

* pair_qc.nf  
  Generates pair-level QC metrics and plots from final deduplicated pairs.

* annotate_pairs.nf  
  Adds per-fragment methylation information to each pair using the reference FASTA.

* methylation_qc.nf  
  Generates methylation bias QC plots from annotated methylation pairs.
  
* validate_pairs.nf  
  Checks that the final annotated pairs file are correctly generated.



### Supporting components

* bin/  
  Contains custom Python scripts for methylation annotation, validation, and QC generation.

* envs/  
  Defines conda environments required for different pipeline steps.

* prebuilt_envs/  
  Contains pre-created environments used to avoid runtime environment conflicts on cluster.

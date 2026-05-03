## Pipeline structure

```text
methyl-microc/
├── main.nf
├── nextflow.config
├── o2_lab.config
├── modules/
│   ├── split_fastq.nf
│   ├── run_methylseq.nf
│   ├── parse_pairs.nf
│   ├── merge_dedup_pairs.nf
│   ├── annotate_pairs.nf
│   └── validate_pairs.nf
├── bin/
│   ├── annotate_pairs_methylation.py
│   └── validate_pairs_methylation.py
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
SPLIT_FASTQ
→ RUN_METHYLSEQ per chunk   (parallel)
→ PARSE_PAIRS per chunk     (parse + sort only, no dedup)
→ MERGE_DEDUP_PAIRS         (global merge + dedup)
→ ANNOTATE_PAIRS
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

* split_fastq.nf  
  Splits paired-end FASTQ files into chunks based on # reads_per_chunk for parallel processing.

* run_methylseq.nf  
  Aligns sequencing reads to the reference genome and generates BAM files.

* parse_pairs.nf  
  Converts BAM files into pairs format to extract chromatin contact information.

* merge_dedup_pairs.nf  
  Merges and deduplicates pairs files generated from multiple FASTQ chunks.

* annotate_pairs.nf  
  Adds per-fragment methylation information to each pair using the reference FASTA.
  
* validate_pairs.nf  
  Checks that the final annotated pairs file are correctly generated.


### Supporting components

* bin/  
  Contains custom Python scripts for methylation annotation and validation.

* envs/  
  Defines conda environments required for different pipeline steps.

* prebuilt_envs/  
  Contains pre-created environments used to avoid runtime environment conflicts on cluster.

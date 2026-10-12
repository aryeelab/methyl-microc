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
│   ├── merge_global_bam.nf
│   ├── parse_pairs.nf
│   ├── merge_dedup_pairs.nf
│   ├── build_matrices.nf
│   ├── pair_qc.nf
│   ├── annotate_pairs.nf
│   ├── methylation_qc.nf
│   ├── methyldackel_extract.nf
│   ├── methyldackel_mbias.nf
│   └── validate_pairs.nf
├── bin/
│   ├── annotate_pairs_methylation.py
│   ├── validate_pairs_methylation.py
│   ├── pair_level_qc.py
│   └── methylation_qc.py
├── envs/
│   ├── methylseq.yml
│   ├── pairtools.yml
│   ├── methyl.yml
│   ├── hictk.yml
│   └── methyldackel.yml
└── prebuilt_envs/
    ├── methylseq/
    ├── pairtools/
    ├── methyl/
    ├── hictk/
    └── methyldackel/
```

## Pipeline components overview

### Workflow

```text
READ_QC                         (parallel with SPLIT_FASTQ)
→ SPLIT_FASTQ
→ RUN_METHYLSEQ per chunk      (parallel)
→ MERGE_GLOBAL_BAM             (sample level merged, globally duplicate marked BAM)
    ├── METHYLDACKEL_EXTRACT
    └── METHYLDACKEL_MBIAS
→ PARSE_PAIRS per chunk        (parse + sort only, no dedup)
→ MERGE_DEDUP_PAIRS            (global merge + dedup)
    ├── BUILD_MATRICES
    ├── PAIR_QC
    └── ANNOTATE_PAIRS
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

* merge_global_bam.nf  
  Merges chunk-level BAM files for each sample, performs global duplicate marking, coordinate sorts the BAM, and generates a BAM index.

* methyldackel_extract.nf  
  Generates MethylDackel intermediate files from the final globally merged and duplicate marked BAM.

* methyldackel_mbias.nf  
  Generates MethylDackel M-bias outputs for CpG, CHG, and CHH contexts.

* parse_pairs.nf  
  Converts BAM files into pairs format to extract chromatin contact information.

* merge_dedup_pairs.nf  
  Merges and deduplicates pairs files generated from multiple FASTQ chunks.

* build_matrices.nf  
  Generates contact matrix intermediates from the final deduplicated pairs using hictk, and creates a temporary pairs file with chromosome headers ordered according to the reference FAI. 

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

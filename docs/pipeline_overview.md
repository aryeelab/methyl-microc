## Pipeline structure

```text
methyl-microc/
├── main.nf
├── nextflow.config
├── modules/
│   ├── run_methylseq.nf
│   ├── parse_pairs.nf
│   ├── annotate_pairs.nf
│   └── validate_pairs.nf
├── bin/
│   ├── annotate_pairs_methylation.py
│   └── validate_pairs_methylation.py
└── envs/
    ├── methylseq.yml
    ├── pairtools.yml
    └── methyl.yml

## Pipeline components overview

* main.nf  
  Defines the overall workflow and connects all pipeline steps.

* nextflow.config  
  Specifies execution settings such as environments and resource configurations.

---

* modules/  

  * run_methylseq.nf  
    Aligns sequencing reads to the reference genome and generates BAM files.

  * parse_pairs.nf  
    Converts BAM files into pairs format to extract chromatin contact information.

  * annotate_pairs.nf  
    Adds per-fragment methylation information to each pair using the reference FASTA.

  * validate_pairs.nf  
    Checks that the final annotated pairs file are correctly generated.

---

* bin/  
  Contains custom scripts for methylation annotation and validation.

* envs/  
  Defines conda environments required for different pipeline steps.

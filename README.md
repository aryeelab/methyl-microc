# Methyl-MicroC Analysis Environment

This repository contains the methyl-microc pipeline using nf-core/methylseq.

## Prerequisites

Before setting up this environment, ensure you have:

- **Conda/Miniconda** available in your shell environment
- **Java 17 or later** (required for Nextflow)

## Quick Setup

### Step 1: Create Conda Environment

```bash
conda create -n methyl-microc -c conda-forge -c bioconda \
  python=3.10 nextflow nf-core -y
```

### Step 2: Activate Environment and Download Pipeline

```bash
conda activate methyl-microc

# Download nf-core/methylseq pipeline
nf-core pipelines download methylseq --revision 4.1.0
```

### Step 3: Install additional tools
```bash
pip install pyBigWig
# Note - we use pip for now to avoid conda conflicts
```

### Step 4: Reference Genome Setup

#### Quick test with chromosome 22 only:
```bash
# For faster testing, download just chromosome 22 (~1.5MB)
mkdir -p references/chr22
cd references/chr22
curl -L -o chr22.fa.gz "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"
gunzip chr22.fa.gz
cd ../..
```

#### Full GRCh38 Reference Genome:
```bash
# Create reference directory
mkdir -p references/GRCh38
cd references/GRCh38

# Download GRCh38 reference genome (245MB compressed, ~800MB uncompressed)
curl -L -C - -o GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz \
  'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz'

# Extract the reference genome
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
cd ../..
```



## Usage

### Running the Pipeline

The pipeline is designed to be run using the provided script:

```bash
# Activate environment
conda activate methyl-microc

# Run the methylseq pipeline on a small test sample
time ./run_methylseq.sh \
    -profile conda \
    --input tests/samplesheet.csv \
    --outdir results \
    --fasta $PWD/references/chr22/chr22.fa \
    --aligner bwameth

# Run the methylseq pipeline on Arsh's HCT116 samples
time ./run_methylseq.sh \
    -profile conda \
    --input samplesheet.csv \
    --outdir results_hct116 \
    --fasta $PWD/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    --aligner bwameth    


# [ADD TO PIPELINE]
# Convert methylation track bedgraphs to bigwig 
samtools faidx references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

cat results_hct116/methyldackel/HCT116_Meth_MicroC.markdup.sorted_CpG.bedGraph | grep -v KI | grep -v GL | grep -v random > tmp.bedGraph
time python bin/bedgraph_to_bigwig.py tmp.bedGraph results_hct116/methyldackel/HCT116_Meth_MicroC.markdup.sorted_CpG.bw references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai

cat results_hct116/methyldackel/HCT116_Meth_MicroC_red_klnw.markdup.sorted_CpG.bedGraph | grep -v KI | grep -v GL | grep -v random > tmp.bedGraph
time python bin/bedgraph_to_bigwig.py tmp.bedGraph results_hct116/methyldackel/HCT116_Meth_MicroC_red_klnw.markdup.sorted_CpG.bw references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai

```


# [ADD TO PIPELINE]
# Parse pairs
```bash
conda create -n pairtools -c conda-forge -c bioconda pairtools -y

conda activate pairtools

CHROM_SIZES="references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai"

BAM="results_hct116/bwameth/deduplicated/HCT116_Meth_MicroC.markdup.sorted.bam"
PAIRS="results_hct116/pairs/HCT116_Meth_MicroC.pairs.gz"
STATS="results_hct116/pairs/HCT116_Meth_MicroC.stats.txt"

#BAM="results_hct116/bwameth/deduplicated/HCT116_Meth_MicroC_red_klnw.markdup.sorted.bam"
#PAIRS="results_hct116/pairs/HCT116_Meth_MicroC_red_klnw.pairs.gz"
#STATS="results_hct116/pairs/HCT116_Meth_MicroC_red_klnw.stats.txt"

mkdir -p results_hct116/pairs
time pairtools parse --min-mapq 30 --walks-policy 5unique \
        --max-inter-align-gap 30 --add-columns pos5,pos3 \
        --drop-sam --nproc-in 8 --nproc-out 8 --chroms-path ${CHROM_SIZES} \
        $BAM | \
        pairtools sort --nproc 4  | \
        pairtools dedup -o $PAIRS --output-stats $STATS
  

```


# [ADD TO PIPELINE]
# QC
```bash
conda create -n multiqc-pairtools pip
conda activate multiqc-pairtools 
pip install git+https://github.com/open2c/MultiQC.git


conda activate multiqc-pairtools 
multiqc -f -o results_hct116/multiqc results_hct116/pairs 


cd stats


```


# [ADD TO PIPELINE]
# Create Cooler and Juicebox hic

```bash
#PAIRS="results_hct116/pairs/HCT116_Meth_MicroC_red_klnw.pairs.gz"
#HIC="results_hct116/hic/HCT116_Meth_MicroC_red_klnw.hic"

PAIRS="results_hct116/pairs/HCT116_Meth_MicroC.pairs.gz"
HIC="results_hct116/hic/HCT116_Meth_MicroC.hic"
COOLER="results_hct116/hic/HCT116_Meth_MicroC.cool"

conda create -n hictk -c conda-forge -c bioconda hictk
conda activate hictk

time hictk load --format 4dn --bin-size 100kbp $PAIRS $HIC
time hictk load --format 4dn --bin-size 100kbp $PAIRS $COOLER


```



### Input Format

Create a samplesheet.csv file with your samples:

```csv
sample,fastq_1,fastq_2
HCT116_Meth_MicroC,data/20250612_250528_HCT116_Meth_MicroC_AK13242_S27_L008_R1_001.fastq.gz,data/20250612_250528_HCT116_Meth_MicroC_AK13242_S27_L008_R2_001.fastq.gz
HCT116_Meth_MicroC_red_klnw,data/20250612_250528_HCT116_Meth_MicroC_red_klnw_AK13242_S28_L008_R1_001.fastq.gz,data/20250612_250528_HCT116_Meth_MicroC_red_klnw_AK13242_S28_L008_R2_001.fastq.gz
```




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
  python=3.10.19 nextflow=25.10.2 nf-core=3.5.1 -y
```

### Step 2: Activate Environment and Download Pipeline

```bash
conda activate methyl-microc

# Download nf-core/methylseq pipeline
nf-core pipelines download methylseq --revision 4.1.0
```

### Step 3: Optional: Install tools for manual steps
```bash
# This is separate from the per-process conda envs created by Nextflow when you run with `-profile conda`.

conda install -c conda-forge -c bioconda \
  bwa=0.7.19 \
  bwameth=0.2.9 \
  samtools=1.23 \
  bedtools=2.31.1 \
  picard=3.4.0 \
  qualimap=2.3 \
  methyldackel=0.6.1 \
  hisat2=2.2.1 \
  minimap2=2.30 \
  trim-galore=0.6.10 \
  fastqc=0.12.1 \
  cutadapt=5.1 \
  pigz=2.8 \
  pbzip2=1.1.13 \
  multiqc=1.30 \
  pairtools=1.1.3 \
  pairix=0.3.9 \
  pysam=0.23.0 \
  -y

# Python dependency for bigWig conversion scripts
python -m pip install pyBigWig==0.3.25
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
    --input data/20250612_hct116/samplesheet_20250612_hct116.csv \
    --outdir results/20250612_hct116 \
    --fasta $PWD/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    --aligner bwameth    


# Convert methylation track bedgraphs to bigwig 
# [ADD TO PIPELINE]
samtools faidx references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

cat results/20250612_hct116/methyldackel/HCT116_Meth_MicroC.markdup.sorted_CpG.bedGraph | grep -v KI | grep -v GL | grep -v random > tmp.bedGraph
time python bin/bedgraph_to_bigwig.py tmp.bedGraph results/20250612_hct116/methyldackel/HCT116_Meth_MicroC.markdup.sorted_CpG.bw references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai

cat results/20250612_hct116/methyldackel/HCT116_Meth_MicroC_red_klnw.markdup.sorted_CpG.bedGraph | grep -v KI | grep -v GL | grep -v random > tmp.bedGraph
time python bin/bedgraph_to_bigwig.py tmp.bedGraph results/20250612_hct116/methyldackel/HCT116_Meth_MicroC_red_klnw.markdup.sorted_CpG.bw references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai

```



# Parse pairs
```bash


conda create -n pairtools -c conda-forge -c bioconda \
  python=3.10.19 \
  pairtools=1.1.3 \
  pysam=0.23.3 \
  -y

# [ADD TO PIPELINE]

# NOTE (macOS Apple Silicon / osx-arm64): `pairtools parse --add-columns ... seq` can crash
# with a Bus error for some newer pysam builds. A known-good pin is `pysam=0.23.0`.

CHROM_SIZES="references/chr22/chr22.fa.fai"
BAM="results/bwameth/deduplicated/test_sample.markdup.sorted.bam"
PAIRSAM="results/pairs/test_sample.pairsam.gz"
STATS="results/pairs/test_sample.stats.txt" 
mkdir -p results/pairs

CHROM_SIZES="references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai"
BAM="results/20250612_hct116/bwameth/deduplicated/HCT116_Meth_MicroC.markdup.sorted.bam"
PAIRSAM="results/20250612_hct116/pairs/HCT116_Meth_MicroC.pairsam.gz"
STATS="results/20250612_hct116/pairs/HCT116_Meth_MicroC.stats.txt"

#BAM="results/20250612_hct116/bwameth/deduplicated/HCT116_Meth_MicroC_red_klnw.markdup.sorted.bam"
#PAIRSAM="results/20250612_hct116/pairs/HCT116_Meth_MicroC_red_klnw.pairsam.gz"
#STATS="results/20250612_hct116/pairs/HCT116_Meth_MicroC_red_klnw.stats.txt"

mkdir -p results/20250612_hct116/pairs

conda activate pairtools
time pairtools parse --min-mapq 30 --walks-policy 5unique \
        --max-inter-align-gap 30 --drop-sam --add-columns pos5,pos3,cigar,seq \
        --nproc-in 8 --nproc-out 8 --chroms-path ${CHROM_SIZES} \
        $BAM | \
        pairtools sort --nproc 4  | \
        pairtools dedup -o $PAIRSAM --output-stats $STATS

# Append per-fragment methylation strings (meth1/meth2) using the reference FASTA
# (requires the FASTA to be indexed: samtools faidx reference.fa)
# NOTE: the annotator requires `cigar{1,2}` and `seq{1,2}` columns.
#
time python bin/annotate_pairsam_methylation.py \
    --fasta references/chr22/chr22.fa \
    --input  $PAIRSAM \
    --output results/pairs/test_sample.meth.pairsam.gz

time python bin/annotate_pairsam_methylation.py \
    --fasta references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    --input  $PAIRSAM \
    --output results/20250612_hct116/pairs/HCT116_Meth_MicroC.meth.pairsam.gz


# Validate methylation annotation on the first pair
#
# This sanity check demonstrates that meth1/meth2 are correctly generated.
# It:
#   1) extracts the first (non-header) pair from results/pairs/test_sample.meth.pairsam.gz
#   2) fetches the reference sequence spanning pos5..pos3 (inclusive) from the FASTA
#   3) checks that len(meth) == abs(pos3-pos5)+1 and that CpG sites in the reference align
#      to CpG calls in the methylation string.
#
# Methylation encoding (per base of the fragment span):
#   '1' = methylated CpG
#   '0' = unmethylated CpG
#   '-' = not a CpG position in the reference
#   '.' = unclear / no-call
#
# Coordinate convention:
#   pairtools reports pos5 and pos3 for each side. The analyzed fragment span is the inclusive
#   region between these endpoints; its length is abs(pos3-pos5)+1.

PAIRS_METH="results/pairs/test_sample.meth.pairsam.gz"
FASTA="references/chr22/chr22.fa"

python bin/validate_pairsam_methylation.py --pairs "${PAIRS_METH}" --fasta "${FASTA}" --record 1

# Or validate a specific record by readID:
READ_ID="LH00547:129:233G5FLT3:8:1101:30945:1689"
python bin/validate_pairsam_methylation.py --pairs "${PAIRS_METH}" --fasta "${FASTA}" --readID "${READ_ID}"

gunzip -c results/pairs/test_sample.meth.pairsam.gz | head -n 10

# HCT116
#python bin/annotate_pairsam_methylation.py \
#   --fasta references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
#   --input  $PAIRSAM \
#   --output results/20250612_hct116/pairs/HCT116_Meth_MicroC.meth.pairsam.gz
  

```



# QC
```bash
# [ADD TO PIPELINE]
conda activate methyl-microc

# (Already included in Quick Setup Step 3, but repeated here for convenience)
conda install -c conda-forge -c bioconda multiqc=1.30 -y

# If you specifically need the open2c fork of MultiQC (e.g. for extra pairtools-related modules),
# install it via pip instead and skip the conda MultiQC install above:
# python -m pip install git+https://github.com/open2c/MultiQC.git

multiqc -f -o results/20250612_hct116/multiqc results/20250612_hct116/pairs 


cd stats


```



# Create Cooler and Juicebox hic

```bash
# Create conda env for hictk
conda create -n hictk -c conda-forge -c bioconda python=3.10.19 hictk=2.2.0 -y

# [ADD TO PIPELINE]
#PAIRS="results/20250612_hct116/pairs/HCT116_Meth_MicroC_red_klnw.pairs.gz"
#HIC="results/20250612_hct116/hic/HCT116_Meth_MicroC_red_klnw.hic"

PAIRS="results/20250612_hct116/pairs/HCT116_Meth_MicroC.pairs.gz"
HIC="results/20250612_hct116/hic/HCT116_Meth_MicroC.hic"
COOLER="results/20250612_hct116/hic/HCT116_Meth_MicroC.cool"

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




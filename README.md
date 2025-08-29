# Methyl-MicroC Analysis Environment

A complete environment setup for running the nf-core/methylseq pipeline on methylation sequencing data using conda-based local execution with bwa-meth aligner.

## Overview

This repository provides a fully configured environment for analyzing bisulfite sequencing data using the nf-core/methylseq pipeline with bwa-meth aligner. The setup includes:

- Conda environment with Nextflow, nf-core tools, and all bioinformatics dependencies
- Local conda-based execution (no Docker required)
- bwa-meth based methylation analysis workflow
- Pre-configured sample sheets and configuration files
- Test data for validation
- Comprehensive documentation and usage examples

## Prerequisites

Before setting up this environment, ensure you have:

- **macOS** (this setup is optimized for macOS with Apple Silicon)
- **Homebrew** package manager
- **Conda/Miniconda** available in your shell environment
- At least 8GB of RAM and 20GB of free disk space
- **Note**: Docker is NOT required for this conda-based setup, but can be used if desired

## Quick Start

### 1. Activate the Environment

```bash
# Activate the conda environment
source ~/.zshrc  # or ~/.bashrc for bash users
conda activate methyl-microc

# Verify installations
nextflow -version
nf-core --version
bwameth.py --version
trim_galore --version
```

### 2. Run Test Analysis

```bash
# Run the pipeline with test data using conda
nextflow run nf-core/methylseq \
    -profile conda,test \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    --aligner bwameth
```

### 3. Monitor Progress

The pipeline will automatically:
- Download required reference genomes
- Create conda environments for each process
- Process the test data through the complete bwa-meth methylation analysis workflow
- Align bisulfite-treated reads using bwa-meth
- Extract methylation calls using MethylDackel
- Generate comprehensive reports and results

## Environment Setup (Detailed)

### Step 1: Create Conda Environment

```bash
# Create the environment with Python 3.11
conda create -n methyl-microc python=3.11 -y

# Activate the environment
conda activate methyl-microc

# Install required packages
conda install -c bioconda nextflow nf-core -y
```

### Step 2: Install Bioinformatics Tools

```bash
# Install all required bioinformatics tools
conda install -c bioconda bwameth bwa methyldackel trim-galore fastqc multiqc samtools bedtools picard qualimap -y

# Verify tool installations
bwameth.py --version
bwa
MethylDackel --version
```

### Step 3: Download and Index GRCh38 Reference Genome (Optional but Recommended)

Pre-downloading the reference genome eliminates network timeouts and speeds up pipeline runs:

```bash
# Create reference directory
mkdir -p references/GRCh38

# Download GRCh38 reference genome (245MB compressed, ~800MB uncompressed)
cd references/GRCh38
curl -L -C - -o GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz \
  'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz'

# Extract the reference genome
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Create bwa-meth index (this takes 1-2 hours and ~15-20GB disk space)
conda activate methyl-microc
bwameth.py index GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

# Return to project directory
cd ../..
```

**Alternative download sources:**
```bash
# Option A: UCSC Genome Browser
curl -L -o hg38.fa.gz "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"

# Option B: Ensembl
curl -L -o Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
  "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
```

**Quick test with chromosome 22 only:**
```bash
# For faster testing, download just chromosome 22 (~1.5MB)
mkdir -p references/chr22
cd references/chr22
curl -L -o chr22.fa.gz "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"
gunzip chr22.fa.gz
bwameth.py index chr22.fa
cd ../..
```

### Step 3: Download nf-core/methylseq Pipeline

```bash
# Pull the latest version of the pipeline
nextflow pull nf-core/methylseq
```

## Configuration Files

### nextflow.config

The main configuration file that enables conda and sets resource limits:

```groovy
// Enable conda by default
conda.enabled = true

// Conda configuration
conda {
    enabled = true
    useMamba = false
    createTimeout = '1 h'
    cacheDir = "$HOME/.nextflow/conda"
}

// Process configuration
process {
    conda = "${projectDir}/environment.yml"
    cpus = 1
    memory = 4.GB
    time = 2.h
    errorStrategy = 'retry'
    maxRetries = 2
}
```

### samplesheet.csv

Sample sheet format for input data:

```csv
sample,fastq_1,fastq_2
test_sample,tests/tiny_r1.fastq.gz,tests/tiny_r2.fastq.gz
```

For your own data, update the paths to point to your FASTQ files:

```csv
sample,fastq_1,fastq_2
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
```

### conf/test.config

Test-specific configuration with resource limits:

```groovy
params {
    config_profile_name = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    max_cpus = 2
    max_memory = '6.GB'
    max_time = '6.h'

    input = 'samplesheet.csv'
    genome = 'GRCh38'

    aligner = 'bwameth'
    save_reference = true
}
```

## Usage Examples

### Basic Analysis

```bash
# Activate environment
conda activate methyl-microc

# Option A: Run with your own data using remote genome (requires internet)
nextflow run nf-core/methylseq \
    -profile conda \
    --input your_samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    --aligner bwameth \
    --max_cpus 8

# Option B: Run with local reference genome (faster, no internet required)
nextflow run nf-core/methylseq \
    -profile conda \
    --input your_samplesheet.csv \
    --outdir results \
    --fasta $PWD/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    --aligner bwameth \
    --max_cpus 8

# Option C: Use the local reference configuration profile
nextflow run nf-core/methylseq \
    -profile conda,local_reference \
    --input your_samplesheet.csv \
    --outdir results \
    --max_cpus 8
```

### Advanced Options

```bash
# Run with custom parameters using remote genome
nextflow run nf-core/methylseq \
    -profile conda \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    --aligner bwameth \
    --save_reference \
    --save_trimmed \
    --cytosine_report \
    --relax_mismatches \
    --num_mismatches 0.6 \
    --clip_r1 10 \
    --clip_r2 10 \
    --max_cpus 8

# Run with custom parameters using local reference
nextflow run nf-core/methylseq \
    -profile conda \
    --input samplesheet.csv \
    --outdir results \
    --fasta $PWD/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    --aligner bwameth \
    --save_reference \
    --save_trimmed \
    --save_align_intermeds \
    --cytosine_report \
    --relax_mismatches \
    --num_mismatches 0.6 \
    --clip_r1 10 \
    --clip_r2 10 \
    --max_cpus 8

# Quick test with chromosome 22 only
nextflow run nf-core/methylseq \
    -profile conda \
    --input samplesheet.csv \
    --outdir results_chr22 \
    --fasta $PWD/references/chr22/chr22.fa \
    --aligner bwameth \
    --max_cpus 8
```

### Resume Failed Runs

```bash
# Resume a previously failed run
nextflow run nf-core/methylseq \
    -profile conda \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    --max_cpus 8 \
    -resume
```

## Output Structure

After successful completion, your results directory will contain:

```
results/
├── bwameth/                    # bwa-meth alignment results
├── fastqc/                     # Quality control reports
├── methyldackel/               # Methylation extraction results
├── multiqc/                    # Comprehensive quality report
├── pipeline_info/              # Pipeline execution information
├── trimgalore/                 # Trimmed reads (if enabled)
└── reference_genome/           # Reference genome files (if saved)
```

## Supported Genomes

The pipeline supports many reference genomes through iGenomes:

- **Human**: GRCh37, GRCh38/hg38
- **Mouse**: GRCm38/mm10, GRCm39/mm39
- **Other**: Many additional organisms available

You can also provide custom reference genomes using `--fasta` parameter.

## bwa-meth Workflow

This setup uses **bwa-meth** as the primary methylation sequence aligner, which offers several advantages:

### **bwa-meth Methodology**
- **Fast alignment**: Uses BWA-MEM algorithm optimized for bisulfite-treated reads
- **Accurate mapping**: Handles C→T and G→A conversions from bisulfite treatment
- **Memory efficient**: Lower memory requirements compared to other aligners
- **Comprehensive output**: Compatible with standard downstream analysis tools

### **Analysis Pipeline**
1. **Quality Control**: FastQC analyzes raw read quality
2. **Read Trimming**: Trim Galore removes adapters and low-quality bases
3. **Reference Preparation**: bwa-meth prepares bisulfite-converted reference genome
4. **Alignment**: bwa-meth aligns reads to the converted reference
5. **Methylation Extraction**: MethylDackel extracts methylation calls from alignments
6. **Quality Assessment**: MultiQC aggregates all quality metrics and results

## Local Reference Configuration

This setup includes a `local_reference` configuration profile that automatically uses your downloaded GRCh38 reference:

### **Benefits of Local Reference**
- ✅ **No network dependency**: Avoid download timeouts and connection issues
- ✅ **Faster startup**: Skip reference download step (saves 10-30 minutes)
- ✅ **Reproducible**: Same reference version every time
- ✅ **Offline capability**: Run without internet connection
- ✅ **Storage efficient**: Reuse index across multiple runs

### **Local Reference Files**
After downloading and indexing, you'll have:
```
references/GRCh38/
├── GCA_000001405.15_GRCh38_no_alt_analysis_set.fna          # Original reference (~800MB)
├── GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwameth.c2t     # C→T converted reference
├── GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwameth.c2t.amb # BWA index files
├── GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwameth.c2t.ann # (~15-20GB total)
├── GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwameth.c2t.bwt
├── GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwameth.c2t.pac
└── GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwameth.c2t.sa
```

### **Checking Index Completion**
```bash
# Verify all index files are present
ls -la references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwameth.c2t.*

# Should show: .amb, .ann, .bwt, .pac, .sa files
```

## Troubleshooting

### Common Issues

#### 1. Conda Environment Issues

```bash
# If conda activate fails, ensure conda is properly initialized
conda init zsh  # or bash
source ~/.zshrc  # or ~/.bashrc

# Recreate environment if needed
conda env remove -n methyl-microc
conda create -n methyl-microc python=3.11 -y
conda activate methyl-microc
conda install -c bioconda nextflow nf-core bwameth bwa methyldackel trim-galore fastqc multiqc samtools -y
```

#### 2. Memory Issues

```bash
# If processes fail due to memory, increase limits in nextflow.config
process {
    memory = 8.GB  # Increase from 4.GB
}

# Or use a profile with higher limits
nextflow run nf-core/methylseq -profile docker,local --max_memory 16.GB
```

#### 3. Conda Environment Issues

```bash
# If conda activate fails, ensure conda is properly initialized
conda init zsh  # or bash
source ~/.zshrc  # or ~/.bashrc

# Recreate environment if needed
conda env remove -n methyl-microc
conda create -n methyl-microc python=3.11 -y
conda activate methyl-microc
conda install -c bioconda nextflow nf-core -y
```

#### 4. Tool Installation Issues

```bash
# If specific tools fail to install, try installing them individually
conda install -c bioconda bwameth -y
conda install -c bioconda bwa -y
conda install -c bioconda methyldackel -y
conda install -c bioconda trim-galore -y
conda install -c bioconda fastqc -y

# Check tool versions and compatibility
bwameth.py --version
bwa
MethylDackel --version
```

#### 5. Reference Genome Download Issues

```bash
# If reference download fails, try alternative sources
# Option A: Resume interrupted download
curl -L -C - -o GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz \
  'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz'

# Option B: Use UCSC source
curl -L -o hg38.fa.gz "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"

# Option C: Use Ensembl source
curl -L -o Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
  "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
```

#### 6. Index Creation Issues

```bash
# If bwa-meth indexing fails or is interrupted
cd references/GRCh38

# Remove incomplete index files
rm -f *.bwameth.c2t*

# Restart indexing
conda activate methyl-microc
bwameth.py index GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

# Monitor progress (indexing takes 1-2 hours)
# Check available disk space (needs ~20GB free)
df -h .
```

### Getting Help

1. **Check the logs**: Look in the `work/` directory for detailed error logs
2. **nf-core community**: Visit [nf-co.re](https://nf-co.re) for documentation
3. **GitHub Issues**: Report bugs at [nf-core/methylseq](https://github.com/nf-core/methylseq/issues)
4. **Slack**: Join the nf-core Slack workspace for community support

## Environment Management

### Deactivating the Environment

```bash
# Deactivate the conda environment
conda deactivate
```

### Updating the Pipeline

```bash
# Update to the latest version
conda activate methyl-microc
nextflow pull nf-core/methylseq

# Check for updates to nf-core tools
conda update -c bioconda nf-core
```

### Cleaning Up

```bash
# Remove old work directories (saves disk space)
rm -rf work/

# Clean up conda caches (optional)
conda clean --all -y

# Remove conda environment (if needed)
conda env remove -n methyl-microc
```

## File Descriptions

- **`nextflow.config`**: Main pipeline configuration with conda settings
- **`environment.yml`**: Conda environment specification with all required tools
- **`samplesheet.csv`**: Input sample sheet defining your data
- **`conf/test.config`**: Test-specific configuration parameters
- **`conf/local_reference.config`**: Configuration for using local reference genome
- **`references/`**: Directory for local reference genomes and indices
  - `GRCh38/`: GRCh38 reference genome and bwa-meth index files
  - `chr22/`: Chromosome 22 for quick testing (optional)
- **`tests/`**: Directory containing test FASTQ files
  - `tiny_r1.fastq.gz`: Forward reads for testing
  - `tiny_r2.fastq.gz`: Reverse reads for testing

## References

- **nf-core/methylseq**: [https://nf-co.re/methylseq](https://nf-co.re/methylseq)
- **Nextflow**: [https://www.nextflow.io](https://www.nextflow.io)
- **Conda**: [https://docs.conda.io](https://docs.conda.io)
- **Bioconda**: [https://bioconda.github.io](https://bioconda.github.io)
- **bwa-meth**: [https://github.com/brentp/bwa-meth](https://github.com/brentp/bwa-meth)
- **MethylDackel**: [https://github.com/dpryan79/MethylDackel](https://github.com/dpryan79/MethylDackel)

## Citation

If you use this environment setup or the nf-core/methylseq pipeline, please cite:

> Ewels, P. A., Peltzer, A., Fillinger, S., Patel, H., Alneberg, J., Wilm, A., Garcia, M. U., Di Tommaso, P., & Nahnsen, S. (2020). The nf-core framework for community-curated bioinformatics pipelines. Nature Biotechnology, 38(3), 276-278. doi: 10.1038/s41587-020-0439-x

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
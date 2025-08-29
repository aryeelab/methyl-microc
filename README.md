# Methyl-MicroC Analysis Environment

A complete environment setup for running the nf-core/methylseq pipeline on methylation sequencing data using conda-based local execution.

## Overview

This repository provides a fully configured environment for analyzing bisulfite sequencing data using the nf-core/methylseq pipeline. The setup includes:

- Conda environment with Nextflow, nf-core tools, and all bioinformatics dependencies
- Local conda-based execution (no Docker required)
- Pre-configured sample sheets and configuration files
- Test data for validation
- Comprehensive documentation and usage examples

## Prerequisites

Before setting up this environment, ensure you have:

- **macOS** (this setup is optimized for macOS with Apple Silicon)
- **Homebrew** package manager
- **Conda/Miniconda** available in your shell environment
- At least 8GB of RAM and 20GB of free disk space
- **Note**: Docker is NOT required for this conda-based setup

## Quick Start

### 1. Activate the Environment

```bash
# Activate the conda environment
source ~/.zshrc  # or ~/.bashrc for bash users
conda activate methylseq-env

# Verify installations
nextflow -version
nf-core --version
bismark --version
trim_galore --version
```

### 2. Run Test Analysis

```bash
# Run the pipeline with test data using conda
nextflow run nf-core/methylseq \
    -profile conda,test \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh37 \
    --max_cpus 8
```

### 3. Monitor Progress

The pipeline will automatically:
- Download required reference genomes
- Pull necessary Docker containers
- Process the test data through the complete methylation analysis workflow
- Generate comprehensive reports and results

## Environment Setup (Detailed)

### Step 1: Create Conda Environment

```bash
# Create the environment with Python 3.11
conda create -n methylseq-env python=3.11 -y

# Activate the environment
conda activate methylseq-env

# Install required packages
conda install -c bioconda nextflow nf-core -y
```

### Step 2: Install Bioinformatics Tools

```bash
# Install all required bioinformatics tools
conda install -c bioconda bismark trim-galore fastqc multiqc samtools bowtie2 bedtools picard qualimap -y

# Verify tool installations
bismark --version
trim_galore --version
fastqc --version
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
    genome = 'GRCh37'

    aligner = 'bismark'
    save_reference = true
}
```

## Usage Examples

### Basic Analysis

```bash
# Activate environment
conda activate methylseq-env

# Run with your own data
nextflow run nf-core/methylseq \
    -profile conda \
    --input your_samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    --aligner bismark \
    --max_cpus 8
```

### Advanced Options

```bash
# Run with custom parameters
nextflow run nf-core/methylseq \
    -profile conda \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    --aligner bismark \
    --save_reference \
    --save_trimmed \
    --cytosine_report \
    --relax_mismatches \
    --num_mismatches 0.6 \
    --clip_r1 10 \
    --clip_r2 10 \
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
├── bismark/                    # Bismark alignment results
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

## Troubleshooting

### Common Issues

#### 1. Conda Environment Issues

```bash
# If conda activate fails, ensure conda is properly initialized
conda init zsh  # or bash
source ~/.zshrc  # or ~/.bashrc

# Recreate environment if needed
conda env remove -n methylseq-env
conda create -n methylseq-env python=3.11 -y
conda activate methylseq-env
conda install -c bioconda nextflow nf-core bismark trim-galore fastqc multiqc samtools -y
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
conda env remove -n methylseq-env
conda create -n methylseq-env python=3.11 -y
conda activate methylseq-env
conda install -c bioconda nextflow nf-core -y
```

#### 4. Tool Installation Issues

```bash
# If specific tools fail to install, try installing them individually
conda install -c bioconda bismark -y
conda install -c bioconda trim-galore -y
conda install -c bioconda fastqc -y

# Check tool versions and compatibility
bismark --version
trim_galore --version
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
conda activate methylseq-env
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
conda env remove -n methylseq-env
```

## File Descriptions

- **`nextflow.config`**: Main pipeline configuration with conda settings
- **`environment.yml`**: Conda environment specification with all required tools
- **`samplesheet.csv`**: Input sample sheet defining your data
- **`conf/test.config`**: Test-specific configuration parameters
- **`tests/`**: Directory containing test FASTQ files
  - `tiny_r1.fastq.gz`: Forward reads for testing
  - `tiny_r2.fastq.gz`: Reverse reads for testing

## References

- **nf-core/methylseq**: [https://nf-co.re/methylseq](https://nf-co.re/methylseq)
- **Nextflow**: [https://www.nextflow.io](https://www.nextflow.io)
- **Conda**: [https://docs.conda.io](https://docs.conda.io)
- **Bioconda**: [https://bioconda.github.io](https://bioconda.github.io)
- **Bismark**: [https://github.com/FelixKrueger/Bismark](https://github.com/FelixKrueger/Bismark)

## Citation

If you use this environment setup or the nf-core/methylseq pipeline, please cite:

> Ewels, P. A., Peltzer, A., Fillinger, S., Patel, H., Alneberg, J., Wilm, A., Garcia, M. U., Di Tommaso, P., & Nahnsen, S. (2020). The nf-core framework for community-curated bioinformatics pipelines. Nature Biotechnology, 38(3), 276-278. doi: 10.1038/s41587-020-0439-x

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
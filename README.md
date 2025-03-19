# RNA-seq Analysis Pipeline

A streamlined, user-friendly pipeline for processing RNA-seq data from raw FASTQ files to aligned BAMs with quality control metrics.

## Overview

This pipeline automates the following steps:
- Quality control using FastQC
- Read trimming with Trimmomatic
- Alignment using STAR
- Duplicate marking with Picard
- Quality metrics collection
- Summary report generation with MultiQC

## Quick Start

### Installation

1. Clone this repository:
```bash
git clone https://github.com/IkramInf/rnaseq-pipeline.git
cd rnaseq-pipeline
```

2. Run the installation script:
```bash
bash install.sh
```

3. Activate the environment:
```bash
source ~/activate_rna_seq_pipeline.sh
```

### Running the Pipeline

Basic usage:
```bash
./rna_seq_pipeline.sh -i /path/to/fastq/files -g hg38 -o results
```

For help and all options:
```bash
./rna_seq_pipeline.sh --help
```

## Requirements

- Linux-based operating system (Ubuntu/Debian or CentOS/RHEL recommended)
- 16GB+ RAM (more for large genomes)
- 100GB+ disk space
- Python 3.8+
- Java Runtime Environment 11+

## Command-line Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i, --input` | Directory containing FASTQ files | Current directory |
| `-g, --genome` | Reference genome (hg19 or hg38) | hg19 |
| `-o, --output` | Output directory | results |
| `-c, --crop` | Crop reads to this length | 50 |
| `-l, --leading` | Remove leading low quality bases | 0 |
| `-t, --trailing` | Remove trailing low quality bases | 0 |
| `-p, --threads` | Number of threads to use | All but one CPU |
| `-h, --help` | Display help message | - |

## Output Structure

```
results/
├── sample1/
│   ├── fastqc_reports/
│   ├── sample1_PE_Aligned.sortedByCoord.out.bam
│   ├── sample1_marked_duplicates.bam
│   ├── sample1_RNA_Metrics.txt
│   ├── sample1_mapping_metrics.txt
│   ├── sample1_duplicate_metrics.txt
│   └── multiqc_report.html
├── sample2/
...
```

## Workflow Details

1. **Quality Control**: FastQC analyzes read quality
2. **Trimming**: Trimmomatic removes adapters and low-quality bases
3. **Alignment**: STAR aligns reads to the reference genome
4. **Metrics Collection**: Picard tools collect RNA-seq and alignment metrics
5. **Duplicate Marking**: Picard identifies PCR duplicates
6. **Summary Report**: MultiQC creates a comprehensive report

## Customization

The pipeline supports both single-end and paired-end data automatically. It will detect FASTQ file pairs based on the standard naming convention: *_R1.fastq.gz/*_R2.fastq.gz.

## Troubleshooting

**Not enough memory for STAR**:
- Reduce the number of threads with `-p` option
- Run on a machine with more RAM

**Missing dependencies**:
- Re-run the installation script
- Check system requirements

## Citation

If you use this pipeline, please cite:
- STAR aligner: Dobin A, et al. (2013) PMID: 23104886
- Trimmomatic: Bolger AM, et al. (2014) PMID: 24695404
- Picard: https://broadinstitute.github.io/picard/
- FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- MultiQC: Ewels P, et al. (2016) PMID: 27312411

## License

This project is licensed under the MIT License - see the LICENSE file for details.

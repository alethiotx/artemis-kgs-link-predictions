# artemis-kgs-link-predictions

Knowledge Graph Link Prediction Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker-enabled-blue.svg)](https://www.docker.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Supported Datasets](#supported-datasets)
- [Pipeline Architecture](#pipeline-architecture)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Configuration](#configuration)
- [Output Files](#output-files)
- [Pipeline Scripts](#pipeline-scripts)
- [Docker Image](#docker-image)
- [Contributing](#contributing)
- [License](#license)

## Overview

artemis-kgs-link-predictions is a scalable Nextflow pipeline that leverages knowledge graph embeddings to predict gene-term associations across multiple biomedical knowledge graphs. The pipeline uses PyKEEN-trained models to generate predictions for diseases, pathways, biological processes, and other biomedical entities.

Built for precision medicine and drug discovery applications, artemis-kgs-link-predictions automates the entire workflow from data preparation to prediction aggregation, supporting parallel processing of thousands of terms across different knowledge graph databases.

## Features

- **Multi-Dataset Support**: Works with Hetionet, BioKG, OpenBioLink, and PrimeKG
- **Scalable Predictions**: Parallel processing of terms for efficient large-scale predictions
- **Knowledge Graph Embeddings**: Utilizes PyKEEN models for accurate association scoring
- **Flexible Sampling**: Optional downsampling to 10,000 terms for faster testing
- **Cloud-Native**: Built-in S3 support for data storage and retrieval
- **Reproducible**: Containerized execution with consistent environments
- **Test Mode**: Quick validation with limited terms before full production runs

## Supported Datasets

### Hetionet
A network of biology and disease knowledge integrating 47,031 nodes (11 types) and 2,250,197 relationships (24 types). Predicts associations for:
- Biological Processes
- Pathways
- Cellular Components
- Molecular Functions
- Anatomies
- Compounds
- Diseases
- Genes

### BioKG
Comprehensive biomedical knowledge graph focusing on proteins, diseases, drugs, and pathways. Extracts:
- Human Proteins
- Diseases
- Pathways
- Drugs
- Protein Complexes
- Genetic Disorders

### OpenBioLink
Large-scale open biomedical knowledge graph containing gene, disease, and phenotype associations. Uses NCBI gene information for protein-coding genes.

### PrimeKG
Precision Medicine Knowledge Graph integrating 20+ biomedical resources with focus on drug discovery. Covers gene/protein entities across diverse biomedical contexts.

## Pipeline Architecture

![](pipeline.png)

The pipeline consists of three main stages:

### 1. Prepare
Extracts and processes knowledge graph data:
- Loads the specified dataset using PyKEEN
- Concatenates training, testing, and validation triples
- Extracts genes and terms based on dataset-specific criteria
- Generates hash tables mapping IDs to human-readable names
- Optional downsampling to 10,000 terms

**Outputs:**
- `triples.npy`: All knowledge graph triples
- `terms.csv`: List of terms for prediction
- `genes_hash_table.csv`: Gene ID to name mappings
- `terms_hash_table.csv`: Term ID to name and type mappings

### 2. Predict
Generates predictions for each term in parallel:
- Uses trained PyKEEN models for target prediction
- Processes each term independently for scalability
- Determines appropriate relation types based on dataset and term type
- Outputs individual prediction files per term

**Outputs:**
- `*_predictions.csv`: Per-term prediction files

### 3. Summarize
Aggregates and annotates all predictions:
- Combines predictions from all terms
- Maps IDs to human-readable gene and term names
- Generates final consolidated prediction matrix

**Outputs:**
- `predictions.csv`: Final aggregated predictions with annotations

## Requirements

- **Nextflow** >= 23.04.0
- **Docker** or **Singularity** (for containerized execution)
- **AWS CLI** (optional, for S3 access)

### Python Dependencies (included in Docker image)
- pandas
- numpy
- pykeen == 1.11.0
- torch
- s3fs

## Installation

1. **Install Nextflow:**

```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

2. **Clone the repository:**

```bash
git clone https://github.com/alethiotx/artemis-kgs-link-predictions.git
cd artemis-kgs-link-predictions
```

3. **Verify installation:**

```bash
nextflow -version
```

## Quick Start

### Test Run (Local)
Run with a small subset of terms to validate setup:

```bash
nextflow run main.nf \
  -profile local \
  --dataset hetionet \
  --model /path/to/trained_model.pkl \
  --env test \
  --sample true
```

### Production Run (Hetionet)
Full prediction run with all terms:

```bash
nextflow run main.nf \
  -profile hetionet \
  --dataset hetionet \
  --model s3://bucket/path/trained_model.pkl \
  --env prod \
  --sample false
```

## Usage

### Basic Command Structure

```bash
nextflow run main.nf \
  -profile <PROFILE> \
  --dataset <DATASET> \
  --model <MODEL_PATH> \
  [--env <ENV>] \
  [--sample <BOOLEAN>] \
  [--outdir <OUTPUT_DIR>]
```

### Parameters

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `--dataset` | Yes | `null` | Dataset name: `hetionet`, `biokg`, `openbiolink`, or `primekg` |
| `--model` | Yes | `null` | Path to trained PyKEEN model (local or S3 URI) |
| `--env` | No | `prod` | Environment mode: `prod` (full run) or `test` (2 terms) |
| `--sample` | No | `true` | Downsample to 10,000 terms: `true` or `false` |
| `--outdir` | No | S3 bucket | Output directory for results |

### Profiles

| Profile | Description | Use Case |
|---------|-------------|----------|
| `local` | Local execution | Development and testing |
| `hetionet` | Hetionet-specific config | Hetionet predictions |
| `biokg` | BioKG-specific config | BioKG predictions |
| `openbiolink` | OpenBioLink-specific config | OpenBioLink predictions |
| `primekg` | PrimeKG-specific config | PrimeKG predictions |

### Examples

#### 1. Test run with BioKG (2 terms only)
```bash
nextflow run main.nf \
  -profile biokg \
  --dataset biokg \
  --model s3://bucket/biokg/trained_model.pkl \
  --env test
```

#### 2. Full production run with OpenBioLink (no sampling)
```bash
nextflow run main.nf \
  -profile openbiolink \
  --dataset openbiolink \
  --model s3://bucket/openbiolink/trained_model.pkl \
  --env prod \
  --sample false
```

#### 3. Sampled production run with PrimeKG (10,000 terms)
```bash
nextflow run main.nf \
  -profile primekg \
  --dataset primekg \
  --model s3://bucket/primekg/trained_model.pkl \
  --env prod \
  --sample true
```

#### 4. Custom output directory
```bash
nextflow run main.nf \
  -profile hetionet \
  --dataset hetionet \
  --model /local/path/trained_model.pkl \
  --outdir /custom/output/path
```

## Configuration

### Base Configuration ([conf/base.config](conf/base.config))

Defines resource allocation for process labels:
- `process_single`: 1 CPU, 4 GB RAM, 4h timeout
- `process_low`: 2 CPUs, 12 GB RAM, 4h timeout
- `process_medium`: 6 CPUs, 36 GB RAM, 8h timeout
- `process_high`: 12 CPUs, 72 GB RAM, 16h timeout
- `process_long`: 20h timeout

All processes have automatic retry (up to 5 times) with exponential resource scaling.

### Dataset-Specific Configurations

Each dataset profile ([conf/hetionet.config](conf/hetionet.config), [conf/biokg.config](conf/biokg.config), etc.) sets:
- Dataset name
- Default model path in S3

### Custom Configuration

Create a custom config file:

```nextflow
// custom.config
params {
    dataset = 'hetionet'
    model = '/my/custom/model.pkl'
    sample = false
}

process {
    withLabel:process_single {
        cpus = 2
        memory = '8.GB'
    }
}
```

Run with custom config:
```bash
nextflow run main.nf -c custom.config -profile hetionet
```

## Output Files

### Directory Structure

```
<outdir>/
└── <dataset>/
    ├── prepare/
    │   ├── terms.csv                    # List of terms for prediction
    │   ├── triples.npy                  # All knowledge graph triples
    │   ├── genes_hash_table.csv         # Gene ID to name mapping
    │   └── terms_hash_table.csv         # Term ID to name/type mapping
    ├── predict/
    │   ├── term1_predictions.csv        # Predictions for term 1
    │   ├── term2_predictions.csv        # Predictions for term 2
    │   └── ...
    └── summarize/
        └── predictions.csv              # Final aggregated predictions
```

### Output File Formats

#### predictions.csv (Final Output)
Matrix format with genes as rows and terms as columns:
```
Gene Name,Disease A,Pathway B,Biological Process C
GENE1,0.95,0.87,0.76
GENE2,0.82,0.91,0.68
...
```

#### terms.csv
One term per line:
```
Disease::Alzheimer's disease
Pathway::MAPK signaling pathway
...
```

#### Hash Tables
Two-column CSV files mapping IDs to names:
```
GENE123,TP53
DISEASE456,Alzheimer's disease
...
```

## Pipeline Scripts

### [bin/prepare.py](bin/prepare.py)
Prepares knowledge graph data by extracting terms, genes, and triples.

**Key Functions:**
- `validate_arguments()`: Validates command line arguments
- `load_dataset()`: Loads PyKEEN dataset
- `save_triples()`: Concatenates and saves all triples
- `process_hetionet()`: Extracts Hetionet-specific terms
- `process_biokg()`: Extracts BioKG-specific terms
- `process_openbiolink()`: Extracts OpenBioLink-specific terms
- `process_primekg()`: Extracts PrimeKG-specific terms

### [bin/predict.py](bin/predict.py)
Generates predictions for individual terms using trained PyKEEN models.

**Key Functions:**
- `pipeline()`: Core prediction function using PyKEEN
- Handles dataset-specific relation types
- Supports both head→tail and tail→head predictions

### [bin/summarize.py](bin/summarize.py)
Aggregates predictions and maps IDs to human-readable names.

**Features:**
- Combines all prediction files
- Maps gene IDs to gene names
- Maps term IDs to term names
- Transposes and sorts final output

## Docker Image

The pipeline uses the Docker image from the [artemis-kgs-embeddings](https://github.com/alethiotx/artemis-kgs-embeddings) repository:

```
public.ecr.aws/alethiotx/artemis-kgs-embeddings:latest
```

**Includes:**
- Python 3.x
- PyKEEN 1.11.0
- PyTorch
- pandas, numpy
- s3fs for S3 access

To use a custom image, modify [nextflow.config](nextflow.config):
```nextflow
process {
    container = 'your-registry/your-image:tag'
}
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

**Main Contributor:**
- Vladimir Kiselev

## License

This project is licensed under the [MIT License](LICENSE).


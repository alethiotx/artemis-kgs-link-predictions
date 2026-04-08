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
- [Data Leakage Prevention](#data-leakage-prevention)
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

Built for precision medicine and drug discovery applications, artemis-kgs-link-predictions automates the entire workflow from data preparation to prediction aggregation. Each pipeline run processes a single dataset across all 4 embedding models (ComplEx, DistMult, RotatE, TransE) with full parallel processing of terms using containerized execution. To run all 4 datasets, launch 4 separate runs via Seqera Platform.

This pipeline consumes embeddings produced by the upstream [artemis-kgs-embeddings](https://github.com/alethiotx/artemis-kgs-embeddings) pipeline, which trains models on triples that have been filtered to prevent data leakage (Drug × Clinical Gene triples are removed before training).

## Features

- **Multi-Dataset Support**: Works with Hetionet, BioKG, OpenBioLink, and PrimeKG
- **Multi-Embedding Support**: Runs ComplEx, DistMult, RotatE, and TransE models per dataset in a single invocation
- **Scalable Predictions**: Parallel processing of terms for efficient large-scale predictions
- **Entity Vocabulary Alignment**: Filters entities to match upstream model vocabulary, preventing silent prediction errors from entity ID mismatches
- **Knowledge Graph Embeddings**: Utilizes PyKEEN 1.11.0 models for accurate association scoring
- **Dual Output Formats**: Generates both CSV and Parquet prediction matrices
- **Flexible Sampling**: Optional downsampling to 10,000 terms for faster testing
- **Cloud-Native**: Built-in S3 support for model loading and data storage
- **Reproducible**: Containerized execution with Docker and deterministic random seeds
- **Test Mode**: Quick validation with 2 terms before full production runs
- **Dynamic Resources**: PrimeKG processes automatically receive increased memory (48 GB) and CPUs
- **Automatic Retry**: Processes automatically retry with exponential resource scaling (up to 5 retries)

## Supported Datasets

### Hetionet
A network of biology and disease knowledge integrating 47,031 nodes (11 types) and 2,250,197 relationships (24 types). Uses dataset-specific relation types (GpBP, GpPW, GpCC, GpMF, AeG, CbG, DaG, GiG) for predictions. Predicts associations for:
- Biological Processes (Gene participates in Biological Process)
- Pathways (Gene participates in Pathway)
- Cellular Components (Gene participates in Cellular Component)
- Molecular Functions (Gene participates in Molecular Function)
- Anatomies (Anatomy expresses Gene)
- Compounds (Compound binds Gene)
- Diseases (Disease associates with Gene)
- Genes (Gene interacts with Gene)

### BioKG
Comprehensive biomedical knowledge graph focusing on proteins, diseases, drugs, and pathways. Uses metadata files (biokg.metadata.protein.tsv, biokg.metadata.disease.tsv, biokg.metadata.pathway.tsv, biokg.metadata.drug.tsv) and links file (biokg.links.tsv) for extraction. Extracts:
- Human Proteins (filtered by SPECIES='HUMAN')
- Diseases (PROTEIN_DISEASE_ASSOCIATION)
- Pathways (PROTEIN_PATHWAY_ASSOCIATION)
- Drugs (Drug-Protein Interaction, reverse direction)
- Protein Complexes (MEMBER_OF_COMPLEX)
- Genetic Disorders (RELATED_GENETIC_DISORDER)

### OpenBioLink
Large-scale open biomedical knowledge graph containing gene, disease, and phenotype associations. Uses NCBI gene information (Homo_sapiens.gene_info) for protein-coding genes and nodes.csv for graph structure. All predictions use Gene → Term direction with relation types:
- GENE_EXPRESSED_ANATOMY (anatomy)
- GENE_GO (GO terms)
- GENE_DRUG (drugs)
- GENE_PHENOTYPE (phenotypes)
- GENE_PATHWAY (pathways)
- GENE_DIS (diseases)
- GENE_GENE (gene interactions)

### PrimeKG
Precision Medicine Knowledge Graph integrating 20+ biomedical resources with focus on drug discovery. Uses nodes.csv for entity extraction. All predictions use Term → Gene (reverse) direction. Covers gene/protein entities across diverse biomedical contexts:
- anatomy (anatomy_protein_present)
- biological_process (bioprocess_protein)
- cellular_component (cellcomp_protein)
- disease (disease_protein)
- drug (drug_protein)
- exposure (exposure_protein)
- molecular_function (molfunc_protein)
- pathway (pathway_protein)
- effect/phenotype (phenotype_protein)
- gene/protein (protein_protein)

## Pipeline Architecture

The pipeline consists of three main stages, executed for each dataset × embedding combination:

### 1. Prepare
Extracts and processes knowledge graph data, filtered to the upstream model's entity vocabulary:
- Loads the specified dataset using PyKEEN (Hetionet, BioKG, OpenBioLink, or PrimeKG)
- Loads the upstream training triples via `TriplesFactory.from_path_binary()` to obtain the model's entity vocabulary
- Filters extracted entities to only those present in the model vocabulary, ensuring entity-to-ID alignment
- Extracts genes and terms based on dataset-specific criteria using metadata files from S3
- Generates hash tables mapping IDs to human-readable names
- Optional downsampling to 10,000 terms using random seed 42 for reproducibility

**Outputs:**
- `terms.csv`: List of terms for prediction (one per line), filtered to model vocabulary
- `genes_hash_table.csv`: Gene ID to name mappings (2 columns: ID, Name), filtered to model vocabulary
- `terms_hash_table.csv`: Term ID to name and type mappings (2-3 columns: ID, Name, [Type])

### 2. Predict
Generates predictions for each term in parallel:
- Loads the upstream training triples via `TriplesFactory.from_path_binary()` to reconstruct the exact entity-to-ID mapping used during model training
- Uses trained PyKEEN models loaded with `torch.load` (weights_only=False)
- Filters genes to only those in the model's entity vocabulary
- If a term is not in the model vocabulary, writes an empty prediction CSV and exits gracefully
- Processes each term independently using Nextflow's parallel processing
- Determines appropriate relation types and prediction direction (head→tail or tail→head) based on dataset and term type
- Uses PyKEEN's `predict.predict_target()` function for scoring
- Handles dataset-specific relation mappings (HETIONET_RELATIONS, BIOKG_RELATIONS, etc.)
- Applies random seeds for reproducibility (seed=42)
- Creates safe filenames by replacing special characters with underscores
- Outputs individual prediction files per term in transposed format (term as row, genes as columns)

**Outputs:**
- `<term>_predictions.csv`: Per-term prediction files with scores rounded to 1 decimal place

### 3. Summarize
Aggregates and annotates all predictions:
- Combines predictions from all terms using `collectFile` to create a metafile
- Transposes data from terms×genes to genes×terms format
- Maps IDs to human-readable gene and term names using hash tables
- Sorts columns alphabetically for consistency
- Removes duplicate genes and terms
- Calculates final output statistics (file size, total predictions)
- Saves in both CSV and Parquet formats

**Outputs:**
- `predictions.csv`: Final aggregated predictions with gene names as rows and term names as columns
- `predictions.parquet`: Same data in Apache Parquet format for efficient downstream processing

## Data Leakage Prevention

This pipeline is designed to work with the upstream [artemis-kgs-embeddings](https://github.com/alethiotx/artemis-kgs-embeddings) pipeline, which filters Drug × Clinical Gene triples (495 clinical genes across 7 disease CSVs) from each knowledge graph before training embeddings. This prevents the model from memorizing known drug-gene associations during training.

To ensure entity alignment between the filtered model and this prediction pipeline:
- **prepare.py** loads the upstream `training_triples/` directory to obtain the entity vocabulary, and filters all extracted genes and terms to only those present in the model's vocabulary
- **predict.py** loads the same `training_triples/` directory via `TriplesFactory.from_path_binary()` to reconstruct the exact entity-to-ID mapping used during training, ensuring predictions use correct entity indices

## Requirements

- **Nextflow** >= 23.04.0
- **Docker** or **Singularity** (for containerized execution)
- **AWS CLI** (optional, for S3 access)

### Python Dependencies (included in Docker image)
As specified in [requirements.txt](requirements.txt):
- pandas
- numpy
- pykeen==1.11.0
- torch
- s3fs
- pyarrow

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
Run Hetionet with ComplEx and 2 terms to validate setup:

```bash
nextflow run main.nf -profile local
```

This uses [conf/local.config](conf/local.config) which defaults to Hetionet/ComplEx, test mode (2 terms), and local output directory.

### Production Run
Full prediction run for a single dataset across all 4 embeddings:

```bash
nextflow run main.nf \
  --dataset hetionet
```

## Usage

### Basic Command Structure

```bash
nextflow run main.nf \
  --dataset <DATASET> \
  [-profile <PROFILE>] \
  [--embeddings <EMBEDDINGS>] \
  [--env <ENV>] \
  [--sample <BOOLEAN>] \
  [--outdir <OUTPUT_DIR>]
```

### Parameters

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `--dataset` | Yes | `null` | Dataset to process: `hetionet`, `biokg`, `openbiolink`, or `primekg` |
| `--embeddings` | No | All 4 models | Embedding models to run (default: ComplEx, DistMult, RotatE, TransE) |
| `--env` | No | `prod` | Environment mode: `prod` (full run) or `test` (2 terms) |
| `--sample` | No | `true` | Downsample to 10,000 terms: `true` or `false` |
| `--outdir` | No | S3 bucket | Output directory for results |
| `--s3_embeddings_base` | No | S3 bucket | Base S3 path for upstream embedding artifacts |

Models and training triples are automatically resolved from the S3 embeddings base path:
- Model: `{s3_embeddings_base}/{dataset}/{embedding}/trained_model.pkl`
- Training triples: `{s3_embeddings_base}/{dataset}/{embedding}/training_triples/`

### Profiles

| Profile | Description | Use Case |
|---------|-------------|----------|
| `local` | Local execution, Docker disabled, Hetionet/ComplEx, test mode | Development and testing |
| `standard` | Default profile with base resource configuration | Production runs |

### Examples

#### 1. Hetionet with all embeddings
```bash
nextflow run main.nf --dataset hetionet
```

#### 2. PrimeKG with all embeddings (no sampling)
```bash
nextflow run main.nf --dataset primekg --sample false
```

#### 3. Single embedding only
```bash
nextflow run main.nf --dataset biokg --embeddings ComplEx
```

#### 4. Test mode (2 terms only)
```bash
nextflow run main.nf --dataset openbiolink --env test
```

#### 5. Custom output directory
```bash
nextflow run main.nf --dataset hetionet --outdir /custom/output/path
```

#### 6. Run all 4 datasets (bash loop)
```bash
for dataset in hetionet biokg openbiolink primekg; do
  nextflow run main.nf --dataset $dataset -resume
done
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

### Dynamic Resource Allocation

The `prepare` and `predict` processes dynamically allocate resources based on the dataset:
- **PrimeKG**: 2 CPUs, 48 GB memory (scales with retry attempts)
- **All other datasets**: 1 CPU, 4 GB memory (scales with retry attempts)

### Local Configuration ([conf/local.config](conf/local.config))

Overrides for local development:
- Dataset: `hetionet`, Embedding: `ComplEx` only
- Output to local `output/` directory
- Test mode (`env = 'test'`)
- No downsampling (`sample = false`)
- Docker disabled
- Terminate on error (no retries)

### Custom Configuration

Create a custom config file:

```nextflow
// custom.config
params {
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
nextflow run main.nf --dataset hetionet -c custom.config
```

## Output Files

### Directory Structure

```
<outdir>/
└── <dataset>/
    └── <embedding>/
        ├── prepare/
        │   ├── terms.csv                    # List of terms for prediction
        │   ├── genes_hash_table.csv         # Gene ID to name mapping
        │   └── terms_hash_table.csv         # Term ID to name/type mapping
        ├── predict/
        │   ├── term1_predictions.csv        # Predictions for term 1
        │   ├── term2_predictions.csv        # Predictions for term 2
        │   └── ...
        └── summarize/
            ├── predictions.csv              # Final aggregated predictions (CSV)
            └── predictions.parquet          # Final aggregated predictions (Parquet)
```

Each pipeline run produces output directories for the specified dataset across all embeddings (4 directories by default).

### Output File Formats

#### predictions.csv / predictions.parquet (Final Output)
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
Prepares knowledge graph data by extracting terms and genes, filtered to the upstream model's entity vocabulary.

**Usage:** `./prepare.py <dataset> <downsample> <training_triples>`

**Key Functions:**
- `set_random_seeds(seed=42)`: Sets random seeds for reproducibility
- `validate_arguments()`: Validates 3 command line arguments (dataset, downsample, training_triples)
- `load_entity_vocabulary()`: Loads entity vocabulary from upstream training triples via `TriplesFactory.from_path_binary()`
- `load_dataset()`: Loads PyKEEN dataset (Hetionet, BioKG, OpenBioLink, PrimeKG)
- `get_entity_list()`: Extracts entities from KG and filters to model vocabulary
- `process_hetionet()`: Extracts Hetionet terms using node TSV with specific prefixes
- `process_biokg()`: Extracts BioKG terms from metadata files (proteins, diseases, pathways, drugs) and links file
- `process_openbiolink()`: Extracts OpenBioLink terms using NCBI gene info and nodes.csv
- `process_primekg()`: Extracts PrimeKG terms from nodes.csv
- `save_results()`: Saves terms, genes hash table, and terms hash table with optional downsampling

**Constants:**
- `RANDOM_SEED = 42`
- `S3_BASE_PATH = 's3://alethiotx-artemis/data/kgs/raw'`
- `DOWNSAMPLE_SIZE = 10000`

### [bin/predict.py](bin/predict.py)
Generates predictions for individual terms using trained PyKEEN models.

**Usage:** `./predict.py <dataset> <term> <terms_hash_table> <genes_hash_table> <training_triples> <model>`

**Key Functions:**
- `set_random_seeds(seed=42)`: Sets random seeds for PyTorch, NumPy, and Python random
- `validate_arguments()`: Validates 6 command line arguments
- `load_data()`: Loads training triples via `TriplesFactory.from_path_binary()`, filters genes to model vocabulary, loads model with `torch.load(weights_only=False)`
- `get_term_type()`: Extracts term type from identifier or hash table
- `generate_predictions()`: Core prediction function using PyKEEN's `predict.predict_target()`
- `process_hetionet()`: Generates Hetionet predictions with dataset-specific relations
- `process_biokg()`: Generates BioKG predictions with relation types
- `process_openbiolink()`: Generates OpenBioLink predictions (all Gene→Term direction)
- `process_primekg()`: Generates PrimeKG predictions (all Term→Gene direction)
- `save_predictions()`: Saves predictions with safe filenames and rounded scores

**Safety Features:**
- Terms not present in the model vocabulary produce an empty CSV and exit successfully
- Genes are filtered to model vocabulary before prediction

**Relation Mappings:**
- `HETIONET_RELATIONS`: 8 relation types with direction (GpBP, GpPW, GpCC, GpMF, AeG, CbG, DaG, GiG)
- `BIOKG_RELATIONS`: 6 relation types (PROTEIN_DISEASE_ASSOCIATION, PROTEIN_PATHWAY_ASSOCIATION, PPI, etc.)
- `OPENBIOLINK_RELATIONS`: 7 relation types (GENE_EXPRESSED_ANATOMY, GENE_GO, etc.)
- `PRIMEKG_RELATIONS`: 10 relation types (anatomy_protein_present, bioprocess_protein, etc.)

### [bin/summarize.py](bin/summarize.py)
Aggregates predictions and maps IDs to human-readable names.

**Usage:** `./summarize.py <metafile> <genes_hash_table> <terms_hash_table>`

**Key Functions:**
- `validate_arguments()`: Validates command line arguments and checks file existence
- `load_data()`: Loads combined predictions from metafile, transposes to genes×terms format
- `process_predictions()`: Maps IDs to names, sorts columns alphabetically, removes duplicates
- `save_predictions()`: Saves final predictions in both CSV and Parquet formats with file size statistics

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
- pyarrow for Parquet output

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

## Acknowledgements

- Public knowledge graph providers (Hetionet, BioKG, OpenBioLink, PrimeKG)
- PyKEEN, scikit-learn, and Nextflow communities
- Portions of this codebase were assisted using GitHub Copilot (Claude Sonnet 4.5) for code generation, refactoring, cleaning and documentation. The authors reviewed, modified, and validated all AI-assisted code. Responsibility for the correctness, performance, and reproducibility of the code rests entirely with the authors. No AI tools were used to generate scientific conclusions or interpretations in this study.
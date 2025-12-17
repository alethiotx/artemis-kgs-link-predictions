#!/usr/bin/env python3
"""
Knowledge Graph Dataset Preparation Script

This script processes multiple knowledge graph datasets (Hetionet, BioKG, OpenBioLink, PrimeKG)
to extract genes, terms, and triples for downstream analysis.

Usage:
    ./prepare.py <dataset> <downsample>

Arguments:
    dataset: Dataset name - 'hetionet', 'biokg', 'openbiolink', or 'primekg'
    downsample: Whether to downsample terms to 10,000 - 'true' or 'false'

Outputs:
    - triples.npy: All concatenated triples from train/test/validation sets
    - terms.csv: Extracted terms (diseases, pathways, etc.)
    - genes_hash_table.csv: Gene IDs with human-readable names
    - terms_hash_table.csv: Term IDs with names and types

Dependencies:
    - pykeen
    - numpy
    - pandas

Example:
    ./prepare.py hetionet false
    ./prepare.py biokg true
"""

import sys
import random
import numpy as np
import pandas as pd
from pykeen.datasets import Hetionet, BioKG, OpenBioLink, PrimeKG
from pykeen.triples import TriplesFactory

# Constants
RANDOM_SEED = 42
S3_BASE_PATH = 's3://alethiotx-artemis/data/kgs/raw'
DOWNSAMPLE_SIZE = 10000
SUPPORTED_DATASETS = ['hetionet', 'biokg', 'openbiolink', 'primekg']


def set_random_seeds(seed: int = RANDOM_SEED):
    """
    Set random seeds for reproducibility across all libraries.
    
    Args:
        seed: Random seed value (default: RANDOM_SEED)
    """
    random.seed(seed)
    np.random.seed(seed)
    
    # Set environment variable for additional reproducibility
    import os
    os.environ['PYTHONHASHSEED'] = str(seed)


def validate_arguments():
    """Validate command line arguments."""
    if len(sys.argv) != 3:
        print("Error: Invalid number of arguments")
        print("Usage: ./prepare.py <dataset> <downsample>")
        print(f"Supported datasets: {', '.join(SUPPORTED_DATASETS)}")
        sys.exit(1)
    
    dataset = sys.argv[1].lower()
    downsample = sys.argv[2].lower()
    
    if dataset not in SUPPORTED_DATASETS:
        print(f"Error: Unsupported dataset '{dataset}'")
        print(f"Supported datasets: {', '.join(SUPPORTED_DATASETS)}")
        sys.exit(1)
    
    if downsample not in ['true', 'false']:
        print("Error: downsample must be 'true' or 'false'")
        sys.exit(1)
    
    return dataset, downsample == 'true'


def load_dataset(dataset_name):
    """Load the specified knowledge graph dataset."""
    dataset_map = {
        'hetionet': Hetionet,
        'biokg': BioKG,
        'openbiolink': OpenBioLink,
        'primekg': PrimeKG
    }
    
    print(f"Loading {dataset_name} dataset...")
    return dataset_map[dataset_name]()


def save_triples(kg):
    """Concatenate and save all triples from the knowledge graph."""
    print("Concatenating triples...")
    triples = np.concatenate([
        kg.training.triples,
        kg.testing.triples,
        kg.validation.triples
    ])
    
    print(f"Saving {len(triples)} triples...")
    np.save('triples.npy', triples)
    
    triples_factory = TriplesFactory.from_labeled_triples(triples)
    return list(triples_factory.entity_to_id.keys())

def process_hetionet(terms_all, should_downsample):
    """
    Process Hetionet dataset to extract genes and terms.
    
    Extracts nodes with prefixes: Biological Process, Pathway, Cellular Component,
    Molecular Function, Anatomy, Compound, Disease, Gene.
    """
    print("Processing Hetionet data...")
    nodes_path = f'{S3_BASE_PATH}/hetionet/hetionet-v1.0-nodes.tsv'
    nodes = pd.read_csv(nodes_path, sep="\t")
    
    # Extract genes
    gene_ids = [x for x in terms_all if 'Gene::' in x]
    
    # Extract terms with specific prefixes
    term_prefixes = [
        'Biological Process::', 'Pathway::', 'Cellular Component::',
        'Molecular Function::', 'Anatomy::', 'Compound::', 'Disease::', 'Gene::'
    ]
    term_ids = [x for x in terms_all if any(prefix in x for prefix in term_prefixes)]
    
    # Create hash tables with names
    gene_hash_table = nodes[nodes['id'].isin(gene_ids)][['id', 'name']]
    term_hash_table = nodes[nodes['id'].isin(term_ids)][['id', 'name']]
    
    # Save results
    save_results(term_ids, gene_hash_table, term_hash_table, should_downsample)

def process_biokg(terms_all, should_downsample):
    """
    Process BioKG dataset to extract genes and terms.
    
    Extracts human proteins, diseases, pathways, drugs, complexes, and genetic disorders.
    """
    print("Processing BioKG data...")
    biokg_path = f'{S3_BASE_PATH}/biokg'
    
    # Load protein metadata and filter for human genes
    protein = pd.read_csv(f'{biokg_path}/biokg.metadata.protein.tsv', sep="\t", header=None)
    genes_hash_table = protein[(protein[0].isin(terms_all)) & (protein[1] == 'NAME')][[0, 2]]
    human_uniprot_ids = protein[
        (protein[1] == 'SPECIES') &
        (protein[0].isin(genes_hash_table[0])) &
        (protein[2] == 'HUMAN')
    ][0].reset_index(drop=True)
    genes_hash_table = genes_hash_table[genes_hash_table[0].isin(human_uniprot_ids)]
    
    # Extract different term types
    terms_list = []
    
    # Human proteins
    terms1 = protein[(protein[0].isin(terms_all)) & (protein[1] == 'SPECIES') & (protein[2] == 'HUMAN')][[0]]
    terms1[1] = 'Protein'
    terms_list.append(terms1)
    
    # Diseases
    disease = pd.read_csv(f'{biokg_path}/biokg.metadata.disease.tsv', sep="\t", header=None)
    terms2 = disease[(disease[0].isin(terms_all)) & (disease[1] == 'NAME')][[0]]
    terms2[1] = 'Disease'
    terms_list.append(terms2)
    
    # Pathways
    pathway = pd.read_csv(f'{biokg_path}/biokg.metadata.pathway.tsv', sep="\t", header=None)
    terms3 = pathway[(pathway[0].isin(terms_all)) & (pathway[1] == 'NAME')][[0]]
    terms3[1] = 'Pathway'
    terms_list.append(terms3)
    
    # Drugs
    drug = pd.read_csv(f'{biokg_path}/biokg.metadata.drug.tsv', sep="\t", header=None)
    terms4 = drug[(drug[0].isin(terms_all)) & (drug[1] == 'NAME')][[0]]
    terms4[1] = 'Drug'
    terms_list.append(terms4)
    
    # Complexes (from links file)
    links = pd.read_csv(f'{biokg_path}/biokg.links.tsv', sep="\t", header=None)
    terms5 = pd.DataFrame(links[links[1] == 'MEMBER_OF_COMPLEX'][2].drop_duplicates())
    terms5.columns = [0]
    terms5[1] = 'Complex'
    terms_list.append(terms5)
    
    # Genetic disorders (from links file)
    terms6 = pd.DataFrame(links[links[1] == 'RELATED_GENETIC_DISORDER'][2].drop_duplicates())
    terms6.columns = [0]
    terms6[1] = 'Genetic Disorder'
    terms_list.append(terms6)
    
    # Combine all terms
    all_terms = pd.concat(terms_list)
    term_hash_table = all_terms
    
    # Save results
    save_results(all_terms[0].tolist(), genes_hash_table, term_hash_table, should_downsample)

def process_openbiolink(terms_all, should_downsample):
    """
    Process OpenBioLink dataset to extract genes and terms.
    
    Uses NCBI gene information to map protein-coding genes.
    """
    print("Processing OpenBioLink data...")
    obl_path = f'{S3_BASE_PATH}/openbiolink'
    
    # Load gene information (protein-coding genes only)
    gene_info = pd.read_csv(f'{obl_path}/Homo_sapiens.gene_info', sep='\t')
    gene_info = gene_info[gene_info['type_of_gene'] == 'protein-coding'][['GeneID', 'Symbol']]
    gene_info = gene_info[~gene_info['Symbol'].duplicated()].drop_duplicates().reset_index(drop=True)
    
    # Load graph nodes
    nodes = pd.read_csv(f'{obl_path}/graph_files/nodes.csv', sep="\t", header=None)
    
    # Create genes hash table
    genes = nodes[nodes[1] == 'GENE'].copy()
    genes.loc[:, 2] = genes[0].apply(lambda x: int(x.split(':')[1]))
    genes.columns = ['id', 'type', 'GeneID']
    genes = genes.merge(gene_info, on='GeneID')
    genes_hash_table = genes[genes['id'].isin(terms_all)][['id', 'Symbol']]
    
    # Filter nodes to terms present in the knowledge graph
    nodes = nodes[nodes[0].isin(terms_all)]
    
    # Create term hash table (node_id, node_id, node_type)
    term_hash_table = nodes[[0, 0, 1]]
    
    # Save results
    save_results(nodes[0].tolist(), genes_hash_table, term_hash_table, should_downsample)

def process_primekg(terms_all, should_downsample):
    """
    Process PrimeKG dataset to extract genes and terms.
    
    Uses node names and types from the PrimeKG graph structure.
    """
    print("Processing PrimeKG data...")
    primekg_path = f'{S3_BASE_PATH}/primekg'
    
    # Load graph nodes
    nodes = pd.read_csv(f'{primekg_path}/nodes.csv')
    
    # Create genes hash table (gene/protein nodes only)
    genes_hash_table = nodes[
        (nodes['node_name'].isin(terms_all)) &
        (nodes['node_type'] == 'gene/protein')
    ][['node_name', 'node_name']]
    
    # Filter nodes to terms present in the knowledge graph
    nodes = nodes[nodes['node_name'].isin(terms_all)]
    
    # Create term hash table (node_name, node_name, node_type)
    term_hash_table = nodes[['node_name', 'node_name', 'node_type']]
    
    # Save results
    save_results(nodes['node_name'].tolist(), genes_hash_table, term_hash_table, should_downsample)


def save_results(term_ids, genes_hash_table, term_hash_table, should_downsample):
    """
    Save extracted data to CSV files.
    
    Args:
        term_ids: List of term IDs to save
        genes_hash_table: DataFrame with gene IDs and names
        term_hash_table: DataFrame with term IDs, names, and types
        should_downsample: Whether to downsample terms to DOWNSAMPLE_SIZE
    """
    print("Saving results...")
    
    # Save genes hash table
    genes_hash_table.to_csv('genes_hash_table.csv', index=False, header=False)
    print(f"  - genes_hash_table.csv: {len(genes_hash_table)} genes")
    
    # Save terms hash table
    term_hash_table.to_csv('terms_hash_table.csv', index=False, header=False)
    print(f"  - terms_hash_table.csv: {len(term_hash_table)} terms")
    
    # Save terms (with optional downsampling)
    terms_df = pd.DataFrame(term_ids)
    if should_downsample and len(terms_df) > DOWNSAMPLE_SIZE:
        terms_df = terms_df.sample(DOWNSAMPLE_SIZE, random_state=RANDOM_SEED)
        print(f"  - terms.csv: {len(terms_df)} terms (downsampled from {len(term_ids)})")
    else:
        print(f"  - terms.csv: {len(terms_df)} terms")
    
    terms_df.to_csv('terms.csv', index=False, header=False)


def main():
    """Main execution function."""
    # Set random seeds for reproducibility
    set_random_seeds()
    
    # Validate arguments
    dataset, should_downsample = validate_arguments()
    
    print(f"Dataset: {dataset}")
    print(f"Downsample: {should_downsample}")
    print()
    
    # Load dataset and save triples
    kg = load_dataset(dataset)
    terms_all = save_triples(kg)
    
    print(f"Total unique entities: {len(terms_all)}")
    print()
    
    # Process dataset-specific extraction
    dataset_processors = {
        'hetionet': process_hetionet,
        'biokg': process_biokg,
        'openbiolink': process_openbiolink,
        'primekg': process_primekg
    }
    
    processor = dataset_processors[dataset]
    processor(terms_all, should_downsample)
    
    print()
    print("Processing complete!")


if __name__ == '__main__':
    main()
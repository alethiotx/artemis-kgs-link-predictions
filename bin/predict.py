#!/usr/bin/env python3
"""
Knowledge Graph Target Prediction Script

This script generates gene-term association predictions using trained PyKEEN models
on knowledge graph embeddings. It handles dataset-specific relation types and
prediction directions for Hetionet, BioKG, OpenBioLink, and PrimeKG.

Usage:
    ./predict.py <dataset> <term> <terms_hash_table> <genes_hash_table> <triples> <model>

Arguments:
    dataset: Dataset name (hetionet, biokg, openbiolink, primekg)
    term: Term to generate predictions for
    terms_hash_table: CSV mapping term IDs to names and types
    genes_hash_table: CSV mapping gene IDs to names
    triples: NumPy array of knowledge graph triples
    model: Trained PyKEEN model file

Output:
    <term>_predictions.csv: Prediction scores for the given term

Dependencies:
    - pykeen
    - numpy
    - pandas
    - torch
"""

import sys
import re
import random
import numpy as np
import pandas as pd
import torch
from typing import List, Tuple, Optional
from pandas import DataFrame, concat
import pykeen
from pykeen.triples import TriplesFactory
from pykeen import predict

# ============================================================================
# Configuration Constants
# ============================================================================

# Random seed for reproducibility
RANDOM_SEED = 42

# ============================================================================
# Relation Mapping Configurations
# ============================================================================

# Hetionet: Maps term types to relation labels and prediction directions
HETIONET_RELATIONS = {
    'Biological Process': ('GpBP', 'forward'),  # Gene participates in Biological Process
    'Gene': ('GiG', 'forward'),                  # Gene interacts with Gene
    'Pathway': ('GpPW', 'forward'),              # Gene participates in Pathway
    'Cellular Component': ('GpCC', 'forward'),   # Gene participates in Cellular Component
    'Molecular Function': ('GpMF', 'forward'),   # Gene participates in Molecular Function
    'Anatomy': ('AeG', 'reverse'),               # Anatomy expresses Gene
    'Compound': ('CbG', 'reverse'),              # Compound binds Gene
    'Disease': ('DaG', 'reverse')                # Disease associates with Gene
}

# BioKG: Maps term types to relation labels and prediction directions
BIOKG_RELATIONS = {
    'Disease': ('PROTEIN_DISEASE_ASSOCIATION', 'forward'),
    'Pathway': ('PROTEIN_PATHWAY_ASSOCIATION', 'forward'),
    'Protein': ('PPI', 'forward'),               # Protein-Protein Interaction
    'Complex': ('MEMBER_OF_COMPLEX', 'forward'),
    'Genetic Disorder': ('RELATED_GENETIC_DISORDER', 'forward'),
    'Drug': ('DPI', 'reverse')                   # Drug-Protein Interaction
}

# OpenBioLink: Maps term types to relation labels
OPENBIOLINK_RELATIONS = {
    'ANATOMY': 'GENE_EXPRESSED_ANATOMY',
    'GO': 'GENE_GO',
    'DRUG': 'GENE_DRUG',
    'PHENOTYPE': 'GENE_PHENOTYPE',
    'PATHWAY': 'GENE_PATHWAY',
    'DIS': 'GENE_DIS',
    'GENE': 'GENE_GENE'
}

# PrimeKG: Maps term types to relation labels
PRIMEKG_RELATIONS = {
    'anatomy': 'anatomy_protein_present',
    'biological_process': 'bioprocess_protein',
    'cellular_component': 'cellcomp_protein',
    'disease': 'disease_protein',
    'drug': 'drug_protein',
    'exposure': 'exposure_protein',
    'molecular_function': 'molfunc_protein',
    'pathway': 'pathway_protein',
    'effect/phenotype': 'phenotype_protein',
    'gene/protein': 'protein_protein'
}


def generate_predictions(
    model: pykeen.models,
    triples_factory: TriplesFactory,
    heads: List[str],
    tails: List[str]
) -> DataFrame:
    """
    Generate predictions for given head-tail entity combinations.
    
    Iterates over all head and tail entity combinations, uses the model to predict
    targets, and aggregates results into a single DataFrame.
    
    Args:
        model: Trained PyKEEN model for predictions
        triples_factory: TriplesFactory containing the knowledge graph
        heads: List of head entity identifiers
        tails: List of tail entity identifiers
    
    Returns:
        DataFrame with predictions including head, tail, relation, and scores
    """
    predictions_all = DataFrame()
    
    for head in heads:
        for tail in tails:
            predictions = predict.predict_target(
                model=model,
                head=head,
                tail=tail,
                triples_factory=triples_factory
            ).df
            predictions['head'] = head
            predictions['tail'] = tail
            predictions_all = concat([predictions_all, predictions], ignore_index=True)
    
    return predictions_all


def set_random_seeds(seed: int = RANDOM_SEED):
    """
    Set random seeds for reproducibility across all libraries.
    
    Args:
        seed: Random seed value (default: RANDOM_SEED)
    """
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    
    # Enable deterministic operations in PyTorch
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    
    # Set environment variable for additional reproducibility
    import os
    os.environ['PYTHONHASHSEED'] = str(seed)


def validate_arguments():
    """Validate command line arguments."""
    if len(sys.argv) != 7:
        print("Error: Invalid number of arguments")
        print("Usage: ./predict.py <dataset> <term> <terms_hash_table> <genes_hash_table> <triples> <model>")
        sys.exit(1)
    
    dataset = sys.argv[1].lower()
    if dataset not in ['hetionet', 'biokg', 'openbiolink', 'primekg']:
        print(f"Error: Unsupported dataset '{dataset}'")
        print("Supported datasets: hetionet, biokg, openbiolink, primekg")
        sys.exit(1)
    
    return dataset


def load_data():
    """
    Load all required data files.
    
    Returns:
        Tuple of (dataset, term, terms_hash_table, genes, triples, model, triples_factory)
    """
    print("Loading data...")
    
    dataset = sys.argv[1].lower()
    term = sys.argv[2]
    
    # Load hash tables
    terms_hash_table = pd.read_csv(sys.argv[3], header=None)
    genes_hash_table = pd.read_csv(sys.argv[4], header=None)
    genes = genes_hash_table[0].tolist()
    
    # Load triples and create factory
    triples = np.load(sys.argv[5])
    triples_factory = TriplesFactory.from_labeled_triples(triples)
    
    # Load model (weights_only=False needed for full PyKEEN models with class definitions)
    model = torch.load(sys.argv[6], map_location=torch.device('cpu'), weights_only=False)
    
    print(f"  Dataset: {dataset}")
    print(f"  Term: {term}")
    print(f"  Genes: {len(genes)}")
    print(f"  Triples: {len(triples)}")
    
    return dataset, term, terms_hash_table, genes, model, triples_factory

def get_term_type(term: str, terms_hash_table: DataFrame, dataset: str) -> str:
    """
    Extract the term type from the term identifier.
    
    Args:
        term: Term identifier
        terms_hash_table: DataFrame mapping terms to types
        dataset: Dataset name
    
    Returns:
        Term type string
    """
    if dataset == 'hetionet':
        # Hetionet terms use prefix format: "Type::Name"
        return term.split('::')[0]
    else:
        # Other datasets store type in hash table
        term_row = terms_hash_table[terms_hash_table[0] == term]
        if len(term_row) == 0:
            raise ValueError(f"Term '{term}' not found in hash table")
        
        # Column index depends on dataset
        type_col = 1 if dataset == 'biokg' else 2
        return term_row[type_col].values[0]


def process_hetionet(
    term: str,
    genes: List[str],
    model: pykeen.models,
    triples_factory: TriplesFactory
) -> DataFrame:
    """
    Generate predictions for Hetionet dataset.
    
    Args:
        term: Term to predict for
        genes: List of gene identifiers
        model: Trained PyKEEN model
        triples_factory: Knowledge graph triples factory
    
    Returns:
        DataFrame with predictions indexed by gene
    """
    term_type = term.split('::')[0]
    
    if term_type not in HETIONET_RELATIONS:
        raise ValueError(f"Unknown Hetionet term type: {term_type}")
    
    relation_label, direction = HETIONET_RELATIONS[term_type]
    
    if direction == 'forward':
        # Gene -> Term relationships
        predictions = generate_predictions(model, triples_factory, heads=genes, tails=[term])
        predictions = predictions[predictions['relation_label'] == relation_label]
        predictions.index = predictions['head']
    else:
        # Term -> Gene relationships
        predictions = generate_predictions(model, triples_factory, heads=[term], tails=genes)
        predictions = predictions[predictions['relation_label'] == relation_label]
        predictions.index = predictions['tail']
    
    return predictions


def process_biokg(
    term: str,
    genes: List[str],
    terms_hash_table: DataFrame,
    model: pykeen.models,
    triples_factory: TriplesFactory
) -> DataFrame:
    """
    Generate predictions for BioKG dataset.
    
    Args:
        term: Term to predict for
        genes: List of gene identifiers
        terms_hash_table: DataFrame mapping terms to types
        model: Trained PyKEEN model
        triples_factory: Knowledge graph triples factory
    
    Returns:
        DataFrame with predictions indexed by gene
    """
    term_type = get_term_type(term, terms_hash_table, 'biokg')
    
    if term_type not in BIOKG_RELATIONS:
        raise ValueError(f"Unknown BioKG term type: {term_type}")
    
    relation_label, direction = BIOKG_RELATIONS[term_type]
    
    if direction == 'forward':
        # Gene -> Term relationships
        predictions = generate_predictions(model, triples_factory, heads=genes, tails=[term])
        predictions = predictions[predictions['relation_label'] == relation_label]
        predictions.index = predictions['head']
    else:
        # Term -> Gene relationships (e.g., Drug-Protein Interaction)
        predictions = generate_predictions(model, triples_factory, heads=[term], tails=genes)
        predictions = predictions[predictions['relation_label'] == relation_label]
        predictions.index = predictions['tail']
    
    return predictions

def process_openbiolink(
    term: str,
    genes: List[str],
    terms_hash_table: DataFrame,
    model: pykeen.models,
    triples_factory: TriplesFactory
) -> DataFrame:
    """
    Generate predictions for OpenBioLink dataset.
    
    All OpenBioLink predictions are Gene -> Term (forward direction).
    
    Args:
        term: Term to predict for
        genes: List of gene identifiers
        terms_hash_table: DataFrame mapping terms to types
        model: Trained PyKEEN model
        triples_factory: Knowledge graph triples factory
    
    Returns:
        DataFrame with predictions indexed by gene
    """
    term_type = get_term_type(term, terms_hash_table, 'openbiolink')
    
    if term_type not in OPENBIOLINK_RELATIONS:
        raise ValueError(f"Unknown OpenBioLink term type: {term_type}")
    
    relation_label = OPENBIOLINK_RELATIONS[term_type]
    
    # All OpenBioLink predictions are Gene -> Term
    predictions = generate_predictions(model, triples_factory, heads=genes, tails=[term])
    predictions = predictions[predictions['relation_label'] == relation_label]
    predictions.index = predictions['head']
    
    return predictions


def process_primekg(
    term: str,
    genes: List[str],
    terms_hash_table: DataFrame,
    model: pykeen.models,
    triples_factory: TriplesFactory
) -> DataFrame:
    """
    Generate predictions for PrimeKG dataset.
    
    All PrimeKG predictions are Term -> Gene (reverse direction).
    
    Args:
        term: Term to predict for
        genes: List of gene identifiers
        terms_hash_table: DataFrame mapping terms to types
        model: Trained PyKEEN model
        triples_factory: Knowledge graph triples factory
    
    Returns:
        DataFrame with predictions indexed by gene
    """
    term_type = get_term_type(term, terms_hash_table, 'primekg')
    
    if term_type not in PRIMEKG_RELATIONS:
        raise ValueError(f"Unknown PrimeKG term type: {term_type}")
    
    relation_label = PRIMEKG_RELATIONS[term_type]
    
    # All PrimeKG predictions are Term -> Gene
    predictions = generate_predictions(model, triples_factory, heads=[term], tails=genes)
    predictions = predictions[predictions['relation_label'] == relation_label]
    predictions.index = predictions['tail']
    
    return predictions


def save_predictions(predictions: DataFrame, term: str):
    """
    Save predictions to a CSV file.
    
    Args:
        predictions: DataFrame with predictions indexed by gene
        term: Term identifier for filename generation
    """
    # Extract scores and prepare output
    output = predictions[['score']]
    output = output.rename(columns={'score': term})
    
    # Create safe filename by replacing special characters
    filename = re.sub(r'[^\w_.-]', '_', term)
    
    # Transpose and save (genes as columns, term as row)
    output.T.round(1).to_csv(f'{filename}_predictions.csv', header=False)
    
    print(f"Saved predictions to {filename}_predictions.csv")
    print(f"  Predictions: {len(predictions)}")


def main():
    """Main execution function."""
    # Set random seeds for reproducibility
    set_random_seeds()
    
    # Validate arguments
    validate_arguments()
    
    # Load data
    dataset, term, terms_hash_table, genes, model, triples_factory = load_data()
    
    print()
    print("Generating predictions...")
    
    # Process dataset-specific predictions
    try:
        if dataset == 'hetionet':
            predictions = process_hetionet(term, genes, model, triples_factory)
        elif dataset == 'biokg':
            predictions = process_biokg(term, genes, terms_hash_table, model, triples_factory)
        elif dataset == 'openbiolink':
            predictions = process_openbiolink(term, genes, terms_hash_table, model, triples_factory)
        elif dataset == 'primekg':
            predictions = process_primekg(term, genes, terms_hash_table, model, triples_factory)
        else:
            raise ValueError(f"Unsupported dataset: {dataset}")
        
        # Save results
        print()
        save_predictions(predictions, term)
        print()
        print("Processing complete!")
        
    except Exception as e:
        print(f"Error during prediction: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()

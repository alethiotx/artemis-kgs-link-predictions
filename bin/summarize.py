#!/usr/bin/env python3
"""
Prediction Aggregation and Annotation Script

This script aggregates predictions from multiple terms, maps gene and term identifiers
to human-readable names using hash tables, and outputs a final prediction matrix.

The script performs the following operations:
1. Loads combined predictions from all terms
2. Transposes data to genes x terms format
3. Maps gene IDs to gene names
4. Maps term IDs to term names/descriptions
5. Sorts columns alphabetically for consistency
6. Saves final annotated prediction matrix

Usage:
    ./summarize.py <metafile> <genes_hash_table> <terms_hash_table>

Arguments:
    metafile: Path to metafile containing paths to all prediction CSVs
    genes_hash_table: CSV mapping gene IDs to gene names (2 columns: ID, Name)
    terms_hash_table: CSV mapping term IDs to term names (2+ columns: ID, Name, [Type])

Output:
    predictions.csv: Final prediction matrix with gene names as rows and term names as columns

Dependencies:
    - pandas
"""

import sys
from pathlib import Path
import pandas as pd
from typing import Dict, Tuple


def validate_arguments() -> None:
    """Validate command line arguments."""
    if len(sys.argv) != 4:
        print("Error: Invalid number of arguments")
        print("Usage: ./summarize.py <metafile> <genes_hash_table> <terms_hash_table>")
        sys.exit(1)
    
    # Check if input files exist
    for i, arg_name in enumerate(['metafile', 'genes_hash_table', 'terms_hash_table'], start=1):
        filepath = Path(sys.argv[i])
        if not filepath.exists():
            print(f"Error: {arg_name} file not found: {filepath}")
            sys.exit(1)


def load_data() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load all required data files.
    
    Returns:
        Tuple of (combined_predictions, genes_hash_table, terms_hash_table)
    """
    print("Loading data...")
    
    # Load hash tables
    genes_hash_table = pd.read_csv(sys.argv[2], index_col=0, header=None)
    terms_hash_table = pd.read_csv(sys.argv[3], index_col=0, header=None)
    
    print(f"  Genes: {len(genes_hash_table)}")
    print(f"  Terms: {len(terms_hash_table)}")
    
    # Load combined predictions
    combined = pd.read_csv(sys.argv[1], index_col=0, header=None)
    print(f"  Predictions shape (terms x genes): {combined.shape}")
    
    # Transpose to genes x terms format
    combined = combined.T
    combined.index = genes_hash_table.index.tolist()
    print(f"  Predictions shape (genes x terms): {combined.shape}")
    
    return combined, genes_hash_table, terms_hash_table


def process_predictions(
    predictions: pd.DataFrame,
    genes_hash_table: pd.DataFrame,
    terms_hash_table: pd.DataFrame
) -> pd.DataFrame:
    """
    Process predictions by mapping IDs to names and finalizing output.
    
    Args:
        predictions: Raw prediction DataFrame with IDs
        genes_hash_table: DataFrame mapping gene IDs to names
        terms_hash_table: DataFrame mapping term IDs to names
    
    Returns:
        Processed DataFrame with human-readable names, sorted and deduplicated
    """
    print("\nProcessing predictions...")
    
    # Create ID to name mappings
    gene_mapping = dict(zip(genes_hash_table.index, genes_hash_table[1]))
    term_mapping = dict(zip(terms_hash_table.index, terms_hash_table[1]))
    
    print(f"  Gene mappings: {len(gene_mapping)}")
    print(f"  Term mappings: {len(term_mapping)}")
    
    # Map IDs to names
    predictions = predictions.rename(index=gene_mapping, columns=term_mapping)
    
    # Sort columns alphabetically
    predictions = predictions.sort_index(axis=1)
    
    # Remove duplicates
    predictions = predictions.loc[:, ~predictions.columns.duplicated()]
    predictions = predictions.loc[~predictions.index.duplicated(), :]
    
    print(f"  Final shape: {predictions.shape} ({len(predictions)} genes × {len(predictions.columns)} terms)")
    
    return predictions


def save_predictions(predictions: pd.DataFrame, output_file: str = 'predictions.csv') -> None:
    """
    Save predictions to CSV file.
    
    Args:
        predictions: Final prediction DataFrame
        output_file: Output filename
    """
    print(f"\nSaving predictions to {output_file}...")
    
    predictions.to_csv(output_file)
    
    # Calculate file size
    file_size = Path(output_file).stat().st_size
    size_mb = file_size / (1024 * 1024)
    
    print(f"  File size: {size_mb:.2f} MB")
    print(f"  Total predictions: {predictions.size:,}")


def main():
    """Main execution function."""
    print("="*60)
    print("Prediction Summarization Pipeline")
    print("="*60)
    print()
    
    # Validate arguments
    validate_arguments()
    
    # Load all data
    combined, genes_hash_table, terms_hash_table = load_data()
    
    # Process and annotate predictions
    predictions = process_predictions(combined, genes_hash_table, terms_hash_table)
    
    # Save results
    save_predictions(predictions)
    
    print()
    print("="*60)
    print("Summarization complete!")
    print("="*60)


if __name__ == '__main__':
    main()

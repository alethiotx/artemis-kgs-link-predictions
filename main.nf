#!/usr/bin/env nextflow

/*
 * Knowledge Graph-Based Target Prediction Pipeline
 * ==============================================================
 * 
 * This workflow performs target prediction using knowledge graph embeddings
 * across multiple biomedical datasets (Hetionet, BioKG, OpenBioLink, PrimeKG).
 * 
 * Workflow steps:
 *   1. Prepare: Extract terms, triples, and metadata from knowledge graphs
 *   2. Predict: Generate predictions for each term using trained models
 *   3. Summarize: Aggregate and annotate all predictions
 * 
 * Parameters:
 *   - params.dataset: Dataset name (hetionet, biokg, openbiolink, primekg)
 *   - params.sample: Whether to downsample terms (true/false)
 *   - params.model: Path to trained prediction model
 *   - params.outdir: Output directory for results
 *   - params.env: Environment mode (test for limited runs)
 */

nextflow.enable.dsl=2

// ============================================================================
// Process Definitions
// ============================================================================

/*
 * Process: prepare
 * ----------------
 * Prepares knowledge graph data by extracting terms, triples, and metadata.
 * 
 * Outputs:
 *   - terms.csv: List of extracted terms for prediction
 *   - triples.npy: Concatenated triples from train/test/validation sets
 *   - genes_hash_table.csv: Gene IDs with human-readable names
 *   - terms_hash_table.csv: Term IDs with names and types
 */
process prepare {
  label 'process_single'

  publishDir "${params.outdir}/${params.dataset}/prepare", mode: 'copy'

  output:
    path 'terms.csv', emit: terms
    path 'triples.npy', emit: triples
    path 'genes_hash_table.csv', emit: genes_hash_table
    path 'terms_hash_table.csv', emit: terms_hash_table
  
  script:
  """
  prepare.py ${params.dataset} ${params.sample}
  """
}

/*
 * Process: predict
 * ----------------
 * Generates predictions for a single term using the trained model.
 * This process runs in parallel for each term extracted in the prepare step.
 * 
 * Inputs:
 *   - dataset: Name of the dataset being processed
 *   - term: Individual term to generate predictions for
 *   - terms_hash_table: Mapping of term IDs to names and types
 *   - genes_hash_table: Mapping of gene IDs to names
 *   - triples: All knowledge graph triples
 *   - model: Trained prediction model
 * 
 * Outputs:
 *   - *_predictions.csv: Predictions for the given term
 */
process predict {
  tag "${term}"
  label 'process_single'

  publishDir "${params.outdir}/${params.dataset}/predict", mode: 'copy'

  input:
    val dataset
    val term
    path terms_hash_table
    path genes_hash_table
    path triples
    path model

  output:
    path '*_predictions.csv', emit: predictions
  
  script:
  """
  predict.py ${dataset} "${term}" ${terms_hash_table} ${genes_hash_table} ${triples} ${model}
  """
}

/*
 * Process: summarize
 * ------------------
 * Aggregates predictions from all terms and annotates with gene/term names.
 * 
 * Inputs:
 *   - predictions: Collected prediction files from all predict processes
 *   - genes_hash_table: Mapping of gene IDs to names
 *   - terms_hash_table: Mapping of term IDs to names and types
 * 
 * Outputs:
 *   - predictions.csv: Final aggregated and annotated predictions
 */
process summarize {
  label 'process_low'

  publishDir "${params.outdir}/${params.dataset}/summarize", mode: 'copy'

  input:
    path predictions
    path genes_hash_table
    path terms_hash_table

  output:
    path 'predictions.csv'

  script:
  """
  summarize.py ${predictions} ${genes_hash_table} ${terms_hash_table}
  """
}

// ============================================================================
// Main Workflow
// ============================================================================

workflow {
  
  // Step 1: Prepare knowledge graph data
  prepare()

  // Step 2: Parse terms and prepare for parallel processing
  terms = prepare
    .out
    .terms
    .splitCsv(quote: '"')
    .map { it[0] }

  // Limit terms in test environment for faster validation
  if (params.env == 'test') {
    terms = terms.take(2)
  }

  // Step 3: Generate predictions for each term in parallel
  predict(
    params.dataset,
    terms, 
    prepare.out.terms_hash_table,
    prepare.out.genes_hash_table,
    prepare.out.triples,
    params.model
  )

  // Step 4: Aggregate and summarize all predictions
  summarize(
    predict.out.predictions.collectFile(name: 'metafile.txt'),
    prepare.out.genes_hash_table,
    prepare.out.terms_hash_table
  )
}

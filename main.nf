#!/usr/bin/env nextflow

/*
 * Knowledge Graph-Based Target Prediction Pipeline
 * ==============================================================
 * 
 * This workflow performs target prediction using knowledge graph embeddings
 * across multiple biomedical datasets (Hetionet, BioKG, OpenBioLink, PrimeKG)
 * and multiple embedding models (ComplEx, DistMult, RotatE, TransE).
 * 
 * By default, runs all 4 datasets × 4 embeddings (16 combinations).
 * Can be restricted to specific datasets/embeddings via parameters.
 * 
 * Workflow steps:
 *   1. Prepare: Extract terms, triples, and metadata from knowledge graphs
 *   2. Predict: Generate predictions for each term using trained models
 *   3. Summarize: Aggregate and annotate all predictions
 * 
 * Parameters:
 *   - params.datasets: List of dataset names (default: all 4)
 *   - params.embeddings: List of embedding model names (default: all 4)
 *   - params.sample: Whether to downsample terms (true/false)
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
 * Prepares knowledge graph data by extracting terms and metadata,
 * filtered to entities present in the upstream model's vocabulary.
 * 
 * Inputs:
 *   - dataset: Name of the dataset being processed
 *   - embedding: Name of the embedding model
 *   - training_triples: Training triples directory from upstream embeddings pipeline
 * 
 * Outputs:
 *   - terms.csv: List of extracted terms for prediction
 *   - genes_hash_table.csv: Gene IDs with human-readable names
 *   - terms_hash_table.csv: Term IDs with names and types
 */
process prepare {
  tag "${dataset}/${embedding}"
  label 'process_single'
  cpus { dataset == 'primekg' ? 2 : 1 }
  memory { dataset == 'primekg' ? 48.GB * task.attempt : 4.GB * task.attempt }

  publishDir "${params.outdir}/${dataset}/${embedding}/prepare", mode: 'copy'

  input:
    val dataset
    val embedding
    path training_triples

  output:
    tuple val(dataset), val(embedding), path('terms.csv'), path('genes_hash_table.csv'), path('terms_hash_table.csv'), emit: results
  
  script:
  """
  prepare.py ${dataset} ${params.sample} ${training_triples}
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
 *   - embedding: Name of the embedding model
 *   - term: Individual term to generate predictions for
 *   - terms_hash_table: Mapping of term IDs to names and types
 *   - genes_hash_table: Mapping of gene IDs to names
 *   - training_triples: Training triples directory from upstream embeddings pipeline
 *   - model: Trained prediction model
 * 
 * Outputs:
 *   - *_predictions.csv: Predictions for the given term
 */
process predict {
  tag "${dataset}/${embedding}/${term}"
  label 'process_single'
  cpus { dataset == 'primekg' ? 2 : 1 }
  memory { dataset == 'primekg' ? 48.GB * task.attempt : 4.GB * task.attempt }

  publishDir "${params.outdir}/${dataset}/${embedding}/predict", mode: 'copy'

  input:
    val dataset
    val embedding
    val term
    path terms_hash_table
    path genes_hash_table
    path training_triples
    path model

  output:
    tuple val(dataset), val(embedding), path('*_predictions.csv'), emit: predictions
  
  script:
  """
  predict.py ${dataset} "${term}" ${terms_hash_table} ${genes_hash_table} ${training_triples} ${model}
  """
}

/*
 * Process: summarize
 * ------------------
 * Aggregates predictions from all terms and annotates with gene/term names.
 * 
 * Inputs:
 *   - dataset: Name of the dataset being processed
 *   - embedding: Name of the embedding model
 *   - predictions: Collected prediction files from all predict processes
 *   - genes_hash_table: Mapping of gene IDs to names
 *   - terms_hash_table: Mapping of term IDs to names and types
 * 
 * Outputs:
 *   - predictions.csv: Final aggregated and annotated predictions
 *   - predictions.parquet: Final aggregated and annotated predictions in Parquet format
 */
process summarize {
  tag "${dataset}/${embedding}"
  label 'process_low'

  publishDir "${params.outdir}/${dataset}/${embedding}/summarize", mode: 'copy'

  input:
    val dataset
    val embedding
    path predictions
    path genes_hash_table
    path terms_hash_table

  output:
    path 'predictions.csv'
    path 'predictions.parquet'

  script:
  """
  summarize.py ${predictions} ${genes_hash_table} ${terms_hash_table}
  """
}

// ============================================================================
// Main Workflow
// ============================================================================

workflow {

  // Build channel of all dataset × embedding combinations
  datasets_ch = Channel.from(params.datasets)
  embeddings_ch = Channel.from(params.embeddings)
  combos = datasets_ch.combine(embeddings_ch)
    .map { dataset, embedding ->
      def s3_base = "${params.s3_embeddings_base}/${dataset}/${embedding}"
      def training_triples = "${s3_base}/training_triples"
      def model = "${s3_base}/trained_model.pkl"
      return [dataset, embedding, training_triples, model]
    }

  // Step 1: Prepare knowledge graph data for each combination
  prepare(
    combos.map { it[0] },  // dataset
    combos.map { it[1] },  // embedding
    combos.map { it[2] }   // training_triples
  )

  // Step 2: Expand prepare results into per-term rows for parallel prediction
  terms_ch = prepare.out.results
    .flatMap { dataset, embedding, terms_csv, genes_hash, terms_hash ->
      def s3_base = "${params.s3_embeddings_base}/${dataset}/${embedding}"
      def training_triples = "${s3_base}/training_triples"
      def model = "${s3_base}/trained_model.pkl"
      def lines = terms_csv.text.trim().split('\n')
      def term_list = lines.collect { it.trim().replaceAll('^"|"$', '') }
      if (params.env == 'test') {
        term_list = term_list.take(2)
      }
      return term_list.collect { term ->
        [dataset, embedding, term, terms_hash, genes_hash, training_triples, model]
      }
    }

  // Step 3: Generate predictions for each term in parallel
  predict(
    terms_ch.map { it[0] },  // dataset
    terms_ch.map { it[1] },  // embedding
    terms_ch.map { it[2] },  // term
    terms_ch.map { it[3] },  // terms_hash_table
    terms_ch.map { it[4] },  // genes_hash_table
    terms_ch.map { it[5] },  // training_triples
    terms_ch.map { it[6] }   // model
  )

  // Step 4: Group predictions by dataset/embedding and summarize
  // Concatenate prediction CSVs into a metafile per dataset/embedding
  metafiles = predict.out.predictions
    .collectFile(sort: true) { dataset, embedding, pred_file ->
      ["${dataset}___${embedding}.txt", pred_file.text]
    }
    .map { metafile ->
      def parts = metafile.baseName.replace('.txt', '').split('___')
      [parts[0], parts[1], metafile]
    }

  // Get genes/terms hash tables from prepare results
  prepare_lookup = prepare.out.results
    .map { dataset, embedding, terms_csv, genes_hash, terms_hash ->
      [dataset, embedding, genes_hash, terms_hash]
    }

  // Join metafiles with hash tables
  summarize_input = metafiles
    .join(prepare_lookup, by: [0, 1])

  summarize(
    summarize_input.map { it[0] },  // dataset
    summarize_input.map { it[1] },  // embedding
    summarize_input.map { it[2] },  // predictions (concatenated metafile)
    summarize_input.map { it[3] },  // genes_hash_table
    summarize_input.map { it[4] }   // terms_hash_table
  )
}

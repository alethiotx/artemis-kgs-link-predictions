#!/usr/bin/env nextflow

/*
 * Knowledge Graph-Based Target Prediction Pipeline
 * ==============================================================
 * 
 * This workflow performs target prediction using knowledge graph embeddings
 * for a single dataset across all 4 embedding models (ComplEx, DistMult,
 * RotatE, TransE).
 * 
 * Supported datasets: Hetionet, BioKG, OpenBioLink, PrimeKG
 * 
 * To run all 4 datasets, launch 4 separate pipeline runs via Seqera Platform.
 * 
 * Workflow steps:
 *   1. Prepare: Extract terms, triples, and metadata from knowledge graphs
 *   2. Predict: Generate predictions for each term using trained models
 *   3. Summarize: Aggregate and annotate all predictions
 * 
 * Parameters:
 *   - params.dataset: Dataset name (required)
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
 *   - batch_start: Starting index into terms.csv for this batch
 *   - batch_size: Number of terms to process in this batch
 *   - terms_csv: Full list of terms from prepare step
 *   - terms_hash_table: Mapping of term IDs to names and types
 *   - genes_hash_table: Mapping of gene IDs to names
 *   - training_triples: Training triples directory from upstream embeddings pipeline
 *   - model: Trained prediction model
 * 
 * Outputs:
 *   - *_predictions.csv: Predictions for each term in the batch
 */
process predict {
  tag "${dataset}/${embedding}/${batch_start}"
  label 'process_single'
  cpus { dataset == 'primekg' ? 2 : 1 }
  memory { dataset == 'primekg' ? 48.GB * task.attempt : 4.GB * task.attempt }

  publishDir "${params.outdir}/${dataset}/${embedding}/predict", mode: 'copy'

  input:
    val dataset
    val embedding
    val batch_start
    val batch_size
    path terms_csv
    path terms_hash_table
    path genes_hash_table
    path training_triples
    path model

  output:
    tuple val(dataset), val(embedding), path('*_predictions.csv'), emit: predictions
  
  script:
  """
  predict.py ${dataset} ${terms_csv} ${batch_start} ${batch_size} ${terms_hash_table} ${genes_hash_table} ${training_triples} ${model}
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

  // Validate required parameters
  if (!params.dataset) {
    error "ERROR: --dataset is required. Options: hetionet, biokg, openbiolink, primekg"
  }

  // Build channel of embeddings for the single dataset
  embeddings_ch = Channel.from(params.embeddings)
    .map { embedding ->
      def s3_base = "${params.s3_embeddings_base}/${params.dataset}/${embedding}"
      def training_triples = "${s3_base}/training_triples"
      def model = "${s3_base}/trained_model.pkl"
      return [params.dataset, embedding, training_triples, model]
    }

  // Step 1: Prepare knowledge graph data for each embedding
  prepare(
    embeddings_ch.map { it[0] },  // dataset
    embeddings_ch.map { it[1] },  // embedding
    embeddings_ch.map { it[2] }   // training_triples
  )

  // Step 2: Expand prepare results into batch indices for parallel prediction
  terms_ch = prepare.out.results
    .flatMap { dataset, embedding, terms_csv, genes_hash, terms_hash ->
      def s3_base = "${params.s3_embeddings_base}/${dataset}/${embedding}"
      def training_triples = "${s3_base}/training_triples"
      def model = "${s3_base}/trained_model.pkl"
      def lines = terms_csv.text.trim().split('\n')
      def term_count = lines.size()
      if (params.env == 'test') {
        term_count = Math.min(2, term_count)
      }
      // Generate batch start indices
      def batch_starts = (0..<term_count).step(params.batch_size)
      return batch_starts.collect { start ->
        [dataset, embedding, start, params.batch_size, terms_csv, terms_hash, genes_hash, training_triples, model]
      }
    }

  // Step 3: Generate predictions for each term batch in parallel
  predict(
    terms_ch.map { it[0] },  // dataset
    terms_ch.map { it[1] },  // embedding
    terms_ch.map { it[2] },  // batch_start
    terms_ch.map { it[3] },  // batch_size
    terms_ch.map { it[4] },  // terms_csv
    terms_ch.map { it[5] },  // terms_hash_table
    terms_ch.map { it[6] },  // genes_hash_table
    terms_ch.map { it[7] },  // training_triples
    terms_ch.map { it[8] }   // model
  )

  // Step 4: Group predictions by embedding and summarize
  metafiles = predict.out.predictions
    .transpose()
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

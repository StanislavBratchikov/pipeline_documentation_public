RNA-seq Analysis Pipeline
=========================

The RNA-seq pipeline provides comprehensive analysis from raw count data to functional enrichment.

Script Location
---------------

``DE_pipelines/RNAseq/run_pipeline_command_line.py``

Basic Usage
-----------

.. code-block:: bash

   # Complete analysis (preprocessing + differential expression)
   python DE_pipelines/RNAseq/run_pipeline_command_line.py \
       --analysis_type both \
       --feature_counts_pattern "*/featureCounts.txt" \
       --sample_metadata_path metadata/samples.csv \
       --preprocessing_output_dir results/preprocessing/ \
       --output_dir results/de_analysis/ \
       --condition_pairs control treatment \
       --verbose

   # Preprocessing only
   python DE_pipelines/RNAseq/run_pipeline_command_line.py \
       --analysis_type preprocessing \
       --feature_counts_pattern "*/featureCounts.txt" \
       --sample_metadata_path metadata/samples.csv \
       --preprocessing_output_dir results/preprocessing/

Key Parameters
--------------

Analysis Type
~~~~~~~~~~~~

* ``--analysis_type {preprocessing,de,both}`` - Type of analysis to run

Input Data
~~~~~~~~~

* ``--feature_counts_pattern`` - Glob pattern for featureCounts files
* ``--sample_metadata_path`` - Sample metadata file path
* ``--gene_info_path`` - Gene annotation file

Output Options
~~~~~~~~~~~~~

* ``--preprocessing_output_dir`` - Preprocessing results directory
* ``--output_dir`` - DE analysis results directory
* ``--h5ad_filename`` - Output AnnData filename

Differential Expression
~~~~~~~~~~~~~~~~~~~~~~

* ``--condition_pairs`` - Conditions to compare
* ``--logfc_threshold`` - Log2 fold change threshold (default: 0.25)
* ``--pval_threshold`` - P-value threshold (default: 0.05)

Enrichment Analysis
~~~~~~~~~~~~~~~~~

* ``--skip_enrichment`` - Skip enrichment analysis
* ``--enrich_databases`` - Custom database directory
* ``--min_size`` - Minimum gene set size (default: 15)
* ``--max_size`` - Maximum gene set size (default: 1000)

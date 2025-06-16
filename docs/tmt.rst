TMT Proteomics Analysis Pipeline
================================

The TMT pipeline provides comprehensive proteomics analysis with advanced normalization and batch correction.

Script Location
---------------

``DE_pipelines/TMT/tmt_command_line_script.py``

Basic Usage
-----------

.. code-block:: bash

   # Complete analysis with batch correction
   python DE_pipelines/TMT/tmt_command_line_script.py \
       --analysis_type both \
       --expression_data_path data/protein_expression.csv \
       --protein_metadata_path data/protein_metadata.csv \
       --sample_metadata_path data/sample_metadata.csv \
       --preprocessing_output_dir results/preprocessing/ \
       --de_output_dir results/de_analysis/ \
       --condition_pairs control treatment \
       --batch_correct \
       --batch_id_column Batch \
       --verbose

Key Parameters
--------------

Input Data
~~~~~~~~~

* ``--expression_data_path`` - Protein expression data file
* ``--protein_metadata_path`` - Protein metadata file  
* ``--sample_metadata_path`` - Sample metadata file

Normalization Options
~~~~~~~~~~~~~~~~~~~

* ``--skip_sum_normalize`` - Skip per-sample sum normalization
* ``--skip_peptide_normalize`` - Skip peptide count normalization
* ``--skip_mad_normalization`` - Skip MAD normalization
* ``--target_sum`` - Target sum for normalization (default: 1,000,000)

Batch Correction
~~~~~~~~~~~~~~

* ``--batch_correct`` - Enable batch correction
* ``--batch_id_column`` - Batch identifier column
* ``--mod_id_column`` - Model covariates column

Protein Mapping
~~~~~~~~~~~~~

* ``--identify_canonical_proteins`` - Identify canonical isoforms
* ``--apply_uniparc_mapping`` - Use UniParc database mapping
* ``--swiss_prot_ref_path`` - Swiss-Prot reference file

Differential Expression
~~~~~~~~~~~~~~~~~~~~~

* ``--logfc_threshold`` - Log2 fold change threshold (default: 0.1)
* ``--pval_threshold`` - P-value threshold (default: 0.05)
* ``--anndata_normalized_data_layer`` - Data layer for analysis

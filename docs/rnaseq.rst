RNA-seq Analysis Pipeline
=========================

The RNA-seq pipeline provides comprehensive analysis from raw count data to functional enrichment using DESeq2.

Script Location
---------------

``DE_pipelines/RNAseq/run_pipeline_command_line.py``

Input Data Requirements
-----------------------

Feature Count Files
~~~~~~~~~~~~~~~~~~

* **Format**: featureCounts output with ``Geneid`` column containing Ensembl IDs
* **Organization**: Files are expected to be stored in directories named the same as corresponding samples in the passed sample metadata
* **Pattern**: Specified by ``--feature_counts_pattern``

Sample Metadata
~~~~~~~~~~~~~~

* **Required column**: ``RNA-seq Ref`` must contain sample ids
* **Additional columns**: All sample related data such as biological conditions, annotation for PCA and batch information should be contained here
* **Format**: CSV file

Complete Parameter Reference
----------------------------

Basic Arguments
~~~~~~~~~~~~~~

``--analysis_type {preprocessing,de,both}``
    Type of analysis to run: preprocessing, differential expression (de), or both

``--preprocessing_config PREPROCESSING_CONFIG``
    Path to JSON preprocessing configuration file

``--de_config DE_CONFIG``
    Path to JSON de configuration file

``--verbose``
    Enable verbose output

Preprocessing Parameters
~~~~~~~~~~~~~~~~~~~~~~~

``--feature_counts_pattern FEATURE_COUNTS_PATTERN``
    Glob pattern to match feature count files with Geneid column containing EnsemblIDs. Files are expected to be stored in directories named the same as corresponding samples in the passed sample metadata.

``--sample_metadata_path SAMPLE_METADATA_PATH``
    Path to the sample metadata file. Pipeline will use this table as sample metadata table in anndata object. All sample related data such as biological conditions, annotation for PCA and batch information should be contained here. Column "RNA-seq Ref" must be contained sample ids.

``--preprocessing_output_dir PREPROCESSING_OUTPUT_DIR``
    Directory where preprocessing results will be saved

``--h5ad_filename H5AD_FILENAME``
    Filename for the h5ad output file (default: rna-seq_data.h5ad)

``--gene_info_path GENE_INFO_PATH``
    Path to the gene metadata file mapping ENSEMBL ids to Gene Symbols (default: ../DE_pipelines/../ref_files/human/ensemble_id_to_gene_name.tab.gz)

``--mitocarta_path MITOCARTA_PATH``
    Path to the MitoCarta data. "MitoCarta3.0_MitoPathways" column must contain OXPHOS complexes names in the format shown in defaults (default: ../DE_pipelines/../ref_files/human/human_genes_mitocarta3.0.csv.gz)

``--lncRNAs_and_prot_coding_genes_path LNCRNAS_AND_PROT_CODING_GENES_PATH``
    Path to file mapping EnsemblID to long non coding RNAs and protein coding transcripts. This mapping will be used for subsetting count data to only contain specified category of transcripts. Used in pair with transcript_types flag (default: ../DE_pipelines/../ref_files/human/lncRNAs_and_prot_coding_genes.csv.gz)

``--transcript_types TRANSCRIPT_TYPES [TRANSCRIPT_TYPES ...]``
    Subset anndata object used for downstream differential expression analysis to only contain specified transcript types. If "all_transcripts", do not subset anndata (default ['protein_coding', 'lncRNA'])

``--mean_threshold_for_expressed_genes MEAN_THRESHOLD_FOR_EXPRESSED_GENES``
    Threshold for mean gene levels used for identification of expressed vs not expressed genes.(default: 1.25)

``--max_mean MAX_MEAN``
    Maximum mean expression threshold for filtering highly variable genes (default: 8)

``--min_disp MIN_DISP``
    Minimum dispersion threshold for filtering highly variable genes (default: 1)

``--apply_sum_norm``
    Flag to apply per sample sum normalization on raw counts of data. Adds layer with normalized data to anndata object (default: False)

``--target_sum TARGET_SUM``
    Target sum per sample. Used if apply_sum_norm is specified (default: 10000000.0)

``--scale_max_value SCALE_MAX_VALUE``
    Maximum scaling value to use for sc.pp.scale function applied to normalized data. This layer will be used for PCA plots (default: 10)

``--pca_color PCA_COLOR [PCA_COLOR ...]``
    List of metadata column names to use as PCA scatter plots color categories(default: ['Condition', 'Replicate'])

``--pca_components PCA_COMPONENTS [PCA_COMPONENTS ...]``
    List of comma-separated pairs of PCA components (default: ['1,2', '2,3', '1,3'])

``--umap_plot``
    Generate UMAP plots (default: False)

``--n_neighbors N_NEIGHBORS``
    Number of nearest neighbors to compute UMAP on (default: 5)

``--n_pcs N_PCS``
    Number of principal components to use for UMAP construction (default: 9)

``--oxphos_complexes OXPHOS_COMPLEXES [OXPHOS_COMPLEXES ...]``
    List of strings of OXPHOS complexes names to generate heatmaps for (default: ['Complex I', 'Complex II', 'Complex III', 'Complex IV', 'Complex V'])

Differential Expression Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``--input_file INPUT_FILE``
    Path to the AnnData file for DE analysis

``--output_dir OUTPUT_DIR``
    Directory where DE results will be saved

``--skip_expressed_filter``
    Flag to skip using only expressed genes for DE analysis

``--design_factors DESIGN_FACTORS [DESIGN_FACTORS ...]``
    Sample metadata columns to include as design factors for DESeq2.

``--condition_pairs CONDITION_PAIRS [CONDITION_PAIRS ...]``
    Condition pairs for comparison (format: cond1 cond2 [cond3 cond4 ...])

``--logfc_threshold LOGFC_THRESHOLD``
    Log2 fold change threshold for DE analysis (default: 0.25)

``--pval_threshold PVAL_THRESHOLD``
    P-value threshold for DE analysis (default: 0.05)

Enrichment Analysis Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``--skip_enrichment``
    Skip enrichment analysis after DE (default: False)

``--enrich_databases ENRICH_DATABASES``
    Directory containing enrichment analysis databases (.gmt files) (default ../DE_pipelines/../ref_files/human/human_enrichr_databases/)

``--logfc_enrich LOGFC_ENRICH``
    Log2 fold change threshold for enrichment analysis. If not specified, logfc_threshold value is used (default: value passed to logfc_threshold)

``--pval_enrich PVAL_ENRICH``
    P-value threshold for enrichment analysis. If not specified, pval_threshold value is used (default: value passed to pval_threshold)

``--pval_enrich_column {padj,pvalue}``
    P-value column to use for enrichment analysis either "padj" or "pvalue" (default: "padj")

``--min_size MIN_SIZE``
    min_size parameter for gseapy.prerank function (default: 15)

``--max_size MAX_SIZE``
    max_size parameter for gseapy.prerank function (default: 1000)

``--permutation_num PERMUTATION_NUM``
    permutation_num parameter for gseapy.prerank function (default: 100)

Usage Examples
--------------

Complete Analysis::

    python DE_pipelines/RNAseq/run_pipeline_command_line.py \
        --analysis_type both \
        --feature_counts_pattern "data/*/featureCounts.txt" \
        --sample_metadata_path metadata/samples.csv \
        --preprocessing_output_dir results/preprocessing/ \
        --output_dir results/de_analysis/ \
        --condition_pairs control treatment \
        --verbose

Preprocessing Only::

    python DE_pipelines/RNAseq/run_pipeline_command_line.py \
        --analysis_type preprocessing \
        --feature_counts_pattern "data/*/featureCounts.txt" \
        --sample_metadata_path metadata/samples.csv \
        --preprocessing_output_dir results/preprocessing/ \
        --transcript_types protein_coding lncRNA \
        --apply_sum_norm \
        --umap_plot

Differential Expression Only::

    python DE_pipelines/RNAseq/run_pipeline_command_line.py \
        --analysis_type de \
        --input_file results/preprocessing/rna-seq_data.h5ad \
        --output_dir results/de_analysis/ \
        --condition_pairs control treatment diseased healthy \
        --logfc_threshold 0.5 \
        --pval_threshold 0.01

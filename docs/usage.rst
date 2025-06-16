Usage Overview
==============

Quick Start
-----------

Both RNA-seq and TMT pipelines support:

* Command-line execution with extensive parameters
* JSON configuration files for reproducible analyses
* Modular execution (preprocessing only, DE only, or both)

Basic Workflow
--------------

1. **Prepare your data**: Organize count files and metadata
2. **Run preprocessing**: Quality control and normalization
3. **Differential expression**: Statistical analysis between conditions
4. **Enrichment analysis**: Functional annotation of results

Command Line Tips
-----------------

* Always run scripts with ``-h`` flag first to see all available options
* Use ``--verbose`` flag for detailed output during analysis
* Reference files are automatically loaded from ``ref_files/`` directory
* Output files are saved in AnnData (.h5ad) format for interoperability

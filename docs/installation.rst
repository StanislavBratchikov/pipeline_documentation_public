Installation
============

Requirements
------------

* Python 3.8+
* R 4.0+ (for DESeq2)
* Required Python packages (see requirements.txt)

Quick Install
-------------

1. **Clone the repository**

   .. code-block:: bash

      git clone [your-repo-url]
      cd omics_downstream_pipeline

2. **Install dependencies**

   .. code-block:: bash

      pip install -r requirements.txt

3. **Required Python packages include:**

   * pandas, numpy, scipy
   * scanpy, anndata
   * matplotlib, seaborn
   * DESeq2 (via rpy2)
   * gseapy

Repository Structure
-------------------

.. code-block:: text

   omics_downstream_pipeline/
   ├── DE_pipelines/
   │   ├── RNAseq/              # RNA-seq differential expression analysis
   │   │   └── run_pipeline_command_line.py
   │   ├── TMT/                 # TMT proteomics analysis  
   │   │   └── tmt_command_line_script.py
   │   └── utilities/           # Shared utility functions
   ├── docs/                    # Documentation and images
   │   └── pipeline_overview.png
   ├── ref_files/               # Reference files and databases
   └── README.md               # This documentation

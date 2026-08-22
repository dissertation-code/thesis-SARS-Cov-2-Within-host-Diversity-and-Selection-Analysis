# Within-Host SARS-CoV-2 Evolution During Sotrovimab Treatment in Immunocompromised Patients: A Longitudinal Genomics Analysis of Viral Diversity, Selection, and Escape Mutations

This repository contains the computational workflows and R scripts used for the analysis presented in my dissertation investigating SARS-CoV-2 within-host genetic diversity, drug-associated selection, and viral load.

The repository includes the upstream `viralrecon` configuration and execution script, together with the R scripts used for downstream quality control, cohort construction, phylogenetic analysis, within-host diversity analysis, sotrovimab selection analysis, and viral load analysis.

## Repository structure

```text
.
├── README.md
├── .gitignore
│
├── viralrecon/
│   ├── viralrecon_parameters.yml
│   ├── viralrecon_config.yml
│   └── run_viralrecon.sh
│
└── scripts/
    ├── 01_QC_and_analysis_cohort_building.R
    ├── 02_baseline_phylogenetic_analysis.R
    ├── 03_within_host_diversity_analysis.R
    ├── 04_sotrovimab_selection_analysis.R
    └── 05_viral_load_analysis.R
```

### `viralrecon/`

Contains the configuration and execution files used for the `nf-core/viralrecon` workflow:

* `viralrecon_parameters.yml` — workflow parameters used for the viralrecon analysis.
* `viralrecon_config.yml` — configuration settings used during variant calling using iVar in the viralrecon workflow.
* `run_viralrecon.sh` — Bash script used to submit and execute the viralrecon workflow on the HPC environment.

### `scripts/`

Contains the R scripts used for downstream analysis. The scripts are numbered according to the general order of the analytical workflow:

1. **`01_QC_and_analysis_cohort_building.R`**
   Quality control and construction of the longitudinal analysis cohort.

2. **`02_baseline_phylogenetic_analysis.R`**
   Phylogenetic analysis of baseline samples.

3. **`03_within_host_diversity_analysis.R`**
   Analysis of within-host genetic diversity and variation.

4. **`04_sotrovimab_selection_analysis.R`**
   Analysis of genetic changes associated with sotrovimab exposure and selection.

5. **`05_viral_load_analysis.R`**
   Analysis of viral load and associated comparisons.
   

## Data availability

The underlying sequencing and clinical/sample-level data used in this study are **confidential and are therefore not included in this repository**.

The data dictionary describing the variables used in the analysis is provided in the **Appendices (Appendix C) of the dissertation**.

Consequently, the scripts in this repository document the computational and statistical analyses performed for the dissertation but cannot be executed from the repository alone without access to the corresponding input data.


## Requirements

The analyses were performed using:

* **R/RStudio** for downstream data processing, statistical analysis, and visualisation.
* **Nextflow** and the **nf-core/viralrecon** workflow for viral genome sequencing analysis.
* **Apptainer/Singularity** for containerised workflow execution.
* An HPC environment for execution of the viralrecon workflow.

The specific software versions and workflow parameters used for viralrecon are documented in the files within the `viralrecon/` directory.


## Reproducibility

The repository is intended to provide a record of the computational workflow and analysis code used in the dissertation.

Because the underlying data are confidential and are not distributed with this repository, complete reproduction of the reported results requires authorised access to the original input data.

Paths relating to local file locations, HPC environments, and input data may need to be modified before running the scripts in another environment.



## Contact

For questions regarding the analysis or repository, please contact the author of the associated dissertation.


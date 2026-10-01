# A clinician-assessed measure of vaginal epithelial health is associated with multiomic signatures in genitourinary syndrome of menopause
Github author: Lauren Mee

Analysis pipeline linking 16S rRNA microbiome and NMR metabolome data to two clinical measures of vaginal atrophy severity — **NGAT** (clinician-scored) and **DIVA** (patient-reported) — in a cohort of post-menopausal women with GSM. Control patients were included for descriptive comparisons.
Manuscript in process of being submitted and this README will be updated with DOI and authorship details once available.

# QuickStart

## Input

Available input here includes the pre-normalised NMR spectra data (`input/nmr/VANS_tampons_master_data_matrix_updated_Sep2024_missing_values_replaced.csv`) and patient metadata given by sample (`input/SampleMetadata.csv`).

Please note that this pipeline also requires **trimmed** FASTQ forward read files (source: [ENA], see manuscript for trimming parameters) and SILVA SSU training sets ([McLaren, M. R., & Callahan, B. J. (2021)](https://doi.org/10.5281/zenodo.4587955)). 
The pipeline will expect the following file naming conventions and input architecture:

```
input/
├── microbiome/
│   ├── trimmed/                         # primer-trimmed FASTQ (forward reads only)
|      ├── VAN[\d+]_R1_Trimmed.fastq.gz
│   ├── taxaDB/                          # SILVA v138.1 train set + species assignment
|      ├──  silva_nr99_v138.1_train_set.fa
|      ├──  silva_species_assignment_v138.1.fa
```

## Running analysis

Scripts **must be run in order** — each consumes checkpoint files written by the
previous ones. All scripts assume the **project root** as the working directory.

```r
# from the project root, in this order
source("scripts/00_Install.R")           # installs all dependencies
source("scripts/01_Micro_Processing.R")  # microbiome: DADA2 + taxonomy assignment
source("scripts/02_Micro_Analysis.R")    # microbiome: analysis
source("scripts/03_NMR_Analysis.R")      # metabolomic: analysis
source("scripts/04_NetworkBuild.R")      # integration: building networks
source("scripts/05_NetworkComparison.R") # integration: assessing networks
```

Note: Scripts 03 and 04 do not depend directly on script 02 and can be run in parallel with it if desired.

## Requirements

- **R >= 4.4** (developed on 4.5.2)
- ~4 GB RAM, ~3 GB free disk
- Dependencies (covered in script 00)





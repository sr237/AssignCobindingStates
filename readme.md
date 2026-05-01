# Assign Co-binding States (R Implementation)

### This repository contains an R implementation of the co-binding state assignment pipeline. The script processes DSMF BEDPE read data and assigns co-binding states based on methylation patterns.The workflow identifies binding events (e.g., transcription factor binding, nucleosome occupancy) and encodes co-binding relationships between primary and secondary peaks.

### Features
Pure R implementation 
Handles large .bedpe.gz genomic datasets
Robust to missing values and boundary issues
Produces identical structured output as the original Python pipeline

### Input Requirements
1. Main Input File (.bedpe.gz)
Contains genomic coordinates and read information
Format: BEDPE-like
Must be gzip compressed
2. Verbose Input File (.bedpe.gz)
Contains per-read detailed information
Includes:
mvec (methylation string: M/F pattern)
bs_seq (base sequence)


### Option 1: Using RStudio (recommended for testing)
```r
args <- c(
  "data/suppressed_merged_demo_S2_to_example_cobinding_spanning_lf_15_rf_15_extended_left_300_right_300.bedpe.gz",
  "data/suppressed_merged_demo_S2_to_example_cobinding_spanning_lf_15_rf_15_extended_left_300_right_300_verbose.bedpe.gz",
  "15", "15", "300", "300",
  "output/output_states.bedpe.gz",
  "output/output_verbose.bedpe.gz",
  "output/output_150bp.bedpe.gz"
)

commandArgs <- function(trailingOnly = TRUE) args

source("R/assign_cobinding_states.R")
```
### Output Interpretation

Each read is assigned a binding state:
0 → Naked DNA
1 → Transcription Factor (TF)
2 → Nucleosome
3 → Discard
Co-binding is encoded as:
(primary_state, secondary_state) → integer label (0–15)

Example:

(1,1) → TF-TF → 4
(2,2) → Nuc-Nuc → 8
Parameters
lflank   → left flank size (bp)
rflank   → right flank size (bp)
lextend  → left extension (bp)
rextend  → right extension (bp)

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

### Command2: For order_footprints_cobinding
```r
# 1. Clean up the environment and any old commandArgs overrides
rm(list = intersect(c("commandArgs", "args"), ls()))

# 2. Source the script to load its functions into RStudio
source("R/order_footprints_cobinding.R")

# 3. Call the function cleanly with your direct paths
order_footprints_cobinding(
  bedpe_gz = "output/output_150bp.bedpe.gz",
  roi_id   = "peak_110_4_and_peak_110_6",
  out_fp   = "output2/final_ordered_footprints.tsv",
  out_meth = "output2/final_ordered_methylation.tsv"
)
```

### Command 3: For plotting
```r
rm(list = intersect(c("commandArgs", "args"), ls()))

commandArgs <- function(trailingOnly = TRUE) {
    c(
        "output2/final_ordered_footprints110R.tsv",
        "output2/final_ordered_methylation110R.tsv",
        "mnase_peaks/peak_229.tsv",
        "dummy4",
        "dummy5",
        "plots/footprints110.pdf",
        "plots/methylation110.pdf",
        "dummy8",
        "dummy9",
        "150",
        "150",
        "35",
        "35",
        "dummy14",
        "dummy15",
        "dummy16",
        "peak_110_4_and_peak_110_6",
        "output/output_150bp.bedpe"
    )
}

source("R/plot_footprint_and_methylation.R")
```

### Command4: For extended_fragments_for_cobiding
```r
args <- c(
    "input_bed/suppressed_merged_demo_S2_to_example_cobinding_spanning_lf_15_rf_15.bedpe.gz",
    "15", "15", "300", "300",
    "ExtendedOutput/extended_fragments.bed.gz",
    "ExtendedOutput/extended_fragments.verbose.bed.gz"
)
commandArgs <- function(trailingOnly = TRUE) args
source("R/extended_fragments_for_cobinding.R")
```

### Command 5: For final method.r
```r
setwd("C:/Users/hp/Downloads/SMTrackR-main/SMTrackR-main")

# Install the package from source
# Step 1: Detach the installed library version
detach("package:SMTrackR", unload = TRUE)

# Step 2: Load directly from source (includes all your new files)
devtools::load_all("C:/Users/hp/Downloads/SMTrackR-main/SMTrackR-main/")

# Step 3: Confirm function exists
exists("plotCobindingFootprints")   # must return TRUE

# Step 4: Now run your call
result <- plotCobindingFootprints(
    chromosome = "chr2L",
    start1     = 19155158,   # 19155173 - 15
    end1       = 19155188,   
    start2     = 19155251,   
    end2       = 19155266,
    lflank     = 15L,
    rflank     = 15L,
    lextend    = 150L,
    rextend    = 150L,
    mnase_file = NULL,
    label      = "peak_110_4_and_peak_110_6",
    target_dir = "results/cobinding",
)
```

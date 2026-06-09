# =============================================================
# m6A Footprint Generator
# =============================================================
# USAGE:
#   Rscript generate_footprints.R <json_file> <bed_file> <fetch_start> <output_file>
#
# ARGUMENTS:
#   json_file    : sequence JSON file downloaded from UCSC API
#   bed_file     : BED12 file with m6A sites
#   fetch_start  : the start coordinate used in UCSC API (0-based)
#   output_file  : name of output file
#
# EXAMPLE:
#   Rscript generate_footprints.R chr2L.json chr2L.bed 33 chr2L_footprints.txt
# =============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  cat("Usage: Rscript generate_footprints.R <json_file> <bed_file> <fetch_start> <output_file>\n")
  cat("Example: Rscript generate_footprints.R chr2L.json chr2L.bed 33 chr2L_footprints.txt\n")
  quit(status = 1)
}

json_file    <- args[1]
bed_file     <- args[2]
fetch_start  <- as.integer(args[3])
output_file  <- args[4]

cat(sprintf("JSON file   : %s\n", json_file))
cat(sprintf("BED file    : %s\n", bed_file))
cat(sprintf("Fetch start : %d\n", fetch_start))
cat(sprintf("Output file : %s\n", output_file))

# ---- Read sequence from JSON --------------------------------
cat("\nReading sequence from JSON...\n")
raw <- paste(readLines(json_file), collapse = "")

# Extract dna field from JSON
dna_seq <- gsub('.*"dna"\\s*:\\s*"([^"]+)".*', '\\1', raw)
dna_seq <- toupper(trimws(dna_seq))

cat(sprintf("Sequence length: %d bp\n", nchar(dna_seq)))
cat(sprintf("Covers genome : %d to %d\n", fetch_start, fetch_start + nchar(dna_seq)))

# ---- Read BED file ------------------------------------------
cat("\nReading BED file...\n")
bed <- read.table(bed_file, header = FALSE, sep = "\t",
                  stringsAsFactors = FALSE, comment.char = "#")

colnames(bed)[1:12] <- c("chrom","start","end","name","score","strand",
                         "thickStart","thickEnd","rgb","blockCount",
                         "blockSizes","blockStarts")
bed <- bed[1:500, ]   # ← add this line after read.table
cat(sprintf("Found %d molecules\n", nrow(bed)))

# ---- Footprint function -------------------------------------
generate_footprint <- function(mol_start, mol_end, offsets, dna_seq, fetch_start) {
  
  idx_start <- mol_start - fetch_start + 1   # R is 1-based
  idx_end   <- mol_end   - fetch_start       # BED end is exclusive
  
  if (idx_start < 1 || idx_end > nchar(dna_seq)) {
    return(NA)
  }
  
  mol_seq <- toupper(substring(dna_seq, idx_start, idx_end))
  bases   <- strsplit(mol_seq, "")[[1]]
  
  # Default: A → F, everything else (T, C, G) → .
  footprint <- ifelse(bases == "A", "F", ".")
  
  # Methylated A (in offset list) → overwrite F with .
  for (off in offsets) {
    pos <- off + 1  # 0-based offset to 1-based R index
    if (pos >= 1 && pos <= length(footprint)) {
      footprint[pos] <- "."
    }
  }
  
  return(paste(footprint, collapse = ""))
}

# ---- Process all molecules ----------------------------------
cat("\nGenerating footprints...\n")

footprints <- character(nrow(bed))

for (i in seq_len(nrow(bed))) {
  offsets <- as.integer(strsplit(bed$blockStarts[i], ",")[[1]])
  offsets <- offsets[!is.na(offsets)]
  
  footprints[i] <- generate_footprint(
    mol_start   = bed$start[i],
    mol_end     = bed$end[i],
    offsets     = offsets,
    dna_seq     = dna_seq,
    fetch_start = fetch_start
  )
  
  if (i %% 1000 == 0) cat(sprintf("  Processed %d / %d\n", i, nrow(bed)))
}

# ---- Write output -------------------------------------------
results <- data.frame(
  chrom     = bed$chrom,
  start     = bed$start,
  end       = bed$end,
  name      = bed$name,
  score     = bed$score,
  strand    = bed$strand,
  n_m6A     = bed$blockCount,
  footprint = footprints,
  stringsAsFactors = FALSE
)

write.table(results, output_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

cat(sprintf("\nDone! Written to: %s\n", output_file))

# ---- Preview first 3 -------------------------------------
cat("\n--- Preview (first 3 molecules) ---\n")
for (i in 1:min(3, nrow(results))) {
  cat(sprintf("\nMolecule %d: %s:%d-%d  (%d m6A sites)\n",
              i, results$chrom[i], results$start[i],
              results$end[i], results$n_m6A[i]))
  fp <- results$footprint[i]
  if (!is.na(fp)) {
    cat(substr(fp, 1, 80), "...\n")
    f_count <- nchar(gsub("[^F]", "", fp))
    d_count <- nchar(gsub("[^.]", "", fp))
    cat(sprintf("  F (unmethylated A): %d  |  . (methylated A or non-A): %d\n",
                f_count, d_count))
  }
}

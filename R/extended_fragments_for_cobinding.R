#!/usr/bin/env Rscript
# R equivalent of extend_fragments_for_cobinding_bedpe.py
#
# ── Command line ─────────────────────────────────────────────────────────────
#   Rscript scripts/extend_fragments_for_cobinding_bedpe.R \
#       "input_bed/file.bedpe.gz" "15" "15" "300" "300" \
#       "ExtendedOutput/out.bed.gz" "ExtendedOutput/out_verbose.bed.gz"
#
# ── RStudio / interactive ─────────────────────────────────────────────────────
#   args <- c(
#     "input_bed/suppressed_merged_demo_S2_to_example_cobinding_spanning_lf_15_rf_15.bedpe.gz",
#     "15", "15", "300", "300",
#     "ExtendedOutput/extended_fragments.bed.gz",
#     "ExtendedOutput/extended_fragments.verbose.bed.gz"
#   )
#   commandArgs <- function(trailingOnly = TRUE) args
#   source("scripts/extend_fragments_for_cobinding_bedpe.R")
# ─────────────────────────────────────────────────────────────────────────────
# Args:
#   1  bedpe_gz      input .bedpe.gz file
#   2  lflank        left  flank used when building the file (e.g. 15)
#   3  rflank        right flank used when building the file (e.g. 15)
#   4  lextend       left  extension around peak center     (e.g. 300)
#   5  rextend       right extension around peak center     (e.g. 300)
#   6  out_file      output .bed.gz
#   7  out_verb_file output verbose .bed.gz

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop(paste(
    "Usage: Rscript extend_fragments_for_cobinding_bedpe.R",
    "<bedpe_gz> <lflank> <rflank> <lextend> <rextend> <out_file> <out_verb_file>"
  ))
}

bedpe_gz      <- args[1]
lflank        <- as.integer(args[2])
rflank        <- as.integer(args[3])
lextend       <- as.integer(args[4])
rextend       <- as.integer(args[5])
out_file      <- args[6]
out_verb_file <- args[7]

# Strip .gz for intermediate uncompressed files (mirrors shell ${6%.gz})
out_file_plain      <- sub("\\.gz$", "", out_file)
out_verb_file_plain <- sub("\\.gz$", "", out_verb_file)

# Ensure output directories exist
for (p in c(out_file_plain, out_verb_file_plain)) {
  d <- dirname(p)
  if (nzchar(d) && !dir.exists(d)) dir.create(d, recursive = TRUE)
}

out_con      <- file(out_file_plain,      open = "w")
out_verb_con <- file(out_verb_file_plain, open = "w")

# Open input: reads gz directly (works both from RStudio and command line)
con <- gzcon(file(bedpe_gz, open = "rb"))

desired_length <- rextend - (0L - lextend) + 1L

# Helper: Python-style 0-based half-open substring s[a:b]
pyslice <- function(s, a, b) {
  n <- nchar(s)
  if (a >= b || a >= n || b <= 0L) return("")
  substr(s, max(a, 0L) + 1L, min(b, n))
}

while (TRUE) {
  line <- readLines(con, n = 1L, warn = FALSE)
  if (length(line) == 0L) break
  
  fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
  
  # Field layout (1-based R index → Python 0-based d_loc column):
  #  1- 8  peak-pair annotation
  #  9      chr of read
  # 10      read start  (complete_footprint_start)
  # 11      read end    (complete_footprint_end)
  # 12      read name
  # 13      "."
  # 14      complete_footprint
  # 15      complete_mvec
  # 16      complete_bsseq
  cobinding_info           <- paste(fields[1:13], collapse = "\t")
  complete_footprint       <- fields[14]
  complete_mvec            <- fields[15]
  complete_bsseq           <- fields[16]
  complete_footprint_start <- as.integer(fields[10])
  complete_footprint_end   <- as.integer(fields[11])
  
  # Field 5 is like "chr2L:19155173-19155251"; peak_center = the start coord
  peak_center <- as.integer(sub(".*:(\\d+)-.*", "\\1", fields[5], perl = TRUE))
  strand      <- "+"   # hardcoded (same as original Python)
  
  extend_start <- peak_center - lextend        # closed
  extend_stop  <- peak_center + rextend + 1L   # open
  
  extended_footprint <- extended_mvec <- extended_bsseq <- ""
  
  # ---- four boundary conditions (identical logic to the Python script) ------
  
  if (extend_start >= complete_footprint_start &&
      extend_stop  <= complete_footprint_end + 1L) {
    # Condition 1: window fully inside read
    a <- extend_start - complete_footprint_start
    b <- extend_stop  - complete_footprint_start
    if (strand == "+") {
      extended_footprint <- pyslice(complete_footprint, a, b)
      extended_mvec      <- pyslice(complete_mvec,      a, b)
      extended_bsseq     <- pyslice(complete_bsseq,     a, b)
    } else {
      extended_footprint <- paste(rev(strsplit(pyslice(complete_footprint, a, b), "")[[1]]), collapse = "")
      extended_mvec      <- paste(rev(strsplit(pyslice(complete_mvec,      a, b), "")[[1]]), collapse = "")
      extended_bsseq     <- paste(rev(strsplit(pyslice(complete_bsseq,     a, b), "")[[1]]), collapse = "")
    }
    
  } else if (extend_start >= complete_footprint_start &&
             extend_stop  >  complete_footprint_end + 1L) {
    # Condition 2: right side spills beyond read
    a <- extend_start - complete_footprint_start
    b <- complete_footprint_end - complete_footprint_start + 1L
    if (strand == "+") {
      extended_footprint <- pyslice(complete_footprint, a, b)
      extended_mvec      <- pyslice(complete_mvec,      a, b)
      extended_bsseq     <- pyslice(complete_bsseq,     a, b)
      pad <- desired_length - nchar(extended_footprint)
      extended_footprint <- paste0(extended_footprint, strrep("M", pad))
      extended_mvec      <- paste0(extended_mvec,      strrep("M", pad))
      extended_bsseq     <- paste0(extended_bsseq,     strrep("M", pad))
    } else {
      fp  <- paste(rev(strsplit(pyslice(complete_footprint, a, b), "")[[1]]), collapse = "")
      mv  <- paste(rev(strsplit(pyslice(complete_mvec,      a, b), "")[[1]]), collapse = "")
      bs  <- paste(rev(strsplit(pyslice(complete_bsseq,     a, b), "")[[1]]), collapse = "")
      pad <- desired_length - nchar(fp)
      extended_footprint <- paste0(strrep("M", pad), fp)
      extended_mvec      <- paste0(strrep("M", pad), mv)
      extended_bsseq     <- paste0(strrep("M", pad), bs)
    }
    
  } else if (extend_start < complete_footprint_start &&
             extend_stop <= complete_footprint_end + 1L) {
    # Condition 3: left side spills before read start
    b <- extend_stop - complete_footprint_start
    if (strand == "+") {
      extended_footprint <- pyslice(complete_footprint, 0L, b)
      extended_mvec      <- pyslice(complete_mvec,      0L, b)
      extended_bsseq     <- pyslice(complete_bsseq,     0L, b)
      pad <- desired_length - nchar(extended_footprint)
      extended_footprint <- paste0(strrep("M", pad), extended_footprint)
      extended_mvec      <- paste0(strrep("M", pad), extended_mvec)
      extended_bsseq     <- paste0(strrep("M", pad), extended_bsseq)
    } else {
      fp  <- paste(rev(strsplit(pyslice(complete_footprint, 0L, b), "")[[1]]), collapse = "")
      mv  <- paste(rev(strsplit(pyslice(complete_mvec,      0L, b), "")[[1]]), collapse = "")
      bs  <- paste(rev(strsplit(pyslice(complete_bsseq,     0L, b), "")[[1]]), collapse = "")
      pad <- desired_length - nchar(fp)
      extended_footprint <- paste0(fp, strrep("M", pad))
      extended_mvec      <- paste0(mv, strrep("M", pad))
      extended_bsseq     <- paste0(bs, strrep("M", pad))
    }
    
  } else {
    # Condition 4: window wider than read on both sides
    left_pad  <- complete_footprint_start - extend_start
    right_pad <- extend_stop - (complete_footprint_end + 1L)
    if (strand == "+") {
      extended_footprint <- paste0(strrep("M", left_pad),  complete_footprint, strrep("M", right_pad))
      extended_mvec      <- paste0(strrep("M", left_pad),  complete_mvec,      strrep("M", right_pad))
      extended_bsseq     <- paste0(strrep("M", left_pad),  complete_bsseq,     strrep("M", right_pad))
    } else {
      rev_fp <- paste(rev(strsplit(complete_footprint, "")[[1]]), collapse = "")
      rev_mv <- paste(rev(strsplit(complete_mvec,      "")[[1]]), collapse = "")
      rev_bs <- paste(rev(strsplit(complete_bsseq,     "")[[1]]), collapse = "")
      extended_footprint <- paste0(strrep("M", right_pad), rev_fp, strrep("M", left_pad))
      extended_mvec      <- paste0(strrep("M", right_pad), rev_mv, strrep("M", left_pad))
      extended_bsseq     <- paste0(strrep("M", right_pad), rev_bs, strrep("M", left_pad))
    }
  }
  
  # Write outputs (mirrors Python: out_fp = footprint only; out_verb_fp = all three)
  writeLines(paste(cobinding_info, extended_footprint, sep = "\t"), con = out_con)
  writeLines(paste(cobinding_info, extended_footprint, sep = "\t"), con = out_verb_con)
  writeLines(paste(cobinding_info, extended_mvec,      sep = "\t"), con = out_verb_con)
  writeLines(paste(cobinding_info, extended_bsseq,     sep = "\t"), con = out_verb_con)
}

close(con)
close(out_con)
close(out_verb_con)

# Gzip the output files (mirrors: cat file | gzip - > file.gz && rm file)
for (pair in list(c(out_file_plain, out_file), c(out_verb_file_plain, out_verb_file))) {
  plain <- pair[1]; gz <- pair[2]
  if (plain != gz) {          # only compress if the names differ
    con_in  <- file(plain, open = "rb")
    con_out <- gzfile(gz,   open = "wb")
    repeat {
      buf <- readBin(con_in, what = "raw", n = 65536L)
      if (length(buf) == 0L) break
      writeBin(buf, con_out)
    }
    close(con_in)
    close(con_out)
    file.remove(plain)
  }
}

message("Done. Output written to:")
message("  ", out_file)
message("  ", out_verb_file)
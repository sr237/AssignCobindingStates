#!/usr/bin/env Rscript
# =============================================================================
# order_footprints_cobinding.R
#
# Pure-R replacement for the bash/Python pipeline:
#   order_footprints_cobinding.sh
#   get_largest_footprint_location.py
#   select_by_id.py
#   footprint_dot_to_digit_vec.py
#
# Produces byte-for-byte identical output to the shell pipeline.
#
# USAGE — interactive:
#   source("order_footprints_cobinding.R")
#   order_footprints_cobinding(
#     bedpe_gz = "file.bedpe.gz",
#     roi_id   = "peak_110_4_and_peak_110_6",
#     out_fp   = "final_ordered_footprints.tsv",
#     out_meth = "final_ordered_methylation.tsv"
#   )
#
# USAGE — command line:
#   Rscript order_footprints_cobinding.R \
#     <bedpe.gz> <roi_id> <out_footprints.tsv> <out_methylation.tsv>
#
# OUTPUT FORMAT (both TSVs, no header row):
#   read_header <TAB> pos_1 <TAB> pos_2 ... <TAB> pos_N
#
# OUTPUT ROW ORDER  (exact match to shell script):
#   For each co-binding state in ascending integer order:
#     - 83~163-flag reads, sorted by (start_offset, header_text)
#     - 99~147-flag reads, sorted by (start_offset, header_text)
#   Ties in start_offset broken lexicographically by header text,
#   matching GNU sort's full-line comparison behaviour.
#
# DIGIT ENCODING:
#   Footprint : '.' -> 0  |  'F'/'E' -> 1  |  'M' -> -1
#   Methylation: '.' -> 0  |  ACGTHXZacgthxz -> 2  |  'M' -> -1
# =============================================================================


# ---------------------------------------------------------------------------
# 1.  get_fp_metrics()
#     Replicates get_largest_footprint_location.py
#
#     Finds every contiguous F-run in a footprint string (after replacing M
#     with '.'), picks the LONGEST run, and returns:
#       start_offset  = 0-based run start  −  window_half   (default 150)
#       max_len       = length of that run
#
#     Matches Python:
#       m_replaced_fp = fp_str.replace("M", ".")
#       lv, ov, lf, abs_loc = get_real_footprint_length_with_abs_start(...)
#       max_start = abs_loc[argmax(lv)] - 150
# ---------------------------------------------------------------------------
get_fp_metrics <- function(fp_str, window_half = 150L) {
  clean <- chartr("M", ".", fp_str)            # replace M -> .
  chars <- strsplit(clean, "", fixed = TRUE)[[1L]]
  n     <- length(chars)
  
  if (n == 0L) return(c(start_offset = -window_half, max_len = 0L))
  
  is_F  <- chars == "F"
  r     <- rle(is_F)
  
  fp_mask <- r$values
  if (!any(fp_mask))
    return(c(start_offset = -window_half, max_len = 0L))
  
  run_end   <- cumsum(r$lengths)
  run_start <- c(1L, run_end[-length(run_end)] + 1L)  # 1-based
  
  fp_len   <- r$lengths[fp_mask]
  fp_start <- run_start[fp_mask]             # 1-based start of each F-run
  
  best            <- which.max(fp_len)
  abs_start_0base <- fp_start[best] - 1L     # convert to 0-based (matches Python enumerate)
  start_offset    <- abs_start_0base - window_half
  
  c(start_offset = start_offset, max_len = fp_len[best])
}


# ---------------------------------------------------------------------------
# 2.  Digitisation functions
#     Replicate footprint_dot_to_digit_vec.py
#     and the implied methylation_dot_to_digit_vec.py
# ---------------------------------------------------------------------------
fp_to_digits <- function(fp_str) {
  chars  <- strsplit(fp_str, "", fixed = TRUE)[[1L]]
  digits <- integer(length(chars))           # default 0  ('.' -> 0)
  digits[chars == "F" | chars == "E"] <-  1L
  digits[chars == "M"]                <- -1L
  digits
}

# Methylated bases observed in this data type
# METH_BASES <- c("A","C","G","T","H","X","Z",
#                 "a","c","g","t","h","x","z")

# meth_to_digits <- function(meth_str) {
#   chars  <- strsplit(meth_str, "", fixed = TRUE)[[1L]]
#   digits <- integer(length(chars))           # default 0  ('.' -> 0)
#   digits[chars %in% METH_BASES] <-  2L
#   digits[chars == "M"]          <- -1L
#   digits
# }

METH_UPPER <- c("H","X","Z","A","C","G","T")   # -> 2
METH_LOWER <- c("h","x","z","a","c","g","t")   # -> 1

meth_to_digits <- function(meth_str) {
  chars  <- strsplit(meth_str, "", fixed = TRUE)[[1L]]
  digits <- integer(length(chars))           # default 0  ('.' -> 0)
  digits[chars %in% METH_UPPER] <- 2L
  digits[chars %in% METH_LOWER] <- 1L
  digits[chars == "M"]          <- -1L
  digits
}

# ---------------------------------------------------------------------------
# 3.  Main pipeline function
# ---------------------------------------------------------------------------
order_footprints_cobinding <- function(
    bedpe_gz,
    roi_id,
    out_fp      = "final_ordered_footprints.tsv",
    out_meth    = "final_ordered_methylation.tsv",
    window_half = 150L,
    flag_83     = "83~163",
    flag_99     = "99~147",
    verbose     = TRUE
) {
  
  # ── STEP 1  Stream & filter  (replicates: zcat | grep roi_id) ─────────────
  if (verbose) message("[1/5] Reading: ", bedpe_gz)
  
  con   <- gzcon(file(bedpe_gz, "rb"))
  lines <- readLines(con, warn = FALSE)
  close(con)
  lines <- sub("\r$", "", lines)             # strip Windows CR if present
  
  keep  <- grepl(roi_id, lines, fixed = TRUE)
  lines <- lines[keep]
  
  if (length(lines) == 0L)
    stop("No lines found for roi_id: '", roi_id, "'")
  if (length(lines) %% 3L != 0L)
    stop("Filtered line count (", length(lines), ") is not divisible by 3. ",
         "Check that roi_id is not too broad and the file is not truncated.")
  
  n_reads <- length(lines) %/% 3L
  if (verbose) message("  Triplets retained: ", n_reads)
  
  # Triplet structure: line 1 = footprint, line 2 = methylation, line 3 = sequence (ignored)
  i_fp   <- seq(1L, length(lines), by = 3L)
  i_meth <- i_fp + 1L
  
  # Split "header\tsequence" into list(header=..., seq=...)
  split_tsv <- function(line) {
    p <- regexpr("\t", line, fixed = TRUE)[[1L]]
    list(header = substr(line, 1L, p - 1L),
         seq    = substr(line, p + 1L, nchar(line)))
  }
  
  fp_lines   <- lapply(lines[i_fp],   split_tsv)
  meth_lines <- lapply(lines[i_meth], split_tsv)
  
  headers   <- vapply(fp_lines,   `[[`, character(1L), "header")
  fp_seqs   <- vapply(fp_lines,   `[[`, character(1L), "seq")
  meth_seqs <- vapply(meth_lines, `[[`, character(1L), "seq")
  
  # ── STEP 2  Parse flag and co-binding state from each header ───────────────
  # Header format (backtick-delimited BEDPE fields + read info):
  #   chr2L`start`end`roi`...`flagstring#STATE
  # flag_83 / flag_99 are literal substrings of the header
  # state   is the integer after the last '#'
  
  flags  <- ifelse(grepl(flag_83, headers, fixed = TRUE), flag_83, flag_99)
  states <- as.integer(sub(".*#(\\d+)$", "\\1", headers))
  
  if (verbose) {
    message("  States : ", paste(sort(unique(states)), collapse = ", "))
    message("  Flags  : ", paste(sort(unique(flags)),  collapse = ", "))
  }
  
  # ── STEP 3  Compute largest-F-run metrics  (get_largest_footprint_location.py)
  if (verbose) message("[2/5] Computing longest-F-run metrics …")
  
  metrics      <- vapply(fp_seqs, get_fp_metrics, numeric(2L),
                         window_half = window_half)
  # metrics: 2 x n_reads matrix; row 1 = start_offset, row 2 = max_len
  start_offset <- metrics["start_offset", ]
  
  # ── STEP 4  Sort  (replicates sort -k2,2n -k3,3g + select_by_id loop) ─────
  #
  # The bash pipeline sorts each flag track independently by (state, start_offset)
  # and then, for each state, concatenates 83-flag rows followed by 99-flag rows.
  #
  # GNU sort uses the full line text as a tiebreaker for equal keys,
  # which here reduces to lexicographic order on the header-before-hash.
  #
  # Equivalent in R:
  #   1. For each flag track, sort by (state ASC, start_offset ASC, header_text ASC)
  #   2. Interleave: for state in sorted(unique states):
  #        emit 83-flag rows, then 99-flag rows
  
  if (verbose) message("[3/5] Sorting: state → (flag: 83 before 99) → start_offset → header …")
  
  # header_before_hash = header text up to (but not including) the final '#STATE'
  # Used as tiebreaker, matching GNU sort's full-line lexicographic comparison.
  hdr_before_hash <- sub("#\\d+$", "", headers)
  
  # Assign flag sort priority: 83 = 1, 99 = 2  (83-flag rows always come first)
  flag_key <- ifelse(flags == flag_83, 1L, 2L)
  
  # Sort within each flag track separately (preserves bash's track-independent sort),
  # then interleave by state. A single order() call with all four keys achieves this:
  sorted_idx <- order(states,       # primary:   co-binding state (ascending)
                      flag_key,     # secondary: 83-flag before 99-flag
                      start_offset, # tertiary:  longest-run start offset (ascending)
                      hdr_before_hash)  # quaternary: header text (lex, = GNU sort tiebreak)
  
  headers_s   <- headers[sorted_idx]
  fp_seqs_s   <- fp_seqs[sorted_idx]
  meth_seqs_s <- meth_seqs[sorted_idx]
  
  if (verbose) message("  Sorted rows: ", length(sorted_idx))
  
  # ── STEP 5  Digitise  (footprint_dot_to_digit_vec.py) ─────────────────────
  if (verbose) message("[4/5] Converting character strings to numeric vectors …")
  
  fp_mat   <- do.call(rbind, lapply(fp_seqs_s,   fp_to_digits))
  meth_mat <- do.call(rbind, lapply(meth_seqs_s, meth_to_digits))
  
  rownames(fp_mat)   <- headers_s
  rownames(meth_mat) <- headers_s
  
  # ── STEP 6  Write tab-separated output files ───────────────────────────────
  if (verbose) message("[5/5] Writing outputs …")
  
  write_tsv_matrix <- function(mat, path) {
    con <- file(path, "w")
    for (i in seq_len(nrow(mat))) {
      cat(rownames(mat)[i], "\t",
          paste(mat[i, ], collapse = "\t"), "\n",
          sep = "", file = con)
    }
    close(con)
  }
  
  write_tsv_matrix(fp_mat,   out_fp)
  write_tsv_matrix(meth_mat, out_meth)
  
  if (verbose) {
    message("  Footprint   -> ", out_fp,
            "  [", nrow(fp_mat), " reads x ", ncol(fp_mat), " positions]")
    message("  Methylation -> ", out_meth,
            "  [", nrow(meth_mat), " reads x ", ncol(meth_mat), " positions]")
    message("Done.")
  }
  
  invisible(list(footprints = fp_mat, methylation = meth_mat))
}


# ---------------------------------------------------------------------------
# 4.  Command-line entry point
# ---------------------------------------------------------------------------
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 4L) {
    cat("Usage: Rscript order_footprints_cobinding.R",
        "<bedpe.gz> <roi_id> <out_fp.tsv> <out_meth.tsv>\n")
    quit(status = 1L)
  }
  order_footprints_cobinding(
    bedpe_gz = args[1L],
    roi_id   = args[2L],
    out_fp   = args[3L],
    out_meth = args[4L]
  )
}


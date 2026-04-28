
# assign_cobinding_states.R
# SMTrackR package — assign cobinding states to DSMF read pairs
#
# SYNOPSIS (called by assign_cobinding_states.sh):
#   Rscript assign_cobinding_states.R \
#     <input.bedpe.gz>          \   # arg 1: main spanning reads (bedpe.gz)
#     <verbose.bedpe.gz>        \   # arg 2: verbose 3-line-per-read file (bedpe.gz)
#     <lflank>                  \   # arg 3: left flank (bp)
#     <rflank>                  \   # arg 4: right flank (bp)
#     <lextend>                 \   # arg 5: left extension (bp)
#     <rextend>                 \   # arg 6: right extension (bp)
#     <out_states.bedpe.gz>     \   # arg 7: main output (one line per read)
#     <out_verbose.bedpe.gz>    \   # arg 8: verbose output (3 lines per read)
#     <out_150bp.bedpe.gz>          # arg 9: 150bp-trimmed verbose output
#
# OUTPUT LABEL ENCODING:
#   Per-peak:  0=Naked  1=TF  2=Nucleosome  3=Discard
#   Combined pair_map (primary-secondary):
#     0-0→0   0-1→1   0-2→2
#     1-0→3   1-1→4   1-2→5
#     2-0→6   2-1→7   2-2→8
#     3-0→9   3-1→10  3-2→11
#     0-3→12  1-3→13  2-3→14
#     3-3→15



# SECTION 1 — helper functions

count_m_at_start_end <- function(fp_str) {
  chars <- strsplit(fp_str, "")[[1L]]
  n <- length(chars)
  
  start_count <- 0L
  for (i in seq_len(n)) {
    if (chars[i] == "M") start_count <- start_count + 1L else break
  }
  end_count <- 0L
  for (i in rev(seq_len(n))) {
    if (chars[i] == "M") end_count <- end_count + 1L else break
  }
  c(start = start_count, end = end_count)
}


# get_read_start_and_end_from_lex_rex_reads(fp_str)
# First and last non-M character positions, returned as 0-based indices
# matching Python list indexing.
# Python equivalent: from modules_for_cobinding import ...

get_read_start_and_end_from_lex_rex_reads <- function(fp_str) {
  chars <- strsplit(fp_str, "")[[1L]]
  non_m <- which(chars != "M")
  if (length(non_m) == 0L) return(c(start = NA_integer_, end = NA_integer_))
  # subtract 1: convert R 1-based → Python 0-based
  c(start = non_m[1L] - 1L, end = non_m[length(non_m)] - 1L)
}


# get_count_and_percentage_methylation(mvec_substr)
# Count upper (methylated) and lower (unmethylated) letters; return fraction.
# Python equivalent: from modules_for_cobinding import ...
# NOTE: computed in the main loop for completeness but not used in label logic
#       (Python script computes and discards these values too).

get_count_and_percentage_methylation <- function(mvec_substr) {
  chars <- strsplit(mvec_substr, "")[[1L]]
  upper <- sum(grepl("[A-LN-Z]", chars))   # methylated (uppercase, excluding M)
  lower <- sum(grepl("[a-z]",    chars))   # unmethylated (lowercase)
  total <- upper + lower
  per_c <- if (total > 0L) upper / total else 0.0
  list(per_c = per_c, total_upper = upper,
       total_lower = lower, total = total)
}


# get_per_orange_for_each_footprint(boundary_str)
# For each contiguous non-'.' run in the boundary substring, return the
# fraction of 'F' (orange) characters as a percentage.
# Python equivalent: from modules_for_cobinding import ...

get_per_orange_for_each_footprint <- function(boundary_str) {
  if (nchar(boundary_str) == 0L) return(numeric(0))
  
  chars    <- strsplit(boundary_str, "")[[1L]]
  result   <- numeric(0)
  in_run   <- FALSE
  run_chars <- character(0)
  
  for (ch in chars) {
    if (ch != ".") {
      in_run    <- TRUE
      run_chars <- c(run_chars, ch)
    } else {
      if (in_run) {
        result    <- c(result, sum(run_chars == "F") / length(run_chars) * 100)
        run_chars <- character(0)
        in_run    <- FALSE
      }
    }
  }
  if (in_run && length(run_chars) > 0L)
    result <- c(result, sum(run_chars == "F") / length(run_chars) * 100)
  
  result
}


# get_real_footprint_length_with_abs_start(boundary_str, lf, rf, full_m_replaced)
# For each contiguous non-'.' run: length, 0-based start within boundary,
# and 0-based start within the full M-replaced string.
# Python equivalent: from length_and_loc_with_absolute import ...

get_real_footprint_length_with_abs_start <- function(boundary_str, lf, rf,
                                                     full_m_replaced) {
  if (nchar(boundary_str) == 0L)
    return(list(flen_list = integer(0), dummy = NULL,
                loc_first = integer(0), abs_start = integer(0)))
  
  chars     <- strsplit(boundary_str, "")[[1L]]
  flen_list <- integer(0)
  loc_first <- integer(0)
  abs_start <- integer(0)
  in_run    <- FALSE
  run_start <- NA_integer_
  run_len   <- 0L
  
  for (i in seq_along(chars)) {
    if (chars[i] != ".") {
      if (!in_run) { in_run <- TRUE; run_start <- i - 1L; run_len <- 1L }
      else         { run_len <- run_len + 1L }
    } else {
      if (in_run) {
        flen_list <- c(flen_list, run_len)
        loc_first <- c(loc_first, run_start)
        abs_start <- c(abs_start, lf + run_start)   # lf is 0-based
        in_run <- FALSE; run_len <- 0L
      }
    }
  }
  if (in_run) {
    flen_list <- c(flen_list, run_len)
    loc_first <- c(loc_first, run_start)
    abs_start <- c(abs_start, lf + run_start)
  }
  
  list(flen_list = flen_list, dummy = NULL,
       loc_first = loc_first, abs_start = abs_start)
}



# SECTION 2 — Line-parsing helpers

# find_tab_positions(line)
# Return 0-based positions of all tab characters (matches Python re.finditer).
# gregexpr is 1-based → subtract 1.

find_tab_positions <- function(line) {
  hits <- gregexpr("\t", line, fixed = TRUE)[[1L]]
  if (hits[1L] == -1L) return(integer(0))
  hits - 1L
}


# build_read_id(line, d_loc)
# Construct the composite key used in Python as:
#   "`".join([ "`".join(line[0:d_loc[7]].split()), line[d_loc[10]+1:d_loc[11]] ])
# d_loc is 0-based; R substr is 1-based and inclusive.

build_read_id <- function(line, d_loc) {
  # First 7 fields (cols 0–6): slice up to the 8th tab (d_loc[8] in 1-based R)
  part1_raw  <- substr(line, 1L, d_loc[8L])          # d_loc[7] 0-based = d_loc[8] 1-based index
  part1_cols <- strsplit(trimws(part1_raw), "\\s+")[[1L]]
  part1      <- paste(part1_cols, collapse = "`")
  
  # Read-name field: between tab[10] and tab[11] (0-based)
  # 0-based d_loc[10]+1  →  1-based start = d_loc[11]+2  (R vector index 11 for 0-based [10])
  part2_start <- d_loc[11L] + 2L
  part2_end   <- d_loc[12L]          # 0-based d_loc[11]  →  1-based end = d_loc[12]
  part2       <- substr(line, part2_start, part2_end)
  
  paste(part1, part2, sep = "`")
}



# extract_last_field(line, d_loc)
# Python: line[d_loc[-1]+1 : len(line)-1]  (readLines already strips \n)

extract_last_field <- function(line, d_loc) {
  start_pos <- d_loc[length(d_loc)] + 2L   # 0-based last tab + 1 char → 1-based: +2
  substr(line, start_pos, nchar(line))
}



# SECTION 2b — Safe substring helper


# -----------------------------------------------------------------------------
# safe_substr(x, start, end)
# Bounds-checked substr: clamps indices, returns "" for NULL/empty/bad range.
# Prevents crashes on short strings or NULL mvec/bs_seq entries.
# -----------------------------------------------------------------------------
safe_substr <- function(x, start, end) {
  if (is.null(x) || nchar(x) == 0L) return("")
  start <- max(1L, start)
  end   <- min(nchar(x), end)
  if (start > end) return("")
  substr(x, start, end)
}



# SECTION 3 — Binding-state assignment

PAIR_MAP <- c(
  "0-0" =  0L, "0-1" =  1L, "0-2" =  2L,
  "1-0" =  3L, "1-1" =  4L, "1-2" =  5L,
  "2-0" =  6L, "2-1" =  7L, "2-2" =  8L,
  "3-0" =  9L, "3-1" = 10L, "3-2" = 11L,
  "0-3" = 12L, "1-3" = 13L, "2-3" = 14L,
  "3-3" = 15L
)



# assign_binding_labels(...)
# Mirrors the large if/elif/else block in the Python script exactly.
# Returns list(primary, secondary, cobinding).
assign_binding_labels <- function(primary_pct_orange,   secondary_pct_orange,
                                  flen_primary,          flen_secondary,
                                  abs_start_primary,     abs_start_secondary,
                                  read_start_abs,        read_end_abs,
                                  s_peak) {
  primary_label   <- NA_character_
  secondary_label <- NA_character_
  cobinding       <- FALSE
  
  has_p <- length(primary_pct_orange)   > 0L
  has_s <- length(secondary_pct_orange) > 0L
  
  # ── CASE A: footprints present in BOTH windows ──────────────────────────────
  if (has_p && has_s) {
    
    idx_p        <- which.max(primary_pct_orange)
    max_p        <- primary_pct_orange[idx_p]
    flen_p       <- flen_primary[idx_p]
    abs_p        <- abs_start_primary[idx_p]
    
    idx_s        <- which.max(secondary_pct_orange)
    max_s        <- secondary_pct_orange[idx_s]
    flen_s       <- flen_secondary[idx_s]
    abs_s        <- abs_start_secondary[idx_s]
    
    if (max_p > 30 && max_s > 30) {
      # Both peaks have orange signal > 30% ───────────────────────────────────
      if (abs_p == abs_s) {
        # Same footprint underlies both peaks
        is_edge <- (abs_p == read_start_abs) ||
          (abs_p + flen_p - 1L == read_end_abs)
        if (!is_edge) {
          if (s_peak <= 100L) {
            if (flen_p <= 100L) {
              cobinding       <- TRUE   # TF-TF co-binding on single footprint
              primary_label   <- "1"
              secondary_label <- "1"
            } else {
              primary_label   <- "2"   # long → Nucleosome
              secondary_label <- "2"
            }
          } else {
            # s_peak > 100 and same abs_start → footprint must span both → NUC
            primary_label   <- "2"
            secondary_label <- "2"
          }
        } else {
          # Footprint touches a read edge → edge artefact
          if (flen_p <= 50L) {
            primary_label   <- "3"   # Discard
            secondary_label <- "3"
          } else {
            primary_label   <- "2"   # NUC
            secondary_label <- "2"
          }
        }
        
      } else {
        # Different footprints for primary and secondary → assign independently
        
        # Primary
        is_edge_p <- (abs_p == read_start_abs) ||
          (abs_p + flen_p - 1L == read_end_abs)
        if (!is_edge_p) {
          primary_label <- if (flen_p <= 50L) "1" else "2"
        } else {
          primary_label <- if (flen_p >  50L) "2" else "3"
        }
        
        # Secondary
        is_edge_s <- (abs_s == read_start_abs) ||
          (abs_s + flen_s - 1L == read_end_abs)
        if (!is_edge_s) {
          secondary_label <- if (flen_s <= 50L) "1" else "2"
        } else {
          secondary_label <- if (flen_s >  50L) "2" else "3"
        }
      }
      
    } else if (max_p > 30 && max_s <= 30) {
      # Primary > 30%, Secondary ≤ 30% ────────────────────────────────────────
      # Secondary: Naked unless BOTH edge conditions are met 
      if ((abs_s != read_start_abs) || (abs_s + flen_s - 1L != read_end_abs)) {
        secondary_label <- "0"   # Naked
      } else {
        secondary_label <- "3"   # Discard
      }
      # Primary: TF / NUC / Discard
      if ((abs_p == read_start_abs) || (abs_p + flen_p - 1L == read_end_abs)) {
        primary_label <- if (flen_p > 50L) "2" else "3"
      } else {
        primary_label <- if (flen_p <= 50L) "1" else "2"
      }
      
    } else if (max_p <= 30 && max_s > 30) {
      # Primary ≤ 30%, Secondary > 30% ────────────────────────────────────────
      if ((abs_p != read_start_abs) || (abs_p + flen_p - 1L != read_end_abs)) {
        primary_label <- "0"
      } else {
        primary_label <- "3"
      }
      if ((abs_s == read_start_abs) || (abs_s + flen_s - 1L == read_end_abs)) {
        secondary_label <- if (flen_s > 50L) "2" else "3"
      } else {
        secondary_label <- if (flen_s <= 50L) "1" else "2"
      }
      
    } else {
      # Both ≤ 30%: Naked or Discard ──────────────────────────────────────────
      primary_label <- if ((abs_p != read_start_abs) ||
                           (abs_p + flen_p - 1L != read_end_abs)) "0" else "3"
      secondary_label <- if ((abs_s != read_start_abs) ||
                             (abs_s + flen_s - 1L != read_end_abs)) "0" else "3"
    }
    
    # ── CASE B: footprint ONLY in primary window ────────────────────────────────
  } else if (has_p && !has_s) {
    secondary_label <- "0"   # No footprint in secondary → Naked
    
    idx_p  <- which.max(primary_pct_orange)
    max_p  <- primary_pct_orange[idx_p]
    flen_p <- flen_primary[idx_p]
    abs_p  <- abs_start_primary[idx_p]
    
    if (max_p > 30) {
      is_edge <- (abs_p == read_start_abs) || (abs_p + flen_p - 1L == read_end_abs)
      if (is_edge) {
        primary_label <- if (flen_p > 50L) "2" else "3"
      } else if (flen_p > 50L) {
        primary_label <- "2"
      } else {
        primary_label <- "1"
      }
    } else {
      primary_label <- "0"
    }
    
    # ── CASE C: footprint ONLY in secondary window ──────────────────────────────
  } else if (!has_p && has_s) {
    primary_label <- "0"   # No footprint in primary → Naked
    
    idx_s  <- which.max(secondary_pct_orange)
    max_s  <- secondary_pct_orange[idx_s]
    flen_s <- flen_secondary[idx_s]
    abs_s  <- abs_start_secondary[idx_s]
    
    if (max_s > 30) {
      is_edge <- (abs_s == read_start_abs) || (abs_s + flen_s - 1L == read_end_abs)
      if (is_edge) {
        secondary_label <- if (flen_s > 50L) "2" else "3"
      } else if (flen_s > 50L) {
        secondary_label <- "2"
      } else {
        secondary_label <- "1"
      }
    } else {
      secondary_label <- "0"
    }
    
    # ── CASE D: no footprints in either window ──────────────────────────────────
  } else {
    primary_label   <- "0"
    secondary_label <- "0"
  }
  
  list(primary   = primary_label,
       secondary = secondary_label,
       cobinding = cobinding)
}


# SECTION 4 — Verbose-file parser

# open_input(path)
# Open a plain or gzip-compressed file for reading via readLines.
# Detects .gz extension automatically.

open_input <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzcon(gzfile(path, open = "rb"))
  } else {
    file(path, open = "r")
  }
}


# read_all_lines(path)
# Read all lines from a plain or .gz file.

read_all_lines <- function(path) {
  con <- open_input(path)
  on.exit(close(con))
  readLines(con, warn = FALSE)
}



# parse_verbose_file(verbose_file)
# Reads the 3-line-per-record verbose bedpe(.gz) into a named list:
#   verb_dict[[read_id]]  →  list(footprint=, mvec=, bs_seq=)
# Corresponds to the first loop in the Python script.

parse_verbose_file <- function(verbose_file) {
  message("  Reading verbose file: ", verbose_file)
  lines     <- read_all_lines(verbose_file)
  n_records <- length(lines) %/% 3L
  if (n_records == 0L) stop("Verbose file appears empty or malformed.")
  
  verb_dict  <- vector("list", n_records)
  names_vec  <- character(n_records)
  type_names <- c("footprint", "mvec", "bs_seq")
  
  for (i in seq_len(n_records)) {
    base <- (i - 1L) * 3L      # 0-based offset into lines vector
    
    # Key is the same for all 3 lines; compute from line 1
    line0  <- lines[base + 1L]
    d_loc0 <- find_tab_positions(line0)
    key    <- build_read_id(line0, d_loc0)
    
    rec        <- vector("list", 3L)
    names(rec) <- type_names
    for (j in 1:3) {
      raw    <- lines[base + j]
      d_loc  <- find_tab_positions(raw)
      rec[[j]] <- extract_last_field(raw, d_loc)
    }
    
    verb_dict[[i]] <- rec
    names_vec[i]   <- key
  }
  
  names(verb_dict) <- names_vec
  verb_dict
}



# SECTION 5 — Main processing function

# assign_cobinding_states_R(...)
# Core logic: reads input bedpe.gz + verbose bedpe.gz, assigns labels, writes
# three gzipped output files.

assign_cobinding_states_R <- function(input_file,
                                      verbose_file,
                                      lflank, rflank,
                                      lextend, rextend,
                                      output_file,
                                      output_verbose,
                                      output_150bp) {
  
  # ── 1. Parse verbose file ───────────────────────────────────────────────────
  verb_dict <- parse_verbose_file(verbose_file)
  
  # ── 2. Read main input ──────────────────────────────────────────────────────
  message("  Reading input file:   ", input_file)
  input_lines <- read_all_lines(input_file)
  input_lines <- input_lines[nchar(trimws(input_lines)) > 0L]   # drop blank lines
  message("  Records to process:   ", length(input_lines))
  
  # ── 3. Open gzip output connections ────────────────────────────────────────
  con_fp      <- gzfile(output_file,    "wb")
  con_verbose <- gzfile(output_verbose, "wb")
  con_150     <- gzfile(output_150bp,   "wb")
  on.exit({
    close(con_fp)
    close(con_verbose)
    close(con_150)
  }, add = TRUE)
  
  # ── 4. Process each read ────────────────────────────────────────────────────
  summary_rows <- vector("list", length(input_lines))
  
  for (line_idx in seq_along(input_lines)) {
    line  <- input_lines[line_idx]
    d_loc <- find_tab_positions(line)
    
    # s_peak: column 6 (0-based) — between tab[6] and tab[7]
    # Python: int(line[d_loc[6]+1 : d_loc[7]])
    # R d_loc is 0-based vector; element [7] = 0-based tab 6 → 1-based index 7
    s_peak <- suppressWarnings(as.integer(substr(line, d_loc[7L] + 2L, d_loc[8L])))
    if (is.na(s_peak)) next
    
    read_id <- build_read_id(line, d_loc)
    fp_str  <- extract_last_field(line, d_loc)
    
    # ── Boundary indices (0-based, matching Python) ──────────────────────────
    primary_lf   <- lextend - lflank
    primary_rf   <- lextend + rflank + 1L
    secondary_lf <- lextend + s_peak - lflank
    secondary_rf <- lextend + s_peak + rflank + 1L
    
    # ── Read extent (0-based) ────────────────────────────────────────────────
    read_bounds    <- get_read_start_and_end_from_lex_rex_reads(fp_str)
    if (is.na(read_bounds["start"]) || is.na(read_bounds["end"])) next
    read_start_abs <- read_bounds["start"]
    read_end_abs   <- read_bounds["end"]
    
    # ── M → '.' replacement & boundary substrings ───────────────────────────
    # Python: m_replaced_str = fp_str.replace("M", ".")
    # Python slice [a:b] → R substr(s, a+1, b)
    m_replaced_str <- gsub("M", ".", fp_str, fixed = TRUE)
    primary_str    <- substr(m_replaced_str, primary_lf   + 1L, primary_rf)
    secondary_str  <- substr(m_replaced_str, secondary_lf + 1L, secondary_rf)
    
    # ── Footprint lengths & absolute starts ──────────────────────────────────
    res_p <- get_real_footprint_length_with_abs_start(
      primary_str,   primary_lf,   primary_rf,   m_replaced_str)
    res_s <- get_real_footprint_length_with_abs_start(
      secondary_str, secondary_lf, secondary_rf, m_replaced_str)
    
    # ── % orange per footprint ───────────────────────────────────────────────
    primary_pct_orange   <- get_per_orange_for_each_footprint(primary_str)
    secondary_pct_orange <- get_per_orange_for_each_footprint(secondary_str)
    
    # ── Methylation (computed to mirror Python; not used in label logic) ─────
    rec <- verb_dict[[read_id]]
    
    if (is.null(rec)) {
      warning("Missing read_id in verbose file: ", read_id, call. = FALSE)
      next
    }
    
    mvec_str <- rec[["mvec"]]
    bs_str   <- rec[["bs_seq"]]
    
    if (!is.null(mvec_str)) {
      get_count_and_percentage_methylation(
        safe_substr(mvec_str, primary_lf - 15L + 1L, primary_rf + 1L))
      get_count_and_percentage_methylation(
        safe_substr(mvec_str, primary_lf + 1L,        primary_rf + 15L + 1L))
      get_count_and_percentage_methylation(
        safe_substr(mvec_str, secondary_lf - 15L + 1L, secondary_rf + 1L))
      get_count_and_percentage_methylation(
        safe_substr(mvec_str, secondary_lf + 1L,        secondary_rf + 15L + 1L))
    }
    
    # ── Assign labels ───
    labels <- assign_binding_labels(
      primary_pct_orange   = primary_pct_orange,
      secondary_pct_orange = secondary_pct_orange,
      flen_primary         = res_p$flen_list,
      flen_secondary       = res_s$flen_list,
      abs_start_primary    = res_p$abs_start,
      abs_start_secondary  = res_s$abs_start,
      read_start_abs       = read_start_abs,
      read_end_abs         = read_end_abs,
      s_peak               = s_peak
    )
    
    label_comb  <- paste0(labels$primary, "-", labels$secondary)
    one_of_nine <- PAIR_MAP[[label_comb]]
    if (is.null(one_of_nine)) {
      warning("Unknown label combination '", label_comb,
              "' for read: ", read_id, call. = FALSE)
      next
    }
    
    # ── Build output strings ──────
    key_tag  <- paste0(read_id, "#", one_of_nine)
    mvec_out <- if (!is.null(rec[["mvec"]]))   rec[["mvec"]]   else ""
    bs_out   <- if (!is.null(rec[["bs_seq"]])) rec[["bs_seq"]] else ""
    
    fp_line   <- paste0(key_tag, "\t", fp_str)
    mvec_line <- paste0(key_tag, "\t", mvec_out)
    bs_line   <- paste0(key_tag, "\t", bs_out)
    
    # ── 150bp centre trim ────────────────────────────────────────────────────
    # Python: int((len(modified_fp_str) - 1) / 2)  → 0-based centre
    fp_len     <- nchar(fp_str)
    centre     <- as.integer((fp_len - 1L) / 2L)   # 0-based
    t_start    <- centre - 150L                      # 0-based inclusive
    t_end      <- centre + 150L                      # 0-based inclusive
    # R substr is 1-based inclusive: add 1 to both
    fp_150_str   <- safe_substr(fp_str,   t_start + 1L, t_end + 1L)
    mvec_150_str <- safe_substr(mvec_out, t_start + 1L, t_end + 1L)
    bs_150_str   <- safe_substr(bs_out,   t_start + 1L, t_end + 1L)
    
    fp_line_150   <- paste0(key_tag, "\t", fp_150_str)
    mvec_line_150 <- paste0(key_tag, "\t", mvec_150_str)
    bs_line_150   <- paste0(key_tag, "\t", bs_150_str)
    
    # ── Write to output connections ───
    writeLines(fp_line,                    con_fp)        # main output
    
    writeLines(c(fp_line, mvec_line, bs_line),            con_verbose)  # verbose
    
    writeLines(c(fp_line_150, mvec_line_150, bs_line_150), con_150)     # 150bp
    
    summary_rows[[line_idx]] <- data.frame(
      read_id         = read_id,
      primary_label   = labels$primary,
      secondary_label = labels$secondary,
      cobinding_label = one_of_nine,
      cobinding       = labels$cobinding,
      stringsAsFactors = FALSE
    )
  }
  
  # ── 5. Report
  message("  Written: ", output_file)
  message("  Written: ", output_verbose)
  message("  Written: ", output_150bp)
  
  invisible(do.call(rbind, Filter(Negate(is.null), summary_rows)))
}


# SECTION 6 — CLI entry point (invoked by assign_cobinding_states.sh)


# commandArgs(trailingOnly = TRUE) receives the 9 positional arguments passed
# by the shell wrapper in this order:
#   1  input.bedpe.gz         main spanning reads
#   2  verbose.bedpe.gz       3-line-per-read verbose file
#   3  lflank                 integer
#   4  rflank                 integer
#   5  lextend                integer
#   6  rextend                integer
#   7  out_states.bedpe.gz    main output
#   8  out_verbose.bedpe.gz   verbose output
#   9  out_150bp.bedpe.gz     150bp-trimmed output

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 9L) {
  cat("Usage: Rscript assign_cobinding_states.R",
      "<input.bedpe.gz> <verbose.bedpe.gz>",
      "<lflank> <rflank> <lextend> <rextend>",
      "<out.bedpe.gz> <out_verbose.bedpe.gz> <out_150bp.bedpe.gz>\n")
  quit(status = 1L)
}

input_file    <- args[1L]
verbose_file  <- args[2L]
lflank        <- as.integer(args[3L])
rflank        <- as.integer(args[4L])
lextend       <- as.integer(args[5L])
rextend       <- as.integer(args[6L])
output_file   <- args[7L]
output_verbose <- args[8L]
output_150bp  <- args[9L]

# Create output directory if needed
for (out_path in c(output_file, output_verbose, output_150bp)) {
  out_dir <- dirname(out_path)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
}

message("=== assign_cobinding_states ===")
message("  input        : ", input_file)
message("  verbose      : ", verbose_file)
message("  lflank/rflank: ", lflank, " / ", rflank)
message("  lextend/rext : ", lextend, " / ", rextend)

assign_cobinding_states_R(
  input_file    = input_file,
  verbose_file  = verbose_file,
  lflank        = lflank,
  rflank        = rflank,
  lextend       = lextend,
  rextend       = rextend,
  output_file   = output_file,
  output_verbose = output_verbose,
  output_150bp  = output_150bp
)

message("=== Done ===")
